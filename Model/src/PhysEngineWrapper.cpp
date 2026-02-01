#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"

#include "../include/PhysEngineWrapper.h"

using namespace std;
using namespace std::chrono;

using StepEnergyTrace = EngineCore::StepEnergyTrace;

static inline double safe_double_from_high_prec(const high_prec& v) {
    try { return v.convert_to<double>(); }
    catch (...) { return std::numeric_limits<double>::infinity(); }
}
static inline double to_double(const high_prec& v) { return safe_double_from_high_prec(v); }

static inline high_prec abs_hp(high_prec v) { return v < high_prec(0) ? -v : v; }

static inline high_prec safe_rel_error(high_prec post, high_prec pre) {
    high_prec diff = post - pre;
    if (pre == high_prec(0)) return diff;
    return diff / pre;
}

// =======================
// OPTIONAL STEP REPORTING
// =======================
static const bool kReportEnergyPerStep = false;
static const int  kStepPauseInterval   = 100;

static inline void pause_for_user() {
    std::cout << "Press enter to continue..." << std::endl;
    std::string _;
    std::getline(std::cin, _);
}

static inline void print_step_report(
    int step_idx,
    int steps_total,
    const StepEnergyTrace& tr,
    double step_seconds)
{
    const int width = 104;
    auto line = [&](char ch) {
        for (int i = 0; i < width; ++i) std::cout << ch;
        std::cout << "\n";
    };

    int pct = 0;
    if (steps_total > 0) {
        pct = (int)(((long long)(step_idx + 1) * 100LL) / (long long)steps_total);
        if (pct > 100) pct = 100;
        if (pct < 0) pct = 0;
    }

    std::cout << std::fixed;

    line('=');
    std::cout << "Step " << step_idx << " / " << (steps_total > 0 ? (steps_total - 1) : 0)
              << "  (" << pct << "%)   "
              << "dt=" << std::setprecision(6) << to_double(EngineCore::dt)
              << "   step_time=" << std::setprecision(4) << step_seconds << " s\n";

    line('-');
    std::cout << std::setprecision(12);
    std::cout << "Energies:\n";
    std::cout << "  TE = " << to_double(tr.TE)
              << "   KE = " << to_double(tr.KE)
              << "   PE = " << to_double(tr.PE)
              << "   HE = " << to_double(tr.HE) << "\n";

    line('-');
    std::cout << "Substep energy drift (TE checkpoints):\n";
    std::cout << "  Overlaps  : TE0 -> TE1   dE = " << to_double(tr.dE_overlap)
              << "   rel = " << to_double(tr.rel_overlap) << "\n";
    std::cout << "  Collisions: TE1 -> TE2   dE = " << to_double(tr.dE_collision)
              << "   rel = " << to_double(tr.rel_collision) << "\n";
    std::cout << "  Verlet    : TE2 -> TE3   dE = " << to_double(tr.dE_verlet)
              << "   rel = " << to_double(tr.rel_verlet) << "\n";

    line('-');
    std::cout << "TE checkpoints:\n";
    std::cout << "  TE0(before overlaps) = " << to_double(tr.TE0) << "\n";
    std::cout << "  TE1(after overlaps ) = " << to_double(tr.TE1) << "\n";
    std::cout << "  TE2(after collisions) = " << to_double(tr.TE2) << "\n";
    std::cout << "  TE3(after verlet   ) = " << to_double(tr.TE3) << "\n";
    line('=');
}


// Equality threshold for comparing CURRENT vs BENCHMARK error (relative difference), default 0.05%.
static const high_prec kBenchmarkEqualRelThreshFrac = high_prec("0.0005"); // 0.05% = 0.0005 fraction
// If benchmark is exactly 0, treat "equal" as absolute current error <= 0.05% (in percent units)
static const high_prec kBenchmarkEqualAbsThreshPct  = high_prec("0.05");

// Benchmarks stored as TOTAL RELATIVE TE ERROR IN PERCENT (|ΔTE|/|TE0| * 100).
// Scenario keys match scenario->name exactly (as before splitting).
static const std::unordered_map<std::string, high_prec> kBenchmarkTotalTEErrorPct = {
    {"Planet + Moon System", high_prec("0.000012043152104")},
    {"Planet + Moon System_short",  high_prec("0.001522")},
    {"Planet + Moon System_shorter", high_prec("0.000000092628923")}
};

static const std::unordered_map<std::string, high_prec> kBenchmarkSimTimeSec = {
    {"Planet + Moon System", high_prec("4.513e+04")},
    {"Planet + Moon System_short",  high_prec("1292.951344")},
    {"Planet + Moon System_shorter", high_prec("184.122661")}
};

static inline bool lookup_benchmark_te_error_pct(const std::string& scenario_name, high_prec& out_benchmark_pct) {
    auto it = kBenchmarkTotalTEErrorPct.find(scenario_name);
    if (it == kBenchmarkTotalTEErrorPct.end()) return false;
    out_benchmark_pct = it->second;
    return true;
}

static inline bool lookup_benchmark_sim_time_sec(const std::string& scenario_name, high_prec& out_benchmark_sec) {
    auto it = kBenchmarkSimTimeSec.find(scenario_name);
    if (it == kBenchmarkSimTimeSec.end()) return false;
    out_benchmark_sec = it->second;
    return true;
}

static inline void print_benchmark_comparison_report(
    const std::string& scenario_name,
    const high_prec& benchmark_err_pct,
    const high_prec& current_err_pct,
    const high_prec& equal_rel_thresh_frac,
    const high_prec& equal_abs_thresh_pct,
    const high_prec& benchmark_time_sec,
    const high_prec& current_time_sec,
    bool have_time_benchmark
) {
    const int width = 104;
    auto line = [&](char ch) {
        for (int i = 0; i < width; ++i) std::cout << ch;
        std::cout << "\n";
    };

    std::string verdict;
    bool have_rel = (benchmark_err_pct != high_prec(0));
    high_prec rel_diff_frac = 0;   // (current - benchmark) / benchmark
    high_prec rel_diff_pct  = 0;   // rel_diff_frac * 100

    if (have_rel) {
        rel_diff_frac = (current_err_pct - benchmark_err_pct) / benchmark_err_pct;
        rel_diff_pct  = rel_diff_frac * high_prec(100);
        if (abs_hp(rel_diff_frac) <= equal_rel_thresh_frac) verdict = "EQUAL ERROR";
        else if (current_err_pct > benchmark_err_pct)       verdict = "WORSE ERROR";
        else                                                verdict = "BETTER ERROR";
    } else {
        if (abs_hp(current_err_pct) <= equal_abs_thresh_pct) verdict = "EQUAL ERROR";
        else if (current_err_pct > high_prec(0))             verdict = "WORSE ERROR";
        else                                                 verdict = "BETTER ERROR";
    }

    line('#');
    std::cout << "ENERGY BENCHMARK\n";
    line('#');

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Scenario               : " << scenario_name << "\n";
    std::cout << "Benchmark error (%)    : " << benchmark_err_pct << "\n";
    std::cout << "Current error (%)      : " << current_err_pct << "\n";

    if (have_rel) {
        std::cout << "Relative diff vs bench : " << rel_diff_pct << " %"
                  << "   (equal if |diff| <= " << (equal_rel_thresh_frac * high_prec(100)) << " %)\n";
    } else {
        std::cout << "Relative diff vs bench : N/A (benchmark == 0)\n";
        std::cout << "Equal threshold (abs)  : " << equal_abs_thresh_pct << " %\n";
    }

    line('-');
    std::cout << "Result                : " << verdict << "\n";
    line('#');

    cout << endl;

    line('#');
    cout << "SIM TIME BENCHMARK\n";
    line('#');
    std::cout << "Scenario               : " << scenario_name << "\n";
    if (have_time_benchmark) {
        std::cout << "Benchmark time (s)     : " << benchmark_time_sec << "\n";
        std::cout << "Current time (s)       : " << current_time_sec << "\n";
        std::cout << "Relative diff vs bench : "
                  << ((current_time_sec - benchmark_time_sec) / benchmark_time_sec * high_prec(100)) << " %\n";


    } else {
        std::cout << "No benchmark entry found for sim time.\n";
        std::cout << "Current time (s)       : " << current_time_sec << "\n";
    }
    line('#');
    cout << endl;
}

// =======================
// RUN
// =======================
shared_ptr<snapshots> Engine::run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles)
{

    //set core physics parameters
    core.set_overlap_beta(scenario);


    // Reset run-level counters (wrapper-side)
    high_prec total_TE_error_overlap   = 0;
    high_prec total_TE_error_collision = 0;
    high_prec total_TE_error_verlet    = 0;
    EngineCore::collisions = 0;

    cout << "Engine (Wrapper) is initialized." << endl;

    shared_ptr<snapshots> particle_states = make_shared<snapshots>();
    cout << scenario->name << "'s initial states loaded" << endl << endl;

    cout << "Press enter to start the simulation." << endl;
    cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    // dt is configured in core (restored single call point)
    core.configure_dt(scenario, particles);

    high_prec steps_db = scenario->time / scenario->dt;
    int steps = static_cast<int>(steps_db);
    cout << "Number of steps to be simulated: " << steps << endl << endl;

    particle_states->metrics.resize(steps, nullptr);

    // Initial TE baseline for final reporting
    high_prec TE_initial = core.calc_TE(particles);
    high_prec TE_final   = TE_initial;

    EngineCore::clock_t::time_point sim_start = EngineCore::clock_t::now();
    int last_reported_pct = -5;

    for (int i = 0; i < steps; i++) {
        if (!particle_states->metrics[i]) {
            particle_states->metrics[i] = make_shared<test_metrics_t>();
        }

        EngineCore::clock_t::time_point step_start = EngineCore::clock_t::now();
        EngineCore::update_iter = i;

        StepEnergyTrace tr;
        core.step(particles, &tr);

        // Accumulate stage-wise TE deltas (signed)
        total_TE_error_overlap   += tr.dE_overlap;
        total_TE_error_collision += tr.dE_collision;
        total_TE_error_verlet    += tr.dE_verlet;
        TE_final = tr.TE3;

        EngineCore::clock_t::time_point step_end = EngineCore::clock_t::now();
        double step_seconds = EngineCore::seconds_between(step_start, step_end);

        if (kReportEnergyPerStep) {
            print_step_report(i, steps, tr, step_seconds);
            if (kStepPauseInterval > 0 && ((i + 1) % kStepPauseInterval == 0)) {
                pause_for_user();
            }
        }

        auto particles_copy = make_unique<Particles>(*particles);
        particle_states->snaps.push_back(move(particles_copy));

        particle_states->metrics[i]->margin_TE_error = core.get_margin_TE_error();
        //particle_states->metrics[i]->margin_TE_error_overlap   = core.get_margin_TE_error_overlap();
        particle_states->metrics[i]->margin_TE_error_collision = core.get_margin_TE_error_collision();
        particle_states->metrics[i]->margin_TE_error_integrate = core.get_margin_TE_error_integrate();
        //particle_states->metrics[i]->overlap_iters_in_step = EngineCore::overlap_iter;
        //particle_states->metrics[i]->margin_TE_error_overlap_ij_transl = core.get_margin_TE_error_overlap_ij_transl();
        //particle_states->metrics[i]->margin_TE_error_overlap_ij_corrected = core.get_margin_TE_error_overlap_ij_corrected();

        // progress reporting: 5% buckets
        if (steps > 0) {
            if (i == 0 && last_reported_pct < 0) {
                cout << "0% of the simulation complete." << endl;
                last_reported_pct = 0;
            }

            int pct_complete = (int)(((long long)(i + 1) * 100LL) / (long long)steps);
            if (pct_complete > 100) pct_complete = 100;
            if (pct_complete < 0) pct_complete = 0;

            int pct_bucket = (pct_complete / 5) * 5;
            if (pct_bucket >= 0 && pct_bucket <= 100 && pct_bucket > last_reported_pct) {
                cout << pct_bucket << "% of the simulation complete." << endl;
                last_reported_pct = pct_bucket;
            }

            if (i == steps - 1 && last_reported_pct < 100) {
                cout << "100% of the simulation complete." << endl;
                last_reported_pct = 100;
            }
        }

        duration<double> time_taken = duration<double>(step_end - step_start);
        particle_states->metrics[i]->fps = 1 / time_taken.count();
    }

    EngineCore::clock_t::time_point sim_end = EngineCore::clock_t::now();
    double total_seconds = EngineCore::seconds_between(sim_start, sim_end);
    high_prec current_time_sec = high_prec(total_seconds);

    // =======================
    // END-OF-SIM LOGGING (RESTORED)
    // =======================
    high_prec total_stage_sum = total_TE_error_overlap + total_TE_error_collision + total_TE_error_verlet;
    high_prec abs_stage_sum = abs_hp(total_TE_error_overlap) + abs_hp(total_TE_error_collision) + abs_hp(total_TE_error_verlet);

    auto pct_share = [&](high_prec part_abs) -> high_prec {
        if (abs_stage_sum == high_prec(0)) return high_prec(0);
        return (part_abs / abs_stage_sum) * high_prec(100);
    };

    high_prec dE_total = TE_final - TE_initial;
    high_prec dE_total_abs = abs_hp(dE_total);
    high_prec dE_total_rel = abs_hp(safe_rel_error(TE_final, TE_initial));
    high_prec dE_total_rel_pct = dE_total_rel * high_prec(100); // total relative TE error (%)

    cout << scenario->name << " simulation completed." << endl << endl;
    cout << "Steps: " << steps << endl;
    cout << "Total collisions: " << EngineCore::collisions << endl;
    cout << "Total sim time (s): " << std::setprecision(4) << std::defaultfloat << total_seconds << endl;

    cout << std::fixed << std::setprecision(15);

    // legacy-style totals
    cout << "Total TE error (sum of parts): " << total_stage_sum << endl;
    if (total_stage_sum != high_prec(0)) {
        cout << "% overlap error:   " << (total_TE_error_overlap   / total_stage_sum) * high_prec(100) << endl;
        cout << "% collision error: " << (total_TE_error_collision / total_stage_sum) * high_prec(100) << endl;
        cout << "% verlet error:    " << (total_TE_error_verlet    / total_stage_sum) * high_prec(100) << endl;
    }

    // stage signed deltas + absolute share
    cout << "Stage TE error (signed delta E):" << endl;
    cout << "  overlaps:   " << total_TE_error_overlap
         << " (abs share " << pct_share(abs_hp(total_TE_error_overlap)) << "%)" << endl;
    cout << "  collisions: " << total_TE_error_collision
         << " (abs share " << pct_share(abs_hp(total_TE_error_collision)) << "%)" << endl;
    cout << "  verlet:     " << total_TE_error_verlet
         << " (abs share " << pct_share(abs_hp(total_TE_error_verlet)) << "%)" << endl;

    cout << "Stage TE error sum (signed): " << total_stage_sum << endl;

    cout << "Total TE error vs initial:" << endl;
    cout << "  initial TE: " << TE_initial << endl;
    cout << "  final TE:   " << TE_final << endl;
    cout << "  delta TE:   " << dE_total << endl;
    cout << "  |delta TE|: " << dE_total_abs << endl;
    cout << "  rel delta:  " << dE_total_rel << " (" << dE_total_rel_pct << "%)" << endl;

    if (dE_total > high_prec(0))      cout << "ENERGY INCREASED" << endl;
    else if (dE_total < high_prec(0)) cout << "ENERGY DECREASED" << endl;
    else                              cout << "ENERGY UNCHANGED" << endl;

    // =======================
    // BENCHMARK COMPARISON (RESTORED)
    // =======================
    {
        high_prec bench_pct = 0;
        high_prec bench_time = 0;
        bool have_time_bench = lookup_benchmark_sim_time_sec(scenario->name, bench_time);

        if (lookup_benchmark_te_error_pct(scenario->name, bench_pct)) {
            print_benchmark_comparison_report(
                scenario->name,
                bench_pct,
                dE_total_rel_pct,
                kBenchmarkEqualRelThreshFrac,
                kBenchmarkEqualAbsThreshPct,
                bench_time,
                current_time_sec,
                have_time_bench
            );
        } else {
            const int width = 104;
            auto line = [&](char ch) {
                for (int i = 0; i < width; ++i) std::cout << ch;
                std::cout << "\n";
            };

            line('#');
            std::cout << "BENCHMARK COMPARISON\n";
            line('#');
            std::cout << "No benchmark entry found for scenario: " << scenario->name << "\n";
            std::cout << "Current total relative TE error (%): " << std::fixed << std::setprecision(6) << dE_total_rel_pct << "\n";
            std::cout << "Current sim time (s): " << std::fixed << std::setprecision(6) << to_double(current_time_sec) << "\n";
            line('#');
        }
    }

    run_to_cache(scenario, particle_states);
    return particle_states;
}

// =======================
// CACHE I/O (WRAPPER)
// =======================
void Engine::run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {
    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    ofstream file(file_name);
    if (!file.is_open()) {
        cout << "Failed to write snapshots to cache!" << endl;
        return;
    }

    file << "step_id, particle_id,r,g,b,x,y,z,vx,vy,vz,m,rad,temp,rest" << endl;

    for (int i = 0; i < (int)particle_states->snaps.size(); i++) {
        for (int j = 0; j < (int)particle_states->snaps[i]->particle_list.size(); j++) {
            file << std::fixed << std::setprecision(15)
                << i << ","
                << particle_states->snaps[i]->particle_list[j]->particle_id << ","
                << (particle_states->snaps[i]->particle_list[j]->r) << ","
                << (particle_states->snaps[i]->particle_list[j]->g) << ","
                << (particle_states->snaps[i]->particle_list[j]->b) << ","
                << (particle_states->snaps[i]->particle_list[j]->x) << ","
                << (particle_states->snaps[i]->particle_list[j]->y) << ","
                << (particle_states->snaps[i]->particle_list[j]->z) << ","
                << (particle_states->snaps[i]->particle_list[j]->vx) << ","
                << (particle_states->snaps[i]->particle_list[j]->vy) << ","
                << (particle_states->snaps[i]->particle_list[j]->vz) << ","
                << (particle_states->snaps[i]->particle_list[j]->m) << ","
                << (particle_states->snaps[i]->particle_list[j]->rad) << ","
                << (particle_states->snaps[i]->particle_list[j]->temp) << ","
                << (particle_states->snaps[i]->particle_list[j]->rest) << endl;
        }
    }

    file.close();
    cout << "Snapshots saved to " << file_name << endl;
}

bool Engine::cache_exists(shared_ptr<scenario> scenario) {
    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";
    ifstream file(file_name);
    return file.is_open();
}

shared_ptr<snapshots> Engine::run_from_cache(shared_ptr<scenario> scenario) {
    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    ifstream file(file_name);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open cache file: " + file_name);
    }

    shared_ptr<snapshots> particle_states = make_shared<snapshots>();

    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;
    const int column_count = 15;
    io::CSVReader<column_count, TrimPolicy, QuotePolicy> in(file_name);

    string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15;
    int step_id = 0;
    shared_ptr<Particles> particles = make_shared<Particles>();

    in.read_header(io::ignore_extra_column, "step_id", "particle_id", "r", "g", "b", "x", "y", "z",
                   "vx", "vy", "vz", "m", "rad", "rest", "temp");

    while (in.read_row(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10,
                       col11, col12, col13, col14, col15)) {
        if (stoi(col1) != step_id) {
            particle_states->snaps.push_back(particles);
            particles = make_shared<Particles>();
            step_id = stoi(col1);
        }

        shared_ptr<Particle> particle = make_shared<Particle>();
        particle->particle_id = stoi(col2);
        particle->r = stod(col3);
        particle->g = stod(col4);
        particle->b = stod(col5);
        particle->x = stod(col6);
        particle->y = stod(col7);
        particle->z = stod(col8);
        particle->vx = stod(col9);
        particle->vy = stod(col10);
        particle->vz = stod(col11);
        particle->m = stod(col12);
        particle->rad = stod(col13);
        particle->rest = stod(col14);
        particle->temp = stod(col15);

        particles->particle_list.push_back(particle);
    }

    particle_states->snaps.push_back(particles);
    return particle_states;
}