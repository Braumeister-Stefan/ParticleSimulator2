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

#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"

#include "../include/PhysEngineWrapper.h"

using namespace std;
using namespace std::chrono;

using StepEnergyTrace = EngineCore::StepEnergyTrace;

static inline double safe_double_from_high_prec(const high_prec& v) {
    return v;
}
static inline double to_double(const high_prec& v) { return v; }

static inline high_prec abs_hp(high_prec v) { return v < 0.0 ? -v : v; }

static inline high_prec safe_rel_error(high_prec post, high_prec pre) {
    high_prec diff = post - pre;
    if (pre == 0.0) return diff;
    return diff / pre;
}

// =======================
// DEBUG & REPORTING CONFIG
// =======================
static const int  kStepPauseInterval   = 100;

// =======================
// CACHE CONFIG
// =======================
static const int  kCacheWriteEveryN    = 1;     // Save every Nth snapshot to cache (1 = all)

bool Engine::debug_mode(std::shared_ptr<scenario> scenario) { return scenario->debug_mode; }

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
static const high_prec kBenchmarkEqualRelThreshFrac = 0.0005; // 0.05% = 0.0005 fraction
// If benchmark is exactly 0, treat "equal" as absolute current error <= 0.05% (in percent units)
static const high_prec kBenchmarkEqualAbsThreshPct  = 0.05;

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
    bool have_rel = (benchmark_err_pct != 0.0);
    high_prec rel_diff_frac = 0;   // (current - benchmark) / benchmark
    high_prec rel_diff_pct  = 0;   // rel_diff_frac * 100

    if (have_rel) {
        rel_diff_frac = (current_err_pct - benchmark_err_pct) / benchmark_err_pct;
        rel_diff_pct  = rel_diff_frac * 100.0;
        if (abs_hp(rel_diff_frac) <= equal_rel_thresh_frac) verdict = "EQUAL ERROR";
        else if (current_err_pct > benchmark_err_pct)       verdict = "WORSE ERROR";
        else                                                verdict = "BETTER ERROR";
    } else {
        if (abs_hp(current_err_pct) <= equal_abs_thresh_pct) verdict = "EQUAL ERROR";
        else if (current_err_pct > 0.0)             verdict = "WORSE ERROR";
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
                  << "   (equal if |diff| <= " << (equal_rel_thresh_frac * 100.0) << " %)\n";
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
                  << ((current_time_sec - benchmark_time_sec) / benchmark_time_sec * 100.0) << " %\n";
    } else {
        std::cout << "No benchmark entry found for sim time.\n";
        std::cout << "Current time (s)       : " << current_time_sec << "\n";
    }
    line('#');
    cout << endl;
}

// =======================
// LIGHTWEIGHT SNAPSHOT (MINIMAL FIELDS COPY)
// =======================
static inline unique_ptr<Particles> make_light_snapshot(const Particles& src) {
    auto dst = make_unique<Particles>();

    const size_t n = src.particle_list.size();
    dst->particle_list.resize(n);

    for (size_t j = 0; j < n; ++j) {
        if (!dst->particle_list[j]) dst->particle_list[j] = make_shared<Particle>();

        const auto& sp = src.particle_list[j];
        if (!sp) continue;

        auto& dp = dst->particle_list[j];

        dp->particle_id = sp->particle_id;
        dp->r    = sp->r;
        dp->g    = sp->g;
        dp->b    = sp->b;
        dp->x    = sp->x;
        dp->y    = sp->y;
        dp->z    = sp->z;
        dp->vx   = sp->vx;
        dp->vy   = sp->vy;
        dp->vz   = sp->vz;
        dp->m    = sp->m;
        dp->rad  = sp->rad;
        dp->temp = sp->temp;
        dp->rest = sp->rest;
    }

    return dst;
}

// =======================
// RUN
// =======================
shared_ptr<snapshots> Engine::run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles)
{
    //set core physics parameters
    core.set_overlap_beta(scenario);
    core.set_collision_distance_tolerance(scenario);

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

    // Match plotter display stride: plotter steps by (int)(1/dt), so only snapshot that often
    const int write_k_frames = std::max(1, static_cast<int>(1.0 / scenario->dt));

    high_prec steps_db = scenario->time / scenario->dt;
    int steps = static_cast<int>(steps_db);
    cout << "Number of steps to be simulated: " << steps << endl << endl;

    // Pre-allocate metrics alongside snapshots (not per-step)
    // Reserve approximate snapshot capacity to reduce reallocations
    if (write_k_frames > 0) {
        const int expected_snaps = (steps + write_k_frames - 1) / write_k_frames;
        particle_states->snaps.reserve(static_cast<size_t>(expected_snaps));
        particle_states->metrics.reserve(static_cast<size_t>(expected_snaps));
    }

    // Initial TE baseline for final reporting
    high_prec TE_initial = core.calc_TE(particles);
    high_prec TE_final   = TE_initial;

    // Debug-mode timing accumulators
    double total_collision_secs = 0.0;
    double total_verlet_secs    = 0.0;
    double total_other_secs     = 0.0;

    EngineCore::clock_t::time_point sim_start = EngineCore::clock_t::now();

    // Progress reporting (5% buckets) with O(#reports) math
    int last_reported_bucket = -5;
    int next_report_bucket   = -1;
    int next_report_step     = std::numeric_limits<int>::max();

    if (steps > 0) {
        cout << "0% of the simulation complete." << endl;
        last_reported_bucket = 0;
        next_report_bucket   = 5;
        next_report_step     = (int)(((long long)steps * (long long)next_report_bucket + 99LL) / 100LL);
        if (next_report_step < 1) next_report_step = 1;
    }

    for (int i = 0; i < steps; i++) {
        EngineCore::update_iter = i;

        StepEnergyTrace tr;  // declared at loop scope so snapshot code can read it

        if (scenario->debug_mode) {

            EngineCore::clock_t::time_point step_start_report;
            EngineCore::clock_t::time_point step_end_report;

            auto iter_t0 = EngineCore::clock_t::now();
            if (scenario->report_energy_per_step) step_start_report = iter_t0;
            core.step(particles, &tr);
            auto iter_t1 = EngineCore::clock_t::now();
            if (scenario->report_energy_per_step) step_end_report = iter_t1;

            // Accumulate stage-wise TE deltas (signed)
            total_TE_error_overlap   += tr.dE_overlap;
            total_TE_error_collision += tr.dE_collision;
            total_TE_error_verlet    += tr.dE_verlet;

            // Accumulate stage timing
            total_collision_secs += tr.collision_seconds;
            total_verlet_secs    += tr.verlet_seconds;
            double iter_total = EngineCore::seconds_between(iter_t0, iter_t1);
            // Exclude plotting time from overhead
            // (Assume plotting is handled outside this timing, or add explicit timing if needed)
            total_other_secs += (iter_total - tr.collision_seconds - tr.verlet_seconds);

            if (scenario->report_energy_per_step) {
                double step_seconds = EngineCore::seconds_between(step_start_report, step_end_report);
                print_step_report(i, steps, tr, step_seconds);
                if (kStepPauseInterval > 0 && ((i + 1) % kStepPauseInterval == 0)) {
                    pause_for_user();
                }
            }
        } else {
            // No debug mode: run without trace (no energy tracking overhead)
            core.step(particles, nullptr);
        }

        if (write_k_frames > 0) {
            const bool should_write = ((i + 1) % write_k_frames == 0) || (i == steps - 1);
            if (should_write) {
                particle_states->snaps.push_back(make_light_snapshot(*particles));

                auto m = make_shared<test_metrics_t>();
                m->margin_TE_error = core.get_margin_TE_error();
                m->margin_TE_error_collision = core.get_margin_TE_error_collision();
                m->margin_TE_error_integrate = core.get_margin_TE_error_integrate();

                if (scenario->debug_mode) {
                    m->KE = tr.KE;
                    m->PE = tr.PE;
                    m->HE = tr.HE;
                    m->TE = tr.TE;
                    m->relative_error = abs_hp(safe_rel_error(tr.TE, TE_initial));

                    // Compute total momentum from current particle state
                    high_prec mx = 0.0, my = 0.0;
                    for (int pi = 0; pi < static_cast<int>(particles->particle_list.size()); ++pi) {
                        const auto& p = particles->particle_list[pi];
                        mx += p->m * p->vx;
                        my += p->m * p->vy;
                    }
                    m->mom_x = mx;
                    m->mom_y = my;
                }

                particle_states->metrics.push_back(m);
            }
        }

        // progress reporting: 5% buckets (compute pct only when needed)
        if (steps > 0 && (i + 1) >= next_report_step) {
            int pct_complete = (int)(((long long)(i + 1) * 100LL) / (long long)steps);
            if (pct_complete > 100) pct_complete = 100;
            if (pct_complete < 0) pct_complete = 0;

            int pct_bucket = (pct_complete / 5) * 5;
            if (pct_bucket >= 0 && pct_bucket <= 100 && pct_bucket > last_reported_bucket) {
                cout << pct_bucket << "% of the simulation complete." << endl;
                last_reported_bucket = pct_bucket;
            }

            next_report_bucket = last_reported_bucket + 5;
            if (next_report_bucket > 100) {
                next_report_step = std::numeric_limits<int>::max();
            } else {
                next_report_step = (int)(((long long)steps * (long long)next_report_bucket + 99LL) / 100LL);
                if (next_report_step < (i + 2)) next_report_step = (i + 2); // keep it forward-moving
            }
        }
    }

    EngineCore::clock_t::time_point sim_end = EngineCore::clock_t::now();
    double total_seconds = EngineCore::seconds_between(sim_start, sim_end);
    high_prec current_time_sec = total_seconds;

    // Calculate final TE once after all steps complete (2nd and final TE calculation)
    TE_final = core.calc_TE(particles);

    // =======================
    // END-OF-SIM LOGGING
    // =======================
    high_prec total_stage_sum = total_TE_error_overlap + total_TE_error_collision + total_TE_error_verlet;
    high_prec abs_stage_sum = abs_hp(total_TE_error_overlap) + abs_hp(total_TE_error_collision) + abs_hp(total_TE_error_verlet);

    auto pct_share = [&](high_prec part_abs) -> high_prec {
        if (abs_stage_sum == 0.0) return 0.0;
        return (part_abs / abs_stage_sum) * 100.0;
    };

    high_prec dE_total = TE_final - TE_initial;
    high_prec dE_total_rel = abs_hp(safe_rel_error(TE_final, TE_initial));
    high_prec dE_total_rel_pct = dE_total_rel * 100.0; // total relative TE error (%)

    cout << endl;
    cout << "========================================" << endl;
    cout << scenario->name << " simulation completed." << endl;
    cout << "========================================" << endl;
    cout << "Steps:            " << steps << endl;
    cout << "Total collisions: " << EngineCore::collisions << endl;
    cout << "Sim wall time:    " << std::setprecision(4) << std::defaultfloat << total_seconds << " s" << endl;

    cout << std::fixed << std::setprecision(15);

    // --- Energy summary (always printed) ---
    cout << endl;
    cout << "--- Energy Summary ---" << endl;
    cout << "  Initial TE: " << TE_initial << endl;
    cout << "  Final TE:   " << TE_final << endl;
    cout << "  Rel error:  " << dE_total_rel << " (" << dE_total_rel_pct << "%)" << endl;

    if (dE_total > 0.0)      cout << "  Verdict:    ENERGY INCREASED" << endl;
    else if (dE_total < 0.0) cout << "  Verdict:    ENERGY DECREASED" << endl;
    else                      cout << "  Verdict:    ENERGY UNCHANGED" << endl;

    // --- Stage error breakdown (debug mode only) ---
    if (scenario->debug_mode) {
           cout << endl;
           cout << "--- Stage TE Error Breakdown (debug) ---" << endl;
           cout << "  Sum of parts (signed): " << total_stage_sum << endl;
           cout << "  Overlaps:   " << total_TE_error_overlap
               << " (abs share " << pct_share(abs_hp(total_TE_error_overlap)) << "%" << endl;
           cout << "  Collisions: " << total_TE_error_collision
               << " (abs share " << pct_share(abs_hp(total_TE_error_collision)) << "%" << endl;
           cout << "  Verlet:     " << total_TE_error_verlet
               << " (abs share " << pct_share(abs_hp(total_TE_error_verlet)) << "%" << endl;
    }

    // --- Stage timing breakdown (debug mode only) ---
    if (scenario->debug_mode) {
        double timing_sum = total_collision_secs + total_verlet_secs + total_other_secs;
        auto timing_pct = [&](double part) -> double {
            return timing_sum > 0.0 ? (part / timing_sum) * 100.0 : 0.0;
        };
        cout << endl;
        cout << "--- Stage Timing Breakdown (debug) ---" << endl;
        cout << std::fixed << std::setprecision(4);
        cout << "  Collision handling: " << total_collision_secs << " s (" << timing_pct(total_collision_secs) << "%)" << endl;
        cout << "  Velocity Verlet:    " << total_verlet_secs    << " s (" << timing_pct(total_verlet_secs)    << "%)" << endl;
        cout << "  Other (overhead):   " << total_other_secs     << " s (" << timing_pct(total_other_secs)     << "%)" << endl;
        cout << "  Total timed:        " << timing_sum           << " s" << endl;
    }

    // =======================
    // BENCHMARK COMPARISON 
    // =======================
    {
        high_prec bench_pct  = scenario->benchmark_te_error_pct;
        high_prec bench_time = scenario->benchmark_sim_time_sec;
        bool have_te_bench   = (bench_pct  >= 0.0);
        bool have_time_bench = (bench_time >= 0.0);

        if (have_te_bench) {
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

    if (scenario->save_scenario) {
        run_to_cache(scenario, particle_states);
    }
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

    const int total_snaps = (int)particle_states->snaps.size();
    const int n_particles = total_snaps > 0 ? (int)particle_states->snaps[0]->particle_list.size() : 0;
    const int stride = std::max(1, kCacheWriteEveryN);
    const int written_snaps = (total_snaps + stride - 1) / stride;
    const double compression = total_snaps > 0 ? (1.0 - (double)written_snaps / total_snaps) * 100.0 : 0.0;
    cout << "Writing " << written_snaps << "/" << total_snaps << " states for " << n_particles
         << " particles" << std::fixed << std::setprecision(1) << "("
         << compression << "% reduction)" << endl;

    file << std::fixed << std::setprecision(15);
    file << "step_id, particle_id,r,g,b,x,y,z,vx,vy,vz,m,rad,temp,rest\n";

    int last_reported_bucket = -5;
    int next_report_step = 0;
    int written = 0;

    for (int i = 0; i < total_snaps; i += stride) {
        // 5% progress reporting (based on written count)
        if (written_snaps > 0 && written >= next_report_step) {
            int pct = (int)(((long long)(written + 1) * 100LL) / (long long)written_snaps);
            int bucket = (pct / 20) * 20;
            if (bucket > last_reported_bucket) {
                cout << "  " << bucket << "%" << endl;
                last_reported_bucket = bucket;
            }
            int nb = last_reported_bucket + 5;
            next_report_step = (nb > 100) ? written_snaps + 1
                : (int)(((long long)written_snaps * (long long)nb + 99LL) / 100LL);
        }
        ++written;

        for (int j = 0; j < (int)particle_states->snaps[i]->particle_list.size(); j++) {
            const auto& p = particle_states->snaps[i]->particle_list[j];
            
            // Check for negative temperature before writing
            // if (p->temp < 0) {
            //     cout << "Negative temperature detected for particle " << p->particle_id
            //          << " at step " << i << ". Pausing for user input." << endl;
            //     pause_for_user();
            // }
            
            file << i << ","
                << p->particle_id << ","
                << p->r << ","
                << p->g << ","
                << p->b << ","
                << p->x << ","
                << p->y << ","
                << p->z << ","
                << p->vx << ","
                << p->vy << ","
                << p->vz << ","
                << p->m << ","
                << p->rad << ","
                << p->temp << ","
                << p->rest << '\n';
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

    unique_ptr<Particles> particles = make_unique<Particles>();

    in.read_header(io::ignore_extra_column, "step_id", "particle_id", "r", "g", "b", "x", "y", "z",
                   "vx", "vy", "vz", "m", "rad", "rest", "temp");

    while (in.read_row(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10,
                       col11, col12, col13, col14, col15)) {
        if (stoi(col1) != step_id) {
            particle_states->snaps.push_back(std::move(particles));
            particles = make_unique<Particles>();
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

    particle_states->snaps.push_back(std::move(particles));
    return particle_states;
}
