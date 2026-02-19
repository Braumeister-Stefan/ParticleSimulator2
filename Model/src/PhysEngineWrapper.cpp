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
#include <filesystem>

#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"

#include "../include/PhysEngineWrapper.h"

using namespace std;
using namespace std::chrono;

using StepEnergyTrace = EngineCore::StepEnergyTrace;

static inline double safe_double_from_high_prec(const high_prec& v) { return v; }
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
static const int  kStepPauseInterval = 100;

// =======================
// CACHE CONFIG
// =======================
static const int  kCacheWriteEveryN     = 20;     // Save every Nth snapshot within a flush chunk
static const int  kCacheFlushEverySteps = 5000;   // flush every N simulation timesteps

// =======================
// INCREMENTAL CACHE STATE (per run, file-scope)
// =======================
static bool      gCacheInitializedThisRun = false;
static long long gCacheStepIdOffset       = 0;
static long long gCacheMaxStepId          = -1;

static inline std::string cache_file_from_name_or_path(const std::string& name_or_path) {
    std::string file_name = name_or_path;

    const bool has_sep =
        (file_name.find('/')  != std::string::npos) ||
        (file_name.find('\\') != std::string::npos);

    const bool ends_with_csv =
        (file_name.size() >= 4) &&
        (file_name.substr(file_name.size() - 4) == ".csv");

    if (!has_sep) {
        if (!ends_with_csv) file_name = "Inputs/rendered_scenarios/" + file_name + ".csv";
        else                file_name = "Inputs/rendered_scenarios/" + file_name;
    }
    return file_name;
}

// =======================
// CACHE SCHEMA HELPER ("dictionary")
// =======================
struct CacheSchema {
    static constexpr int kColumnCount = 15;

    static inline void write_header(std::ostream& os) {
        os << "step_id, particle_id,r,g,b,x,y,z,vx,vy,vz,m,rad,temp,rest\n";
    }

    template <typename TrimPolicy, typename QuotePolicy>
    static inline void read_header(io::CSVReader<kColumnCount, TrimPolicy, QuotePolicy>& in) {
        in.read_header(io::ignore_extra_column,
            "step_id",
            "particle_id",
            "r", "g", "b",
            "x", "y", "z",
            "vx", "vy", "vz",
            "m",
            "rad",
            "temp",
            "rest"
        );
    }

    struct Row {
        std::string step_id;
        std::string particle_id;
        std::string r, g, b;
        std::string x, y, z;
        std::string vx, vy, vz;
        std::string m;
        std::string rad;
        std::string temp;
        std::string rest;
    };

    static inline void write_row(std::ostream& os, long long step_id, const Particle& p) {
        os << step_id << ","
           << p.particle_id << ","
           << p.r << "," << p.g << "," << p.b << ","
           << p.x << "," << p.y << "," << p.z << ","
           << p.vx << "," << p.vy << "," << p.vz << ","
           << p.m << ","
           << p.rad << ","
           << p.temp << ","
           << p.rest << '\n';
    }

    static inline long long parse_step_id(const Row& row) {
        return std::stoll(row.step_id);
    }

    static inline void parse_particle_fields(const Row& row, Particle& p) {
        p.particle_id = std::stoi(row.particle_id);
        p.r = std::stod(row.r);
        p.g = std::stod(row.g);
        p.b = std::stod(row.b);
        p.x = std::stod(row.x);
        p.y = std::stod(row.y);
        p.z = std::stod(row.z);
        p.vx = std::stod(row.vx);
        p.vy = std::stod(row.vy);
        p.vz = std::stod(row.vz);
        p.m = std::stod(row.m);
        p.rad = std::stod(row.rad);
        p.temp = std::stod(row.temp);
        p.rest = std::stod(row.rest);
    }
};

static bool load_particles_from_cache_step(const std::string& name_or_path,
                                           long long target_step_id,
                                           Particles& out)
{
    const std::string file_name = cache_file_from_name_or_path(name_or_path);

    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;

    io::CSVReader<CacheSchema::kColumnCount, TrimPolicy, QuotePolicy> in(file_name);
    CacheSchema::read_header(in);

    CacheSchema::Row row;

    out.particle_list.clear();
    bool collecting = false;

    while (in.read_row(
        row.step_id,
        row.particle_id,
        row.r, row.g, row.b,
        row.x, row.y, row.z,
        row.vx, row.vy, row.vz,
        row.m,
        row.rad,
        row.temp,
        row.rest))
    {
        const long long sid = CacheSchema::parse_step_id(row);

        if (sid == target_step_id) {
            collecting = true;
            auto particle = std::make_shared<Particle>();
            CacheSchema::parse_particle_fields(row, *particle);
            out.particle_list.push_back(particle);
        } else if (collecting) {
            break;
        }
    }

    return !out.particle_list.empty();
}

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
        pct = std::clamp(pct, 0, 100);
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
    std::cout << "Work counters (trace):\n";
    std::cout << "  Collision checks = " << tr.collision_checks
              << "   Collisions = " << tr.collision_hits
              << " (" << std::setprecision(4) << tr.collision_hit_pct << "%)\n";
    std::cout << std::setprecision(0);
    std::cout << "  Gravity pair-forces (i<j) evaluated = " << tr.gravity_pair_calcs << "\n";

    line('-');
    std::cout << std::setprecision(12);
    std::cout << "TE checkpoints:\n";
    std::cout << "  TE0(before overlaps) = " << to_double(tr.TE0) << "\n";
    std::cout << "  TE1(after overlaps ) = " << to_double(tr.TE1) << "\n";
    std::cout << "  TE2(after collisions) = " << to_double(tr.TE2) << "\n";
    std::cout << "  TE3(after verlet   ) = " << to_double(tr.TE3) << "\n";
    line('=');
}

static const high_prec kBenchmarkEqualRelThreshFrac = 0.0005;
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
    high_prec rel_diff_frac = 0;
    high_prec rel_diff_pct  = 0;

    if (have_rel) {
        rel_diff_frac = (current_err_pct - benchmark_err_pct) / benchmark_err_pct;
        rel_diff_pct  = rel_diff_frac * 100.0;
        if (abs_hp(rel_diff_frac) <= equal_rel_thresh_frac) verdict = "EQUAL ERROR";
        else if (current_err_pct > benchmark_err_pct)       verdict = "WORSE ERROR";
        else                                                verdict = "BETTER ERROR";
    } else {
        if (abs_hp(current_err_pct) <= equal_abs_thresh_pct) verdict = "EQUAL ERROR";
        else if (current_err_pct > 0.0)                      verdict = "WORSE ERROR";
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

static inline void clear_particle_states_storage(shared_ptr<snapshots> particle_states) {
    if (!particle_states) return;
    particle_states->snaps.clear();
    particle_states->metrics.clear();
}

double Engine::check_mb(shared_ptr<scenario> scenario) {
    try {
        namespace fs = std::filesystem;
        const std::string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";
        if (!fs::exists(file_name)) return 0.0;
        const auto bytes = static_cast<double>(fs::file_size(file_name));
        return bytes / (1024.0 * 1024.0);
    } catch (...) {
        return 0.0;
    }
}

bool Engine::initiate_cache(shared_ptr<scenario> scenario) {
    const bool existed = cache_exists(scenario);

    cout << "cache_exists: " << existed << endl;

    const string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";
    ofstream file(file_name, std::ios::out | std::ios::trunc);
    if (!file.is_open()) {
        cout << "Failed to create a clean cache file!" << endl;
        return false;
    }

    cout << "Initializing cache file: " << file_name << endl;

    file << std::fixed << std::setprecision(15);
    CacheSchema::write_header(file);
    file.close();

    gCacheInitializedThisRun = true;
    gCacheStepIdOffset       = 0;

    if (existed) cout << "Cache existed and was overwritten: " << file_name << endl;
    else         cout << "Cache initialized: " << file_name << endl;

    return true;
}

unique_ptr<Particles> Engine::make_light_snapshot(const Particles& src) {
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

void Engine::run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles)
{
    core.set_overlap_beta(scenario);
    core.set_collision_distance_tolerance(scenario);

    size_t initial_n_particles = particles ? particles->particle_list.size() : 0;

    high_prec total_TE_error_overlap   = 0;
    high_prec total_TE_error_collision = 0;
    high_prec total_TE_error_verlet    = 0;
    EngineCore::collisions = 0;

    long long total_collision_checks = 0;
    long long total_collision_hits   = 0;
    long long total_gravity_pairs    = 0;

    cout << "Engine (Wrapper) is initialized." << endl;

    shared_ptr<snapshots> particle_states = make_shared<snapshots>();
    cout << scenario->name << "'s initial states loaded" << endl << endl;

    cout << "Press enter to start the simulation." << endl;
    cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    core.configure_dt(scenario, particles);

    const int write_k_frames = std::max(1, static_cast<int>(1.0 / scenario->dt));

    high_prec steps_db = scenario->time / scenario->dt;
    int steps = static_cast<int>(steps_db);
    cout << "Number of steps to be simulated: " << steps << endl << endl;


    int start_iter = 0;

    bool cache_exists_flag = cache_exists(scenario);

    bool cache_enabled = true; // WIP. User should have choice to ignore cached files and start fresh, or look for it.

    if (cache_enabled && !cache_exists_flag) {
        cache_enabled = initiate_cache(scenario);
        if (!cache_enabled) {
            cout << "Cache initialization failed. Proceeding without cache." << endl;
        }
    } else if (cache_enabled && cache_exists_flag) {
        cout << "Cache exists. Inspecting contents..." << endl;
        inspect_cache(scenario->name);

        if (gCacheMaxStepId >= 0) {
            if (!particles) particles = std::make_shared<Particles>();

            if (load_particles_from_cache_step(scenario->name, gCacheMaxStepId, *particles)) {
                start_iter = static_cast<int>(gCacheMaxStepId + 1);

                gCacheInitializedThisRun = true;
                gCacheStepIdOffset       = static_cast<long long>(start_iter);

                cout << "Restored up to " << gCacheMaxStepId << endl;
                     
            } else {
                cout << "Failed to load cached step_id " << gCacheMaxStepId
                     << ". Starting from step 0." << endl;
                gCacheMaxStepId = -1;
            }
        }
    }

    if (write_k_frames > 0) {
        int window_steps = steps;
        if (cache_enabled && kCacheFlushEverySteps > 0) {
            window_steps = std::min(steps, kCacheFlushEverySteps);
        }
        const int expected_snaps = (window_steps + write_k_frames - 1) / write_k_frames;
        particle_states->snaps.reserve(static_cast<size_t>(expected_snaps));
        particle_states->metrics.reserve(static_cast<size_t>(expected_snaps));
    }

    high_prec TE_initial = core.calc_TE(particles);
    high_prec TE_final   = TE_initial;

    double total_collision_secs     = 0.0;
    double total_verlet_secs        = 0.0;
    double total_gravity_force_secs = 0.0;
    double total_integration_secs   = 0.0;
    double total_other_secs         = 0.0;

    EngineCore::clock_t::time_point sim_start = EngineCore::clock_t::now();

    int last_reported_bucket = -5;
    int next_report_bucket   = -1;
    int next_report_step     = std::numeric_limits<int>::max();

    if (steps > 0) {
        //cout << "0% of the simulation complete." << endl;
        last_reported_bucket = 0;
        next_report_bucket   = 5;
        next_report_step     = (int)(((long long)steps * (long long)next_report_bucket + 99LL) / 100LL);
        if (next_report_step < 1) next_report_step = 1;

        if (cache_enabled) {
            auto old_flags = cout.flags();
            auto old_prec  = cout.precision();
            cout << "  Cached: " << std::fixed << std::setprecision(2) << check_mb(scenario) << " MB" << endl;
            cout.flags(old_flags);
            cout.precision(old_prec);
        }
    }

    for (int i = start_iter; i < steps; i++) {
        EngineCore::update_iter = i;

        StepEnergyTrace tr;

        if (scenario->debug_mode) {
            EngineCore::clock_t::time_point step_start_report;
            EngineCore::clock_t::time_point step_end_report;

            auto iter_t0 = EngineCore::clock_t::now();
            if (scenario->report_energy_per_step) step_start_report = iter_t0;
            core.step(particles, &tr);
            auto iter_t1 = EngineCore::clock_t::now();
            if (scenario->report_energy_per_step) step_end_report = iter_t1;

            total_TE_error_overlap   += tr.dE_overlap;
            total_TE_error_collision += tr.dE_collision;
            total_TE_error_verlet    += tr.dE_verlet;

            total_collision_secs     += tr.collision_seconds;
            total_verlet_secs        += tr.verlet_seconds;
            total_gravity_force_secs += tr.gravity_force_seconds;
            total_integration_secs   += tr.integration_seconds;

            double iter_total = EngineCore::seconds_between(iter_t0, iter_t1);
            total_other_secs += (iter_total - tr.collision_seconds - tr.verlet_seconds);

            total_collision_checks += tr.collision_checks;
            total_collision_hits   += tr.collision_hits;
            total_gravity_pairs    += tr.gravity_pair_calcs;

            if (scenario->report_energy_per_step) {
                double step_seconds = EngineCore::seconds_between(step_start_report, step_end_report);
                print_step_report(i, steps, tr, step_seconds);
                if (kStepPauseInterval > 0 && ((i + 1) % kStepPauseInterval == 0)) {
                    pause_for_user();
                }
            }
        } else {
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

        if (steps > 0 && (i + 1) >= next_report_step) {
            int pct_complete = (int)(((long long)(i + 1) * 100LL) / (long long)steps);
            pct_complete = std::clamp(pct_complete, 0, 100);

            int pct_bucket = (pct_complete / 5) * 5;
            if (pct_bucket >= 0 && pct_bucket <= 100 && pct_bucket > last_reported_bucket) {
                cout << pct_bucket << "% of the simulation complete." << endl;
                last_reported_bucket = pct_bucket;

                if (cache_enabled) {
                    auto old_flags = cout.flags();
                    auto old_prec  = cout.precision();
                    cout << "  Cache written: " << std::fixed << std::setprecision(2) << check_mb(scenario) << " MB" << endl;
                    cout.flags(old_flags);
                    cout.precision(old_prec);
                }
            }

            next_report_bucket = last_reported_bucket + 5;
            if (next_report_bucket > 100) {
                next_report_step = std::numeric_limits<int>::max();
            } else {
                next_report_step = (int)(((long long)steps * (long long)next_report_bucket + 99LL) / 100LL);
                if (next_report_step < (i + 2)) next_report_step = (i + 2);
            }
        }

        if (cache_enabled && kCacheFlushEverySteps > 0 && ((i + 1) % kCacheFlushEverySteps == 0)) {
            if (!particle_states->snaps.empty()) {
                run_to_cache(scenario, particle_states);
                clear_particle_states_storage(particle_states);
            }
        }
    }

    EngineCore::clock_t::time_point sim_end = EngineCore::clock_t::now();
    double total_seconds = EngineCore::seconds_between(sim_start, sim_end);
    high_prec current_time_sec = total_seconds;

    TE_final = core.calc_TE(particles);

    high_prec total_stage_sum = total_TE_error_overlap + total_TE_error_collision + total_TE_error_verlet;
    high_prec abs_stage_sum = abs_hp(total_TE_error_overlap) + abs_hp(total_TE_error_collision) + abs_hp(total_TE_error_verlet);

    auto pct_share = [&](high_prec part_abs) -> high_prec {
        if (abs_stage_sum == 0.0) return 0.0;
        return (part_abs / abs_stage_sum) * 100.0;
    };

    high_prec dE_total = TE_final - TE_initial;
    high_prec dE_total_rel = abs_hp(safe_rel_error(TE_final, TE_initial));
    high_prec dE_total_rel_pct = dE_total_rel * 100.0;

    cout << endl;
    cout << "========================================" << endl;
    cout << scenario->name << " simulation completed." << endl;
    cout << "========================================" << endl;
    cout << "Steps:            " << steps << endl;
    cout << "Total collisions: " << EngineCore::collisions << endl;
    cout << "Sim wall time:    " << std::setprecision(4) << std::defaultfloat << total_seconds << " s" << endl;

    cout << std::fixed << std::setprecision(15);

    cout << endl;
    cout << "--- Energy Summary ---" << endl;
    cout << "  Initial TE: " << TE_initial << endl;
    cout << "  Final TE:   " << TE_final << endl;
    cout << "  Rel error:  " << dE_total_rel << " (" << dE_total_rel_pct << "%)" << endl;

    if (dE_total > 0.0)      cout << "  Verdict:    ENERGY INCREASED" << endl;
    else if (dE_total < 0.0) cout << "  Verdict:    ENERGY DECREASED" << endl;
    else                      cout << "  Verdict:    ENERGY UNCHANGED" << endl;

    if (scenario->debug_mode) {
        cout << endl;
        cout << "--- Stage TE Error Breakdown (debug) ---" << endl;
        cout << "  Sum of parts: " << total_stage_sum << endl;
        cout << "  Collisions: " << total_TE_error_collision
             << " (abs share " << pct_share(abs_hp(total_TE_error_collision)) << "%)" << endl;
        cout << "  Verlet:     " << total_TE_error_verlet
             << " (abs share " << pct_share(abs_hp(total_TE_error_verlet)) << "%)" << endl;
    }

    if (scenario->debug_mode) {
        const double hit_pct = (total_collision_checks > 0)
            ? (100.0 * (double)total_collision_hits / (double)total_collision_checks)
            : 0.0;

        const double avg_checks = (steps > 0) ? ((double)total_collision_checks / (double)steps) : 0.0;
        const double avg_hits   = (steps > 0) ? ((double)total_collision_hits   / (double)steps) : 0.0;
        const double avg_gpairs = (steps > 0) ? ((double)total_gravity_pairs    / (double)steps) : 0.0;

        long long n2_pairs_per_step = (long long)initial_n_particles * (initial_n_particles - 1) / 2;
        long long n2_pairs = n2_pairs_per_step * steps;
        double bh_reduction = (n2_pairs > 0)
            ? 100.0 * (1.0 - (double)total_gravity_pairs / (double)n2_pairs)
            : 0.0;

        cout << endl;
        cout << "--- Work Counters (debug) ---" << endl;
        cout << std::fixed << std::setprecision(2);

        long long n2_collisions_per_step = n2_pairs_per_step;
        long long n2_collisions = n2_collisions_per_step * steps;
        double bh_collision_reduction = (n2_collisions > 0)
            ? 100.0 * (1.0 - (double)total_collision_checks / (double)n2_collisions)
            : 0.0;

        cout << "  Collision checks (total): " << total_collision_checks
             << " [O(n^2) = " << n2_collisions << ", BH reduction " << std::setprecision(1) << bh_collision_reduction << "%]" << endl;
        cout << std::setprecision(2);
        cout << "  Collisions (total):       " << total_collision_hits << " (" << hit_pct << "% of checks)" << endl;
        cout << std::setprecision(1);
        cout << "  Collision checks/step:    " << avg_checks << endl;
        cout << "  Collisions/step:          " << avg_hits << endl;
        cout << "  Gravity pair-forces/step: " << avg_gpairs << endl;
        cout << std::setprecision(0);
        cout << "  Gravity pair-forces total: " << total_gravity_pairs
             << " [O(n^2) = " << n2_pairs << ", BH reduction " << std::setprecision(1) << bh_reduction << "%]" << endl;
        cout << std::setprecision(0);
    }

    if (scenario->debug_mode) {
        double timing_sum = total_collision_secs
                          + total_gravity_force_secs
                          + total_integration_secs
                          + total_other_secs;

        auto timing_pct = [&](double part) -> double {
            return timing_sum > 0.0 ? (part / timing_sum) * 100.0 : 0.0;
        };

        cout << endl;
        cout << "--- Stage Timing Breakdown (debug) ---" << endl;
        cout << std::fixed << std::setprecision(4);

        cout << "  Collision handling: " << total_collision_secs << " s (" << timing_pct(total_collision_secs) << "%)" << endl;
        cout << "  Gravity forces: " << total_gravity_force_secs << " s (" << timing_pct(total_gravity_force_secs) << "%)" << endl;
        cout << "  Integration:    " << total_integration_secs   << " s (" << timing_pct(total_integration_secs)   << "%)" << endl;
        cout << "  Other (overhead):   " << total_other_secs     << " s (" << timing_pct(total_other_secs)     << "%)" << endl;
        cout << "  Total timed:        " << timing_sum           << " s (" << 100.0 << "%)" << endl;
    }

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

    if (cache_enabled) {
        if (!particle_states->snaps.empty()) {
            run_to_cache(scenario, particle_states);
            clear_particle_states_storage(particle_states);
        }
        gCacheInitializedThisRun = false;
        gCacheStepIdOffset       = 0;
        gCacheMaxStepId          = -1;
    }

    
}

void Engine::run_to_cache(std::shared_ptr<scenario> scenario, std::shared_ptr<snapshots> particle_states) {
    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    std::ios::openmode mode = std::ios::out;
    if (gCacheInitializedThisRun) mode |= std::ios::app;
    else                          mode |= std::ios::trunc;

    ofstream file(file_name, mode);
    if (!file.is_open()) {
        cout << "Failed to write snapshots to cache!" << endl;
        return;
    }

    const int total_snaps = (int)particle_states->snaps.size();
    if (total_snaps <= 0) return;

    const int n_particles = (int)particle_states->snaps[0]->particle_list.size();
    const int stride = std::max(1, kCacheWriteEveryN);
    const int written_snaps = (total_snaps + stride - 1) / stride;
    const double compression = total_snaps > 0 ? (1.0 - (double)written_snaps / total_snaps) * 100.0 : 0.0;

    cout << "Writing " << written_snaps << "/" << total_snaps << " states for " << n_particles
         << " particles" << std::fixed << std::setprecision(1)
         << "(" << compression << "% reduction)" << endl;

    file << std::fixed << std::setprecision(15);

    if (!gCacheInitializedThisRun) {
        CacheSchema::write_header(file);
    }

    const long long base_step_id = gCacheInitializedThisRun ? gCacheStepIdOffset : 0LL;

    int last_reported_bucket = -5;
    int next_report_step = 0;
    int written = 0;

    for (int i = 0; i < total_snaps; i += stride) {
        if (written_snaps > 0 && written >= next_report_step) {
            int pct = (int)(((long long)(written + 1) * 100LL) / (long long)written_snaps);
            int bucket = (pct / 20) * 20;
            if (bucket > last_reported_bucket) {
                //cout << "  " << bucket << "%" << endl;
                last_reported_bucket = bucket;
            }
            int nb = last_reported_bucket + 5;
            next_report_step = (nb > 100) ? written_snaps + 1
                : (int)(((long long)written_snaps * (long long)nb + 99LL) / 100LL);
        }
        ++written;

        const long long step_id = base_step_id + (long long)i;

        for (int j = 0; j < (int)particle_states->snaps[i]->particle_list.size(); j++) {
            const auto& sp = particle_states->snaps[i]->particle_list[j];
            if (!sp) continue;
            CacheSchema::write_row(file, step_id, *sp);
        }
    }

    file.close();
    cout << "Snapshots saved to " << file_name << endl;

    if (gCacheInitializedThisRun) {
        gCacheStepIdOffset += (long long)total_snaps;
    }
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
    typedef io::double_quote_escape<',', '"'> QuotePolicy;

    io::CSVReader<CacheSchema::kColumnCount, TrimPolicy, QuotePolicy> in(file_name);
    CacheSchema::read_header(in);

    CacheSchema::Row row;

    bool have_any = false;
    long long current_step_id = 0;
    unique_ptr<Particles> particles = make_unique<Particles>();

    long long step_count = 0;
    long long progressStep = 0;
    long long total_steps = 0;

    // First pass: count total steps
    {
        io::CSVReader<CacheSchema::kColumnCount, TrimPolicy, QuotePolicy> count_in(file_name);
        CacheSchema::read_header(count_in);
        CacheSchema::Row count_row;
        long long last_step_id = -1;
        while (count_in.read_row(
            count_row.step_id,
            count_row.particle_id,
            count_row.r, count_row.g, count_row.b,
            count_row.x, count_row.y, count_row.z,
            count_row.vx, count_row.vy, count_row.vz,
            count_row.m,
            count_row.rad,
            count_row.temp,
            count_row.rest)) {
            long long sid = CacheSchema::parse_step_id(count_row);
            if (sid != last_step_id) {
                ++total_steps;
                last_step_id = sid;
            }
        }
        progressStep = total_steps / 4; // 5% increments
        if (progressStep < 1) progressStep = 1;
    }

    // Second pass: actual reading
    long long current_progress = 0;
    while (in.read_row(
        row.step_id,
        row.particle_id,
        row.r, row.g, row.b,
        row.x, row.y, row.z,
        row.vx, row.vy, row.vz,
        row.m,
        row.rad,
        row.temp,
        row.rest))
    {
        const long long sid = CacheSchema::parse_step_id(row);

        if (!have_any) {
            have_any = true;
            current_step_id = sid;
            ++step_count;
            current_progress = step_count;
            if (progressStep > 0 && (current_progress % progressStep == 0)) {
                std::cout << "Cache loading " << (current_progress * 100) / total_steps << "% complete." << std::endl;
            }
        } else if (sid != current_step_id) {
            particle_states->snaps.push_back(std::move(particles));
            particles = make_unique<Particles>();
            current_step_id = sid;
            ++step_count;
            current_progress = step_count;
            if (progressStep > 0 && (current_progress % progressStep == 0)) {
                std::cout << "Cache loading " << (current_progress * 100) / total_steps << "% complete." << std::endl;
            }
        }

        shared_ptr<Particle> particle = make_shared<Particle>();
        CacheSchema::parse_particle_fields(row, *particle);
        particles->particle_list.push_back(particle);

    }

    if (have_any) {
        particle_states->snaps.push_back(std::move(particles));
    }

    return particle_states;
}

void Engine::inspect_cache(const std::string& name_or_path) {
    const std::string file_name = cache_file_from_name_or_path(name_or_path);

    {
        std::ifstream f(file_name);
        if (!f.is_open()) {
            std::cout << "Cache file not found/unreadable: " << file_name << std::endl;
            gCacheMaxStepId = -1;
            return;
        }
    }

    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;

    io::CSVReader<CacheSchema::kColumnCount, TrimPolicy, QuotePolicy> in(file_name);
    CacheSchema::read_header(in);

    CacheSchema::Row row;
    long long max_value = std::numeric_limits<long long>::min();

    while (in.read_row(
        row.step_id,
        row.particle_id,
        row.r, row.g, row.b,
        row.x, row.y, row.z,
        row.vx, row.vy, row.vz,
        row.m,
        row.rad,
        row.temp,
        row.rest))
    {
        long long val = CacheSchema::parse_step_id(row);
        if (val > max_value) max_value = val;
    }

    if (max_value == std::numeric_limits<long long>::min()) {
        std::cout << "Cache file is empty or unreadable: " << file_name << std::endl;
        gCacheMaxStepId = -1;
    } else {
        
        gCacheMaxStepId = max_value;
    }
}
