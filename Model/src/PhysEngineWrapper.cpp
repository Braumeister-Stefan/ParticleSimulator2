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
static const int  kCacheWriteEveryN     = 50;     // Save every Nth snapshot within a flush chunk
static const int  kCacheFlushEverySteps = 5000;   // flush every N simulation timesteps
// Heat-brightness (glow mode) is now controlled by scenario->glow_mode

// =======================
// INCREMENTAL CACHE STATE (per run, file-scope)
// =======================
static bool      gCacheInitializedThisRun = false;
static long long gCacheStepIdOffset       = 0;     // (kept as-is; no longer used to generate step_id values)
static long long gCacheMaxStepId          = -1;

// Per-snapshot true simulation step ids (kept in sync with particle_states->snaps)
static std::vector<long long> gCacheSnapStepIds;

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

    static inline void write_row(std::ostream& os, long long step_id, const Particle& p,
                                  bool zero_temp = false) {
        os << step_id << ","
           << p.particle_id << ","
           << p.r << "," << p.g << "," << p.b << ","
           << p.x << "," << p.y << "," << p.z << ","
           << p.vx << "," << p.vy << "," << p.vz << ","
           << p.m << ","
           << p.rad << ","
           << (zero_temp ? 0.0 : static_cast<double>(p.temp)) << ","
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
        p.rgb = (static_cast<int>(p.r * 255.0) << 16)
              | (static_cast<int>(p.g * 255.0) << 8)
              |  static_cast<int>(p.b * 255.0);
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
    gCacheSnapStepIds.clear();
}

double Engine::check_mb(shared_ptr<scenario> scenario) {
    try {
        namespace fs = std::filesystem;
        const std::string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

        
        if (!fs::exists(file_name)) return 0.0;
            double bytes = static_cast<double>(fs::file_size(file_name));
        return bytes / (1024.0 * 1024.0);
    } catch (...) {
        return 0.0;
    }
}

bool Engine::initiate_cache(shared_ptr<scenario> scenario) {
    if (gCacheInitializedThisRun) {
        cout << "Cache already initialized for this run. Skipping re-initialization." << endl;
        return true;
    }

    const bool existed = cache_exists(scenario);

    

    const string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";
    ofstream file(file_name, std::ios::out | std::ios::trunc);
    if (!file.is_open()) {
        cout << "Failed to create a clean cache file!" << endl;
        return false;
    }

    cout << "Initializing cache..." << endl;

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


    shared_ptr<snapshots> particle_states = make_shared<snapshots>();
    gCacheSnapStepIds.clear();

    cout << scenario->name << "'s initial states loaded" << endl << endl;

    cout << "Press enter to start the simulation." << endl;
    cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    core.configure_dt(scenario, particles);

    const int write_k_frames = std::max(1, static_cast<int>(1.0 / scenario->dt));

    high_prec steps_db = scenario->time / scenario->dt;
    int steps = static_cast<int>(steps_db);
    cout << "Steps: " << steps  << ", Particles: " << initial_n_particles << endl << endl;

    int start_iter = 0;

    bool cache_enabled = scenario->try_cache;

    //cout << "cache_enabled:" << cache_enabled << endl;
    bool cache_exists_flag = cache_exists(scenario);

    // ============================================================
    // NEW: If try_cache == false, delete any existing cache and
    //      reinitialize from scratch.
    // ============================================================
    if (!cache_enabled) {
        if (cache_exists_flag) {
            const std::string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";
            try {
                std::filesystem::remove(file_name);
                cout << "try_cache=false: Deleted existing cache file: " << file_name << endl;
                
            } catch (...) {
                // best-effort; initiate_cache will still attempt trunc/create
                cout << "try_cache=false: Failed to delete cache file (will reinit anyway): " << file_name << endl;
            }
            

            cache_enabled = initiate_cache(scenario); 
            
            cache_exists_flag = true; // since we just created it, probally redundant

            

            
        } else {
            cache_enabled = initiate_cache(scenario);

            cache_exists_flag = true;
        }

        // Ensure we start from scratch regardless
        start_iter = 0;
        gCacheMaxStepId = -1;
        gCacheInitializedThisRun = false;
        gCacheStepIdOffset = 0;
        gCacheSnapStepIds.clear();
    }

    if (scenario->try_cache && !cache_exists_flag) {
        cache_enabled = initiate_cache(scenario);
        if (!cache_enabled) {
            cout << "Cache initialization failed. Proceeding without cache." << endl;
        }
    } else if (cache_enabled && cache_exists_flag) {
        cout << "Cache exists. Inspecting contents..." << endl;
        if (!inspect_cache(scenario->name)) {
            cout << "Cache invalid. Reinitializing..." << endl;
            initiate_cache(scenario);
        } else if (gCacheMaxStepId >= 0) {
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
        gCacheSnapStepIds.reserve(static_cast<size_t>(expected_snaps));
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
        last_reported_bucket = 0;
        next_report_bucket   = 5;
        next_report_step     = (int)(((long long)steps * (long long)next_report_bucket + 99LL) / 100LL);
        if (next_report_step < 1) next_report_step = 1;

        if (cache_enabled) {
            auto old_flags = cout.flags();
            auto old_prec  = cout.precision();
            cout << "  Cached: "  << fixed << std::setprecision(2) << check_mb(scenario) << " MB" << endl;
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
                gCacheSnapStepIds.push_back((long long)i);

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

                //cout << "cache_enabled2:" << cache_enabled << endl;

                if (cache_enabled) {
                    auto old_flags = cout.flags();
                    auto old_prec  = cout.precision();
                    cout << "  Cache written: " << fixed << std::setprecision(2) << check_mb(scenario) << " MB" << endl;
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

        if (kCacheFlushEverySteps > 0 && ((i + 1) % kCacheFlushEverySteps == 0)) {
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
        gCacheSnapStepIds.clear();

        // Post-simulation: apply heat brightness if glow mode is on.
        // When glow_mode is off, temps were already written as 0 during run_to_cache,
        // so no post-processing is needed — skip the expensive load+rewrite cycle.
        if (scenario->glow_mode) {
            cout << "Post-processing cached states (glow mode)..." << endl;
            shared_ptr<snapshots> cached_states = run_from_cache(scenario);
            apply_heat_brightness_and_pack_rgb(cached_states, scenario);
            rewrite_cache(scenario, cached_states);
        } else {
            cout << "Glow mode disabled. Temps already zeroed during cache write." << endl;
        }
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

    // Base stride writes + optional forced final snapshot (if not already on stride)
    const bool need_force_last = (((total_snaps - 1) % stride) != 0);
    const int written_snaps_base = (total_snaps + stride - 1) / stride;
    const int written_snaps = written_snaps_base + (need_force_last ? 1 : 0);

    const double compression = total_snaps > 0 ? (1.0 - (double)written_snaps / total_snaps) * 100.0 : 0.0;


    file << std::fixed << std::setprecision(15);

    if (!gCacheInitializedThisRun) {
        CacheSchema::write_header(file);
    }

    int last_reported_bucket = -5;
    int next_report_step = 0;
    int written = 0;

    const bool zero_temp = !scenario->glow_mode;

    auto write_snapshot_index = [&](int idx) {
        const long long step_id =
            (idx >= 0 && idx < (int)gCacheSnapStepIds.size())
                ? gCacheSnapStepIds[(size_t)idx]
                : (long long)idx;

        for (int j = 0; j < (int)particle_states->snaps[idx]->particle_list.size(); j++) {
            const auto& sp = particle_states->snaps[idx]->particle_list[j];
            if (!sp) continue;
            CacheSchema::write_row(file, step_id, *sp, zero_temp);
        }
    };

    int last_written_i = -1;

    for (int i = 0; i < total_snaps; i += stride) {
        if (written_snaps > 0 && written >= next_report_step) {
            int pct = (int)(((long long)(written + 1) * 100LL) / (long long)written_snaps);
            int bucket = (pct / 20) * 20;
            if (bucket > last_reported_bucket) {
                last_reported_bucket = bucket;
            }
            int nb = last_reported_bucket + 5;
            next_report_step = (nb > 100) ? written_snaps + 1
                : (int)(((long long)written_snaps * (long long)nb + 99LL) / 100LL);
        }
        ++written;

        write_snapshot_index(i);
        last_written_i = i;
    }

    // ---- ensure final buffered snapshot is always cached ----
    const int last_idx = total_snaps - 1;
    if (last_idx >= 0 && last_written_i != last_idx) {
        if (written_snaps > 0 && written >= next_report_step) {
            int pct = (int)(((long long)(written + 1) * 100LL) / (long long)written_snaps);
            int bucket = (pct / 20) * 20;
            if (bucket > last_reported_bucket) {
                last_reported_bucket = bucket;
            }
            int nb = last_reported_bucket + 5;
            next_report_step = (nb > 100) ? written_snaps + 1
                : (int)(((long long)written_snaps * (long long)nb + 99LL) / 100LL);
        }
        ++written;

        write_snapshot_index(last_idx);
        last_written_i = last_idx;
    }

    file.close();

    if (!gCacheInitializedThisRun) {
        gCacheInitializedThisRun = true;
    }
    gCacheStepIdOffset += (long long)total_snaps;
}

void Engine::rewrite_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {
    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    // First, collect the step_ids from the existing cache file
    std::vector<long long> step_ids;
    {
        typedef io::trim_chars<' ', '\t'> TrimPolicy;
        typedef io::double_quote_escape<',', '\"'> QuotePolicy;

        io::CSVReader<CacheSchema::kColumnCount, TrimPolicy, QuotePolicy> in(file_name);
        CacheSchema::read_header(in);

        CacheSchema::Row row;
        long long last_sid = -1;
        while (in.read_row(
            row.step_id, row.particle_id,
            row.r, row.g, row.b,
            row.x, row.y, row.z,
            row.vx, row.vy, row.vz,
            row.m, row.rad, row.temp, row.rest))
        {
            long long sid = CacheSchema::parse_step_id(row);
            if (sid != last_sid) {
                step_ids.push_back(sid);
                last_sid = sid;
            }
        }
    }

    // Verify step_id count matches snapshot count
    if ((int)step_ids.size() != (int)particle_states->snaps.size()) {
        cout << "WARNING: step_id count (" << step_ids.size()
             << ") != snapshot count (" << particle_states->snaps.size()
             << "). Cache rewrite aborted." << endl;
        return;
    }

    // Rewrite the cache file with updated particle data
    ofstream file(file_name, std::ios::out | std::ios::trunc);
    if (!file.is_open()) {
        cout << "Failed to rewrite cache file: " << file_name << endl;
        return;
    }

    file << std::fixed << std::setprecision(15);
    CacheSchema::write_header(file);

    for (int i = 0; i < (int)particle_states->snaps.size(); ++i) {
        const long long sid = step_ids[i];
        const auto& snap = particle_states->snaps[i];
        for (const auto& p : snap->particle_list) {
            if (!p) continue;
            CacheSchema::write_row(file, sid, *p);
        }
    }

    file.close();
    cout << "Cache rewritten with post-processed data: " << file_name << endl;
    cout << endl;
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

bool Engine::inspect_cache(const std::string& name_or_path) {
    const std::string file_name = cache_file_from_name_or_path(name_or_path);

    // Basic open check
    std::ifstream fin(file_name, std::ios::binary);
    if (!fin.is_open()) {
        std::cout << "Cache file not found/unreadable: " << file_name << std::endl;
        gCacheMaxStepId = -1;
        return false;
    }

    long long old_bytes = 0;
    try { old_bytes = (long long)std::filesystem::file_size(file_name); }
    catch (...) { old_bytes = 0; }

    // Read header
    std::string line;
    if (!std::getline(fin, line)) {
        std::cout << "Cache file is empty or unreadable: " << file_name << std::endl;
        gCacheMaxStepId = -1;
        fin.close();
        filesystem::remove(file_name);

        cout <<  "Deleted empty or unreadable cache file: " << file_name << endl;       // Delete cache file
        try {
            std::filesystem::remove(file_name);
            std::cout << "Deleted empty or unreadable cache file: " << file_name << std::endl;
        } catch (...) {
            std::cout << "Failed to delete cache file: " << file_name << std::endl;
        }
        return false;
    }

    auto strip_cr = [&](std::string& s) {
        if (!s.empty() && s.back() == '\r') s.pop_back();
    };
    strip_cr(line);

    // Parse step_id from prefix up to first comma (robust, trims spaces/tabs)
    auto parse_step_prefix = [&](const std::string& s, long long& out_sid) -> bool {
        const size_t comma = s.find(',');
        if (comma == std::string::npos) return false;

        size_t a = 0;
        while (a < comma && (s[a] == ' ' || s[a] == '\t')) ++a;

        size_t b = comma;
        while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t')) --b;

        if (a >= b) return false;

        try {
            out_sid = std::stoll(s.substr(a, b - a));
            return true;
        } catch (...) {
            return false;
        }
    };

    const int expected_cols = CacheSchema::kColumnCount;

    bool corrupt_found = false;
    long long corrupt_step_id = -1;
    long long corrupt_line_no = -1;

    bool any_good = false;
    long long last_good_step_id = -1;
    long long max_step_id = std::numeric_limits<long long>::min();

    long long line_no = 1; // header line
    while (std::getline(fin, line)) {
        ++line_no;
        strip_cr(line);

        // Ignore blank lines (don’t treat trailing blank as corruption)
        if (line.empty()) continue;

        const int cols = 1 + (int)std::count(line.begin(), line.end(), ',');
        long long sid = 0;
        const bool sid_ok = parse_step_prefix(line, sid);

        if (cols != expected_cols) {
            corrupt_found = true;
            corrupt_line_no = line_no;

            if (sid_ok) {
                corrupt_step_id = sid;
            } else if (last_good_step_id >= 0) {
                corrupt_step_id = last_good_step_id + 1;
            } else {
                corrupt_step_id = -1;
            }
            break;
        }

        if (sid_ok) {
            any_good = true;
            last_good_step_id = sid;
            if (sid > max_step_id) max_step_id = sid;
        }
    }
    fin.close();

    // If no corruption: just set max and return (same semantics as before)
    if (!corrupt_found) {
        if (!any_good || max_step_id == std::numeric_limits<long long>::min()) {
            //std::cout << "Cache file is empty: " << file_name << std::endl;
            gCacheMaxStepId = -1;
        } else {
            gCacheMaxStepId = max_step_id;
        }
        return true;
    }

    // If we cannot determine the corrupt timestep safely, DO NOT overwrite.
    if (corrupt_step_id < 0) {
        std::cout << "CACHE CORRUPTION DETECTED but could not determine cutoff step_id. "
                  << "File left unchanged. first_bad_line=" << corrupt_line_no
                  << " file=" << file_name << std::endl;
        // best effort: keep existing max (may be unreliable if file corrupt)
        gCacheMaxStepId = (any_good && max_step_id != std::numeric_limits<long long>::min()) ? max_step_id : -1;
        return false;
    }

    std::cout << "CACHE CORRUPTION DETECTED: file=" << file_name
              << " first_bad_line=" << corrupt_line_no
              << " cutoff_step_id(exclusive)=" << corrupt_step_id << std::endl;

    // Rebuild using the normal grouped-by-step semantics up to (excluding) corrupt_step_id
    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;

    const std::string tmp_name = file_name + ".tmp_repair";

    std::ofstream out(tmp_name, std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        std::cout << "Failed to create temp repair cache: " << tmp_name << std::endl;
        gCacheMaxStepId = -1;
        return true;
    }
    out << std::fixed << std::setprecision(15);
    CacheSchema::write_header(out);

    std::shared_ptr<snapshots> repaired_states = std::make_shared<snapshots>(); // as requested
    CacheSchema::Row row;

    long long last_written_step = -1;
    bool have_step = false;
    long long current_step = -1;
    std::unique_ptr<Particles> current_particles = std::make_unique<Particles>();

    auto flush_step = [&](long long sid) {
        if (!have_step) return;

        // Store in a snapshots object "the usual way" (Particles per step)
        repaired_states->snaps.push_back(std::move(current_particles));

        // Write that step out
        const auto& snap = repaired_states->snaps.back();
        for (int j = 0; j < (int)snap->particle_list.size(); ++j) {
            const auto& sp = snap->particle_list[j];
            if (!sp) continue;
            CacheSchema::write_row(out, sid, *sp);
        }

        last_written_step = sid;

        // Keep memory bounded while still using snapshots structure
        repaired_states->snaps.clear();
        current_particles = std::make_unique<Particles>();
    };

    try {
        io::CSVReader<CacheSchema::kColumnCount, TrimPolicy, QuotePolicy> in(file_name);
        CacheSchema::read_header(in);

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
            if (sid >= corrupt_step_id) break;

            if (!have_step) {
                have_step = true;
                current_step = sid;
            } else if (sid != current_step) {
                flush_step(current_step);
                current_step = sid;
            }

            auto p = std::make_shared<Particle>();
            CacheSchema::parse_particle_fields(row, *p);
            current_particles->particle_list.push_back(p);
        }
    } catch (...) {
        // If the CSV reader trips on the corrupted line, that's OK: we stop where we are.
        // The requirement is "up to but excluding this timestep".
    }

    if (have_step && current_step >= 0 && current_step < corrupt_step_id) {
        flush_step(current_step);
    }

    out.close();


    // Replace original with repaired
    bool replaced = false;
    try {
        std::filesystem::remove(file_name);
    } catch (...) {}

    try {
        std::filesystem::rename(tmp_name, file_name);
        replaced = true;
    } catch (...) {
        try {
            std::filesystem::copy_file(tmp_name, file_name,
                                       std::filesystem::copy_options::overwrite_existing);
            std::filesystem::remove(tmp_name);
            replaced = true;
        } catch (...) {
            replaced = false;
        }
    }

    if (!replaced) {
        std::cout << "FAILED to replace cache with repaired version. Temp remains: " << tmp_name << std::endl;
        gCacheMaxStepId = -1;
        return false;
    }

    long long new_bytes = 0;
    try { new_bytes = (long long)std::filesystem::file_size(file_name); }
    catch (...) { new_bytes = 0; }

    double pct_deleted = 0.0;
    if (old_bytes > 0 && new_bytes >= 0 && new_bytes <= old_bytes) {
        pct_deleted = (1.0 - (double)new_bytes / (double)old_bytes) * 100.0;
    }

    gCacheMaxStepId = last_written_step;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "CACHE REPAIRED: deleted ~" << pct_deleted << "% due to corruption." << std::endl;
    std::cout << "Rewritten cache saved: " << file_name << std::endl;

    std::cout << "Press Enter to exit..." << std::endl;
    std::string _;
    std::getline(std::cin >> std::ws, _);

    throw std::runtime_error("Cache repaired due to corruption. Please restart the program.");
}

void Engine::apply_heat_brightness_and_pack_rgb(shared_ptr<snapshots> snaps, shared_ptr<scenario> scenario) {

    // Config: cutoff, gamma (cast once to double)
    const double HEAT_CUTOFF_FRAC   = static_cast<double>(scenario->heat_cutoff_frac);
    const double HEAT_GAMMA         = static_cast<double>(scenario->heat_gamma);
    const double LOWER_OUTLIER_FRAC = 0.005;// cut bottom 2% coldest samples for normalization
    const double UPPER_OUTLIER_FRAC = 0.005;//cut top 2% hottest samples for normalization

    auto& frames = snaps->snaps;
    const int nf = static_cast<int>(frames.size());

    double max_temp = 0.0;
    double min_temp = std::numeric_limits<double>::max();
    std::vector<double> temp_samples;

    for (int i = 0; i < nf; ++i) {
        auto& plist = frames[i]->particle_list;
        const int np = static_cast<int>(plist.size());
        temp_samples.reserve(temp_samples.size() + static_cast<size_t>(np));
        for (int j = 0; j < np; ++j) {
            const auto& p = plist[j];
            const double t = static_cast<double>(p->temp);
            temp_samples.push_back(t);
            if (t > max_temp) max_temp = t;
            if (t < min_temp) min_temp = t;
        }
    }

    if (frames.empty() || frames[0]->particle_list.empty()) min_temp = 0.0;

    double effective_min_temp = min_temp;
    double effective_max_temp = max_temp;

    if (!temp_samples.empty()) {
        const size_t n = temp_samples.size();

        if (LOWER_OUTLIER_FRAC > 0.0 && LOWER_OUTLIER_FRAC < 1.0) {
            const size_t idx_low = static_cast<size_t>(std::floor(LOWER_OUTLIER_FRAC * static_cast<double>(n - 1)));
            std::nth_element(temp_samples.begin(), temp_samples.begin() + idx_low, temp_samples.end());
            effective_min_temp = temp_samples[idx_low];
        }

        if (UPPER_OUTLIER_FRAC > 0.0 && UPPER_OUTLIER_FRAC < 1.0) {
            const size_t idx_high = static_cast<size_t>(std::floor((1.0 - UPPER_OUTLIER_FRAC) * static_cast<double>(n - 1)));
            std::nth_element(temp_samples.begin(), temp_samples.begin() + idx_high, temp_samples.end());
            effective_max_temp = temp_samples[idx_high];
        }

        if (effective_min_temp < min_temp) effective_min_temp = min_temp;
        if (effective_max_temp > max_temp) effective_max_temp = max_temp;
        if (effective_max_temp < effective_min_temp) {
            effective_min_temp = min_temp;
            effective_max_temp = max_temp;
        }
    }

    const double temp_range = effective_max_temp - effective_min_temp;

    // Fast path: pack base rgb only (no heat tinting)
    auto pack_base_rgb = [&]() {
        for (int i = 0; i < nf; ++i) {
            auto& plist = frames[i]->particle_list;
            const int np = static_cast<int>(plist.size());
            for (int j = 0; j < np; ++j) {
                auto& p = plist[j];
                const double r = static_cast<double>(p->r);
                const double g = static_cast<double>(p->g);
                const double b = static_cast<double>(p->b);

                const int r255 = static_cast<int>(r * 255.0);
                const int g255 = static_cast<int>(g * 255.0);
                const int b255 = static_cast<int>(b * 255.0);
                p->rgb = (r255 << 16) | (g255 << 8) | (b255);
            }
        }
    };
        const double raw_temp_span = max_temp - min_temp;


        // Profile log: 11 timeline points (0%, 10%, ..., 100%) of mass-weighted mean temperature.
        {
            cout << "Brightness profile (mass-weighted mean temp, 0%..100%):" << endl;
            for (int k = 0; k <= 10; ++k) {
                int idx = 0;
                if (nf > 1) {
                    idx = static_cast<int>(std::llround((static_cast<double>(nf - 1) * static_cast<double>(k)) / 10.0));
                }

                double sum_m = 0.0;
                double sum_mt = 0.0;
                double sum_t = 0.0;
                int count = 0;

                if (nf > 0) {
                    auto& plist = frames[idx]->particle_list;
                    const int np = static_cast<int>(plist.size());
                    for (int j = 0; j < np; ++j) {
                        const auto& p = plist[j];
                        const double m = static_cast<double>(p->m);
                        const double t = static_cast<double>(p->temp);

                        if (m > 0.0) {
                            sum_m += m;
                            sum_mt += m * t;
                        }
                        sum_t += t;
                        ++count;
                    }
                }

                double mean_t = 0.0;
                if (sum_m > 0.0) {
                    mean_t = sum_mt / sum_m;
                } else if (count > 0) {
                    mean_t = sum_t / static_cast<double>(count);
                }

                cout << "  " << (k * 10) << "%: "
                     << std::scientific << std::setprecision(6) << mean_t
                     << std::fixed << endl;
            }
        }

    if (max_temp < 1e-5) {
        cout << "RGB values calculated (max temp < 0.00001, no heating applied)." << endl;
        pack_base_rgb();
        return;
    }
    if (temp_range <= 0.0) {
        cout << "RGB values calculated (uniform temp: " << max_temp << ")." << endl;
        pack_base_rgb();
        return;
    }

        cout << "Temperature range (raw): "
            << std::scientific << std::setprecision(6) << raw_temp_span
            << ", (clipped): " << temp_range
            << std::fixed << endl;

    const double cutoff_temp = effective_min_temp + HEAT_CUTOFF_FRAC * temp_range;
    const double denom = (effective_max_temp - cutoff_temp);
    if (denom <= 0.0) {
        pack_base_rgb();
        return;
    }

    const bool gamma_is_1 = (std::abs(HEAT_GAMMA - 1.0) < 1e-12);

    int next_pct_threshold = 10;
    for (int i = 0; i < nf; ++i) {
        int pct = (nf > 1) ? static_cast<int>((static_cast<double>(i) / (nf - 1)) * 100.0) : 100;
        if (pct >= next_pct_threshold) {
            //cout << "  Brightening progress: " << next_pct_threshold << "%" << endl;
            next_pct_threshold += 10;
        }
        auto& plist = frames[i]->particle_list;
        const int np = static_cast<int>(plist.size());
        for (int j = 0; j < np; ++j) {
            auto& p = plist[j];

            const double t = static_cast<double>(p->temp);

            double fraction = 0.0;
            if (t > cutoff_temp) {
                fraction = (t - cutoff_temp) / denom;
                if (fraction < 0.0) fraction = 0.0;
                if (fraction > 1.0) fraction = 1.0;

                if (!gamma_is_1) {
                    fraction = std::pow(fraction, HEAT_GAMMA);
                }
            }

            double r = static_cast<double>(p->r);
            double g = static_cast<double>(p->g);
            double b = static_cast<double>(p->b);

            const double brighten = fraction * 0.25;

            r = r + (1.0 - r) * fraction + (1.0 - r) * brighten;
            g = g * (1.0 - fraction) + (1.0 - g) * brighten;
            b = b * (1.0 - fraction) + (1.0 - b) * brighten;

            if (r > 1.0) r = 1.0;
            if (g > 1.0) g = 1.0;
            if (b > 1.0) b = 1.0;

            p->r = r;
            p->g = g;
            p->b = b;

            const int r255 = static_cast<int>(r * 255.0);
            const int g255 = static_cast<int>(g * 255.0);
            const int b255 = static_cast<int>(b * 255.0);
            p->rgb = (r255 << 16) | (g255 << 8) | (b255);
        }
    }
    if (next_pct_threshold <= 100) {
        cout << "  Brightening progress: 100%" << endl;
    }
}