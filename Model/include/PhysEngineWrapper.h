// =======================
// PhysEngineWrapper.h (UPDATED)
//   - wrapper owns reporting + caching; core owns physics stepping
// =======================

#ifndef PHYSENGINEWRAPPER_H
#define PHYSENGINEWRAPPER_H

#include <memory>
#include <string>

#define CSV_IO_NO_THREAD
#include "3party/csv.h"

#include "InitStructs.h"
#include "PhysEngineCore.h"

using namespace std;

class Engine {
public:
    Engine() = default;
    ~Engine() = default;

    // main
    void run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles);

    // cache I/O
    void run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states);
    bool cache_exists(shared_ptr<scenario> scenario);
    shared_ptr<snapshots> run_from_cache(shared_ptr<scenario> scenario);

    // optional dt access
    static high_prec dt() { return EngineCore::dt; }

    // expose debug mode flag for other components (e.g. plotter)
    static bool debug_mode(std::shared_ptr<scenario> scenario);

    // make light snapshot
    static std::unique_ptr<Particles> make_light_snapshot(const Particles& src);

private:
    // initializes/truncates cache and writes schema header; uses cache_exists()
    bool initiate_cache(shared_ptr<scenario> scenario);

    bool inspect_cache(const string& scenario_name);

    // NEW: helper for progress logging (cache file size in MB)
    double check_mb(shared_ptr<scenario> scenario);

    // Heat-brightness RGB packing (moved from Plotter)
    void apply_heat_brightness_and_pack_rgb(shared_ptr<snapshots> snaps, shared_ptr<scenario> scenario);

    // Rewrite entire cache file with current snapshot data (e.g. after RGB post-processing)
    void rewrite_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states);

    EngineCore core;
};

#endif // PHYSENGINEWRAPPER_H
