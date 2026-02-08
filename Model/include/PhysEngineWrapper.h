// =======================
// PhysEngineWrapper.h (UPDATED)
//   - add cache_exists / run_from_cache / run_to_cache back onto Engine (wrapper)
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
    shared_ptr<snapshots> run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles);

    void run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states);
    bool cache_exists(shared_ptr<scenario> scenario);
    shared_ptr<snapshots> run_from_cache(shared_ptr<scenario> scenario);

    // optional dt access
    static high_prec dt() { return EngineCore::dt; }

    // expose debug mode flag for other components (e.g. plotter)
    static bool debug_mode();

private:
    EngineCore core;
};

#endif // PHYSENGINEWRAPPER_H
