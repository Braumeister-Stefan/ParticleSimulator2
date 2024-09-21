#include "Model.h"

// Constructor
ParticleModel::ParticleModel() {
    

    unique_ptr<Interfacer> interfacer;
    unique_ptr<Engine> engine;
    unique_ptr<Metrics> metrics;
    unique_ptr<Plotter> plotter;
    unique_ptr<ObjHandler> obj_handler;
    unique_ptr<ParticleDumper> dumper;

}


// Destructor
ParticleModel::~ParticleModel() {
    // The unique_ptr will automatically clean up

}
