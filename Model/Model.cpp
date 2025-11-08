

//Internal libraries
#include "Model.h"

// Constructor

ParticleModel::ParticleModel() {
    // Initialize the MEMBER variables using make_unique
    interfacer = make_unique<Interfacer>();
    engine = make_unique<Engine>();
    metrics = make_unique<Metrics>();
    plotter = make_unique<Plotter>();
    obj_handler = make_unique<ObjHandler>();
    cout << "DEBUG: ParticleModel constructed, engine address: " << engine.get() << endl;
}

// Destructor
ParticleModel::~ParticleModel() {
    // The unique_ptr will automatically clean up

}
