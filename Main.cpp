// This project allows the user to simulate collissions and gravitational forces between spheres in 2D space


#include "Model/Model.h"
#include "Model/include/InitStructs.h"


#include <memory>



using namespace std;

int main() {

    //initialize the model
    ParticleModel model;

    //select the scenario
    shared_ptr<scenario> selected_scenario = model.interfacer->select_scenario();

    //render and store the objects at t0

    shared_ptr<Particles> particles = model.obj_handler-> process_objs(selected_scenario);

    //run the model and save the snapshots
    shared_ptr<snapshots> particle_states = model.engine->run(selected_scenario, particles);

    //validate the simulation and generate metrics

    shared_ptr<test_metrics> metrics = model.metrics->compute_metrics(selected_scenario, particle_states);

    //visualize the model results using the snapshots
    model.plotter->plot_run(selected_scenario, particle_states, metrics);

    return 0;
    
}