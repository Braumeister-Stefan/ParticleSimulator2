// This project allows the user to simulate collissions and gravitational forces between spheres in 2D space


#include "Model/Model.h"
#include "Model/include/InitStructs.h"


#include <memory>



using namespace std;

int main() {

    //1. initialize the model
    ParticleModel model;

    //2. select the scenario
    shared_ptr<scenario> selected_scenario = model.interfacer->select_scenario();

    //3. run the model
    shared_ptr<snapshots> particle_states;

    //3a. Run from cache if allowed
    if(selected_scenario->try_cache && model.engine->cache_exists(selected_scenario)) {

       particle_states = model.engine->run_from_cache(selected_scenario);

    } else {
    //3b. Run the model from scratch and save the snapshots

        //3b.1. render and store the objects at t0

        shared_ptr<Particles> particles = model.obj_handler-> process_objs(selected_scenario);

        //3b.2. run the model
        particle_states = model.engine->run(selected_scenario, particles);
        
    }

    //4. validate the simulation and generate metrics

    shared_ptr<test_metrics> metrics = model.metrics->compute_metrics(selected_scenario, particle_states);

    //5. visualize the model results using the snapshots
    model.plotter->plot_run(selected_scenario, particle_states, metrics);

    return 0;
    
}