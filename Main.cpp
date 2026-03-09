//Author: Stephane Mertens de Wilmars

// This project allows the user to simulate collissions and gravitational forces between spheres in 2D space
//Project to be built from the Main.cpp file. For more information, please refer to the README.md file.


//Standard libraries
#include <memory>
#include <iostream>
#include <random>
#include <cstdlib>

//Internal libraries
#include "Model/Model.h"
#include "Model/include/InitStructs.h"

//namespaces
using namespace std;

int main() {

    //1. initialize the model
    ParticleModel model;
    
    //2. select the scenario
    shared_ptr<scenario> selected_scenario = model.interfacer->select_scenario();

    //3. run the model
    shared_ptr<snapshots> particle_states;

    // render and store the objects at t0
    shared_ptr<Particles> particles = model.obj_handler->process_objs(selected_scenario);

    //create a snapshots object and store the initial state before any simulation steps
    shared_ptr<snapshots> initial_states = make_shared<snapshots>();
    initial_states->snaps.push_back(Engine::make_light_snapshot(*particles));

    model.engine->run(selected_scenario, particles);

    //4. retrieve states from cache

    if (model.engine->cache_exists(selected_scenario)) {
        particle_states = model.engine->run_from_cache(selected_scenario);

        particle_states  = model.metrics->compute_metrics(selected_scenario, particle_states);
    } else {
        cout << "No cache found for this scenario. Unable to load particle states." << endl;
        return 1;
    }



    if (selected_scenario->save_obj != "FALSE") {
        //save the last particle state as a complex object
        shared_ptr<Particles> final_state = particle_states->snaps.back();
        model.obj_handler->state_to_cache(final_state, selected_scenario->save_obj);
        
    }

    
    //5. visualize the model results using the snapshots
    model.plotter->plot_run(selected_scenario, particle_states);

    return 0;
    
}