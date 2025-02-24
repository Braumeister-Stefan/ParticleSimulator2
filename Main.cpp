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

    particle_states = model.metrics->compute_metrics(selected_scenario, particle_states);

    int returnCode = system(
        
        R"(C:\\Users\\smdw1\\anaconda3\\python.exe C:\\Users\\smdw1\\OneDrive\\Bureaublad\\Development\\Projects\\cpp\\ParticleSimulator2\\Misc\\Energy_plots\\main.py)"

    );

    if (returnCode == 0) {
        cout << "Python script executed successfully." << endl;
    } else {
        cout << "Python script failed to execute." << endl;
    }

    //5. visualize the model results using the snapshots
    model.plotter->plot_run(selected_scenario, particle_states);

    return 0;
    
}