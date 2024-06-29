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

    //run the model
    model.engine->run(selected_scenario, particles);

    //visualize the model results
    model.plotter->plotcast();

    return 0;
    


}