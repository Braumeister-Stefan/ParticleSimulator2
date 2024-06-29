#include "../include/PhysEngine.h"
#include "../include/InitStructs.h"

#include <memory>

//namespaces
using namespace std;

// Constructor
Engine::Engine() {
    cout << "Engine initialized." << endl;
}

// Destructor
Engine::~Engine() {
    cout << "Engine destroyed." << endl;
}

// Example method
void Engine::run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles) {
    cout << "Engine is initialized." << endl;

    //confirm that the run should start

    //print the list of particle ids

    cout << "The following particles are available:" << endl;

    for (int i = 0; i < particles->particle_list.size(); i++) {
        cout << particles->particle_list[i]->particle_id << endl;
        cout <<"values should follow"<< endl;
        cout << particles->particle_list[i]->x << endl;
        cout << particles->particle_list[i]->vx << endl;
        cout << particles->particle_list[i]->vy << endl;


    }

    cout << scenario->name << " is ready to start" << endl << endl;
    cout << "Press enter to continue." << endl;
    cin.ignore();
    cin.get();



    

}