#ifndef PHYS_ENGINE_H
#define PHYS_ENGINE_H

//Standard libraries
#include <iostream>
#include <memory>

//Internal libraries
#include "InitStructs.h"


using namespace std;

class Engine {
public:

    //declare static member variables
    static double dt;


    // Constructor
    Engine();

    // Destructor
    ~Engine();

    //Function to run the simulation and return snapshots of each time step
    shared_ptr<snapshots> run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles);

    //Function to update the particles
    void update_particles(shared_ptr<Particles> particles);

    //Functions to resolve overlap between particles
    bool resolve_overlap(shared_ptr<Particles> particles);
    
    //functions to resolve overlap between a pair of particles
    double backtrack_pair(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    double impulse_pair(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);

    //Functions to resolve collissions between particles
    void resolve_collisions(shared_ptr<Particles> particles);

    



    void resolve_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    bool check_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);



    //Functions to resolve gravitational attraction between particles
    void resolve_gravity_euler(shared_ptr<Particles> particles);
    void resolve_gravity_verlet(shared_ptr<Particles> particles);

    //Function to update the locations of the particles
    void update_locations(shared_ptr<Particles> particles, shared_ptr<backed_scaler> scaler);

    //Function to interact with cache
    void run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> snapshots);
    bool cache_exists(shared_ptr<scenario> scenario);
    shared_ptr<snapshots> run_from_cache(shared_ptr<scenario> scenario);

    //Function to validate TE
    double calc_TE(shared_ptr<Particles> particles);
    
};

#endif
