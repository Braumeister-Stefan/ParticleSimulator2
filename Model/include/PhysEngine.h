#ifndef PHYS_ENGINE_H
#define PHYS_ENGINE_H

#include "InitStructs.h"

#include <iostream>
#include <memory>

using namespace std;

class Engine {
public:
    // Constructor
    Engine();

    // Destructor
    ~Engine();

    //Function to run the simulation and return snapshots of each time step
    shared_ptr<snapshots> run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles);

    //Function to update the particles
    void update_particles(shared_ptr<Particles> particles);

    //Functions to resolve collissions between particles
    shared_ptr<backed_scaler> resolve_collisions(shared_ptr<Particles> particles);

    double backtrack_pair(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    void resolve_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    bool check_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);



    //Function to resolve gravitational attraction between particles
    void resolve_gravity(shared_ptr<Particles> particles);

    //Function to update the locations of the particles
    void update_locations(shared_ptr<Particles> particles, shared_ptr<backed_scaler> scaler);

    //Function to interact with cache
    void run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> snapshots);
    bool cache_exists(shared_ptr<scenario> scenario);
    shared_ptr<snapshots> run_from_cache(shared_ptr<scenario> scenario);
    
};

#endif
