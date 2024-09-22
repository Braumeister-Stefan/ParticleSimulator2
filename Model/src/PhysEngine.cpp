#include "../include/PhysEngine.h"
#include "../include/InitStructs.h"

#include <memory>
#include <iostream>
#include <cmath>
#include <vector>



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

// Run the simulation and return snapshots of each time step
shared_ptr<snapshots> Engine::run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles) {
    cout << "Engine is initialized." << endl;

    //0. initialize the snapshots object
    shared_ptr<snapshots> particle_states = make_shared<snapshots>();


    //confirm that the run should start. If a dump file is available, ask the user if they want to use it, if yes save to snapshots object and return it.
    //to be implemented


    cout << scenario->name << " is ready to start" << endl << endl;
    cout << "Press enter to continue." << endl;
    cin.ignore();
    cin.get();




    //loop through the steps as defined in the scenario, do nothing for now

    for (int i = 0; i < scenario->steps; i++) {

        cout << "Starting step " << i << "..." << endl;

        //1.add current state of particles to the snapshots object
        auto particles_copy = make_unique<Particles>(*particles);

        cout << "Continuing step " << i << "..." << endl;

        particle_states->snaps.push_back(move(particles_copy));
        //if instead we want to add it to the first element of the snapshots object, we can use the following line

        //2. Update the state of the particles

        cout << "Updating particles..." << endl;
        
        update_particles(particles);

        cout << "Step " << i << " completed." << endl;

    }

    //confirm that the run has ended

    cout << scenario->name << " simulation completed." << endl;

    //save down snapshots to the csv using the dumper
    //dumper->dump_snapshots(snapshots); TO BE IMPLEMENTED

    return particle_states;


}

//update particles steps to be implemented

//2) debug the code

//3) make complex object (circle)


void Engine::update_particles(shared_ptr<Particles> particles) {
    //this function will update the particles in the particles object

    
    //1 resolve collissions

    cout << "Resolving collissions..." << endl;
    shared_ptr<backed_scaler> scaler = resolve_collisions(particles);

    //2. resolve gravitational attraction

    //cout << "Resolving gravity..." << endl;
    //resolve_gravity(particles);

    //3. update locations with velocities

    cout << "Updating locations..." << endl;
    update_locations(particles, scaler);


}

shared_ptr<backed_scaler> Engine::resolve_collisions(shared_ptr<Particles> particles) {
    //this function will resolve collissions between particles

    //1.initialize the scaler object
    shared_ptr<backed_scaler> scaler = make_shared<backed_scaler>(); 
    scaler->scaler.reserve(particles->particle_list.size()); //reserve the amount of scalers needed
    for (int i = 0; i < particles->particle_list.size(); i++) {
        scaler->scaler.push_back(1);
    }

    cout <<"Scalers initialized." << endl;

    //2.check for collissions between particles and resolve them



    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            //2a. check for collission between particle i and particle j

            bool collission = check_collission(particles->particle_list[i], particles->particle_list[j]);
            //if collission detected, backtrack the particles to the point of collission and store the amount of time of the timestep that was left
            
            
            if (collission) {
                //2b. backtrack the particles
                 int scaler_i = backtrack_pair(particles->particle_list[i], particles->particle_list[j]);

                //2c. store the scaler for the particle

                scaler->scaler[i] = scaler_i;

                //2d. resolve the collission

                resolve_collission(particles->particle_list[i], particles->particle_list[j]);

                
            }
            

        }
    }

    
    return scaler;
}


bool Engine::check_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    //this function will check if two particles are colliding

    //1. calculate the distance between the particles
    double distance = sqrt(pow(particle1->x - particle2->x, 2) + pow(particle1->y - particle2->y, 2) + pow(particle1->z - particle2->z, 2));

    //2. check if the distance is smaller than the sum of the radii of the particles
    if (distance < particle1->rad + particle2->rad) {
        return true;
    }
    else {
        return false;
    }
}

int Engine::backtrack_pair(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    //this function will backtrack two particles to the point of collission

    //1. calculate the distance between the particles
    double distance = sqrt(pow(particle1->x - particle2->x, 2) + pow(particle1->y - particle2->y, 2) + pow(particle1->z - particle2->z, 2));

    //2. calculate the time of the timestep that was left
    double time_left = (particle1->rad + particle2->rad - distance) / (particle1->vx + particle2->vx + particle1->vy + particle2->vy + particle1->vz + particle2->vz);

    //3. backtrack the particles
    particle1->x -= particle1->vx * time_left;
    particle1->y -= particle1->vy * time_left;
    particle1->z -= particle1->vz * time_left;

    particle2->x -= particle2->vx * time_left;
    particle2->y -= particle2->vy * time_left;
    particle2->z -= particle2->vz * time_left;

    return time_left;
}

void Engine::resolve_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    //this function will resolve the collission between two particles

    //1. calculate the distance between the particles
    double distance = sqrt(pow(particle1->x - particle2->x, 2) + pow(particle1->y - particle2->y, 2) + pow(particle1->z - particle2->z, 2));

    //2. calculate the normal vector of the collission
    double nx = (particle2->x - particle1->x) / distance;
    double ny = (particle2->y - particle1->y) / distance;
    double nz = (particle2->z - particle1->z) / distance;

    //3. calculate the relative velocity of the particles
    double vrelx = particle2->vx - particle1->vx;
    double vrely = particle2->vy - particle1->vy;
    double vrelz = particle2->vz - particle1->vz;

    //4. calculate the dot product of the relative velocity and the normal vector
    double dot = vrelx * nx + vrely * ny + vrelz * nz;

    //5. calculate the impulse
    double impulse = 2 * dot / (1 / particle1->m + 1 / particle2->m);

    //6. update the velocities of the particles
    particle1->vx += impulse / particle1->m * nx;
    particle1->vy += impulse / particle1->m * ny;
    particle1->vz += impulse / particle1->m * nz;

    particle2->vx -= impulse / particle2->m * nx;
    particle2->vy -= impulse / particle2->m * ny;
    particle2->vz -= impulse / particle2->m * nz;

}


void Engine::resolve_gravity(shared_ptr<Particles> particles) {
    //this function will resolve the gravitational attraction between particles

    //1. loop through the particles
    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = 0; j < particles->particle_list.size(); j++) {
            if (i != j) {
                //1a. calculate the distance between the particles
                double distance = sqrt(pow(particles->particle_list[j]->x - particles->particle_list[i]->x, 2) + pow(particles->particle_list[j]->y - particles->particle_list[i]->y, 2) + pow(particles->particle_list[j]->z - particles->particle_list[i]->z, 2));

                //1b. calculate the force of gravity
                double force = 6.674 * pow(10, -11) * particles->particle_list[i]->m * particles->particle_list[j]->m / pow(distance, 2);

                //1c. calculate the direction of the force
                double fx = force * (particles->particle_list[j]->x - particles->particle_list[i]->x) / distance;
                double fy = force * (particles->particle_list[j]->y - particles->particle_list[i]->y) / distance;
                double fz = force * (particles->particle_list[j]->z - particles->particle_list[i]->z) / distance;

                //1d. update the velocity of the particle
                particles->particle_list[i]->vx += fx / particles->particle_list[i]->m;
                particles->particle_list[i]->vy += fy / particles->particle_list[i]->m;
                //particles->particle_list[i]->vz += fz / particles->particle_list[i]->m;
            }
        }
    }
}

void Engine::update_locations(shared_ptr<Particles> particles, shared_ptr<backed_scaler> scaler) {
    //this function will update the locations of the particles

    //1. loop through the particles
    for (int i = 0; i < particles->particle_list.size(); i++) {
        
        particles->particle_list[i]->x += particles->particle_list[i]->vx * scaler->scaler[i];
        particles->particle_list[i]->y += particles->particle_list[i]->vy * scaler->scaler[i];
        //particles->particle_list[i]->z += particles->particle_list[i]->vz * scaler->scaler[i];
    }
}

