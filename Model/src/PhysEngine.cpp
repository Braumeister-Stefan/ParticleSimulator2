#include "../include/PhysEngine.h"
#include "../include/InitStructs.h"
#include "../include/MathUtils.h"

#include <memory>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>



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

    cout << "Resolving gravity..." << endl;
    resolve_gravity(particles);

    

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
                 double scaler_i = backtrack_pair(particles->particle_list[i], particles->particle_list[j]);

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


double Engine::backtrack_pair(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    // Define the dot product function locally within this function
    

    // This function will backtrack two particles to the point of collision
    // and return the amount of time of the timestep that was left.

    // 1. Store particle positions and velocities in 2D
    Vector2D pos1 = { particle1->x, particle1->y };
    Vector2D pos2 = { particle2->x, particle2->y };
    Vector2D vel1 = { particle1->vx, particle1->vy };
    Vector2D vel2 = { particle2->vx, particle2->vy };

    // 2. Calculate relative velocity and position difference
    Vector2D rel_pos = pos2 - pos1;
    Vector2D rel_vel = vel2 - vel1;

    // 3. Calculate quadratic terms for collision detection
    double a = dot(rel_vel, rel_vel);  // Coefficient of t^2
    double b = 2 * dot(rel_pos, rel_vel);  // Coefficient of t
    double c = dot(rel_pos, rel_pos) - pow(particle1->rad + particle2->rad, 2);  // Constant term

    // 4. Solve quadratic equation for time of collision
    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        // No collision occurs
        return 1.0;  // Full timestep remains
    }

    double sqrt_discriminant = sqrt(discriminant);
    double t_collision1 = (-b - sqrt_discriminant) / (2 * a);
    double t_collision2 = (-b + sqrt_discriminant) / (2 * a);

    // Use the root that falls within the current timestep 
    double t_collision = 0;
    if (t_collision1 >= -1 && t_collision1 <= 0) {
        t_collision = t_collision1;
    } else if (t_collision2 >= -1 && t_collision2 <= 0) {
        t_collision = t_collision2;
    }



    // 5. Backtrack positions to collision point
    pos1 = pos1 + t_collision * vel1;
    pos2 = pos2 + t_collision * vel2;

    // 6. Calculate the time scaler
    double time_scaler = 1.0 - t_collision;  // Remaining time in the timestep

    return time_scaler;
}





void Engine::resolve_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    // This function will resolve the collision between two particles

    // Calculate the normal vector
    Vector2D normal = { particle2->x - particle1->x, particle2->y - particle1->y };
    double distance = sqrt(normal.x * normal.x + normal.y * normal.y);
    normal.x /= distance;
    normal.y /= distance;

    // Calculate the tangential vector
    Vector2D tangent = { -normal.y, normal.x };

    // Decompose velocities into normal and tangential components
    double v1n = normal.x * particle1->vx + normal.y * particle1->vy;
    double v1t = tangent.x * particle1->vx + tangent.y * particle1->vy;
    double v2n = normal.x * particle2->vx + normal.y * particle2->vy;
    double v2t = tangent.x * particle2->vx + tangent.y * particle2->vy;

    // Calculate new normal velocities using conservation of momentum and kinetic energy
    double v1n_new = (v1n * (particle1->m - particle2->m) + (1 + particle1->rest) * particle2->m * v2n) / (particle1->m + particle2->m);
    double v2n_new = (v2n * (particle2->m - particle1->m) + (1 + particle2->rest) * particle1->m * v1n) / (particle1->m + particle2->m);

    // The tangential components remain unchanged
    double v1t_new = v1t;
    double v2t_new = v2t;

    // Convert the scalar normal and tangential velocities into vectors
    Vector2D v1n_vec = { v1n_new * normal.x, v1n_new * normal.y };
    Vector2D v1t_vec = { v1t_new * tangent.x, v1t_new * tangent.y };
    Vector2D v2n_vec = { v2n_new * normal.x, v2n_new * normal.y };
    Vector2D v2t_vec = { v2t_new * tangent.x, v2t_new * tangent.y };

    // Update the velocities of the particles
    particle1->vx = v1n_vec.x + v1t_vec.x;
    particle1->vy = v1n_vec.y + v1t_vec.y;
    particle2->vx = v2n_vec.x + v2t_vec.x;
    particle2->vy = v2n_vec.y + v2t_vec.y;
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

