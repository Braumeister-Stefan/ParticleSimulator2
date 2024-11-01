//Standard libraries
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <chrono>
#include <iomanip>

//External libraries
#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"
#include <boost/multiprecision/cpp_dec_float.hpp>

//Internal libraries
#include "../include/PhysEngine.h"
#include "../include/InitStructs.h"
#include "../include/MathUtils.h"


//namespaces
using namespace std;
using namespace std::chrono;
using namespace boost::multiprecision;

//constants
//Gravitational constant
const double G = 6.674 * pow(10, -11); //m^3 kg^-1 s^-2

//define static member variables
double Engine::dt = 1; //time step

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

    
    

    //1. initialize the snapshots object
    shared_ptr<snapshots> particle_states = make_shared<snapshots>();


    //2.confirm that the run should start. If a dump file is available, ask the user if they want to use it, if yes save to snapshots object and return it.


    cout << scenario->name << "'s initial states loaded" << endl << endl;
    cout << "Press enter to start the simulation." << endl;
    cin.ignore();
    cin.get();




    //3. Initialize the time step globally and calculate the number of steps

    Engine::dt = scenario->dt;
    //cout << "Total time of the simulation: " << scenario->time << " s" << endl;
    //cout << "Time step length: " << scenario->dt << " s" << endl;

    double steps_db = scenario->time / scenario->dt;
    int steps = static_cast<int>(steps_db);

    cout << "Number of steps to be simulated: " << steps << endl << endl;

    //4. loop through the steps as defined in the scenario

    for (int i = 0; i < steps; i++) {


        //save start time
        auto start_time = high_resolution_clock::now();

        //1.add current state of particles to the snapshots object
        auto particles_copy = make_unique<Particles>(*particles);

        particle_states->snaps.push_back(move(particles_copy));
        //if instead we want to add it to the first element of the snapshots object, we can use the following line

        //2. Update the state of the particles

        //cout << "Updating particles..." << endl;
        double error_threshold = 0.00001;
        
        //print the current step
        //cout << "Step " << i << " of the simulation." << endl;

        double te_pre_update = calc_TE(particles);
        update_particles(particles);
        double te_post_update = calc_TE(particles);

        double te_error_update = (te_post_update - te_pre_update) / te_pre_update;

        if (te_error_update > error_threshold) {
            cout << "TE error in run(): " << te_error_update << endl;
        }

        //print every 5% of the simulation
        if (i % (steps / 20) == 0) {
            cout << i / (steps / 100) << "% of the simulation complete." << endl;
        }

        //save end time
        auto end_time = high_resolution_clock::now();

        //calculate the time taken to complete the step
        duration<double> time_taken = end_time - start_time;

        //save fps to the metrics object in the snapshots object
        particle_states->metrics.push_back(make_shared<test_metrics_t>());

        particle_states->metrics[i]->fps = 1 / time_taken.count();

    }

    //5.confirm that the run has ended

    cout << scenario->name << " simulation completed." << endl << endl;

    //6. save down snapshots to the csv 

    run_to_cache(scenario, particle_states);

    return particle_states;


}




void Engine::update_particles(shared_ptr<Particles> particles) {
    //this function will update the particles in the particles object

    double error_threshold = 0.00001;
    
    //1 resolve overlap
    //double te_pre_overlap = calc_TE(particles);


    //shared_ptr<backed_scaler> scaler;

    //iterate 8 times through the overlap resolution

    bool no_overlap = false;

    double te_pre_overlap = calc_TE(particles);

    for (int i = 0; i < 4; i++) {

        double te_pre_overlap_iter = calc_TE(particles);

        if (no_overlap != true) {
            no_overlap = resolve_overlap(particles);
        } else {
            break;
        }

        double te_post_overlap_iter = calc_TE(particles);

        double te_error_overlap_iter = (te_post_overlap_iter - te_pre_overlap_iter) / te_pre_overlap_iter;

        if (te_error_overlap_iter > error_threshold) {
            cout << "TE error overlap - for iteration in update_particles()" << i << " :" << te_error_overlap_iter << endl;
        }
        
    }
    
    //scaler = resolve_overlap(particles);


    double te_post_overlap = calc_TE(particles);

    double te_error_overlap = (te_post_overlap - te_pre_overlap) / te_pre_overlap;

    if (te_error_overlap > error_threshold) {
        cout << "TE error overlap - update func: " << te_error_overlap << endl;
    }
    
    //2 resolve collissions

    //cout << "Resolving collissions..." << endl;
    double te_pre_collission = calc_TE(particles);
    resolve_collisions(particles);
    double te_post_collission = calc_TE(particles);

    double te_error_collission = (te_post_collission - te_pre_collission) / te_pre_collission;

    if (te_error_collission > error_threshold) {
        cout << "TE error collission: " << te_error_collission << endl;
    }

    //2. resolve gravitational attraction

    //cout << "Resolving gravity..." << endl;
    //resolve_gravity_euler(particles);

    double te_pre_verlet = calc_TE(particles);
    resolve_gravity_verlet(particles); //use velocity verlet method
    double te_post_verlet = calc_TE(particles);
    double te_error_verlet = (te_post_verlet - te_pre_verlet) / te_pre_verlet;

    if (te_error_verlet > error_threshold) {
        cout << "TE error verlet: " << te_error_verlet << endl;
    }
    //cout << "TE error: " << te_error_verlet << endl;

    //3. update the locations of the particles (done inside velocity verlet method, but not in euler method)
    //update_locations(particles, scaler);
    
    

    
}

bool Engine::resolve_overlap(shared_ptr<Particles> particles) {
    bool no_overlap = true;
    high_prec tTE_pre = calc_TE(particles);  // Accumulates total energy error for all particle pairs

    high_prec error_threshold = 0.00001;  // Energy error threshold

    int pair_count = 0;  // To track the number of overlapping particle pairs

    // Loop through particle pairs
    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            if (check_collission(particles->particle_list[i], particles->particle_list[j])) {

                no_overlap = false;
                pair_count++;

<<<<<<< HEAD
                // Call the core function to resolve overlap between particles i and j
                resolve_overlap_ij(particles->particle_list[i], particles->particle_list[j]);
=======
                // Store pre-resolution positions and velocities using high_prec
                high_prec x1_before = particles->particle_list[i]->x;
                high_prec y1_before = particles->particle_list[i]->y;
                high_prec x2_before = particles->particle_list[j]->x;
                high_prec y2_before = particles->particle_list[j]->y;
                high_prec vx1_before = particles->particle_list[i]->vx;
                high_prec vy1_before = particles->particle_list[i]->vy;
                high_prec vx2_before = particles->particle_list[j]->vx;
                high_prec vy2_before = particles->particle_list[j]->vy;

                high_prec TE_pre_ij = calc_TE(particles);

                // Compute masses
                high_prec m1 = particles->particle_list[i]->m;
                high_prec m2 = particles->particle_list[j]->m;

                // Compute overlap
                high_prec dx = x2_before - x1_before;
                high_prec dy = y2_before - y1_before;
                high_prec distance = hypot(dx, dy);

                // Avoid division by zero
                if (distance == 0.0) {
                    dx = 0.001 * (rand() / (double)RAND_MAX - 0.5);
                    dy = 0.001 * (rand() / (double)RAND_MAX - 0.5);
                    distance = hypot(dx, dy);
                }

                high_prec overlap = (particles->particle_list[i]->rad + particles->particle_list[j]->rad) - distance;

                // Compute collision normal
                high_prec nx = dx / distance;
                high_prec ny = dy / distance;

                // Compute the displacement for each particle
                high_prec total_mass = m1 + m2;
                high_prec d1 = (overlap * m2) / total_mass;
                high_prec d2 = (overlap * m1) / total_mass;

                // Move particles proportionally to their masses to resolve overlap
                particles->particle_list[i]->x = (x1_before - nx * d1).convert_to<double>();
                particles->particle_list[i]->y = (y1_before - ny * d1).convert_to<double>();
                particles->particle_list[j]->x = (x2_before + nx * d2).convert_to<double>();
                particles->particle_list[j]->y = (y2_before + ny * d2).convert_to<double>();

                // **Adjust velocities to conserve momentum and energy**

                // Compute relative velocity along the normal
                high_prec rel_vel = (vx2_before - vx1_before) * nx + (vy2_before - vy1_before) * ny;

                // Compute impulse scalar
                high_prec e = (particles->particle_list[i]->rest + particles->particle_list[j]->rest) / 2.0; // Average restitution
                high_prec impulse = (-(1 + e) * rel_vel) / (1 / m1 + 1 / m2);

                // Apply impulse to particles
                high_prec vx1_new = vx1_before - (impulse / m1) * nx;
                high_prec vy1_new = vy1_before - (impulse / m1) * ny;
                high_prec vx2_new = vx2_before + (impulse / m2) * nx;
                high_prec vy2_new = vy2_before + (impulse / m2) * ny;
          


                // Update particle velocities
                particles->particle_list[i]->vx = vx1_new.convert_to<double>();
                particles->particle_list[i]->vy = vy1_new.convert_to<double>();
                particles->particle_list[j]->vx = vx2_new.convert_to<double>();
                particles->particle_list[j]->vy = vy2_new.convert_to<double>();

                // **Compute energies after**

                high_prec TE_post_ij = calc_TE(particles);

                if (TE_post_ij != TE_pre_ij) {
                    // Momentum should be conserved, so Pn_new = Pn
                    // Energy should be adjusted to TE_pre_ij

                    // Compute the energy difference
                    high_prec delta_E = TE_pre_ij - TE_post_ij;

                    // Compute velocities along the collision normal
                    high_prec v1n = vx1_new * nx + vy1_new * ny;
                    high_prec v2n = vx2_new * nx + vy2_new * ny;

                    // Compute the total momentum along the normal
                    high_prec Pn = m1 * v1n + m2 * v2n;

                    // delta_v2n = -(m1 / m2) * delta_v1n (from conservation of momentum)

                    // Set up the quadratic equation:
                    // m1 * (v1n - v2n) * delta_v1n + 0.5 * delta_v1n^2 * m1 * (1 + m1 / m2) = delta_E

                    // Coefficients for the quadratic equation: a * delta_v1n^2 + b * delta_v1n + c = 0
                    high_prec a = 0.5 * m1 * (1 + m1 / m2);
                    high_prec b = m1 * (v1n - v2n);
                    high_prec c = -delta_E;

                    // Solve the quadratic equation for delta_v1n
                    high_prec discriminant = b * b - 4 * a * c;

                    if (discriminant >= 0) {
                        high_prec sqrt_discriminant = sqrt(discriminant);
                        high_prec delta_v1n1 = (-b + sqrt_discriminant) / (2 * a);
                        high_prec delta_v1n2 = (-b - sqrt_discriminant) / (2 * a);

                        // Choose the solution that results in smaller adjustment
                        high_prec delta_v1n = abs(delta_v1n1) < abs(delta_v1n2) ? delta_v1n1 : delta_v1n2;

                        // Compute delta_v2n
                        high_prec delta_v2n = -(m1 / m2) * delta_v1n;

                        // Adjust the normal components of velocities
                        high_prec v1n_new = v1n + delta_v1n;
                        high_prec v2n_new = v2n + delta_v2n;

                        // Compute the tangential components (unchanged)
                        high_prec v1t = vx1_new * (-ny) + vy1_new * nx;
                        high_prec v2t = vx2_new * (-ny) + vy2_new * nx;

                        // Reconstruct the velocities
                        particles->particle_list[i]->vx = (v1n_new * nx - v1t * ny).convert_to<double>();
                        particles->particle_list[i]->vy = (v1n_new * ny + v1t * nx).convert_to<double>();
                        particles->particle_list[j]->vx = (v2n_new * nx - v2t * ny).convert_to<double>();
                        particles->particle_list[j]->vy = (v2n_new * ny + v2t * nx).convert_to<double>();

                        // Calculate the new total energy

                        
                        high_prec TE_post_ij_corrected = calc_TE(particles);

                        high_prec TE_error_ij_corrected = (TE_post_ij_corrected - TE_pre_ij) / TE_pre_ij;

                        if (abs(TE_error_ij_corrected) > error_threshold) {
                            cout << "TE error after energy correction: " << TE_error_ij_corrected << endl;
                        }

                    } else {
                        cout << "Could not apply energy correction!" << endl;
                        cout << "TE error before energy correction: " << (TE_post_ij - TE_pre_ij) / TE_pre_ij << endl;
                    }
                }


>>>>>>> 1c9555868a5bb893f0b409901a228c4c1d973fd9
            }
        }
    }

    // Calculate the final total energy error and print it
    high_prec tTE_post = calc_TE(particles);

    high_prec tTE_error = (tTE_post - tTE_pre) / tTE_pre;

    if (abs(tTE_error) > error_threshold) {
        cout << "Total TE error after overlap resolution - in resolve_overlap(): " << tTE_error << endl;
    }

    return no_overlap;
}

void Engine::resolve_overlap_ij(shared_ptr<Particle> particle_i,shared_ptr<Particle> particle_j) {

    //define error threshold
    high_prec error_threshold = 0.00000001;

    // Store pre-resolution positions and velocities using high_prec
    high_prec x1_before = particle_i->x;
    high_prec y1_before = particle_i->y;
    high_prec x2_before = particle_j->x;
    high_prec y2_before = particle_j->y;
    high_prec vx1_before = particle_i->vx;
    high_prec vy1_before = particle_i->vy;
    high_prec vx2_before = particle_j->vx;
    high_prec vy2_before = particle_j->vy;

    high_prec TE_pre_ij = calc_TE_ij(particle_i, particle_j);

    // Compute masses
    high_prec m1 = particle_i->m;
    high_prec m2 = particle_j->m;

    // Compute overlap
    high_prec dx = x2_before - x1_before;
    high_prec dy = y2_before - y1_before;
    high_prec distance = hypot(dx, dy);

    // Avoid division by zero
    if (distance == 0.0) {
        dx = 0.001 * (rand() / (double)RAND_MAX - 0.5);
        dy = 0.001 * (rand() / (double)RAND_MAX - 0.5);
        distance = hypot(dx, dy);
    }

    high_prec overlap = (particle_i->rad + particle_j->rad) - distance;

    // Compute collision normal
    high_prec nx = dx / distance;
    high_prec ny = dy / distance;

    // Compute the displacement for each particle
    high_prec total_mass = m1 + m2;
    high_prec d1 = (overlap * m2) / total_mass;
    high_prec d2 = (overlap * m1) / total_mass;

    // Move particles proportionally to their masses to resolve overlap
    particle_i->x = (x1_before - nx * d1).convert_to<double>();
    particle_i->y = (y1_before - ny * d1).convert_to<double>();
    particle_j->x = (x2_before + nx * d2).convert_to<double>();
    particle_j->y = (y2_before + ny * d2).convert_to<double>();

    // **Adjust velocities to conserve momentum and energy**

    // Compute relative velocity along the normal
    high_prec rel_vel = (vx2_before - vx1_before) * nx + (vy2_before - vy1_before) * ny;

    // Compute impulse scalar
    high_prec e = (particle_i->rest + particle_j->rest) / 2.0; // Average restitution
    high_prec impulse = (-(1 + e) * rel_vel) / (1 / m1 + 1 / m2);

    // Apply impulse to particles
    high_prec vx1_new = vx1_before - (impulse / m1) * nx;
    high_prec vy1_new = vy1_before - (impulse / m1) * ny;
    high_prec vx2_new = vx2_before + (impulse / m2) * nx;
    high_prec vy2_new = vy2_before + (impulse / m2) * ny;

    // Update particle velocities
    particle_i->vx = vx1_new.convert_to<double>();
    particle_i->vy = vy1_new.convert_to<double>();
    particle_j->vx = vx2_new.convert_to<double>();
    particle_j->vy = vy2_new.convert_to<double>();

    // **Compute energies after**
    high_prec TE_post_ij = calc_TE_ij(particle_i, particle_j);

    cout << "TE error after overlap resolution, before correction: " << (TE_post_ij - TE_pre_ij) / TE_pre_ij << endl;

    if (TE_post_ij != TE_pre_ij) {
        // Momentum should be conserved, so Pn_new = Pn
        // Energy should be adjusted to TE_pre_ij

        // Compute the energy difference
        high_prec delta_E = TE_pre_ij - TE_post_ij;

        // Compute velocities along the collision normal
        high_prec v1n = vx1_new * nx + vy1_new * ny;
        high_prec v2n = vx2_new * nx + vy2_new * ny;

        // Compute the total momentum along the normal
        high_prec Pn = m1 * v1n + m2 * v2n;

        // delta_v2n = -(m1 / m2) * delta_v1n (from conservation of momentum)

        // Set up the quadratic equation:
        // m1 * (v1n - v2n) * delta_v1n + 0.5 * delta_v1n^2 * m1 * (1 + m1 / m2) = delta_E

        // Coefficients for the quadratic equation: a * delta_v1n^2 + b * delta_v1n + c = 0
        high_prec a = 0.5 * m1 * (1 + m1 / m2);
        high_prec b = m1 * (v1n - v2n);
        high_prec c = -delta_E;

        // Solve the quadratic equation for delta_v1n
        high_prec discriminant = b * b - 4 * a * c;

        if (discriminant >= 0) {
            high_prec sqrt_discriminant = sqrt(discriminant);
            high_prec delta_v1n1 = (-b + sqrt_discriminant) / (2 * a);
            high_prec delta_v1n2 = (-b - sqrt_discriminant) / (2 * a);

            // Choose the solution that results in smaller adjustment
            high_prec delta_v1n = abs(delta_v1n1) < abs(delta_v1n2) ? delta_v1n1 : delta_v1n2;

            // Compute delta_v2n
            high_prec delta_v2n = -(m1 / m2) * delta_v1n;

            // Adjust the normal components of velocities
            high_prec v1n_new = v1n + delta_v1n;
            high_prec v2n_new = v2n + delta_v2n;

            // Compute the tangential components (unchanged)
            high_prec v1t = vx1_new * (-ny) + vy1_new * nx;
            high_prec v2t = vx2_new * (-ny) + vy2_new * nx;

            // Reconstruct the velocities
            particle_i->vx = (v1n_new * nx - v1t * ny).convert_to<double>();
            particle_i->vy = (v1n_new * ny + v1t * nx).convert_to<double>();
            particle_j->vx = (v2n_new * nx - v2t * ny).convert_to<double>();
            particle_j->vy = (v2n_new * ny + v2t * nx).convert_to<double>();

            // Calculate the new total energy
            high_prec TE_post_ij_corrected = calc_TE_ij(particle_i, particle_j);

            high_prec TE_error_ij_corrected = (TE_post_ij_corrected - TE_pre_ij) / TE_pre_ij;

            //if (abs(TE_error_ij_corrected) > error_threshold) {
            cout << "TE error after energy correction: " << TE_error_ij_corrected << endl;
            //}

        } else {
            cout << "Could not apply energy correction!" << endl;
            cout << "TE error before energy correction: " << (TE_post_ij - TE_pre_ij) / TE_pre_ij << endl;
        }
    }
}




void Engine::resolve_collisions(shared_ptr<Particles> particles) {
    //this function will resolve collissions between particles

  

    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            //2a. check for collission between particle i and particle j

            bool collission = check_collission(particles->particle_list[i], particles->particle_list[j]);
            //if collission detected, backtrack the particles to the point of collission and store the amount of time of the timestep that was left
            
            

            if (collission) {
                //resolve the collission

                resolve_collission(particles->particle_list[i], particles->particle_list[j]);

                
            }
            

        }
    }

    //cout << "Particles after collission resolution:" << endl;
    //for (int i = 0; i < particles->particle_list.size(); i++) {
    //    cout << "Particle " << i << ": x=" << particles->particle_list[i]->x << ", y=" << particles->particle_list[i]->y << ", z=" << particles->particle_list[i]->z << endl;

    //}

    
    
}


bool Engine::check_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    //this function will check if two particles are colliding

    //1. calculate the distance between the particles
    double distance = hypot(particle1->x - particle2->x, particle1->y - particle2->y);

    //2. check if the distance is smaller than the sum of the radii of the particles
    if (distance < particle1->rad + particle2->rad) {
        return true;
    }
    else {
        return false;
    }
}


using high_prec = cpp_dec_float_50;



void Engine::resolve_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    // This function will resolve the collision between two particles

    // Calculate the normal vector
    Vector2D normal = { high_prec(particle2->x) - high_prec(particle1->x), high_prec(particle2->y) - high_prec(particle1->y) };
    high_prec distance = sqrt(normal.x * normal.x + normal.y * normal.y);
    normal.x /= distance;
    normal.y /= distance;

    // Calculate the tangential vector
    Vector2D tangent = { -normal.y, normal.x };

    // Decompose velocities into normal and tangential components
    high_prec v1n = normal.x * high_prec(particle1->vx) + normal.y * high_prec(particle1->vy);
    high_prec v1t = tangent.x * high_prec(particle1->vx) + tangent.y * high_prec(particle1->vy);
    high_prec v2n = normal.x * high_prec(particle2->vx) + normal.y * high_prec(particle2->vy);
    high_prec v2t = tangent.x * high_prec(particle2->vx) + tangent.y * high_prec(particle2->vy);

    // Calculate new normal velocities using conservation of momentum and kinetic energy
    high_prec v1n_new = (v1n * (high_prec(particle1->m) - high_prec(particle2->m)) + (1 + high_prec(particle1->rest)) * high_prec(particle2->m) * v2n) / (high_prec(particle1->m) + high_prec(particle2->m));
    high_prec v2n_new = (v2n * (high_prec(particle2->m) - high_prec(particle1->m)) + (1 + high_prec(particle2->rest)) * high_prec(particle1->m) * v1n) / (high_prec(particle1->m) + high_prec(particle2->m));

    // The tangential components remain unchanged
    high_prec v1t_new = v1t;
    high_prec v2t_new = v2t;

    // Convert the scalar normal and tangential velocities into vectors
    Vector2D v1n_vec = { v1n_new * normal.x, v1n_new * normal.y };
    Vector2D v1t_vec = { v1t_new * tangent.x, v1t_new * tangent.y };
    Vector2D v2n_vec = { v2n_new * normal.x, v2n_new * normal.y };
    Vector2D v2t_vec = { v2t_new * tangent.x, v2t_new * tangent.y };

    // Update the velocities of the particles
    particle1->vx = (v1n_vec.x + v1t_vec.x).convert_to<double>();
    particle1->vy = (v1n_vec.y + v1t_vec.y).convert_to<double>();
    particle2->vx = (v2n_vec.x + v2t_vec.x).convert_to<double>();
    particle2->vy = (v2n_vec.y + v2t_vec.y).convert_to<double>();
}


void Engine::resolve_gravity_euler(shared_ptr<Particles> particles) {//abandoned technique
    // This function resolves the gravitational attraction between particles using Euler's method. Not sufficiently accurate for energy conservation.

    // Define a small value to avoid division by zero

    high_prec epsilon = 0.0000001;

    // Loop through the particles
    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            // Skip self-interaction
            if (i == j) continue;

            // Calculate the distance between particles i and j
            high_prec dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            high_prec dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            high_prec distance = hypot(dx, dy);

            if (distance < epsilon) {
                
                distance = epsilon;
            }

            // Calculate the gravitational force magnitude
            high_prec force = G * particles->particle_list[i]->m * particles->particle_list[j]->m / (distance * distance);

            // Calculate the components of the force
            high_prec fx = force * (dx / distance);
            high_prec fy = force * (dy / distance);

            // Calculate KE, PE, and TE for the system of i and j before the force update, all are magnitudes
            high_prec KE_pre = 0.5 * particles->particle_list[i]->m * (pow(particles->particle_list[i]->vx, 2) + pow(particles->particle_list[i]->vy, 2)) +
                               0.5 * particles->particle_list[j]->m * (pow(particles->particle_list[j]->vx, 2) + pow(particles->particle_list[j]->vy, 2));
            high_prec PE_pre = -G * particles->particle_list[i]->m * particles->particle_list[j]->m / distance;
            high_prec TE_pre = KE_pre + PE_pre;



            // Update the velocities of both particles to conserve momentum
            particles->particle_list[i]->vx += (fx / particles->particle_list[i]->m).convert_to<double>() * dt;
            particles->particle_list[i]->vy += (fy / particles->particle_list[i]->m).convert_to<double>() * dt;

            particles->particle_list[j]->vx -= (fx / particles->particle_list[j]->m).convert_to<double>() * dt;
            particles->particle_list[j]->vy -= (fy / particles->particle_list[j]->m).convert_to<double>() * dt;

            // Recalculate distance based on new positions
            high_prec new_dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            high_prec new_dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            high_prec new_distance = hypot(new_dx, new_dy);

            
            // Calculate KE, PE, and TE for the system of i and j after the force update using new distance
            high_prec KE_post = 0.5 * particles->particle_list[i]->m * (pow(particles->particle_list[i]->vx, 2) + pow(particles->particle_list[i]->vy, 2)) +
                                0.5 * particles->particle_list[j]->m * (pow(particles->particle_list[j]->vx, 2) + pow(particles->particle_list[j]->vy, 2));

            high_prec PE_post = -G * particles->particle_list[i]->m * particles->particle_list[j]->m / new_distance;

            high_prec TE_post = KE_post + PE_post;

            
            // Calculate the difference in total energy before and after the force update
            high_prec TE_diff = TE_post - TE_pre;
            high_prec TE_error = TE_diff / TE_pre;

            // Threshold for the error
            high_prec TE_error_threshold = 0.1; // Lower threshold for higher precision

            // If the error is greater than the threshold, print the error and distance
            if (abs(TE_error) > TE_error_threshold) {
                cout << "Distance: " << distance << endl;
                cout << "TE error %: " << TE_error * 100 << "%" << endl;
            }
        }
    }
}

void Engine::resolve_gravity_verlet(shared_ptr<Particles> particles) {
    // This function resolves the gravitational attraction between particles using the full Velocity Verlet integration scheme
    double epsilon = 0.001; // To avoid division by zero
    int n = particles->particle_list.size();

    double TE_pre = calc_TE(particles); //calc total energy before the update


    // Step 1: Initialize net forces for all particles
    vector<Vector2D_double> net_force(n, {0.0, 0.0});


    
    // Step 2: Calculate and accumulate gravitational forces for all pairs
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // Calculate distance components between particles i and j
            double dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            double dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            double distance = hypot(dx, dy);
            //if (distance < epsilon) distance = epsilon;
            distance = sqrt(distance*distance + epsilon*epsilon);

            // Calculate gravitational force magnitude
            double force = G * particles->particle_list[i]->m * particles->particle_list[j]->m / (distance * distance);

            // Force components
            double fx = force * (dx / distance);
            double fy = force * (dy / distance);

            // Accumulate forces for particles i and j (Newton's third law)
            net_force[i].x += fx;
            net_force[i].y += fy;
            net_force[j].x -= fx;
            net_force[j].y -= fy;
        }
    }

    // Step 3: Update velocities using the accumulated forces (first half-step)
    for (int i = 0; i < n; i++) {
        double dt_i = dt * 1;
        particles->particle_list[i]->vx += 0.5 * (net_force[i].x / particles->particle_list[i]->m) * dt_i;
        particles->particle_list[i]->vy += 0.5 * (net_force[i].y / particles->particle_list[i]->m) * dt_i;
    }

    // Step 4: Update positions using the updated velocities
    for (int i = 0; i < n; i++) {
        double dt_i = dt * 1;
        particles->particle_list[i]->x += particles->particle_list[i]->vx * dt_i;
        particles->particle_list[i]->y += particles->particle_list[i]->vy * dt_i;
    }


    // Step 5: Recalculate gravitational forces for the updated positions (second half-step)
    fill(net_force.begin(), net_force.end(), Vector2D_double{0.0, 0.0});
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // Recalculate distance components between particles i and j
            double dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            double dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            double distance = hypot(dx, dy);
            if (distance < epsilon) distance = epsilon;

            // Recalculate gravitational force magnitude
            double force = G * particles->particle_list[i]->m * particles->particle_list[j]->m / (distance * distance);

            // Force components
            double fx = force * (dx / distance);
            double fy = force * (dy / distance);

            // Accumulate forces for particles i and j
            net_force[i].x += fx;
            net_force[i].y += fy;
            net_force[j].x -= fx;
            net_force[j].y -= fy;
        }
    }

    // Step 6: Update velocities again using the new forces (second half-step)
    for (int i = 0; i < n; i++) {
        double dt_i = dt * 1;
        particles->particle_list[i]->vx += 0.5 * (net_force[i].x / particles->particle_list[i]->m) * dt_i;
        particles->particle_list[i]->vy += 0.5 * (net_force[i].y / particles->particle_list[i]->m) * dt_i;
    }

    // Step 7: validate conservation of energy
    double TE_post = calc_TE(particles); //calc total energy after the update

    // Calculate the difference in total energy
    double TE_diff = TE_post - TE_pre;
    double TE_error = TE_diff / TE_pre;
    double TE_error_threshold = 0.001; // Lower threshold for higher precision

    // Inform the user if the energy error is above the threshold
    if (abs(TE_error) > TE_error_threshold) {
        cout << "Warning: Energy loss: " << TE_error * 100 << "%" << endl;


    }
}



void Engine::update_locations(shared_ptr<Particles> particles, shared_ptr<backed_scaler> scaler) {
    //this function will update the locations of the particles

    //1. loop through the particles
    for (int i = 0; i < particles->particle_list.size(); i++) {
        
        particles->particle_list[i]->x += particles->particle_list[i]->vx * scaler->scaler[i] * dt;
        particles->particle_list[i]->y += particles->particle_list[i]->vy * scaler->scaler[i] * dt;
        //particles->particle_list[i]->z += particles->particle_list[i]->vz * scaler->scaler[i] * dt;
    }
}

void Engine::run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {
    //this function will save the snapshots to a csv file

    //1. define the file name, located on Inputs\rendered_scenarios

    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    //2. open the file

    ofstream file(file_name);

    //check if the file is open

    if (!file.is_open()) {
        cout << "Failed to write snapshots to cache!" << endl;
        return;
    }

    //3. write headers

    file << "step_id, particle_id,r,g,b,x,y,z,vx,vy,vz,m,rad,rest" << endl;

    //4. loop through the snapshots and write the data to the file, adding the step_id

    for (int i = 0; i < particle_states->snaps.size(); i++) {
    for (int j = 0; j < particle_states->snaps[i]->particle_list.size(); j++) {
        file << std::fixed << std::setprecision(15) // Set precision to 15 decimal places
             << i << ","
             << particle_states->snaps[i]->particle_list[j]->particle_id << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->r) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->g) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->b) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->x) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->y) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->z) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->vx) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->vy) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->vz) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->m) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->rad) << ","
             << static_cast<double>(particle_states->snaps[i]->particle_list[j]->rest) << std::endl;
    }
}

    //5. save and close the file

    cout << "Saving simulation..." << endl;

    file.close();

    cout << "Snapshots saved to " << file_name << endl << endl;
    

}

bool Engine::cache_exists(shared_ptr<scenario> scenario) {
    //this function will check if a cache file exists for the scenario

    //1. define the file name, located on Inputs\rendered_scenarios

    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    //2. open the file

    ifstream file(file_name);

    //check if the file is open

    if (!file.is_open()) {
        return false;
    }

    cout << "Cache file found for " << scenario->name << endl;

    return true;
}

shared_ptr<snapshots> Engine::run_from_cache(shared_ptr<scenario> scenario) {

    //this function will read the snapshots from a csv file

    //1. define the file name, located on Inputs\rendered_scenarios

    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    //2. open the file

    ifstream file(file_name);

    //3. read the file, using the first column to determine the step_id, store every step in a new snapshot

    shared_ptr<snapshots> particle_states = make_shared<snapshots>();

    // Create a CSV reader object
    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;
    const int column_count = 14; // Adjust the column count based on your CSV structure
    io::CSVReader<column_count, TrimPolicy, QuotePolicy> in(file_name);

    // 4. Read the contents of the CSV file and store each row in a particle object. Store all particle objects in a vector.
    string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14;
    int step_id = 0;
    shared_ptr<Particles> particles = make_shared<Particles>();

    // Discard the header row. This list is a guide to the columns in the CSV file.
    in.read_header(io::ignore_extra_column, "step_id", "particle_id", "r", "g", "b", "x", "y", "z", "vx", "vy", "vz", "m", "rad", "rest");


    while (in.read_row(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14)) {
        // Check if the step_id has changed
        if (stoi(col1) != step_id) {
            // Store the current snapshot
            particle_states->snaps.push_back(particles);
            // Create a new snapshot
            particles = make_shared<Particles>();
            step_id = stoi(col1);
        }

        // Create a new particle
        shared_ptr<Particle> particle = make_shared<Particle>();
        particle->particle_id = stoi(col2);
        particle->r = stod(col3);
        particle->g = stod(col4);
        particle->b = stod(col5);
        particle->x = stod(col6);
        particle->y = stod(col7);
        particle->z = stod(col8);
        particle->vx = stod(col9);
        particle->vy = stod(col10);
        particle->vz = stod(col11);
        particle->m = stod(col12);
        particle->rad = stod(col13);
        particle->rest = stod(col14);

        // Add the particle to the snapshot
        particles->particle_list.push_back(particle);

    }

    // Store the last snapshot
    particle_states->snaps.push_back(particles);

    cout << "Snapshots successfully loaded from cache: " << file_name << std::endl;

    return particle_states;
}


// Function to calculate the total kinetic energy (KE) and potential energy (PE) for a system of particles
double Engine::calc_TE(shared_ptr<Particles> particles) {
    int n = particles->particle_list.size();
    double KE = 0.0;
    double PE = 0.0;
    double epsilon = 0.001; // To avoid division by zero

    // Calculate kinetic energy
    for (int i = 0; i < n; i++) {
        KE += 0.5 * particles->particle_list[i]->m * 
              (pow(particles->particle_list[i]->vx, 2) + pow(particles->particle_list[i]->vy, 2));
    }

    // Calculate potential energy
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            double dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            double distance = hypot(dx, dy);
            distance = sqrt(distance * distance + epsilon * epsilon);
            PE += -G * particles->particle_list[i]->m * particles->particle_list[j]->m / distance;
        }
    }

    // Return the total energy (KE + PE)
    return KE + PE;
}

double Engine::calc_TE_ij(shared_ptr<Particle> p1, shared_ptr<Particle> p2) {
    double KE = 0.0;
    double PE = 0.0;
    double epsilon = 0.001; // To avoid division by zero

    // Calculate kinetic energy for both particles
    KE += 0.5 * p1->m * (pow(p1->vx, 2) + pow(p1->vy, 2));
    KE += 0.5 * p2->m * (pow(p2->vx, 2) + pow(p2->vy, 2));

    // Calculate potential energy between the two particles
    double dx = p2->x - p1->x;
    double dy = p2->y - p1->y;
    double distance = hypot(dx, dy);
    distance = sqrt(distance * distance + epsilon * epsilon);
    PE += -G * p1->m * p2->m / distance;

    // Return the total energy (KE + PE)
    return KE + PE;
}

double Engine::calc_mom(shared_ptr<Particles> particles) {
    int n = particles->particle_list.size();
    double total_momentum = 0.0;

    // Calculate momentum for each particle and sum up
    for (int i = 0; i < n; i++) {
        double m = particles->particle_list[i]->m;
        double vx = particles->particle_list[i]->vx;
        double vy = particles->particle_list[i]->vy;
        double momentum_magnitude = m * hypot(vx, vy);
        total_momentum += momentum_magnitude;
    }

    return total_momentum;
}

double Engine::calc_mom_ij(shared_ptr<Particle> p1, shared_ptr<Particle> p2) {
    double total_momentum = 0.0;

    // Calculate momentum magnitude for the first particle
    double momentum1 = p1->m * hypot(p1->vx, p1->vy);
    total_momentum += momentum1;

    // Calculate momentum magnitude for the second particle
    double momentum2 = p2->m * hypot(p2->vx, p2->vy);
    total_momentum += momentum2;

    return total_momentum;
}