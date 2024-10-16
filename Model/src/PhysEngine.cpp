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

        //if i=714, print a message to the console
        //if (i == 714) {

        //    cout << "Step 714 reached." << endl;
        //}

        //save start time
        auto start_time = high_resolution_clock::now();

        //1.add current state of particles to the snapshots object
        auto particles_copy = make_unique<Particles>(*particles);

        particle_states->snaps.push_back(move(particles_copy));
        //if instead we want to add it to the first element of the snapshots object, we can use the following line

        //2. Update the state of the particles

        //cout << "Updating particles..." << endl;
        
        update_particles(particles);

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

    
    

    
    //1 resolve collissions

    //cout << "Resolving collissions..." << endl;
    shared_ptr<backed_scaler> scaler = resolve_collisions(particles);

    //2. resolve gravitational attraction

    //cout << "Resolving gravity..." << endl;
    resolve_gravity(particles);

    //make scaler with all scalers set to 1
    //  shared_ptr<backed_scaler> scaler = make_shared<backed_scaler>();
    //  scaler->scaler.reserve(particles->particle_list.size());
    //  for (int i = 0; i < particles->particle_list.size(); i++) {
    //      scaler->scaler.push_back(1);
    //  }
    

    //3. update locations with velocities

    //cout << "Updating locations..." << endl;
    update_locations(particles, scaler);


}

shared_ptr<backed_scaler> Engine::resolve_collisions(shared_ptr<Particles> particles) {
    //this function will resolve collissions between particles

    //print all particles (1 line per particle) before the collission resolution
    //cout << "Particles before collission resolution:" << endl;
    //for (int i = 0; i < particles->particle_list.size(); i++) {
    //    cout << "Particle " << i << ": x=" << particles->particle_list[i]->x << ", y=" << particles->particle_list[i]->y << ", z=" << particles->particle_list[i]->z << endl;

    //}

    //1.initialize the scaler object
    shared_ptr<backed_scaler> scaler = make_shared<backed_scaler>(); 
    scaler->scaler.reserve(particles->particle_list.size()); //reserve the amount of scalers needed
    for (int i = 0; i < particles->particle_list.size(); i++) {
        scaler->scaler.push_back(1);
    }

    //cout <<"Scalers initialized." << endl;

    //2.check for collissions between particles and resolve them


    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            //2a. check for collission between particle i and particle j

            bool collission = check_collission(particles->particle_list[i], particles->particle_list[j]);
            //if collission detected, backtrack the particles to the point of collission and store the amount of time of the timestep that was left
            
            

            if (collission) {
                //2b. backtrack the particles

                //cout << "Collission detected between particles " << i << " and " << j << endl;
                
                //cout << "pre-backtrack positions: " << endl;
                //cout << "Particle " << i << ": x=" << particles->particle_list[i]->x << ", y=" << particles->particle_list[i]->y << ", z=" << particles->particle_list[i]->z << endl;
                //cout << "Particle " << j << ": x=" << particles->particle_list[j]->x << ", y=" << particles->particle_list[j]->y << ", z=" << particles->particle_list[j]->z << endl;


                 double scaler_i = backtrack_pair(particles->particle_list[i], particles->particle_list[j]);

                //cout << "post-backtrack positions: " << endl;
                //cout << "Particle " << i << ": x=" << particles->particle_list[i]->x << ", y=" << particles->particle_list[i]->y << ", z=" << particles->particle_list[i]->z << endl;
                //cout << "Particle " << j << ": x=" << particles->particle_list[j]->x << ", y=" << particles->particle_list[j]->y << ", z=" << particles->particle_list[j]->z << endl;

                //2c. store the scaler for the particle

                scaler->scaler[i] = scaler_i;

                //2d. resolve the collission

                resolve_collission(particles->particle_list[i], particles->particle_list[j]);

                
            }
            

        }
    }

    //cout << "Particles after collission resolution:" << endl;
    //for (int i = 0; i < particles->particle_list.size(); i++) {
    //    cout << "Particle " << i << ": x=" << particles->particle_list[i]->x << ", y=" << particles->particle_list[i]->y << ", z=" << particles->particle_list[i]->z << endl;

    //}

    
    return scaler;
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

double Engine::backtrack_pair(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    // Define the dot product function locally within this function
    

    // This function will backtrack two particles to the point of collision
    // and return the amount of time of the timestep that was left.

    // 1. Store particle positions and velocities in 2D
    // Vector2D pos1 = { particle1->x, particle1->y };
    // Vector2D pos2 = { particle2->x, particle2->y };
    // Vector2D vel1 = { particle1->vx, particle1->vy };
    // Vector2D vel2 = { particle2->vx, particle2->vy };

    Vector2D pos1 = { high_prec(particle1->x), high_prec(particle1->y) };
    Vector2D pos2 = { high_prec(particle2->x), high_prec(particle2->y) };
    Vector2D vel1 = { high_prec(particle1->vx), high_prec(particle1->vy) };
    Vector2D vel2 = { high_prec(particle2->vx), high_prec(particle2->vy) };

    //cout << "vel1: " << vel1.x << ", " << vel1.y << endl;
    //cout << "vel2: " << vel2.x << ", " << vel2.y << endl;




    //double distance_pre = hypot(particle1->x - particle2->x, particle1->y - particle2->y);
    high_prec distance_pre = hypot(particle1->x - particle2->x, particle1->y - particle2->y);
    //cout << "Distance before collission: " << distance_pre << endl;


   

    // 2. Calculate relative velocity and position difference
    Vector2D rel_pos = pos2 - pos1;

    //check if rel_pos is zero
    //if (rel_pos.x == 0 && rel_pos.y == 0) {
    //    cout << "rel_pos is zero" << endl;
        
    //}

    //cout << "rel_pos: " << rel_pos.x << ", " << rel_pos.y << endl;
    Vector2D rel_vel = vel2 - vel1;
    //cout << "rel_vel: " << rel_vel.x << ", " << rel_vel.y << endl;

    //check if rel_vel is zero
    if (rel_vel.x == 0 && rel_vel.y == 0) {
        cout << "rel_vel is zero" << endl;
        //add a very small positive value to the first particle's x velocity
        particle1->vx += 0.0000000001;
        
    }

    // 3. Calculate quadratic terms for collision detection with high precision
    high_prec a = dot(rel_vel, rel_vel);  // Coefficient of t^2

    //cout << "a: " << a << endl;

    high_prec b = 2 * dot(rel_pos, rel_vel);  // Coefficient of t

    //cout << "b: " << b << endl;

    high_prec c = dot(rel_pos, rel_pos) - pow(particle1->rad + particle2->rad, 2);  // Constant term

    //cout << "c: " << c << endl;

    // 4. Solve quadratic equation for time of collision
    high_prec discriminant = b * b - 4 * a * c;

    //cout << "Discriminant: " << discriminant << endl;

    if (discriminant < 0) {
        // No collision occurs
        return 1.0;  // Full timestep remains
    }

    high_prec sqrt_discriminant = sqrt(discriminant);
    high_prec t_collision1 = (-b - sqrt_discriminant) / (2 * a);
    high_prec t_collision2 = (-b + sqrt_discriminant) / (2 * a);

    //cout << "t_collision1: " << t_collision1 << endl;
    //cout << "t_collision2: " << t_collision2 << endl;

    // Use the value for t_collission closest to 0

    high_prec t_collission1_distfrom0 = abs(t_collision1);
    high_prec t_collission2_distfrom0 = abs(t_collision2);


    high_prec t_collision;
    if (t_collission1_distfrom0 < t_collission2_distfrom0) {
        // Use t_collision1

        t_collision = t_collision1;
        
    }
    else {

        t_collision = t_collision2;
  
    }

    
      



    // 5. Backtrack positions to collision point
    pos1 = pos1 + t_collision * vel1;
    pos2 = pos2 + t_collision * vel2;

    particle1->x = pos1.x.convert_to<double>();  // Convert back to double for storage
    particle1->y = pos1.y.convert_to<double>();
    particle2->x = pos2.x.convert_to<double>();
    particle2->y = pos2.y.convert_to<double>();

    //debug: validate that distance is equal to the sum of the radii if not, print the values
    //double distance = sqrt(pow(pos1.x - pos2.x, 2) + pow(pos1.y - pos2.y, 2)); old formula for distance which is less precise

    high_prec distance = hypot(pos1.x - pos2.x, pos1.y - pos2.y);
    high_prec distance_error = abs(distance - (particle1->rad + particle2->rad));

    double tolerance_distance = 0.0000000001;
        //print in scientific notation

    if (distance_error > tolerance_distance) {
        cout << "Distance error: " << scientific << distance_error << endl;
        
    }

    //cout << "Distance after collission: " << distance << endl;

    // 6. Calculate the time scaler
    double time_scaler = 1.0 - t_collision.convert_to<double>();  // Remaining time in the timestep

    //cout << "Time scaler: " << time_scaler << endl;

    return time_scaler;
}



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


void Engine::resolve_gravity(shared_ptr<Particles> particles) {
    // This function resolves the gravitational attraction between particles

    // Define a small value to avoid division by zero

    high_prec epsilon = 0.0000000001;

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