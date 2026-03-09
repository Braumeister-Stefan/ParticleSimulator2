#include "../include/ObjHandler.h"
#include "../include/InitStructs.h"
#include "../include/MathUtils.h"
#include "../include/PhysEngineCore.h"
#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"
#include "../include/Particles.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <limits>

//to use pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// =======================
// CONFIGURATION
// =======================
// Toggle overlap relaxation post-processing when saving objects to cache.
// Set to false to skip relaxation (objects may blow up on reload if overlaps exist).
static bool RELAX_OVERLAPS = true;


//namespaces
using namespace std;

// Constructor
ObjHandler::ObjHandler() {
    cout << "Interfacer initialized." << endl;
}

// Destructor
ObjHandler::~ObjHandler() {
    cout << "Interfacer destroyed." << endl;
}



// Methods for interfacing with the code
shared_ptr<Particles> ObjHandler::process_objs(shared_ptr<scenario> scenario) {
    //This function will load the list of objects for the selected scenario and flatten them into a single struct.

    // Set runtime post-processing toggle from selected scenario input.
    RELAX_OVERLAPS = scenario->relax_overlaps;


    //1. Retrieve the object_inputs from csv file
    string object_input_path = "Inputs/object_inputs.csv";


    //2. Create csv reader object
    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;

    const int column_count = 18;

    io::CSVReader<column_count, TrimPolicy, QuotePolicy> in(object_input_path);

    //3. Read the contents of the csv file and store each row in a object object. Store all object objects in a vector.
    objects object_list;

    string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18;

    int objects_loaded = 1;

    //discard the header row. This list is a guide to the columns in the csv file.
    in.read_header(io::ignore_extra_column, "OBJECT_NAME", "R", "G", "B", "X", "Y", "Z", "VX", "VY", "VZ", "M", "RAD", "REST", "TEMP","COMPLEXITY", "COMPLEXITY_SIZE", "COMPLEXITY_N", "OMEGA");


    while (in.read_row(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17, col18)) {
    

        unique_ptr<object> new_object(new object);

        new_object->object_id = objects_loaded;

        new_object->name = col1;

        new_object->r = safe_stod(col2);   
        new_object->g = safe_stod(col3);   
        new_object->b = safe_stod(col4);   
        new_object->x = safe_stod(col5);  
        //cout << "x: " << new_object->x << endl; 
        new_object->y = safe_stod(col6);   
        new_object->z = safe_stod(col7);
        new_object->vx = safe_stod(col8);  
        new_object->vy = safe_stod(col9);  

        //cout << "vx: " << col8 << endl;
        //cout << "vy: " << col9 << endl;



        new_object->vz = safe_stod(col10);
        new_object->m = safe_stod(col11);  
        new_object->rad = safe_stod(col12); 
        new_object->rest = safe_stod(col13);
        new_object->temp = safe_stod(col14);
        

        if (col15.empty()) {
            new_object->complexity = "SIMPLE";
        } else {
            new_object->complexity = col15;
        }

        new_object->complexity_size = safe_stod(col16);

        new_object->complexity_n = safe_stoi(col17);

        new_object->omega = safe_stod(col18);

        
        //add the new scenario to the list of scenarios
        object_list.object_list.push_back(move(new_object));

        objects_loaded++;

    }

    //print object id and names to the user

    // cout << "The following objects were found in object_inputs.csv:" << endl;

    // for (int i = 0; i < object_list.object_list.size(); i++) {
    //     cout << object_list.object_list[i]->object_id << ". " << object_list.object_list[i]->name << endl;
    // }


    //4a. Retrieve the objects from the object_list that are in the obj_list of the selected scenario

    shared_ptr<objects> requested_objects(new objects);
    

    for (auto &object : object_list.object_list) {
        if (scenario->obj_list.find(object->name) != string::npos) {
            requested_objects->object_list.push_back(move(object));
        }
    }

    //4b print the list of objects in the scenario obj_list, which could not be found in the object_inputs.csv file
    string cache_found_names;

    shared_ptr<Particles> sim_states(new Particles);

    auto trim = [](string& s) {
        while (!s.empty() && (s.front()==' ' || s.front()=='\t')) s.erase(s.begin());
        while (!s.empty() && (s.back() ==' ' || s.back() =='\t')) s.pop_back();
    };

    istringstream ss(scenario->obj_list);
    for (string name; getline(ss, name, ','); ) {
        trim(name);
        bool found = false;


        for (auto &object : requested_objects->object_list) {
            if (object && object->name == name) { found = true; break; }
        }

        if (!found) {
            

            if (auto cached = obj_from_cache(name)) {
                //cout << "Object found in " << name << ".csv\n";

                if (!cache_found_names.empty()) cache_found_names += ", ";
                cache_found_names += name;


                for (auto &p : cached->particle_list) sim_states->particle_list.push_back(std::move(p));
            } else {
                cout << "Object not found in object_inputs.csv or rendered cache!\n" << endl;
            }
        }
    }
    //5 communicate the found objects to the user



    cout << "The following objects have been succesfully selected:" << endl;

    for (int i = 0; i < requested_objects->object_list.size(); i++) {
        cout << requested_objects->object_list[i]->object_id << ". "
            << requested_objects->object_list[i]->name << endl;
    }

    // append cached/simulated objects as extra list entries, marked with '*'
    if (!cache_found_names.empty()) {
        std::istringstream css(cache_found_names);
        std::string n;
        int id = static_cast<int>(requested_objects->object_list.size()) + 1;

        while (std::getline(css, n, ',')) {
            while (!n.empty() && (n.front()==' ' || n.front()=='\t')) n.erase(n.begin());
            while (!n.empty() && (n.back() ==' ' || n.back() =='\t')) n.pop_back();
            if (n.empty()) continue;

            cout << "*" << id++ << ". " << n << endl;
        }
    }
    

    //6a. Flatten the objects and store all objects in a single struct

    shared_ptr<Particles> particles = flatten_objs(requested_objects, scenario);


    //6b. Add sim_states particles to particles
    for (auto &p : sim_states->particle_list) particles->particle_list.push_back(std::move(p));



    //7. Remove overlapping particles with slop 

    remove_overlaps2(particles, scenario->collision_distance_tolerance);
    
    //8. Return the particles struct



    return particles;



}



shared_ptr<Particles> ObjHandler::flatten_objs(shared_ptr<objects> requested_objects, shared_ptr<scenario> scenario) {
    //This function will flatten the simple & complex objects in the requested_objects list and store all objects in a single struct.

    //create a Particles struct to store the flattened objects as particles

    shared_ptr<Particles> particles(new Particles);

    //loop through the requested_objects list and if the object is complex, flatten it and store the particles in the particles struct, else store the object in the particles struct


    int particles_loaded = 0;
    for (auto &object : requested_objects->object_list) {

        shared_ptr<Particles> complex_particles; // Ensure complex_particles is always declared

        if (object->complexity == "SIMPLE") {
            
            

            unique_ptr<Particle> particle = flatten_simple_obj(particles_loaded, object);
            

            particles -> particle_list.push_back(move(particle));

            particles_loaded++;

        } else if (object->complexity == "CIRCLE") {
            //flatten the complex object and store the particles in the particles struct
            //if the refresh_obj flag is false, attempt to retrieve the particles from the cache
            if (!scenario->refresh_obj) {
                complex_particles = obj_from_cache(object->name);

                if (complex_particles == nullptr) {
                    complex_particles = flatten_complex_obj(object);
                } else {
                    cout << "Loaded " << object->name << " from cache." << endl;
                }
            } else {
                complex_particles = flatten_complex_obj(object);
            }

        } else if (obj_from_cache(object->complexity) != nullptr) {

            
            complex_particles = obj_from_cache(object->complexity);
            //manipulate the positions, velocities and rotation of the cached object.
            overwrite_rendered_obj(complex_particles, object);

        } else {
            cout << "Complex object is a " << object->complexity << ". This shape is not supported." << endl;
            continue;
        }

        //add the particles to the particles struct
        if (complex_particles) {
            particles->particle_list.insert(particles->particle_list.end(), complex_particles->particle_list.begin(), complex_particles->particle_list.end());
            particles_loaded += static_cast<int>(complex_particles->particle_list.size());
        }

        //cout << "particle_list has size " << particles->particle_list.size() << endl;
    }

    // Reset particle_id for all particles
    for (int i = 0; i < static_cast<int>(particles->particle_list.size()); ++i) {
        if (particles->particle_list[i]) {
            particles->particle_list[i]->particle_id = i;
        }
    }

    return particles;

}




unique_ptr<Particle> ObjHandler::flatten_simple_obj(int particles_loaded, shared_ptr<object> simple_object) {
    //This function will flatten a simple object and store it in a Particles struct

    //store the object in the particles struct
            unique_ptr<Particle> particle(new Particle);

            particle -> particle_id = particles_loaded;
            particle->r = simple_object->r;
            particle->g = simple_object->g;
            particle->b = simple_object->b;
            particle->x = simple_object->x;
            particle->y = simple_object->y;
            particle->z = simple_object->z;
            particle->vx = simple_object->vx;
            particle->vy = simple_object->vy;
            particle->vz = simple_object->vz;
            particle->m = simple_object->m;
            particle->rad = simple_object->rad;
            particle->rest = simple_object->rest;
            particle->temp = simple_object->temp;

            return particle;

}


shared_ptr<Particles> ObjHandler::flatten_complex_obj(shared_ptr<object> requested_object) {
    //This function will flatten a complex object and store it in a Particles struct

    //1. create a Particles struct to store the flattened object as particles
    shared_ptr<Particles> unfolded_particles(new Particles);

    //2 If complex flag is"circle", call flatten_complex_circle


    if (requested_object->complexity == "CIRCLE") {
            
        unfolded_particles = flatten_complex_circle(requested_object);
      
    } else {

        cout << "Complex object is a " << requested_object->complexity << ". This shape is not supported." << endl;
        
    }

        

    cout << "Object " << requested_object->name << " flattened into " << unfolded_particles->particle_list.size() << " particles." << endl;

    //3. Save the complex object to the cache

    obj_to_cache(requested_object, unfolded_particles);


        
    

    return unfolded_particles;

}


void add_rotation(Particle& particle, const Vector2D& center, const high_prec& omega) {
    high_prec dx = particle.x - center.x;
    high_prec dy = particle.y - center.y;

    high_prec dvx = -omega * dy;
    high_prec dvy =  omega * dx;

    particle.vx += dvx;
    particle.vy += dvy;
}


shared_ptr<Particles> ObjHandler::flatten_complex_circle(shared_ptr<object> complex_object) {
    shared_ptr<Particles> particles(new Particles);

    high_prec circle_rad = complex_object->complexity_size;
    int complexity_n = static_cast<int>(complex_object->complexity_n);
    Vector2D center = { complex_object->x, complex_object->y, complex_object->z };

    // assume omega exists on complex_object
    high_prec omega = complex_object->omega;

    int max_attempts = 10 * std::max(1, complexity_n);
    int last_pct = -1;

    cout << "Sampling " << complexity_n << " particles (batch mode, " << max_attempts << " candidates): ";
    cout.flush();

    // --- Batch sampling: generate all candidates up-front, then greedily accept ---

    // 1. Generate the full batch of candidate particles
    vector<shared_ptr<Particle>> candidates;
    candidates.reserve(max_attempts);

    for (int b = 0; b < max_attempts; ++b) {
        auto particle = make_shared<Particle>();

        Vector2D sample_point = sample_in_circle(center, circle_rad);

        particle->particle_id = b;
        particle->r    = complex_object->r;
        particle->g    = complex_object->g;
        particle->b    = complex_object->b;
        particle->x    = sample_point.x;
        particle->y    = sample_point.y;
        particle->z    = complex_object->z;
        particle->vx   = complex_object->vx;
        particle->vy   = complex_object->vy;
        particle->vz   = complex_object->vz;
        particle->rad  = complex_object->rad;
        particle->rest = complex_object->rest;
        particle->temp = complex_object->temp;

        add_rotation(*particle, center, omega);

        candidates.push_back(std::move(particle));

        // progress bar
        int pct = ((b + 1) * 100) / max_attempts;
        if (pct > last_pct) {
            for (int p = last_pct + 1; p <= pct; ++p) cout << "-";
            cout.flush();
            last_pct = pct;
        }
    }

    // 2. Greedy accept: walk through candidates, keep those that don't overlap any accepted particle
    const high_prec two_rad = complex_object->rad * 2.0;
    const high_prec two_rad_sq = two_rad * two_rad;

    particles->particle_list.reserve(complexity_n);

    for (int c = 0; c < static_cast<int>(candidates.size()) &&
             static_cast<int>(particles->particle_list.size()) < complexity_n; ++c) {
        const auto& cand = candidates[c];
        bool overlaps = false;

        for (int j = static_cast<int>(particles->particle_list.size()) - 1; j >= 0; --j) {
            const auto& placed = particles->particle_list[j];
            high_prec dx = cand->x - placed->x;
            high_prec dy = cand->y - placed->y;
            if (dx * dx + dy * dy < two_rad_sq) {
                overlaps = true;
                break;
            }
        }

        if (!overlaps) {
            particles->particle_list.push_back(cand);
        }
    }

    cout << " done (" << particles->particle_list.size() << "/" << complexity_n
         << " placed from " << max_attempts << " candidates)" << endl;

    // set mass equal to complex_object mass divided by final particle count
    if (!particles->particle_list.empty()) {
        high_prec m_i = complex_object->m / particles->particle_list.size();
        for (int i = 0; i < static_cast<int>(particles->particle_list.size()); i++) {
            particles->particle_list[i]->m = m_i;
        }
    }

    return particles;
}

shared_ptr<Particles> ObjHandler::obj_from_cache(string obj_name){
    //This function will attempt to retrieve the particles from the cache

    //1. check if the cache file exists at Inputs\rendered_objects\obj_name.csv

    bool cache_exists = false;

    string cache_path = "Inputs/rendered_objects/" + obj_name + ".csv";

    ifstream cache_file(cache_path);

    if (cache_file.good()) {
        cache_exists = true;
        //cout << "Cache file found for object " << obj_name << "." << endl;
    }

    //2. if the cache file exists, read the particles from the cache file and store them in a particles struct

    if (cache_exists){

        shared_ptr<Particles> particles(new Particles);

        //create a csv reader object
        typedef io::trim_chars<' ', '\t'> TrimPolicy;
        typedef io::double_quote_escape<',', '\"'> QuotePolicy;

        const int column_count = 17;

        

        io::CSVReader<column_count, TrimPolicy, QuotePolicy> in(cache_path);

        //3. Read the contents of the csv file and store each row in a particle object. Store all particle objects in a vector.
        string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17;

        int particles_loaded = 1;
        
        //discard the header row. This list is a guide to the columns in the csv file.
        in.read_header(io::ignore_extra_column, "PARTICLE_ID", "R", "G", "B", "X", "Y", "Z", "VX", "VY", "VZ", "M", "RAD", "REST", "TEMP","COMPLEXITY", "COMPLEXITY_SIZE", "COMPLEXITY_N");
        
        while (in.read_row(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17)) {
    
            //cout << "Loading particle " << particles_loaded << " from cache." << endl;
            unique_ptr<Particle> new_particle(new Particle);

            new_particle->particle_id = particles_loaded;

            new_particle->r = safe_stod(col2);   
            new_particle->g = safe_stod(col3);   
            new_particle->b = safe_stod(col4);   
            new_particle->x = safe_stod(col5);  
            new_particle->y = safe_stod(col6);   
            new_particle->z = safe_stod(col7);
            new_particle->vx = safe_stod(col8);  
            new_particle->vy = safe_stod(col9);  
            new_particle->vz = safe_stod(col10);
            new_particle->m = safe_stod(col11);  
            new_particle->rad = safe_stod(col12); 
            new_particle->rest = safe_stod(col13);
            new_particle->temp = safe_stod(col14);


            //add the new particle to the list of particles
            particles->particle_list.push_back(move(new_particle));

            particles_loaded++;

            
        }

    
        return particles;





    } else {
        return nullptr;
    }


}

void ObjHandler::obj_to_cache(shared_ptr<object> complex_object, shared_ptr<Particles> particles) {
    // This function will attempt to save the particles to the cache

    // 1. Define the cache file path
    string cache_path = "Inputs/rendered_objects/" + complex_object->name + ".csv";

    // 2. Open the cache file for writing
    ofstream cache_file(cache_path);

    if (!cache_file.is_open()) {
        cout << "Object " << complex_object->name << " could not be saved to object cache." << endl;
        return;
    }

    // 3. Write the header row
    cache_file << "PARTICLE_ID,R,G,B,X,Y,Z,VX,VY,VZ,M,RAD,REST,TEMP,COMPLEXITY,COMPLEXITY_SIZE,COMPLEXITY_N\n";

    // 4. Write the particles to the cache file

    //define complex object parameters
    

    for (const auto& particle : particles->particle_list) {
        cache_file << particle->particle_id << ","
                << particle->r << ","
                << particle->g << ","
                << particle->b << ","
                << particle->x << ","
                << particle->y << ","
                << particle->z << ","
                << particle->vx << ","
                << particle->vy << ","
                << particle->vz << ","
                << particle->m << ","
                << particle->rad << ","
                << particle->rest << ","
                << particle->temp << ","
                << complex_object->complexity << ","
                << complex_object->complexity_size << ","
                << int(1) << "\n";


  
    }

    // 5. Close the cache file
    cache_file.close();

    cout << "Object " << complex_object->name << " saved to object cache." << endl;
}

void ObjHandler::state_to_cache(shared_ptr<Particles> final_state, string obj_name) {
    // This function will save the final particle state as a complex object in the cache

    // 1. Define the cache file path
    string cache_path = "Inputs/rendered_objects/" + obj_name + ".csv";

    // 2. Open the cache file for writing
    ofstream cache_file(cache_path);

    


    // if obj_name is empty/FALSE, return
    if (obj_name.empty() || obj_name == "FALSE") {
        return;
    }
    if (!cache_file.is_open()) {
        cout << "Final state could not be saved to object cache as " << obj_name << "." << endl;
        cout <<"Attempting to force close the file." << endl;
        //return;
        cache_file.close();

        //try to open again
        cache_file.open(cache_path);
    }
 




    // 3. Post-process: resolve any remaining particle overlaps before writing.
    //    Strictly positional — no velocity changes.
    if (RELAX_OVERLAPS) {
        cout << "Running overlap relaxation before saving..." << endl;
        int relax_iters = relax_overlaps(final_state);
        if (relax_iters == 0) {
            cout << "No overlaps detected - object is clean." << endl;
        }
    } else {
        cout << "Overlap relaxation disabled - saving without post-processing." << endl;
    }

    // 4. Write the header row
    cache_file << "PARTICLE_ID,R,G,B,X,Y,Z,VX,VY,VZ,M,RAD,REST,TEMP,COMPLEXITY,COMPLEXITY_SIZE,COMPLEXITY_N\n";

    // 5. Write the particles to the cache file
    for (const auto& particle : final_state->particle_list) {
        cache_file << particle->particle_id << ","
                << particle->r << ","
                << particle->g << ","
                << particle->b << ","
                << particle->x << ","
                << particle->y << ","
                << particle->z << ","
                << particle->vx << ","
                << particle->vy << ","
                << particle->vz << ","
                << particle->m << ","
                << particle->rad << ","
                << particle->rest << ","
                << particle->temp << ","
                << obj_name << ","
                << int(1) << ","
                << int(1) << "\n";

        


    }

    // 6. Close the cache file
    cache_file.close();

    cout << "Simulation result saved as object: " <<  obj_name << "." << endl;
}

// =======================
// ITERATIVE POSITIONAL RELAXATION  (Jacobi + under-relaxation)
// =======================
// Pairwise center-to-center direction. Each iteration:
//   1. Scan all pairs — for each overlapping pair, compute the displacement
//      needed to separate them along the center-to-center axis, split by
//      inverse mass (heavier particle moves less, CoM preserved).
//   2. Accumulate displacements into per-particle buffers (Jacobi: no
//      position is written during the scan).
//   3. Apply all accumulated displacements at once, scaled by a relaxation
//      factor to prevent overcorrection cascades in dense packings.
//
// Strictly positional — velocities are never touched.
// Terminates when:
//   - Zero overlaps remain, OR
//   - max_iterations reached, OR
//   - cumulative overlap has not improved for 100 consecutive iterations.
int ObjHandler::relax_overlaps(shared_ptr<Particles> particles, int max_iterations, high_prec separation_margin) {
    auto& list = particles->particle_list;
    const int n = static_cast<int>(list.size());
    if (n < 2) return 0;

    // Default: 100 * particle count
    if (max_iterations < 0) max_iterations = 0.1 * n;

    // --- Helper: count overlaps and cumulative overlap ---
    auto measure_overlaps = [&](high_prec margin) -> std::pair<int, high_prec> {
        int count = 0;
        high_prec cumulative = 0.0;
        for (int i = 0; i < n; ++i) {
            const Particle* pi = list[i].get();
            for (int j = i + 1; j < n; ++j) {
                const Particle* pj = list[j].get();
                high_prec dx = pj->x - pi->x;
                high_prec dy = pj->y - pi->y;
                high_prec dz = pj->z - pi->z;
                high_prec dist2 = dx * dx + dy * dy + dz * dz;
                high_prec target = pi->rad + pj->rad + margin;
                if (dist2 < target * target) {
                    high_prec dist = sqrt(dist2);
                    cumulative += (target - dist);
                    ++count;
                }
            }
        }
        return {count, cumulative};
    };

    // --- Initial measurements ---
    auto [initial_overlaps, initial_cumulative] = measure_overlaps(separation_margin);

    cout << "[relax_overlaps] " << n << " particles, "
         << initial_overlaps << " overlapping pair(s), "
         << "cumulative overlap: " << std::scientific << std::setprecision(4)
         << initial_cumulative << std::fixed << endl;
    cout << "[relax_overlaps] COLLISION_DIST_TOL = " << std::scientific << std::setprecision(4)
         << EngineCore::collision_distance_tolerance_ << std::fixed
         << " (will stop if worst overlap falls below this)" << endl;

    if (initial_overlaps == 0) return 0;

    // --- Jacobi buffers ---
    struct Disp { high_prec dx = 0, dy = 0, dz = 0; };
    std::vector<Disp> accum(n);

    // Under-relaxation factor
    const high_prec relaxation_factor = 0.4;

    int iter = 0;
    int total_corrections = 0;
    int last_reported_bucket = -5;

    // Stall detection: stop after 10 iterations with no improvement
    high_prec best_cumulative = initial_cumulative;
    int no_improvement_count = 0;
    const int no_improvement_limit = 10;

    // Relative-improvement stall: stop if improvement stays below 0.1%
    // for more than 10 consecutive iterations.
    high_prec prev_cumulative = initial_cumulative;
    int low_improvement_count = 0;
    const high_prec min_relative_improvement = 0.001; // 0.1%
    const int low_improvement_limit = 30;

    for (; iter < max_iterations; ++iter) {
        int corrections_this_pass = 0;
        high_prec worst_overlap = 0.0;
        high_prec cumulative_this_pass = 0.0;

        // Reset accumulation buffers
        for (int k = 0; k < n; ++k) { accum[k].dx = 0; accum[k].dy = 0; accum[k].dz = 0; }

        // --- Accumulate corrections (read current positions, don't write yet) ---
        for (int i = 0; i < n; ++i) {
            const Particle* pi = list[i].get();
            for (int j = i + 1; j < n; ++j) {
                const Particle* pj = list[j].get();

                high_prec dx = pj->x - pi->x;
                high_prec dy = pj->y - pi->y;
                high_prec dz = pj->z - pi->z;
                high_prec dist2 = dx * dx + dy * dy + dz * dz;
                high_prec sum_radii = pi->rad + pj->rad;
                high_prec target = sum_radii + separation_margin;

                if (dist2 >= target * target) continue;

                high_prec dist = sqrt(dist2);
                if (dist < 1e-30) {
                    dx = 1.0;
                    dy = 0.0;
                    dz = 0.0;
                    dist = 1e-30;
                }

                high_prec overlap = target - dist;
                cumulative_this_pass += overlap;
                if (overlap > worst_overlap) worst_overlap = overlap;

                high_prec nx = dx / dist;
                high_prec ny = dy / dist;
                high_prec nz = dz / dist;

                // Mass-weighted split (heavier particle moves less)
                high_prec total_mass = pi->m + pj->m;
                high_prec wi = pj->m / total_mass;
                high_prec wj = pi->m / total_mass;

                accum[i].dx -= wi * overlap * nx;
                accum[i].dy -= wi * overlap * ny;
                accum[i].dz -= wi * overlap * nz;
                accum[j].dx += wj * overlap * nx;
                accum[j].dy += wj * overlap * ny;
                accum[j].dz += wj * overlap * nz;

                ++corrections_this_pass;
            }
        }

        if (corrections_this_pass == 0) break; // fully relaxed

        // --- Apply accumulated displacements with under-relaxation ---
        for (int k = 0; k < n; ++k) {
            list[k]->x += relaxation_factor * accum[k].dx;
            list[k]->y += relaxation_factor * accum[k].dy;
            list[k]->z += relaxation_factor * accum[k].dz;
        }

        total_corrections += corrections_this_pass;

        // --- Relative improvement check on cumulative overlap ---
        high_prec improvement = prev_cumulative - cumulative_this_pass;
        high_prec denom = (prev_cumulative > 1e-30) ? prev_cumulative : 1e-30;
        high_prec relative_improvement = improvement / denom;
        if (relative_improvement < min_relative_improvement) {
            ++low_improvement_count;
        } else {
            low_improvement_count = 0;
        }
        prev_cumulative = cumulative_this_pass;

        if (low_improvement_count > low_improvement_limit) {
            cout << "[relax_overlaps] Relative improvement below "
                 << (min_relative_improvement * 100.0) << "% for more than "
                 << low_improvement_limit << " iterations (latest: "
                 << std::scientific << std::setprecision(4)
                 << (relative_improvement * 100.0) << "%"
                 << std::fixed << ") - stopping early." << endl;
            ++iter;
            break;
        }

        // --- Stall detection on cumulative overlap ---
        if (cumulative_this_pass < best_cumulative) {
            best_cumulative = cumulative_this_pass;
            no_improvement_count = 0;
        } else {
            ++no_improvement_count;
        }
        if (no_improvement_count >= no_improvement_limit) {
            cout << "[relax_overlaps] Cumulative overlap stalled at "
                 << std::scientific << std::setprecision(4) << cumulative_this_pass
                 << std::fixed << " for " << no_improvement_limit
                 << " iterations - stopping early." << endl;
            ++iter;
            break;
        }

        // --- Tolerance check: stop if worst overlap is below collision distance tolerance ---
        if (worst_overlap < EngineCore::collision_distance_tolerance_) {
            cout << "[relax_overlaps] Worst overlap " << std::scientific << std::setprecision(4)
                 << worst_overlap << " is below collision distance tolerance "
                 << EngineCore::collision_distance_tolerance_ << std::fixed
                 << " - stopping." << endl;
            ++iter;
            break;
        }

        // --- Progress logging (every 5%) ---
        if (max_iterations > 0) {
            int pct = (int)(((long long)(iter + 1) * 100LL) / (long long)max_iterations);
            int bucket = (pct / 0.001) * 0.001;
            if (bucket > last_reported_bucket && bucket <= 100) {
                cout << "[relax_overlaps] " << bucket << "% (" << (iter + 1) << "/" << max_iterations
                     << " iters) " << corrections_this_pass << " overlap(s), worst: "
                     << std::scientific << std::setprecision(3) << worst_overlap
                     << " (tol: " << std::setprecision(3) << EngineCore::collision_distance_tolerance_ << ")"
                     << ", cumulative: " << std::setprecision(4) << cumulative_this_pass
                     << ", rel_impr: " << std::setprecision(3) << (relative_improvement * 100.0) << "%"
                     << ", low-impr streak: " << low_improvement_count << "/" << low_improvement_limit
                     << std::fixed << endl;
                last_reported_bucket = bucket;
            }
        }
    }

    // --- Final measurements ---
    auto [final_overlaps, final_cumulative] = measure_overlaps(separation_margin);

    // --- Summary ---
    cout << "[relax_overlaps] Done. "
         << iter << " iteration(s), "
         << total_corrections << " positional correction(s)." << endl;
    cout << "[relax_overlaps] Overlaps: " << initial_overlaps << " -> " << final_overlaps
         << ". Cumulative overlap: " << std::scientific << std::setprecision(4)
         << initial_cumulative << " -> " << final_cumulative
         << std::fixed << endl;
    if (final_overlaps > 0) {
        cout << "[relax_overlaps] WARNING: " << final_overlaps
             << " overlap(s) remain after " << iter << " iterations." << endl;
    }

    return iter;
}

void ObjHandler::remove_overlaps(shared_ptr<Particles> particles) {
    // Iterate from the back so that newly-appended particles (at the end)
    // are the ones removed when they overlap an already-placed particle.
    for (int i = static_cast<int>(particles->particle_list.size()) - 1; i >= 0; --i) {
        bool overlap_found = false;
        for (int j = 0; j < i; ++j) {
            if (remove_overlap(particles->particle_list[i], particles->particle_list[j])) {
                overlap_found = true;
                break;
            }
        }
        if (overlap_found) {
            particles->particle_list.erase(particles->particle_list.begin() + i);
        }
    }
}

void ObjHandler::remove_overlaps2(shared_ptr<Particles> particles, high_prec slop) {
    //this function uses calculate_overlap_amount from PhysEngineCore to remove overlapping particles if overlap is greater than slop
    vector<shared_ptr<Particle>> non_overlapping_particles;
    int removed_overlaps = 0;
    EngineCore engine_core;

    // Iterate over each particle
    for (int i = 0; i < particles->particle_list.size(); ++i) {
        bool overlap_found = false;
        // Compare the current particle with all remaining particles
        for (int j = i + 1; j < particles->particle_list.size(); ++j) {
            high_prec overlap_amount = engine_core.calculate_overlap_amount(particles->particle_list[i], particles->particle_list[j]);
            if (overlap_amount > slop) {
                overlap_found = true;
                ++removed_overlaps;
                break; // Stop further checks for this particle since it's overlapped
            }
        }
        
        // If no overlap was found, keep the particle
        if (!overlap_found) {
            non_overlapping_particles.push_back(particles->particle_list[i]);
        }
    }
    // Replace the original particle list with the non-overlapping particles
    particles->particle_list = std::move(non_overlapping_particles);
    cout << "Removed " << removed_overlaps << " overlapping particles." << endl;
}


bool ObjHandler::remove_overlap(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    // Calculate the distance between the two particles
    high_prec distance = sqrt(pow(particle1->x - particle2->x, 2) +
                           pow(particle1->y - particle2->y, 2) +
                           pow(particle1->z - particle2->z, 2));

    // If particles overlap, determine which to keep based on mass
    if (distance < particle1->rad + particle2->rad) {
        if (particle1->m < particle2->m) {
            //cout << "Particle " << particle1->particle_id << " removed." << endl;
        } else {
            //cout << "Particle " << particle2->particle_id << " removed." << endl;
        }
        return true; // Indicate an overlap was found
    }
    return false; // No overlap found
}

double ObjHandler::safe_stod(string str) {

    // Function to safely convert string to double, defaulting to 0.0 if empty or non-convertible
    if (str.empty()) {
        return 0.0; // Return 0.0 if string is empty
    } 
    //else, if abs(stod(str)) < e-10, return 0.0
    else if (abs(stod(str)) < 1e-10) {
        return 0.0;
    }

    return stod(str); // Convert string to double
    
}

int ObjHandler::safe_stoi(string str) {

    // Function to safely convert string to int, defaulting to 0 if empty or non-convertible
    if (str.empty()) {
        return 0; // Return 0 if string is empty
    }

    return stoi(str); // Convert string to int
    
}



void ObjHandler::overwrite_rendered_obj(shared_ptr<Particles> complex_particles, shared_ptr<object> object) {
    // This function will overwrite the positions, velocities and rotation of the cached object with the parameters of the complex object

    Vector2D new_center = { object->x, object->y, object->z };
    high_prec boost_omega = object->omega;

    // 1. Calculate the current center of mass of the complex_particles

    Vector2D current_center = {0, 0, 0};
    for (const auto& particle : complex_particles->particle_list) {
        current_center.x += particle->x;
        current_center.y += particle->y;
        current_center.z += particle->z;
    }
    current_center.x /= complex_particles->particle_list.size();
    current_center.y /= complex_particles->particle_list.size();
    current_center.z /= complex_particles->particle_list.size();

    // if x,y,vx,vy, omega are empty, set value to zero

    // 2. Calculate the translation vector from the current center to the new center

    Vector2D translation = { new_center.x - current_center.x, new_center.y - current_center.y, new_center.z - current_center.z };

    // 3. Apply the translation and rotation to each particle in complex_particles
    for (const auto& particle : complex_particles->particle_list) {
        // Translate
        particle->x += translation.x;
        particle->y += translation.y;
        particle->z += translation.z;

        //update velocity 
        particle->vx += object->vx;
        particle->vy += object->vy;
        particle->vz += object->vz;


        // Rotate (in XY plane about z-axis)
        add_rotation(*particle, new_center, boost_omega);
    }
 

    cout << "Loaded " << object->name << " from pre-rendered template: " << object->complexity << endl;
}