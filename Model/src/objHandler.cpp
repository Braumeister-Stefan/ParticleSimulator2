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

//to use pi
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


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

    //0. List of objects to be loaded
    cout << "Loading objects for the selected scenario:" << endl;

    for (const auto& name : scenario->obj_list) {
        cout << name;
    }


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
        new_object->z = 0;  
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
            new_object->complexity = "simple";
        } else {
            new_object->complexity = col15;
        }

        new_object->complexity_size = safe_stod(col16);
        new_object->complexity_n = stoi(col17);
        new_object->omega = safe_stod(col18);

        
        //add the new scenario to the list of scenarios
        object_list.object_list.push_back(move(new_object));

        objects_loaded++;

    }

    //print object id and names to the user

    //for (int i = 0; i < object_list.object_list.size(); i++) {
    //    cout << object_list.object_list[i]->object_id << ". " << object_list.object_list[i]->name << endl;
    //}


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

    remove_overlaps2(particles, 0.1);
    

    //8 Inform user of the number of particles loaded

    cout << "Loaded " << particles->particle_list.size() << " particles." << endl;

    //9. Return the particles struct



    return particles;



}



shared_ptr<Particles> ObjHandler::flatten_objs(shared_ptr<objects> requested_objects, shared_ptr<scenario> scenario) {
    //This function will flatten the simple & complex objects in the requested_objects list and store all objects in a single struct.

    //create a Particles struct to store the flattened objects as particles

    shared_ptr<Particles> particles(new Particles);

    //loop through the requested_objects list and if the object is complex, flatten it and store the particles in the particles struct, else store the object in the particles struct


    int particles_loaded = 0;
    for (auto &object : requested_objects->object_list) {

        
        if (object->complexity == "simple") {
            
            

            unique_ptr<Particle> particle = flatten_simple_obj(particles_loaded, object);
            

            particles -> particle_list.push_back(move(particle));

            particles_loaded++;

        } else {
            //flatten the complex object and store the particles in the particles struct
            shared_ptr<Particles> complex_particles(new Particles);
           
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

            //add the particles to the particles struct
            particles->particle_list.insert(particles->particle_list.end(), complex_particles->particle_list.begin(), complex_particles->particle_list.end());
            

           particles_loaded += particles->particle_list.size();

            cout << "particle_list has size " << particles->particle_list.size() << endl;

        }

        
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
            particle->z = 0;
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
    Vector2D center = { complex_object->x, complex_object->y };

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
        particle->z    = 0;
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
    cout << "Checking for cache file at " << cache_path << endl;

    ifstream cache_file(cache_path);

    if (cache_file.good()) {
        cache_exists = true;
        cout << "Cache file found for object " << obj_name << "." << endl;
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
            new_particle->z = 0;   
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
 




    // 3. Write the header row
    cache_file << "PARTICLE_ID,R,G,B,X,Y,Z,VX,VY,VZ,M,RAD,REST,TEMP,COMPLEXITY,COMPLEXITY_SIZE,COMPLEXITY_N\n";

    // 4. Write the particles to the cache file

    

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

    // 5. Close the cache file
    cache_file.close();

    cout << "Simulation result saved as object: " <<  obj_name << "." << endl;
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