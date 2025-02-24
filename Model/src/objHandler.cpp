#include "../include/ObjHandler.h"
#include "../include/InitStructs.h"
#include "../include/MathUtils.h"
#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"
#include "../include/Particles.h"

#include <iostream>
#include <fstream>
#include <memory>
#include <cmath>

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


    //1. Retrieve the object_inputs from csv file
    string object_input_path = "Inputs/object_inputs.csv";


    //2. Create csv reader object
    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;

    const int column_count = 17;

    io::CSVReader<column_count, TrimPolicy, QuotePolicy> in(object_input_path);

    //3. Read the contents of the csv file and store each row in a object object. Store all object objects in a vector.
    objects object_list;

    string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17;

    int objects_loaded = 1;

    //discard the header row. This list is a guide to the columns in the csv file.
    in.read_header(io::ignore_extra_column, "OBJECT_NAME", "R", "G", "B", "X", "Y", "Z", "VX", "VY", "VZ", "M", "RAD", "REST", "TEMP","COMPLEXITY", "COMPLEXITY_SIZE", "COMPLEXITY_N");


    while (in.read_row(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16, col17)) {
    

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

        
        //add the new scenario to the list of scenarios
        object_list.object_list.push_back(move(new_object));

        objects_loaded++;

    }

    //print object id and names to the user

    //for (int i = 0; i < object_list.object_list.size(); i++) {
    //    cout << object_list.object_list[i]->object_id << ". " << object_list.object_list[i]->name << endl;
    //}


    //4. Retrieve the objects from the object_list that are in the obj_list of the selected scenario

    shared_ptr<objects> requested_objects(new objects);
    

    for (auto &object : object_list.object_list) {
        if (scenario->obj_list.find(object->name) != string::npos) {
            requested_objects->object_list.push_back(move(object));
        }
    }

    //5 communicate the found objects to the user

    cout << "The following objects have been succesfully selected:" << endl;

    for (int i = 0; i < requested_objects->object_list.size(); i++) {
        cout << requested_objects->object_list[i]->object_id << ". " << requested_objects->object_list[i]->name << endl;
    }

    

    //6. Flatten the objects and store all objects in a single struct

    shared_ptr<Particles> particles = flatten_objs(requested_objects, scenario);

    //7. Remove overlapping particles

    remove_overlaps(particles);

    //8 Inform user of the number of particles loaded

    cout << "Loaded " << particles->particle_list.size() << " particles." << endl;

    //9. Return the particles struct
    
    //print the xv and yv of the first particle
    //cout << "vx: " << particles->particle_list[0]->vx << endl;
    //cout << "vy: " << particles->particle_list[0]->vy << endl;


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

    //return the particles struct

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

shared_ptr<Particles> ObjHandler::flatten_complex_circle(shared_ptr<object> complex_object) {
    //This function will flatten a complex circle object and store it in a Particles struct

    //create a Particles struct to store the flattened objects as particles
    shared_ptr<Particles> particles(new Particles);

    //1. retrieve complexity parameters

    double circle_rad = complex_object->complexity_size;
    int complexity_n = complex_object->complexity_n;
    Vector2D center = { complex_object->x, complex_object->y };


    //loop through n*3 times the complexity_n and sample random points within the circle radius (defined by complexity_size)
    for (int i = 0; i < 4*complexity_n; i++) {

        //check if the number of particles loaded is equal to the complexity_n
        if (particles->particle_list.size() == complexity_n) {
            break;
        }
        
        //store the sampled point in a particle struct and add it to the particles struct
        unique_ptr<Particle> particle(new Particle);
        
        //sample a random point within the circle radius using the above formula

        Vector2D sample_point = sample_in_circle(center, circle_rad);

        


        particle->particle_id = i;
        particle->r = complex_object->r;
        particle->g = complex_object->g;
        particle->b = complex_object->b;
        particle->x = sample_point.x.convert_to<double>();
        particle->y = sample_point.y.convert_to<double>();
        particle->z = 0;
        particle->vx = complex_object->vx;
        particle->vy = complex_object->vy;

        //print vx and vy
        //cout << "vx: " << particle->vx << " vy: " << particle->vy << endl;


        particle->vz = complex_object->vz;
        particle->rad = complex_object->rad;
        particle->rest = complex_object->rest;
        particle->temp = complex_object->temp;


        particles->particle_list.push_back(move(particle));
    }

    //remove overlapping particles
    
    remove_overlaps(particles);

    //set mass equal to complex_object mass divided by size of particles
    double m_i = complex_object->m / particles->particle_list.size();

    for (int i = 0; i < particles->particle_list.size(); i++) {
        particles->particle_list[i]->m = m_i;
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
                << complex_object->complexity << ","
                << complex_object->complexity_size << ","
                << complex_object->complexity_n << "\n";
    }

    // 5. Close the cache file
    cache_file.close();

    cout << "Object " << complex_object->name << " saved to object cache." << endl;
}


void ObjHandler::remove_overlaps(shared_ptr<Particles> particles) {
    // Sort particles by particle_id to ensure a consistent processing order
    sort(particles->particle_list.begin(), particles->particle_list.end(),
        [](shared_ptr<Particle> a, shared_ptr<Particle> b) {
            return a->particle_id < b->particle_id;
        });

    vector<shared_ptr<Particle>> non_overlapping_particles;
    int removed_overlaps = 0;

    // Iterate over each particle
    for (int i = 0; i < particles->particle_list.size(); ++i) {
        bool overlap_found = false;

        // Compare the current particle with all remaining particles
        for (int j = i + 1; j < particles->particle_list.size(); ++j) {
            if (remove_overlap(particles->particle_list[i], particles->particle_list[j])) {
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
    double distance = sqrt(pow(particle1->x - particle2->x, 2) +
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