#include "../include/ObjHandler.h"
#include "../include/InitStructs.h"
#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"
#include "../include/Particles.h"

#include <iostream>
#include <fstream>
#include <memory>

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

    //1. Load the objects from the csv file 

    cout << "The following objects are available:" << endl;

    //1. Retrieve the object_inputs from csv file
    string object_input_path = "Inputs/object_inputs.csv";


    //2. Create csv reader object
    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;

    const int column_count = 16;

    io::CSVReader<column_count, TrimPolicy, QuotePolicy> in(object_input_path);

    //3. Read the contents of the csv file and store each row in a object object. Store all object objects in a vector.
    objects object_list;

    string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16;

    int objects_loaded = 1;

    //discard the header row. This list is a guide to the columns in the csv file.
    in.read_header(io::ignore_extra_column, "OBJECT_NAME", "R", "G", "B", "X", "Y", "Z", "VX", "VY", "VZ", "M", "RAD", "REST", "COMPLEXITY", "COMPLEXITY_SIZE", "COMPLEXITY_N");


    while (in.read_row(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15, col16)) {
    

        unique_ptr<object> new_object(new object);

        new_object->object_id = objects_loaded;

        new_object->name = col1;

        new_object->r = safe_stod(col2);   
        new_object->g = safe_stod(col3);   
        new_object->b = safe_stod(col4);   
        new_object->x = safe_stod(col5);  
        new_object->y = safe_stod(col6);   
        new_object->z = safe_stod(col7);   
        new_object->vx = safe_stod(col8);  
        new_object->vy = safe_stod(col9);  
        new_object->vz = safe_stod(col10);
        new_object->m = safe_stod(col11);  
        new_object->rad = safe_stod(col12); 
        new_object->rest = safe_stod(col13);

        if (col14.empty()) {
            new_object->complexity = "simple";
        } else {
            new_object->complexity = col14;
        }

        new_object->complexity_size = safe_stod(col15);
        new_object->complexity_n = stoi(col16);

        
        //add the new scenario to the list of scenarios
        object_list.object_list.push_back(move(new_object));

        objects_loaded++;

    }

    //print object id and names to the user

    for (int i = 0; i < object_list.object_list.size(); i++) {
        cout << object_list.object_list[i]->object_id << ". " << object_list.object_list[i]->name << endl;
    }


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

    

    //6. Flatten the complex objects and store all objects in a single struct

    shared_ptr<Particles> particles = flatten_objs(requested_objects);

    //to implement, remove overlapping particles

    //7. Return the particles struct
    

    return particles;



}



shared_ptr<Particles> ObjHandler::flatten_objs(shared_ptr<objects> requested_objects) {
    //This function will flatten the complex objects in the requested_objects list and store all objects in a single struct.

    //create a Particles struct to store the flattened objects as particles

    shared_ptr<Particles> particles(new Particles);

    //loop through the requested_objects list and if the object is complex, flatten it and store the particles in the particles struct, else store the object in the particles struct


    int particles_loaded = 1;
    for (auto &object : requested_objects->object_list) {

        
        if (object->complexity == "simple") {
            
            

            unique_ptr<Particle> particle = flatten_simple_obj(particles_loaded, object);
            

            particles -> particle_list.push_back(move(particle));
        } else {
            //flatten the complex object and store the particles in the particles struct
            
            //TO BE IMPLEMENTED
        }

        particles_loaded++;
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
            particle->z = simple_object->z;
            particle->vx = simple_object->vx;
            particle->vy = simple_object->vy;
            particle->vz = simple_object->vz;
            particle->m = simple_object->m;
            particle->rad = simple_object->rad;
            particle->rest = simple_object->rest;

            return particle;

}


shared_ptr<Particles> ObjHandler::flatten_complex_obj(shared_ptr<objects> requested_objects) {
    //This function will flatten a complex object and store it in a Particles struct

    //create a Particles struct to store the flattened objects as particles

    //loop through the requested_objects list and if the object is complex, flatten it and store the particles in the particles struct, else store the object in the particles struct

    return nullptr;

}

double ObjHandler::safe_stod(string str) {
    // Function to safely convert string to double, defaulting to 0.0 if empty or non-convertible
    if (str.empty()) {
        return 0.0; // Return 0.0 if string is empty
    }

    return std::stod(str); // Convert string to double
    
}