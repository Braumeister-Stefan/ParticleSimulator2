#ifndef OBJ_HANDLER_H
#define OBJ_HANDLER_H

#include <iostream>
#include <vector>
#include <memory>
#include "../include/InitStructs.h"

using namespace std;

class ObjHandler {
public:
    // Constructor
    ObjHandler();

    // Destructor
    ~ObjHandler();

    // Declare the process object function
    shared_ptr<Particles> process_objs(shared_ptr<scenario> scenario);

    //Declare the flatten objects function

    shared_ptr<Particles> flatten_objs(shared_ptr<objects> requested_objects, shared_ptr<scenario> scenario);

    //Declare the flatten simple objects function

    unique_ptr<Particle> flatten_simple_obj(int particles_loaded, shared_ptr<object> requested_objects);

    // Declare the flatten complex objects functions
    shared_ptr<Particles> flatten_complex_obj(shared_ptr<object> requested_object);
    shared_ptr<Particles> flatten_complex_circle(shared_ptr<object> complex_object);


    //Declare the remove overlap functions
    void remove_overlaps(shared_ptr<Particles> particles);
    bool remove_overlap(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    void remove_overlaps2(shared_ptr<Particles> particles, high_prec slop);

    // Post-processing: Jacobi positional relaxation to eliminate overlaps
    // before saving to cache. Strictly positional — no velocity changes.
    // max_iterations defaults to 100 * particle count when not specified.
    // Returns the number of iterations used (0 = no overlaps found).
    int relax_overlaps(shared_ptr<Particles> particles, int max_iterations = -1, high_prec separation_margin = 1e-10);
    

    //cache functions
    shared_ptr<Particles> obj_from_cache(string obj_name);
    void obj_to_cache(shared_ptr<object> complex_object, shared_ptr<Particles> particles);
    void state_to_cache(shared_ptr<Particles> particles, string obj_name);

    void overwrite_rendered_obj(shared_ptr<Particles> complex_particles, shared_ptr<object> object);


    //helper functions
    double safe_stod(string str);
    int safe_stoi(string str);


};

#endif // OBJ_HANDLER_H