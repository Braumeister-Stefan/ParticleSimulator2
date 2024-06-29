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

    // Declare the flatten complex objects function
    shared_ptr<Particles> flatten_complex_objs(shared_ptr<objects> requested_objects);

    //helper functions
    double safe_stod(string str);


};

#endif // OBJ_HANDLER_H