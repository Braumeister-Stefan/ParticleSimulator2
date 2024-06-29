#ifndef INTERFACER_H
#define INTERFACER_H

#include <iostream>
#include <vector>
#include <memory>
#include "../include/InitStructs.h"

using namespace std;

class Interfacer {
public:
    // Constructor
    Interfacer();

    // Destructor
    ~Interfacer();

    // Declare the select_scenario function
    shared_ptr<scenario> select_scenario();

};

#endif // INTERFACER_H