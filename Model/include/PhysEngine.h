#ifndef PHYS_ENGINE_H
#define PHYS_ENGINE_H

#include "InitStructs.h"

#include <iostream>
#include <memory>

using namespace std;

class Engine {
public:
    // Constructor
    Engine();

    // Destructor
    ~Engine();

    // Example method
    void run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles);

private:
    // Add private member variables if needed
};

#endif
