
#ifndef SCENARIO_H
#define SCENARIO_H

#include <string>
#include <memory>
#include "Particles.h"  // Assuming Particles class is used in scenarios

using namespace std;

class scenario {
public:
    string name; // Name of the scenario (used for cache filenames)

    // Constructor
    scenario(const string& scenario_name = "") : name(scenario_name) {}

    // Destructor
    ~scenario() {}
};

#endif
