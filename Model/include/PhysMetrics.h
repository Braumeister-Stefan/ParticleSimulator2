#ifndef PHYS_METRICS_H
#define PHYS_METRICS_H

#include "InitStructs.h"

#include <iostream>
#include <memory>

using namespace std;

class Metrics {
public:
    // Constructor
    Metrics();

    // Destructor
    ~Metrics();

    //Function to compute (validation) metrics of the simulation
    shared_ptr<test_metrics> compute_metrics(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states);

    

private:
    // Add private member variables if needed
};

#endif
