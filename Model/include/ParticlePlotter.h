#ifndef PARTICLE_PLOTTER_H
#define PARTICLE_PLOTTER_H

#include <iostream>
#include <memory>
#include "InitStructs.h"

//namespaces
using namespace std;


class Plotter {
public:
    // Constructor
    Plotter();

    // Destructor
    ~Plotter();

    // Define the plot_run method
    void plot_run(shared_ptr<scenario> scenario, shared_ptr<snapshots> snapshots, shared_ptr<test_metrics> metrics);

    //Define the GNUplot functions
    void init_GNU(shared_ptr<scenario> scenario);
    void plot_GNU(shared_ptr<Particles> particles, shared_ptr<test_metrics_t> metrics_t);
    void close_GNU();




    //convert r,g,b to hex for all particles
    //functions to convert rgb values to hex code
    int intensity_to_rgb(double r, double g, double b);

    shared_ptr<Particles> convert_intensity_to_rgb(shared_ptr<Particles> particles);

private:
    // Add private member variables if needed
};

#endif // PARTICLE_PLOTTER_H