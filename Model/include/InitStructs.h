// Description: This file will contain all structures needed for setting up the scenario.

#ifndef INIT_STRUCTS_H
#define INIT_STRUCTS_H


#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include "Particles.h"

using namespace std;

struct scenario {
    int scenario_id;
    string name;
    string obj_list;
    double time;
    string interaction_func;
    bool try_cache; //if true, check if the scenario has been run before and use the results 
    bool refresh_obj; //if true, rebuild the complex objects
    double dt; //time step
    bool three_d; //if false, the simulation will be in 2D
};

struct scenarios {
    vector<unique_ptr<scenario>> scenario_list;
};


struct object {
    int object_id;
    string name;
    double r; //rgb values
    double g;
    double b;
    double x; //position values
    double y;
    double z;
    double vx; //velocity values
    double vy;
    double vz;
    double m; //mass
    double rad; //radius of sphere
    double rest; //restitution parameter of sphere
    string complexity = "simple"; //if the object is complex, what is its shape, e.g simple, circle, square, etc.
    double complexity_size; //if the object is complex, what is its radius. empty for simple objects
    double complexity_n; //how many particles make up the complex object

};

struct objects {
    vector<shared_ptr<object>> object_list;
};



struct test_metrics_t {
    //this struct will contain the metrics of the test
    double fps;
    double memory;
    double cpu;
    double gpu;
    double KE;
    double PE;
    double TE;
    double mom_x;
    double mom_y;
    double mom_x_change;
    double mom_y_change;

    double TE_change;
    double TE_error;
    double relative_error;
};


struct snapshots {
    //this struct will contain all the snapshots of the particles at different times
    vector<shared_ptr<Particles>> snaps;
    vector<shared_ptr<test_metrics_t>> metrics;

    
};

#endif // INIT_STRUCTS_H

