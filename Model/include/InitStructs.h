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

    //[simulation parameters]
    high_prec time;
    high_prec dt; //time step
    high_prec beta = 0.1; //contact bias beta

    //[run options]
    bool try_cache; //if true, check if the scenario has been run before and use the results 
    bool refresh_obj; //if true, rebuild the complex objects
    string save_obj = "FALSE";

    //[physics engine parameters]
    high_prec collision_distance_tolerance = 1e-5;

    //[debug/output parameters]
    bool debug_mode = false;
    bool report_energy_per_step = false;
    bool save_scenario = false;

    //[plotter parameters]
    high_prec heat_gamma = 0.7;
    high_prec heat_cutoff_frac = 0.0;
    int plot_speed_multiplier = 20;

    //[benchmark parameters] (-1 means no benchmark available)
    high_prec benchmark_te_error_pct = -1.0;
    high_prec benchmark_sim_time_sec = -1.0;

};

struct scenarios {
    vector<unique_ptr<scenario>> scenario_list;
};


struct object {
    int object_id;
    string name;
    high_prec r; //rgb values
    high_prec g;
    high_prec b;
    high_prec x; //position values
    high_prec y;
    high_prec z;
    high_prec vx; //velocity values
    high_prec vy;
    high_prec vz;
    high_prec m; //mass
    high_prec rad; //radius of sphere
    high_prec rest; //restitution parameter of sphere
    high_prec temp; //temperature of the particle
    string complexity = "simple"; //if the object is complex, what is its shape, e.g simple, circle, square, etc.
    high_prec complexity_size; //if the object is complex, what is its radius. empty for simple objects
    high_prec complexity_n = int(1);
    high_prec omega;

};

struct objects {
    vector<shared_ptr<object>> object_list;
};



struct test_metrics_t {
    //this struct will contain the metrics of the test
    high_prec fps;
    high_prec memory;
    high_prec cpu;
    high_prec gpu;
    high_prec KE;
    high_prec PE;
    high_prec TE;
    high_prec mom_x;
    high_prec mom_y;
    high_prec mom_x_change;
    high_prec mom_y_change;
    high_prec HE;

    high_prec TE_change;
    high_prec TE_error;
    high_prec relative_error;

    //marginal errors for plotting
    high_prec margin_TE_error;
    high_prec margin_TE_error_overlap;
    high_prec margin_TE_error_collision;
    high_prec margin_TE_error_integrate;
    int overlap_iters_in_step;

    high_prec margin_TE_error_overlap_ij_transl;
    high_prec margin_TE_error_overlap_ij_corrected;
    
};


struct snapshots {
    //this struct will contain all the snapshots of the particles at different times
    vector<shared_ptr<Particles>> snaps;
    vector<shared_ptr<test_metrics_t>> metrics;

    
};

#endif // INIT_STRUCTS_H

