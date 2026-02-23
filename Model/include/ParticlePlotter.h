#ifndef PARTICLE_PLOTTER_H
#define PARTICLE_PLOTTER_H

#include <iostream>
#include <memory>
#include <string>
#include "InitStructs.h"

//namespaces
using namespace std;

class Plotter {
public:
    Plotter();
    ~Plotter();

    // Main entry point
    void plot_run(shared_ptr<scenario> scenario, shared_ptr<snapshots> snapshots);

    // GNUplot functions (now used for OFFLINE rendering)
    void init_GNU(shared_ptr<scenario> scenario);
    void plot_GNU(shared_ptr<Particles> particles, shared_ptr<test_metrics_t> metrics_t);
    void close_GNU();

    // Playback after offline render completes
    void playback_offline(const std::string& rendered_file);

    // RGB/heat
    void apply_heat_brightness_and_pack_rgb(shared_ptr<snapshots> snaps, shared_ptr<scenario> scenario);
    int intensity_to_rgb(high_prec r, high_prec g, high_prec b);

    shared_ptr<snapshots> heat_to_rgb(shared_ptr<snapshots> snapshots);
    shared_ptr<Particles> convert_intensity_to_rgb(shared_ptr<Particles> particles);

private:
    bool debug_mode_ = false;

    // Used so plot_GNU can embed the correct step label BEFORE each plot frame
    int current_step_label_ = 0;

    // Output for offline render (gif)
    std::string offline_output_path_;
};

#endif // PARTICLE_PLOTTER_H