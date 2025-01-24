//Standard libraries
#include <memory>
#include <iostream>
#include <cmath>

//Internal libraries
#include "../include/PhysMetrics.h"
#include "../include/InitStructs.h"

//namespaces
using namespace std;

// Constructor
Metrics::Metrics() {
    cout << "Metrics computor initialized." << endl;
}

// Destructor
Metrics::~Metrics() {
    cout << "Metrics computor destroyed." << endl;
}

//define the compute_metrics function

shared_ptr<snapshots> Metrics::compute_metrics(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {
    cout << "Computing metrics for " << particle_states->snaps.size() << " snapshots." << endl;

    // Initialize metrics if empty
    if (particle_states->metrics.empty()) {
        particle_states->metrics.resize(particle_states->snaps.size(), make_shared<test_metrics_t>());
    }

    double dt = scenario->dt;  // Time step length
    const double TE_error_threshold = 0.01;
    const double mom_change_threshold = 0.01;

    for (int i = 0; i < particle_states->snaps.size(); i++) {
        double kinetic_energy = 0.0;
        double potential_energy = 0.0;
        double momentum_x = 0.0;
        double momentum_y = 0.0;
        double heating_energy = 0.0;


        // Get the particles for the current snapshot
        shared_ptr<Particles> particles = particle_states->snaps[i];
        int particle_count = particles->particle_list.size();

        // Calculate kinetic energy,momentum and temperature for the snapshot
        for (int j = 0; j < particle_count; j++) {
            double mass = particles->particle_list[j]->m;
            double vx = particles->particle_list[j]->vx;
            double vy = particles->particle_list[j]->vy;
            double vz = particles->particle_list[j]->vz;
            double temp = particles->particle_list[j]->temp;

            kinetic_energy += 0.5 * mass * (vx * vx + vy * vy + vz * vz);
            momentum_x += mass * vx;
            momentum_y += mass * vy;
            heating_energy += temp;

        }

        // Calculate potential energy for each unique pair of particles
        for (int j = 0; j < particle_count; j++) {
            for (int k = j + 1; k < particle_count; k++) {
                double dx = particles->particle_list[j]->x - particles->particle_list[k]->x;
                double dy = particles->particle_list[j]->y - particles->particle_list[k]->y;
                double dz = particles->particle_list[j]->z - particles->particle_list[k]->z;
                double distance = sqrt(dx * dx + dy * dy + dz * dz);
                potential_energy += -6.674e-11 * particles->particle_list[j]->m * particles->particle_list[k]->m / distance;
            }
        }

        // Store calculated values in metrics
        shared_ptr<test_metrics_t> metrics = particle_states->metrics[i];
        metrics->KE = kinetic_energy;
        metrics->PE = potential_energy;
        metrics->TE = kinetic_energy + potential_energy + heating_energy;
        metrics->mom_x = momentum_x;
        metrics->mom_y = momentum_y;
        metrics->HE = heating_energy;

        // Calculate relative changes from the initial snapshot
        if (i == 0) {
            metrics->TE_change = 0.0;
            metrics->mom_x_change = 0.0;
            metrics->mom_y_change = 0.0;
        } else {
            double initial_TE = particle_states->metrics[0]->TE;
            double initial_mom_x = particle_states->metrics[0]->mom_x;
            double initial_mom_y = particle_states->metrics[0]->mom_y;

            metrics->TE_change = (metrics->TE - initial_TE) / initial_TE;
            metrics->mom_x_change = (initial_mom_x != 0) ? (momentum_x - initial_mom_x) / initial_mom_x : 0.0;
            metrics->mom_y_change = (initial_mom_y != 0) ? (momentum_y - initial_mom_y) / initial_mom_y : 0.0;
        }

        // Cumulative TE error over time
        if (i > 0) {
            double prev_TE_error = particle_states->metrics[i - 1]->TE_error;
            double te_error = metrics->TE - particle_states->metrics[i - 1]->TE;
            metrics->TE_error = prev_TE_error + te_error;
            metrics->relative_error = metrics->TE_error / particle_states->metrics[0]->TE;
        } else {
            metrics->TE_error = 0.0;
            metrics->relative_error = 0.0;
        }
    }

    // Check thresholds for final snapshot
    shared_ptr<test_metrics_t> final_metrics = particle_states->metrics.back();
    if (abs(final_metrics->relative_error) > TE_error_threshold) {
        cout << "Total Energy Error (" << final_metrics->relative_error * 100 << "%) exceeds threshold of "
                  << TE_error_threshold * 100 << "%" << endl;
    } else {
        cout << "Total Energy Error (" << final_metrics->relative_error * 100 << "%) within threshold of "
                  << TE_error_threshold * 100 << "%" << endl;
    }

    if (abs(final_metrics->mom_x_change) > mom_change_threshold) {
        cout << "Momentum change (x) (" << final_metrics->mom_x_change * 100 << "%) exceeds threshold of "
                  << mom_change_threshold * 100 << "%" << endl;
    } else {
        cout << "Momentum change (x) (" << final_metrics->mom_x_change * 100 << "%) within threshold of "
                  << mom_change_threshold * 100 << "%" << endl;
    }

    if (abs(final_metrics->mom_y_change) > mom_change_threshold) {
        cout << "Momentum change (y) (" << final_metrics->mom_y_change * 100 << "%) exceeds threshold of "
                  << mom_change_threshold * 100 << "%" << endl;
    } else {
        cout << "Momentum change (y) (" << final_metrics->mom_y_change * 100 << "%) within threshold of "
                  << mom_change_threshold * 100 << "%" << endl;
    }



    

    // Calculate the rolling average by considering the last `fps_smoothing_window` frames
    int fps_smoothing_window = 30/ dt;  // Number of frames to consider for smoothing
    for (int i = 0; i < particle_states->metrics.size(); i++) {
        double avg_fps = 0;
        int count = 0;
        
        for (int j = i; j >= 0 && j > i - fps_smoothing_window; j--) {
            avg_fps += particle_states->metrics[j]->fps;
            count++;
        }
        // Update the FPS with the averaged value for smoother data
        particle_states->metrics[i]->fps = avg_fps / count;
    }

    // Placeholder: If metrics is a null pointer, create a new metrics object and fill with zeros
    if (particle_states->metrics.size() == 0) {
        for (int i = 0; i < particle_states->snaps.size(); i++) {
            particle_states->metrics.push_back(make_shared<test_metrics_t>());
            particle_states->metrics[i]->fps = 0;
            particle_states->metrics[i]->memory = 0;
            particle_states->metrics[i]->cpu = 0;
            particle_states->metrics[i]->gpu = 0;
        }
    }

    cout << "Metrics computed." << endl << endl;
    return particle_states;
}