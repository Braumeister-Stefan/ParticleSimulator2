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
    // This function computes the metrics for the simulation, used for validation purposes.
    cout << "Computing metrics for " << particle_states->snaps.size() << " snapshots." << endl;

    //check if metrics object is empty, if so, create a new one

    if (particle_states->metrics.size() == 0) {
        for (int i = 0; i < particle_states->snaps.size(); i++) {
            particle_states->metrics.push_back(make_shared<test_metrics_t>());
        }
    }

    double dt = scenario->dt;  // Time step length

    double TE_error_threshold = 0.01;  // Threshold for total energy error
    double mom_change_threshold = 0.01;

    

    // Loop through the snapshots and compute the metrics
    for (int i = 0; i < particle_states->snaps.size(); i++) {
        // Compute the kinetic energy of the system in a time step
        double kinetic_energy = 0;
        for (int j = 0; j < particle_states->snaps[i]->particle_list.size(); j++) {
            kinetic_energy += 0.5 * particle_states->snaps[i]->particle_list[j]->m * 
                              (pow(particle_states->snaps[i]->particle_list[j]->vx, 2) +
                               pow(particle_states->snaps[i]->particle_list[j]->vy, 2) +
                               pow(particle_states->snaps[i]->particle_list[j]->vz, 2));
        }

        // Store the kinetic energy in the metrics object
        particle_states->metrics[i]->KE = kinetic_energy;

        // Compute the potential energy of the system in a time step
        double potential_energy = 0;
        for (int j = 0; j < particle_states->snaps[i]->particle_list.size(); j++) {
            for (int k = j + 1; k < particle_states->snaps[i]->particle_list.size(); k++) {
                double distance = sqrt(pow(particle_states->snaps[i]->particle_list[j]->x - particle_states->snaps[i]->particle_list[k]->x, 2) + 
                                       pow(particle_states->snaps[i]->particle_list[j]->y - particle_states->snaps[i]->particle_list[k]->y, 2) + 
                                       pow(particle_states->snaps[i]->particle_list[j]->z - particle_states->snaps[i]->particle_list[k]->z, 2));
                double pe = -6.674 * pow(10, -11) * particle_states->snaps[i]->particle_list[j]->m * 
                            particle_states->snaps[i]->particle_list[k]->m / distance;
                potential_energy += pe;
            }
        }

        //Compute the momentum of the system
        double momentum = 0;

        for (int j = 0; j < particle_states->snaps[i]->particle_list.size(); j++) {
            momentum += particle_states->snaps[i]->particle_list[j]->m * sqrt(pow(particle_states->snaps[i]->particle_list[j]->vx, 2) + pow(particle_states->snaps[i]->particle_list[j]->vy, 2) + pow(particle_states->snaps[i]->particle_list[j]->vz, 2));
        }

        // Store the momentum in the metrics object
        particle_states->metrics[i]->mom = momentum;

        // Calculate the change in momentum relative to the first snapshot
        if (i == 0) {
            particle_states->metrics[i]->mom_change = 0;
        } else {
            double initial_momentum = particle_states->metrics[0]->mom;
        if (initial_momentum != 0) {
                particle_states->metrics[i]->mom_change = (momentum - initial_momentum) / initial_momentum;
            } else {
                particle_states->metrics[i]->mom_change = 0;  // Handle division by zero by setting to 0
                //cout << "Warning: Initial momentum is zero. Momentum change set to 0 to avoid division by zero." << endl;
            }
            //cout << "Momentum change: " << particle_states->metrics[i]->mom_change << endl;
        }
        

        // Store the potential energy in the metrics object
        particle_states->metrics[i]->PE = potential_energy;

        // Calculate total energy (TE) as the sum of kinetic and potential energy
        double total_energy = kinetic_energy + potential_energy;
        particle_states->metrics[i]->TE = total_energy;

        double total_TE_error = 0.0;  // Variable to store aggregated TE_error
        // Calculate TE_change relative to the first snapshot
        if (i == 0) {
            particle_states->metrics[i]->TE_change = 0;
        } else {
            particle_states->metrics[i]->TE_change = 
                (particle_states->metrics[i]->TE - particle_states->metrics[0]->TE) / particle_states->metrics[0]->TE;
        }

        // Calculate TE_error as the difference between the current and previous timestep
        if (i > 0) {
            double te_error = particle_states->metrics[i]->TE - particle_states->metrics[i - 1]->TE;
            particle_states->metrics[i]->TE_error = particle_states->metrics[i-1]->TE_error + fabs(te_error);


            // Calculate relative error as the ratio of TE_error to the initial TE
            particle_states->metrics[i]->relative_error = particle_states->metrics[i]->TE_error / particle_states->metrics[0]->TE;


        } else {
            particle_states->metrics[i]->TE_error = 0;  // No previous timestep for the first snapshot

            particle_states->metrics[i]->relative_error = 0;
        }

        
    }

    double relative_final_error = particle_states->metrics[particle_states->snaps.size() - 1]->relative_error;


    if (abs(relative_final_error) > TE_error_threshold) {
        cout << "Total Energy Error (" << relative_final_error*100 << "%) exceeds threshold of " << TE_error_threshold*100 <<"%"<< endl;
    } else {
        cout << "Total Energy Error (" << relative_final_error*100 << "%) within threshold of " << TE_error_threshold *100 <<"%"<< endl;
    }

    //check if relative momentum change is within threshold
    double relative_mom_change = particle_states->metrics[particle_states->snaps.size() - 1]->mom_change;


    if (abs(relative_mom_change) > mom_change_threshold) {
        cout << "Momentum change (" << relative_mom_change*100 << "%) exceeds threshold of " << mom_change_threshold*100 <<"%"<< endl;
    } else {
        cout << "Momentum change (" << relative_mom_change*100 << "%) within threshold of " << mom_change_threshold *100 <<"%"<< endl;
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