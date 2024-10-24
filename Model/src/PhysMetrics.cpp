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

    double total_TE_error = 0.0;  // Variable to store aggregated TE_error

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

        // Store the potential energy in the metrics object
        particle_states->metrics[i]->PE = potential_energy;

        // Calculate total energy (TE) as the sum of kinetic and potential energy
        double total_energy = kinetic_energy + potential_energy;
        particle_states->metrics[i]->TE = total_energy;

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
            //total_TE_error += fabs(te_error);  // Aggregate absolute value of TE_error
            //cout << "Total aggregated TE_error: " << particle_states->metrics[i]->TE_error << endl;

            // Calculate relative error as the ratio of TE_error to the initial TE
            particle_states->metrics[i]->relative_error = particle_states->metrics[i]->TE_error / particle_states->metrics[0]->TE;


        } else {
            particle_states->metrics[i]->TE_error = 0;  // No previous timestep for the first snapshot

            particle_states->metrics[i]->relative_error = 0;
        }

        //cout << "TE_error (relative): " << particle_states->metrics[i]->relative_error << endl;
    }

    double relative_final_error = particle_states->metrics[particle_states->snaps.size() - 1]->relative_error;


    if (abs(relative_final_error) > TE_error_threshold) {
        cout << "Total Energy Error (" << relative_final_error*100 << "%) exceeds threshold of " << TE_error_threshold*100 <<"%"<< endl;
    } else {
        cout << "Total Energy Error (" << relative_final_error*100 << "%) within threshold of " << TE_error_threshold *100 <<"%"<< endl;
    }
    

    

    // Average the fps metrics every 10 snapshots to smooth the curve
    for (int i = 0; i < particle_states->metrics.size(); i++) {
        if (i % 10 == 0) {
            double avg_fps = 0;
            for (int j = 0; j < 10 && (i + j) < particle_states->metrics.size(); j++) {
                avg_fps += particle_states->metrics[i + j]->fps;
            }
            avg_fps /= 10;
            particle_states->metrics[i]->fps = avg_fps;
        }
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