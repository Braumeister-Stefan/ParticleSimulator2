//Standard libraries
#include <memory>
#include <iostream>
#include <fstream>
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

    high_prec dt = scenario->dt;  // Time step length
    const high_prec TE_error_threshold = 0.01;
    const high_prec mom_change_threshold = 0.01;

    for (int i = 0; i < particle_states->snaps.size(); i++) {
        high_prec kinetic_energy = 0.0;
        high_prec potential_energy = 0.0;
        high_prec momentum_x = 0.0;
        high_prec momentum_y = 0.0;
        high_prec heating_energy = 0.0;


        // Get the particles for the current snapshot
        shared_ptr<Particles> particles = particle_states->snaps[i];
        int particle_count = particles->particle_list.size();

        // Calculate kinetic energy,momentum and temperature for the snapshot
        for (int j = 0; j < particle_count; j++) {
            high_prec mass = particles->particle_list[j]->m;
            high_prec vx = particles->particle_list[j]->vx;
            high_prec vy = particles->particle_list[j]->vy;
            high_prec vz = particles->particle_list[j]->vz;
            high_prec temp = particles->particle_list[j]->temp;

            kinetic_energy += 0.5 * mass * (vx * vx + vy * vy + vz * vz);
            momentum_x += mass * vx;
            momentum_y += mass * vy;
            heating_energy += temp;

        }

        // Calculate potential energy for each unique pair of particles
        for (int j = 0; j < particle_count; j++) {
            for (int k = j + 1; k < particle_count; k++) {
                high_prec dx = particles->particle_list[j]->x - particles->particle_list[k]->x;
                high_prec dy = particles->particle_list[j]->y - particles->particle_list[k]->y;
                high_prec dz = particles->particle_list[j]->z - particles->particle_list[k]->z;
                high_prec distance = sqrt(dx * dx + dy * dy + dz * dz);
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
            high_prec initial_TE = particle_states->metrics[0]->TE;
            high_prec initial_mom_x = particle_states->metrics[0]->mom_x;
            high_prec initial_mom_y = particle_states->metrics[0]->mom_y;

            metrics->TE_change = (metrics->TE - initial_TE) / initial_TE;
            // metrics->mom_x_change = (initial_mom_x != 0) ? (momentum_x - initial_mom_x) / initial_mom_x : high_prec(0.0);
            // metrics->mom_y_change = (initial_mom_y != 0) ? (momentum_y - initial_mom_y) / initial_mom_y : high_prec(0.0);

            // NEW (FIXED):
            metrics->mom_x_change = (fabs(initial_mom_x) > 1e-12) ? (momentum_x - initial_mom_x) / initial_mom_x : high_prec(0.0);
            metrics->mom_y_change = (fabs(initial_mom_y) > 1e-12) ? (momentum_y - initial_mom_y) / initial_mom_y : high_prec(0.0);

            //calculate the current TE versus the TE of the previous step. Call it margin_TE_error
            //high_prec prev_TE = particle_states->metrics[i - 1]->TE;
            //high_prec margin_TE_error = metrics->TE - prev_TE;
            //metrics->margin_TE_error = margin_TE_error;


        }

        // Cumulative TE error over time
        if (i > 0) {
            high_prec prev_TE_error = particle_states->metrics[i - 1]->TE_error;
            high_prec te_error = metrics->TE - particle_states->metrics[i - 1]->TE;
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
    int fps_smoothing_window = static_cast<int>(30 / dt);  // Number of frames to consider for smoothing
    for (int i = 0; i < particle_states->metrics.size(); i++) {
        high_prec avg_fps = 0;
        int count = 0;
        
        for (int j = i; j >= 0 && j > i - fps_smoothing_window; j--) {
            avg_fps += particle_states->metrics[j]->fps;
            count++;
        }
        // Update the FPS with the averaged value for smoother data
        particle_states->metrics[i]->fps = avg_fps / count;

        //print the value of margin_TE_error for step i

        //cout << "Margin TE error for step " << i << " : " << particle_states->metrics[i]->margin_TE_error << endl;
    }

    // Placeholder: If metrics is a null pointer, create a new metrics object and fill with zeros
    if (particle_states->metrics.size() == 0) {
        for (int i = 0; i < particle_states->snaps.size(); i++) {
            particle_states->metrics.push_back(make_shared<test_metrics_t>());
            particle_states->metrics[i]->fps = 0;
            particle_states->metrics[i]->memory = 0;
            particle_states->metrics[i]->cpu = 0;
            particle_states->metrics[i]->gpu = 0;
            particle_states->metrics[i]->margin_TE_error = 0;
            particle_states->metrics[i]->margin_TE_error_overlap = 0;
            particle_states->metrics[i]->margin_TE_error_collision = 0;
            particle_states->metrics[i]->margin_TE_error_integrate = 0;
            particle_states->metrics[i]->overlap_iters_in_step = 0;
            particle_states->metrics[i]->margin_TE_error_overlap_ij_transl = 0;
            particle_states->metrics[i]->margin_TE_error_overlap_ij_corrected = 0;
        }
    }

    cout << "Metrics computed." << endl << endl;

    metrics_to_cache(scenario, particle_states);

    return particle_states;
}

void Metrics::metrics_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {
    // 1. Define file name
    std::string file_name = "Inputs/rendered_scenarios/" + scenario->name + "_metrics.csv";

    // 2. Open the file
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cout << "Failed to write metrics to cache!" << std::endl;
        return;
    }

    // 3. Write headers (one row per step)
    file << "step_id,KE,PE,TE,HE,mom_x,mom_y,TE_change,mom_x_change,mom_y_change,TE_error,fps,margin_TE_error,margin_TE_error_overlap,margin_TE_error_collision,margin_TE_error_integrate, overlap_iters_in_step,margin_TE_error_overlap_ij_transl,margin_TE_error_overlap_ij_corrected\n";

    // 4. Loop over steps only (NOT over particles).
    //    metrics[i] is presumably already the aggregated or averaged data for step i.
    for (int i = 0; i < static_cast<int>(particle_states->metrics.size()); i++) 
    {
        // Because metrics[i] is the step-level data, we just write one line per step.
        file << std::fixed << std::setprecision(15)
             << i << ","
             << particle_states->metrics[i]->KE << ","
             << particle_states->metrics[i]->PE << ","
             << particle_states->metrics[i]->TE << ","
             << particle_states->metrics[i]->HE << ","
             << particle_states->metrics[i]->mom_x << ","
             << particle_states->metrics[i]->mom_y << ","
             << particle_states->metrics[i]->TE_change << ","
             << particle_states->metrics[i]->mom_x_change << ","
             << particle_states->metrics[i]->mom_y_change << ","
             << particle_states->metrics[i]->relative_error << ","
             << particle_states->metrics[i]->fps << ","
             << particle_states->metrics[i]->margin_TE_error<<  ","
            << particle_states->metrics[i]->margin_TE_error_overlap << ","
            << particle_states->metrics[i]->margin_TE_error_collision << ","
            << particle_states->metrics[i]->margin_TE_error_integrate << ","
            << particle_states->metrics[i]->overlap_iters_in_step << ","
            << particle_states->metrics[i]->margin_TE_error_overlap_ij_transl << ","
            << particle_states->metrics[i]->margin_TE_error_overlap_ij_corrected << endl;
            
            

    }

    // 5. Close the file
    std::cout << "Saving simulation metrics..." << std::endl;
    file.close();
    std::cout << "Metrics saved to " << file_name << std::endl << std::endl;
}

    


// void Engine::run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {
//     //this function will save the snapshots to a csv file

//     //1. define the file name, located on Inputs\rendered_scenarios

//     string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

//     //2. open the file

//     ofstream file(file_name);

//     //check if the file is open

//     if (!file.is_open()) {
//         cout << "Failed to write snapshots to cache!" << endl;
//         return;
//     }

//     //3. write headers

//     file << "step_id, particle_id,r,g,b,x,y,z,vx,vy,vz,m,rad,temp,rest" << endl;

//     //4. loop through the snapshots and write the data to the file, adding the step_id

//     for (int i = 0; i < particle_states->snaps.size(); i++) {
//         for (int j = 0; j < particle_states->snaps[i]->particle_list.size(); j++) {
//             file << std::fixed << std::setprecision(15) // Set precision to 15 decimal places
//                 << i << ","
//                 << particle_states->snaps[i]->particle_list[j]->particle_id << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->r) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->g) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->b) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->x) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->y) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->z) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->vx) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->vy) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->vz) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->m) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->rad) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->temp) << ","
//                 << static_cast<high_prec>(particle_states->snaps[i]->particle_list[j]->rest) << endl;
//         }
//     }

//     //5. save and close the file

//     cout << "Saving simulation..." << endl;

//     file.close();

//     cout << "Snapshots saved to " << file_name << endl << endl;
    

// }