//Standard libraries
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <chrono>
#include <iomanip>

//External libraries
#define CSV_IO_NO_THREAD
#include "../include/3party/csv.h"
#include <boost/multiprecision/cpp_dec_float.hpp>

//Internal libraries
#include "../include/PhysEngine.h"
#include "../include/InitStructs.h"
#include "../include/MathUtils.h"


//namespaces
using namespace std;
using namespace std::chrono;
using namespace boost::multiprecision;
using high_prec = cpp_dec_float_50;

//use momentum from the engine struct
using momentum = Engine::momentum;

//constants
//Gravitational constant
const high_prec G = 6.674 * pow(10, -11); //m^3 kg^-1 s^-2

//variables storing total error of overlap, collision and verlet
high_prec total_TE_error_overlap = 0;
high_prec total_TE_error_collision = 0;
high_prec total_TE_error_verlet = 0;

//variable storing counter of energy gap corrections
int energy_gap_corrections = 0;
int energy_gap_corrections_incomplete = 0;

//variables storing the index of overlap, overlap_ij
int overlap_iter = 0;
int overlap_ij_iter = 0;
int update_iter = 0;

//variables storing marginal error metrics
static high_prec s_prev_TE = -1;
static high_prec s_prev_step = -1;
static high_prec s_margin_TE_error = 0.0;

static high_prec s_margin_TE_error_overlap  = 0.0;
static high_prec s_margin_TE_error_collision = 0.0;
static high_prec s_margin_TE_error_integrate = 0.0;

static high_prec s_margin_TE_error_overlap_ij_transl = 0.0;
static high_prec s_margin_TE_error_overlap_ij_corrected = 0.0;

static high_prec collissions = 0; //counter for collissions






//define static member variables
high_prec Engine::dt = 1; //time step

// Constructor
Engine::Engine() {
    cout << "Engine initialized." << endl;
}

// Destructor
Engine::~Engine() {
    cout << "Engine destroyed." << endl;
}

high_prec Engine::get_margin_TE_error() {
    return s_margin_TE_error;
}

high_prec Engine::get_margin_TE_error_overlap()   { return s_margin_TE_error_overlap; }
high_prec Engine::get_margin_TE_error_collision() { return s_margin_TE_error_collision; }
high_prec Engine::get_margin_TE_error_integrate() { return s_margin_TE_error_integrate; }
high_prec Engine::get_margin_TE_error_overlap_ij_transl() { return s_margin_TE_error_overlap_ij_transl; }
high_prec Engine::get_margin_TE_error_overlap_ij_corrected() { return s_margin_TE_error_overlap_ij_corrected; }

// Run the simulation and return snapshots of each time step
shared_ptr<snapshots> Engine::run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles)
{
    cout << "Engine is initialized." << endl;

    // 1. Initialize the snapshots object
    shared_ptr<snapshots> particle_states = make_shared<snapshots>();

    cout << scenario->name << "'s initial states loaded" << endl << endl;
    cout << "Press enter to start the simulation." << endl;
    cin.ignore();
    cin.get();

    // 3. Initialize time step
    Engine::dt = scenario->dt;
    high_prec steps_db = scenario->time / scenario->dt;
    int steps = static_cast<int>(steps_db);
    cout << "Number of steps to be simulated: " << steps << endl << endl;

    // Ensure metrics vector is sized
    particle_states->metrics.resize(steps, nullptr);

    // 4. Main loop
    for (int i = 0; i < steps; i++) {

        if (!particle_states->metrics[i]) {
            particle_states->metrics[i] = make_shared<test_metrics_t>();
        }

        if (i >= 6401) {
            // Potential debug pause
        }

        auto start_time = high_resolution_clock::now();
        update_iter = i;

        high_prec te_pre_update = calc_TE(particles);
        momentum mom_pre_update = calc_mom(particles);

        // Call update_particles
        //cout << "Updating particles at step " << i + 1 << " of " << steps << endl;
        update_particles(particles);

        high_prec te_post_update = calc_TE(particles);
        momentum mom_post_update = calc_mom(particles);

        // Copy final state for this step
        auto particles_copy = make_unique<Particles>(*particles);
        particle_states->snaps.push_back(move(particles_copy));

        // Store the overall margin error
        high_prec margin_TE_error = get_margin_TE_error();
        particle_states->metrics[i]->margin_TE_error = margin_TE_error;

        // Also store the three substep marginal errors
        high_prec margin_overlap   = get_margin_TE_error_overlap();
        high_prec margin_collision = get_margin_TE_error_collision();
        high_prec margin_integrate = get_margin_TE_error_integrate();
        particle_states->metrics[i]->margin_TE_error_overlap   = margin_overlap;
        particle_states->metrics[i]->margin_TE_error_collision = margin_collision;
        particle_states->metrics[i]->margin_TE_error_integrate = margin_integrate;
        particle_states->metrics[i]->overlap_iters_in_step = overlap_iter;
        particle_states->metrics[i]->margin_TE_error_overlap_ij_transl = get_margin_TE_error_overlap_ij_transl();
        particle_states->metrics[i]->margin_TE_error_overlap_ij_corrected = get_margin_TE_error_overlap_ij_corrected();


        high_prec te_error_update = (te_post_update - te_pre_update) / te_pre_update;
        if (te_error_update > 0.05) {
            cout << "TE error in run(): " << te_error_update << endl;
        }

        if (i % (steps / 20) == 0) {
            cout << i / (steps / 100) << "% of the simulation complete." << endl;
        }

        auto end_time = high_resolution_clock::now();
        duration<double> time_taken = end_time - start_time;
        particle_states->metrics[i]->fps = 1 / time_taken.count();
    }

    cout << scenario->name << " simulation completed." << endl << endl;
    cout << "Total TE error (sum of parts): "
         << total_TE_error_overlap + total_TE_error_collision + total_TE_error_verlet << endl;
    cout << "% overlap error: "
         << total_TE_error_overlap
            / (total_TE_error_overlap + total_TE_error_collision + total_TE_error_verlet)
            * 100 << endl;
    cout << "% collision error: "
         << total_TE_error_collision
            / (total_TE_error_overlap + total_TE_error_collision + total_TE_error_verlet)
            * 100 << endl;
    cout << "% verlet error: "
         << total_TE_error_verlet
            / (total_TE_error_overlap + total_TE_error_collision + total_TE_error_verlet)
            * 100 << endl;

    cout << "Total energy gap corrections: " << energy_gap_corrections << endl;
    cout << "Total energy gap corrections incomplete: " << energy_gap_corrections_incomplete << endl;

    cout << "Total collissions: " << collissions << endl;

    run_to_cache(scenario, particle_states);
    return particle_states;
}




void Engine::update_particles(shared_ptr<Particles> particles)
{   
    // Overlap resolution
    bool no_overlap = false;
    high_prec te_pre_overlap = calc_TE(particles);

    for (int i = 0; i < 8; i++) {
        if (no_overlap) break;
        no_overlap = true;
        overlap_iter = i;
        no_overlap = resolve_overlaps(particles);
    }
    high_prec te_post_overlap = calc_TE(particles);
    high_prec te_error_overlap = (te_post_overlap - te_pre_overlap) / te_pre_overlap;
    if (te_error_overlap > 0.05) {
        cout << "TE error overlap: " << te_error_overlap << endl;
    }
    // Store substep margin for overlap
    s_margin_TE_error_overlap = te_post_overlap - te_pre_overlap;

    // Collisions
    high_prec te_pre_collision = calc_TE(particles);
    resolve_collisions(particles);
    high_prec te_post_collision = calc_TE(particles);
    high_prec te_error_collision = (te_post_collision - te_pre_collision) / te_pre_collision;
    if (te_error_collision > 0.05) {
        cout << "TE error collision: " << te_error_collision << endl;
    }
    // Store substep margin for collision
    s_margin_TE_error_collision = te_post_collision - te_pre_collision;

    // Gravity (Velocity Verlet)
    high_prec te_pre_verlet = calc_TE(particles);
    resolve_gravity_verlet(particles);
    high_prec te_post_verlet = calc_TE(particles);
    high_prec te_error_verlet = (te_post_verlet - te_pre_verlet) / te_pre_verlet;
    if (te_error_verlet > 0.05) {
        cout << "TE error verlet: " << te_error_verlet << endl;
    }
    // Store substep margin for integration
    s_margin_TE_error_integrate = te_post_verlet - te_pre_verlet;

    // Now final TE for overall margin
    high_prec current_TE_post = calc_TE(particles);
    if (s_prev_step >= 0 && update_iter == s_prev_step + 1) {
        s_margin_TE_error = current_TE_post - s_prev_TE;
    } else {
        s_margin_TE_error = 0.0;
    }
    s_prev_TE   = current_TE_post;
    s_prev_step = update_iter;
}

bool Engine::resolve_overlaps(shared_ptr<Particles> particles) {

    //cout << "Resolving overlap between particles triggered" << endl;
    bool no_overlap = true;
    high_prec tTE_pre = calc_TE(particles);  // Accumulates total energy error for all particle pairs
    //momentum pre_mom = calc_mom(particles);  // Accumulates total momentum error for all particle pairs

    high_prec error_threshold = 0.001;  // Energy and momentum error threshold

    int pair_count = 0;  // To track the number of overlapping particle pairs

    //reset marginal error metrics
    s_margin_TE_error_overlap_ij_transl = 0.0;
    s_margin_TE_error_overlap_ij_corrected = 0.0;

    // Loop through particle pairs
    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            if (check_overlap(particles->particle_list[i], particles->particle_list[j])) {

                //no_overlap = false;

                pair_count++;

                // Call the core function to resolve overlap between particles i and j

                high_prec TE_pre_overlap = calc_TE(particles);
                
                resolve_overlap_ij(particles->particle_list[i], particles->particle_list[j]);

                high_prec TE_post_overlap = calc_TE(particles);
                high_prec TE_error_overlap = (TE_post_overlap - TE_pre_overlap) / TE_pre_overlap;

                //cout << "------------------------------------------" << endl;
                //cout << "For particle pair: " << particles->particle_list[i]->particle_id << " and " << particles->particle_list[j]->particle_id << endl;
                //cout << "TE error overlap: " << TE_post_overlap - TE_pre_overlap << endl;
                //cout << "TE error overlap: " << TE_error_overlap*100 << "%" << endl;
                //cout << "------------------------------------------" << endl;


            }
        }
    }

    // Calculate the final total energy error and print it
    high_prec tTE_post = calc_TE(particles);
    //momentum tmom_post = calc_mom(particles);

    high_prec tTE_error = (tTE_post - tTE_pre) / tTE_pre;
    //momentum tmom_error = (tmom_post - pre_mom) / pre_mom;

    
    //cout << "Total TE error after overlap resolution - in resolve_overlap(): " << tTE_error*100 << "%" << endl;

    
    

    //check if momentum is conserved, taking into account that there is an x and y component
    //if (abs(tmom_error.x) > error_threshold || abs(tmom_error.y) > error_threshold) {
        //cout << "Momentum error after overlap resolution - in resolve_overlap(): " << tmom_error.x << " " << tmom_error.y << endl;
    //}

    //cout << "returning no_overlap: " << no_overlap << endl;

    return no_overlap;
}

void Engine::resolve_overlap_ij(shared_ptr<Particle> particle_i, shared_ptr<Particle> particle_j) 
{
    //cout << "Resolving overlap between particles_ij triggered" << endl;

    // Increment the overlap_ij_iter
    overlap_ij_iter++;

    // Define thresholds
    high_prec error_threshold = 0.0001;

    // Store pre-resolution positions and velocities
    high_prec x1_before = particle_i->x;
    high_prec y1_before = particle_i->y;
    high_prec x2_before = particle_j->x;
    high_prec y2_before = particle_j->y;
    high_prec vx1_before = particle_i->vx;
    high_prec vy1_before = particle_i->vy;
    high_prec vx2_before = particle_j->vx;
    high_prec vy2_before = particle_j->vy;

    // Store distance, overlap, and relative velocity pre-collision
    high_prec distance_pre = hypot(x2_before - x1_before, y2_before - y1_before);
    high_prec overlap_pre = (particle_i->rad + particle_j->rad) - distance_pre;
    high_prec rel_vel_pre = (vx2_before - vx1_before) * (x2_before - x1_before) +
                            (vy2_before - vy1_before) * (y2_before - y1_before);

    //the above gives the dot product of the relative velocity and the distance vector, which is the relative velocity in the direction of the distance vector
    //instead the magnitude of rel_vel can be calculated as:
    high_prec rel_vel_pre_mag = hypot(vx2_before - vx1_before, vy2_before - vy1_before) * distance_pre / hypot(x2_before - x1_before, y2_before - y1_before);
    //cout << "---" << endl;
    //cout << "pre resolve overlap: " << endl;
    
    //cout << "rel_vel_mag: " << rel_vel_pre_mag << endl;

    // Calculate initial total kinetic energy using calc_TE_ij
    high_prec TE_pre_ij = calc_TE_ij(particle_i, particle_j, true);
    

    // Calculate initial kinetic energies in the normal and tangential directions
    high_prec normal_x = (x2_before - x1_before) / distance_pre;
    high_prec normal_y = (y2_before - y1_before) / distance_pre;
    high_prec tangential_x = -normal_y;
    high_prec tangential_y = normal_x;

    high_prec v1n = vx1_before * normal_x + vy1_before * normal_y;
    high_prec v2n = vx2_before * normal_x + vy2_before * normal_y;
    high_prec v1t = vx1_before * tangential_x + vy1_before * tangential_y;
    high_prec v2t = vx2_before * tangential_x + vy2_before * tangential_y;

    high_prec m1 = particle_i->m;
    high_prec m2 = particle_j->m;
    high_prec reduced_mass = (m1 * m2) / (m1 + m2);

    high_prec KE_rel_n_initial = 0.5 * reduced_mass * (v2n - v1n) * (v2n - v1n);
    high_prec KE_rel_t_initial = 0.5 * reduced_mass * (v2t - v1t) * (v2t - v1t);

    //cout << "KE_rel_n_initial: " << KE_rel_n_initial << "and KE_rel_t_initial: " << KE_rel_t_initial << endl;

    // Print the initial kinetic energy information in a structured way
    // std::cout << "=== Initial Kinetic Energy Information ===" << std::endl;
    // std::cout << "Total Kinetic Energy (TE_pre_ij): " << TE_pre_ij << std::endl;
    // std::cout << "Normal Kinetic Energy: " << KE_rel_n_initial << std::endl;
    // std::cout << "Tangential Kinetic Energy: " << KE_rel_t_initial << std::endl;

    //check that normal+ tangential KE = total KE
    //cout << "NKE + TKE " << KE_rel_n_initial + KE_rel_t_initial << endl;

    momentum mom_pre_ij = calc_mom_ij(particle_i, particle_j);

    // Compute overlap
    high_prec dx = x2_before - x1_before;
    high_prec dy = y2_before - y1_before;
    high_prec distance = hypot(dx, dy);

    // Avoid division by zero
    if (distance == 0.0) {
        dx = 0.001 * (rand() / (high_prec)RAND_MAX - 0.5);
        dy = 0.001 * (rand() / (high_prec)RAND_MAX - 0.5);
        distance = hypot(dx, dy);
    }

    high_prec overlap = (particle_i->rad + particle_j->rad) - distance;
    //cout << "overlap_pre in resolve_overlap_ij(): " << overlap << endl;

    // Collision normal
    high_prec nx = dx / distance;
    high_prec ny = dy / distance;

    // Displacements
    high_prec total_mass = m1 + m2;
    high_prec d1 = (overlap * m2) / total_mass;
    high_prec d2 = (overlap * m1) / total_mass;

    // Forcibly separate
    particle_i->x = (x1_before - nx * d1);
    particle_i->y = (y1_before - ny * d1);
    particle_j->x = (x2_before + nx * d2);
    particle_j->y = (y2_before + ny * d2);

    // Re-check distances
    high_prec dx_new = particle_j->x - particle_i->x;
    high_prec dy_new = particle_j->y - particle_i->y;
    high_prec distance_new = hypot(dx_new, dy_new);

    // Overlap post separation (debug only)
    high_prec overlap_post = (particle_i->rad + particle_j->rad) - distance_new;
    //cout << "overlap_post in resolve_overlap_ij(): " << overlap_post << endl;

    high_prec rel_vel_post = (particle_j->vx - particle_i->vx) * dx_new +
                             (particle_j->vy - particle_i->vy) * dy_new;

    high_prec rel_vel_post_mag = hypot(particle_j->vx - particle_i->vx, particle_j->vy - particle_i->vy) * distance_new / hypot(dx_new, dy_new);
    //cout << "---" << endl;
    //cout << "post resolve overlap: " << endl;
    
    //cout << "rel_vel_post: " << rel_vel_post << endl;
    //cout << "rel_vel_mag: " << rel_vel_post_mag << endl;


    // Compute TE, momentum again
    high_prec TE_post_ij = calc_TE_ij(particle_i, particle_j, false); 
    //momentum mom_post_ij = calc_mom_ij(particle_i, particle_j);

    
    //output the kinetic energy in normal and tangential directions
    high_prec v1n_new = particle_i->vx * nx + particle_i->vy * ny;
    high_prec v2n_new = particle_j->vx * nx + particle_j->vy * ny;
    high_prec v1t_new = particle_i->vx * (-ny) + particle_i->vy * nx;
    high_prec v2t_new = particle_j->vx * (-ny) + particle_j->vy * nx;
    high_prec KE_rel_n_post = 0.5 * reduced_mass * (v2n_new - v1n_new) * (v2n_new - v1n_new);
    high_prec KE_rel_t_post = 0.5 * reduced_mass * (v2t_new - v1t_new) * (v2t_new - v1t_new);

    //cout << "KE_rel_n_post: " << KE_rel_n_post << " and KE_rel_t_post: " << KE_rel_t_post << endl;

    // Print the kinetic energy right before calling resolve_energy_gap
    // std::cout << "=== Kinetic Energy Before Energy Gap Resolution ===" << std::endl;
    // std::cout << "Total Kinetic Energy (TE_post_ij): " << TE_post_ij << std::endl;
    // //breakdown in normal and tangential components
    // high_prec KE_rel_n_post = 0.5 * reduced_mass * (v2n - v1n) * (v2n - v1n);
    // high_prec KE_rel_t_post = 0.5 * reduced_mass * (v2t - v1t) * (v2t - v1t);
    // std::cout << "Normal Kinetic Energy: " << KE_rel_n_post << std::endl;
    // std::cout << "Tangential Kinetic Energy: " << KE_rel_t_post << std::endl;
    // cout << "NKE + TKE " << KE_rel_n_post + KE_rel_t_post << endl;

    // Store post-collision data, positions and velocities
    high_prec distance_post = distance_new;

    
    high_prec x1_post = particle_i->x;
    high_prec y1_post = particle_i->y;
    high_prec x2_post = particle_j->x;
    high_prec y2_post = particle_j->y;
    high_prec vx1_post = particle_i->vx;
    high_prec vy1_post = particle_i->vy;
    high_prec vx2_post = particle_j->vx;
    high_prec vy2_post = particle_j->vy;


    if (TE_post_ij != TE_pre_ij) {
    high_prec delta_E = TE_post_ij - TE_pre_ij;
    cout << "---" << endl;
    //cout << "delta_E to be removed: " << fixed << setprecision(10) << delta_E << defaultfloat << endl;

    //calculate the KE_post_ij based on the post-collision velocities
    high_prec KE_post_ij = 0.5 * m1 * (vx1_post * vx1_post + vy1_post * vy1_post) +
                            0.5 * m2 * (vx2_post * vx2_post + vy2_post * vy2_post);
    // Compute target KE to restore original TE
    high_prec KE_target = KE_post_ij - delta_E;

    // Safety check: don't take sqrt of negative
    if (KE_target < 0) {
        cerr << "ERROR: Target KE would be negative. Aborting correction." << endl;
    } else {
        // Calculate correct scaling factor for velocity
        high_prec scale = sqrt(KE_target / KE_post_ij);

        //cout << "scaling factor: " << fixed << setprecision(10) << scale << defaultfloat << endl;

        // Scale velocities
        particle_i->vx = vx1_post * scale;
        particle_i->vy = vy1_post * scale;
        particle_j->vx = vx2_post * scale;
        particle_j->vy = vy2_post * scale;

        // Recompute energy and log
        //cout << "---" << endl;
        //cout << "post kinetic removal correction:" << endl;
        high_prec TE_postgap_ij = calc_TE_ij(particle_i, particle_j, false);
        //cout << "TE post correction: " << TE_postgap_ij << endl;

        // Overlap geometry
        high_prec dx_after = particle_j->x - particle_i->x;
        high_prec dy_after = particle_j->y - particle_i->y;
        high_prec distance_after = hypot(dx_after, dy_after);
        high_prec overlap_after = (particle_i->rad + particle_j->rad) - distance_after;

        //cout << "Distance after correction: " << distance_after << endl;
        //cout << "Overlap after correction: " << overlap_after << endl;
    }
}

    // //if condition1 and condition2 is true
    // if (TE_post_ij != TE_pre_ij && False){
    //     high_prec delta_E = TE_post_ij - TE_pre_ij;
    //     cout << "---"   << endl;
    //     cout << "delta_E to be corrected:" << delta_E << endl;

    //     //subtract the delta_E from the heating energy
    //     //high_prec heat_energy = heat_ij(-delta_E, particle_i, particle_j);

    //     //high_prec TE_postheat_ij = calc_TE_ij(particle_i, particle_j, true);

    //     //high_prec delta_E_postheat = TE_postheat_ij - TE_pre_ij;

    //     //cout << "delta_E_postheat: " << delta_E_postheat << endl;

    //     //instead you can subtract the delta_E from the tangential kinetic energy
    //     cout << "Subtracting delta_E from tangential kinetic energy." << endl;

    //     //1 calculate the tangential kinetic energy
    //     high_prec KE_rel_t = KE_rel_t_initial;
    //     high_prec KE_rel_n = KE_rel_n_initial;

    //     KE_rel_t -= delta_E;

    //     //change the particle velocities in the tangential direction with the new KE_rel_t

    //     //assume there is sufficient tangential kinetic energy

    //     //compute impulse and apply adjustment
    //     high_prec new_rel_vel_t = sqrt(high_prec(2 * KE_rel_t / reduced_mass));
    //     if (v2t_new - v1t_new < 0) new_rel_vel_t = -new_rel_vel_t;
    //     high_prec delta_rel_vel_t = new_rel_vel_t - (v2t_new - v1t_new);    
    //     high_prec impulse_t = delta_rel_vel_t * reduced_mass;
    //     high_prec dvx1_t = (impulse_t / m1) * (-ny);
    //     high_prec dvy1_t = (impulse_t / m1) * nx;
    //     high_prec dvx2_t = (impulse_t / m2) * (-ny);
    //     high_prec dvy2_t = (impulse_t / m2) * nx;
    //     particle_i->vx -= dvx1_t.convert_to<double>();  
    //     particle_i->vy -= dvy1_t.convert_to<double>();
    //     particle_j->vx += dvx2_t.convert_to<double>();
    //     particle_j->vy += dvy2_t.convert_to<double>();

    //     //print TE post kinetic energy fix
    //     high_prec TE_post_ij_corrected = calc_TE_ij(particle_i, particle_j, true);

       
        
        


    // }

    //calculate the relative velocity in the normal direction. recalculate based on the post-collision positions and velocities
    high_prec distance_post2 = hypot(particle_j->x - particle_i->x, particle_j->y - particle_i->y);
    high_prec nx_post = (particle_j->x - particle_i->x) / distance_post2;
    high_prec ny_post = (particle_j->y - particle_i->y) / distance_post2;

    high_prec rel_vel_n = (particle_j->vx - particle_i->vx) * nx_post + 
                          (particle_j->vy - particle_i->vy) * ny_post;


    // cout << "--------Post overlap correction breakdown----------" << endl;
    // cout << "Particle i ID: " << particle_i->particle_id << endl;
    // cout << "Particle j ID: " << particle_j->particle_id << endl;

    // cout << "nx (post overlap correction): " << nx_post << endl;
    // cout << "ny (post overlap correction): " << ny_post << endl;

    cout << "Rel_vel_normal (post overlap correction): " << rel_vel_n << endl;

}


high_prec Engine::heat_ij(high_prec E, shared_ptr<Particle> particle_i, shared_ptr<Particle> particle_j){

    //convert the given energy to heat
    high_prec heat = E * (particle_i->m + particle_j->m) / 2;

    //cout << "Heat: " << heat << endl;

    //convert the heat to temperature
    high_prec temp = heat / (particle_i->m + particle_j->m);

    particle_i->temp += temp.convert_to<double>();
    particle_j->temp += temp.convert_to<double>();

    return E;


}


void Engine::resolve_collisions(shared_ptr<Particles> particles) {
    //this function will resolve collissions between particles

  

    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            //2a. check for collission between particle i and particle j

            bool collission = check_collission(particles->particle_list[i], particles->particle_list[j]);
            
            //output distance



            if (collission) {
                //resolve the collission

                //add collission to the counter
                collissions++;

                resolve_collission(particles->particle_list[i], particles->particle_list[j]);

                
            }
            

        }
    }

    
    
    
}


bool Engine::check_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    //this function will check if two particles are colliding


    //1. calculate the distance between the particles
    high_prec distance = hypot(particle2->x - particle1->x, particle2->y - particle1->y);
    
    
    //calculate relative velocity 
    high_prec relative_velocity = hypot(particle2->vx - particle1->vx, particle2->vy - particle1->vy);
    // Calculate rel_velocity in the normal direction

    high_prec dx = particle2->x - particle1->x;
    high_prec dy = particle2->y - particle1->y;
    high_prec distance_vector = hypot(dx, dy);

    //collision normal 
    high_prec nx = dx / distance_vector;
    high_prec ny = dy / distance_vector;

    //project relative velocity onto the normal vector
    high_prec rel_vel_normal = (particle2->vx - particle1->vx) * nx + (particle2->vy - particle1->vy) * ny;

    //cout << "rel_vel_normal (collision check): " << rel_vel_normal << endl;



    //calculate the direction of the relative velocity. If away from each other, then the direction is positive, if towards each other, then the direction is negative
    high_prec direction = (particle1->x - particle2->x) * (particle1->vx - particle2->vx) + 
                        (particle1->y - particle2->y) * (particle1->vy - particle2->vy);

    //if (direction < 0) {
        //cout << "Relative velocity is towards each other." << endl;
    //} else if (direction > 0) {
        //cout << "Relative velocity is away from each other." << endl;
    //}


    //cout << "Relative velocity (colission check - note this is magnitude): " << fixed << setprecision(10) << relative_velocity << defaultfloat << endl;
    //cout << "Distance (collission check): " << fixed << setprecision(10) << distance << defaultfloat << endl;


    //2. check if the distance+threshold is smaller than the sum of the radii of the particles

    //allow for a threshold of 0.0001
    high_prec threshold = 0.0000000001;

    if (abs(distance - (particle1->rad + particle2->rad)) < threshold) {
        if (rel_vel_normal < 0) {


            return true; 
        }
        //cout << "Particles are colliding." << endl;
        //cout << "Distance: " << fixed << setprecision(10) << distance << defaultfloat << endl;
        //cout << "Sum of radii: " << fixed << setprecision(10) << (particle1->rad + particle2->rad) << defaultfloat << endl;
        //cout << "Relative velocity: " << fixed << setprecision(10) << relative_velocity << defaultfloat << endl;

    //if (distance == (particle1->rad + particle2->rad) && (relative_velocity > 0)) {
        //cout << "Particles are colliding." << endl;
        //cout << "Distance: " << fixed << setprecision(10) << distance << defaultfloat << endl;
        //cout << "Sum of radii: " << fixed << setprecision(10) << (particle1->rad + particle2->rad) << defaultfloat << endl;
        //cout << "Relative velocity at collission: " << fixed << setprecision(10) << relative_velocity << defaultfloat << endl;

        //cin.get();


      



        return false;
    } else {

        //cout << "Particles are NOT colliding." << endl;
        //cout << "Distance: " << fixed << setprecision(10) << distance << defaultfloat << endl;
        //cout << "Sum of radii: " << fixed << setprecision(10) << (particle1->rad + particle2->rad) << defaultfloat << endl;
        //out << "Relative velocity: " << fixed << setprecision(10) << relative_velocity << defaultfloat << endl;
        
        return false;
    }
}

bool Engine::check_overlap(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    //this function will check if two particles are colliding

 

    //1. calculate the distance between the particles
    high_prec distance = hypot(particle1->x - particle2->x, particle1->y - particle2->y);
    
    
    //2. check if the distance+threshold is smaller than the sum of the radii of the particles
    if ((distance)< particle1->rad + particle2->rad ) {
        //cout << "Particles are overlapping." << endl;

        return true;
    } else {
        
        return false;
    }
}


void Engine::resolve_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    // This function will resolve the collision between two particles
    //cout << "particle collission resolution started" << endl;
    

    //calculate TE before the collision
    high_prec TE_pre = calc_TE_ij(particle1, particle2);



    //calculate the relative velocity, adjusted for mass
    //high_prec relative_velocity = sqrt(pow(particle1->vx - particle2->vx, 2) + pow(particle1->vy - particle2->vy, 2)) * sqrt(particle1->m / particle2->m);
    //cout << "Relative velocity: " << relative_velocity << endl;

    


    // Calculate the normal vector
    Vector2D normal = { high_prec(particle2->x) - high_prec(particle1->x), high_prec(particle2->y) - high_prec(particle1->y) };
    high_prec distance = sqrt(normal.x * normal.x + normal.y * normal.y);
    normal.x /= distance;
    normal.y /= distance;

    // Calculate the tangential vector
    Vector2D tangent = { -normal.y, normal.x };

    // Decompose velocities into normal and tangential components
    high_prec v1n = normal.x * high_prec(particle1->vx) + normal.y * high_prec(particle1->vy);
    high_prec v1t = tangent.x * high_prec(particle1->vx) + tangent.y * high_prec(particle1->vy);
    high_prec v2n = normal.x * high_prec(particle2->vx) + normal.y * high_prec(particle2->vy);
    high_prec v2t = tangent.x * high_prec(particle2->vx) + tangent.y * high_prec(particle2->vy);

    // Calculate the initial kinetic energy in the normal direction
    high_prec initial_kinetic_energy = 0.5 * particle1->m * v1n * v1n + 0.5 * particle2->m * v2n * v2n;

    // Use combined restitution coefficient
    high_prec combined_rest = min(particle1->rest, particle2->rest);  // Changed this line
 
    // Calculate new normal velocities using conservation of momentum and combined restitution
    high_prec total_mass = high_prec(particle1->m) + high_prec(particle2->m);
    high_prec reduced_mass = (high_prec(particle1->m) * high_prec(particle2->m)) / total_mass;
    high_prec normal_velocity = v1n - v2n;

    // For perfectly inelastic collision (rest = 0), both particles should have same final velocity

    high_prec v1n_new = 0;
    high_prec v2n_new = 0;

    //cout << "Combined restitution: " << combined_rest << endl;

    

    if (combined_rest == 0) {

        
        //calculate rel velocities in x and y

        high_prec rel_vel_x = (particle2->vx - particle1->vx);
        high_prec rel_vel_y = (particle2->vy - particle1->vy);

        //recalculate the normal
        high_prec distance_pre_collission = hypot(particle2->x - particle1->x, particle2->y - particle1->y);
        high_prec normal_x = (particle2->x - particle1->x) / distance_pre_collission;
        high_prec normal_y = (particle2->y - particle1->y) / distance_pre_collission;
        //calculate the relative velocity in the normal direction
        high_prec rel_vel_n_pre = rel_vel_x * normal_x + rel_vel_y * normal_y;

        cout << "rel_vel_normal (pre collission resolve) " << rel_vel_n_pre << endl;


        // Calculate common velocity after collision (conservation of momentum)
        high_prec v_final = (particle1->m * v1n + particle2->m * v2n) / total_mass;
        v1n_new = v_final;
        v2n_new = v_final;

       
        
        //cout << "Rel vel x after inelastic collision: " << rel_vel_x << endl;
        //cout << "Rel vel y after inelastic collision: " << rel_vel_y << endl;

        //calculate the relative velocity in the normal direction
        

    } else {
        // For elastic or partially elastic collisions
        v1n_new = v1n - (1 + combined_rest) * (particle2->m / total_mass) * normal_velocity;
        v2n_new = v2n + (1 + combined_rest) * (particle1->m / total_mass) * normal_velocity;
    }

    // The tangential components remain unchanged
    high_prec v1t_new = v1t;
    high_prec v2t_new = v2t;

    // Calculate the final kinetic energy in the normal direction
    high_prec final_kinetic_energy = 0.5 * particle1->m * v1n_new * v1n_new + 0.5 * particle2->m * v2n_new * v2n_new;

    

    // Convert the scalar normal and tangential velocities into vectors
    Vector2D v1n_vec = { v1n_new * normal.x, v1n_new * normal.y };
    Vector2D v1t_vec = { v1t_new * tangent.x, v1t_new * tangent.y };
    Vector2D v2n_vec = { v2n_new * normal.x, v2n_new * normal.y };
    Vector2D v2t_vec = { v2t_new * tangent.x, v2t_new * tangent.y };

    // Update the velocities of the particles
    particle1->vx = (v1n_vec.x + v1t_vec.x);
    particle1->vy = (v1n_vec.y + v1t_vec.y);
    particle2->vx = (v2n_vec.x + v2t_vec.x);
    particle2->vy = (v2n_vec.y + v2t_vec.y);

    //convert the energy loss to heat

    // Calculate the energy loss due to restitution
    high_prec energy_loss = initial_kinetic_energy - final_kinetic_energy;
    //cout << "Energy loss: " << energy_loss << endl;
    //cout << "Collission Energy loss: " << fixed << setprecision(10) << energy_loss << defaultfloat << endl;
    



    


    //convert the energy loss to heat
    heat_ij(energy_loss, particle1, particle2);

        

    

    //calculate TE after the collision
    high_prec TE_post = calc_TE_ij(particle1, particle2);

    //print the TE before and after the collision
    //cout << "Collision TE error: " << TE_post - TE_pre << endl;

    //calculate all inputs needed for calculating the rel_vel_normal
    high_prec dx_post = particle2->x - particle1->x;
    high_prec dy_post = particle2->y - particle1->y;
    high_prec distance_post = hypot(dx_post, dy_post);
    high_prec nx_post = dx_post / distance_post;
    high_prec ny_post = dy_post / distance_post;
    high_prec rel_vel_x_post = particle2->vx - particle1->vx;
    high_prec rel_vel_y_post = particle2->vy - particle1->vy;
    //calculate the relative velocity in the normal direction
    high_prec rel_vel_n_post = rel_vel_x_post * nx_post + rel_vel_y_post * ny_post;

    cout << "rel_vel_normal (post collission resolve) " << rel_vel_n_post << endl;
}


void Engine::resolve_gravity_euler(shared_ptr<Particles> particles) {//abandoned technique
    // This function resolves the gravitational attraction between particles using Euler's method. Not sufficiently accurate for energy conservation.
    

    // Define a small value to avoid division by zero

    high_prec epsilon = 0.0000000;

    // Loop through the particles
    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {


            ///GRAVINHIB ATTEMPT 1///

            // Skip self-interaction
            if (i == j) continue;
            // Skip additionally, if the particles are AND(zero relative velocity, and their distance is the same as the sum of their radii)
            //1) make a variable to check if the particles are at rest,relative to each other
            high_prec rel_vel_x = particles->particle_list[j]->vx - particles->particle_list[i]->vx;
            high_prec rel_vel_y = particles->particle_list[j]->vy - particles->particle_list[i]->vy;
            high_prec rel_vel = hypot(rel_vel_x, rel_vel_y);
            //2) make a variable to check if the particles touching, e.g. their distance is the same as the sum of their radii
            high_prec distance = hypot(particles->particle_list[j]->x - particles->particle_list[i]->x,
                                       particles->particle_list[j]->y - particles->particle_list[i]->y);

            //2) store these particles in 2 boolean variables, make a third that is the AND condition of the first two
            bool at_rest = (abs(rel_vel_x) < epsilon && abs(particles->particle_list[j]->vy - particles->particle_list[i]->vy) < epsilon);
            bool touching = (abs(distance - (particles->particle_list[i]->rad + particles->particle_list[j]->rad)) < epsilon);
            bool should_skip = at_rest && touching;

            //3) if should_skip is true, then skip the iteration
            if (should_skip){
                cout << " zero relative velocity and touching ->Skipping interaction"<< endl;
                continue;
            }
            

            ///GRAVINHIB ATTEMPT 1///        


            //calculate TE before the update
            //high_prec TE_pre = calc_TE_ij(particles->particle_list[i], particles->particle_list[j]);
            
     
            if (distance < epsilon) {
                
                distance = epsilon;
            }

            // Calculate the gravitational force magnitude
            high_prec force = G * particles->particle_list[i]->m * particles->particle_list[j]->m / (distance * distance);

            // Calculate the components of the force
            high_prec dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            high_prec dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            high_prec fx = force * (dx / distance);
            high_prec fy = force * (dy / distance);

            // Calculate KE, PE, and TE for the system of i and j before the force update, all are magnitudes
            high_prec KE_pre = 0.5 * particles->particle_list[i]->m * (pow(particles->particle_list[i]->vx, 2) + pow(particles->particle_list[i]->vy, 2)) +
                               0.5 * particles->particle_list[j]->m * (pow(particles->particle_list[j]->vx, 2) + pow(particles->particle_list[j]->vy, 2));
            high_prec PE_pre = -G * particles->particle_list[i]->m * particles->particle_list[j]->m / distance;
            high_prec TE_pre = KE_pre + PE_pre;



            // Update the velocities of both particles to conserve momentum
            particles->particle_list[i]->vx += (fx / particles->particle_list[i]->m) * dt;
            particles->particle_list[i]->vy += (fy / particles->particle_list[i]->m) * dt;

            particles->particle_list[j]->vx -= (fx / particles->particle_list[j]->m) * dt;
            particles->particle_list[j]->vy -= (fy / particles->particle_list[j]->m) * dt;

            // Recalculate distance based on new positions
            high_prec new_dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            high_prec new_dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            high_prec new_distance = hypot(new_dx, new_dy);

            
            // Calculate KE, PE, and TE for the system of i and j after the force update using new distance
            high_prec KE_post = 0.5 * particles->particle_list[i]->m * (pow(particles->particle_list[i]->vx, 2) + pow(particles->particle_list[i]->vy, 2)) +
                                0.5 * particles->particle_list[j]->m * (pow(particles->particle_list[j]->vx, 2) + pow(particles->particle_list[j]->vy, 2));

            high_prec PE_post = -G * particles->particle_list[i]->m * particles->particle_list[j]->m / new_distance;

            high_prec TE_post = KE_post + PE_post;

            
            // Calculate the difference in total energy before and after the force update
            high_prec TE_diff = TE_post - TE_pre;
            high_prec TE_error = TE_diff / TE_pre;

            // Threshold for the error
            high_prec TE_error_threshold = 0.01; // Lower threshold for higher precision

            // If the error is greater than the threshold, print the error and distance
            if (abs(TE_error) > TE_error_threshold) {
                cout << "Distance: " << distance << endl;
                cout << "TE error %: " << TE_error * 100 << "%" << endl;
            }
        }
    }
}

void Engine::resolve_gravity_verlet(shared_ptr<Particles> particles) {
    // This function resolves the gravitational attraction between particles using the full Velocity Verlet integration scheme
    high_prec epsilon = 0.000000001; 
    int n = particles->particle_list.size();

    high_prec TE_pre = calc_TE(particles); //calc total energy before the update


    // Step 1: Initialize net forces for all particles
    vector<Vector2D> net_force(n, {0.0, 0.0});

    bool debug = false; // Set to true for detailed output
    
    // Step 2: Calculate and accumulate gravitational forces for all pairs
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {

            //Skipping functions
            // Skip self-interaction
            if (i == j) continue;
            // Skip additionally, if the particles are AND(zero relative velocity, and their distance is the same as the sum of their radii)
            //1) make a variable to check if the particles are at rest,relative to each other
            high_prec rel_vel_x = particles->particle_list[j]->vx - particles->particle_list[i]->vx;
            high_prec rel_vel_y = particles->particle_list[j]->vy - particles->particle_list[i]->vy;
            high_prec rel_vel = hypot(rel_vel_x, rel_vel_y);
            //2) make a variable to check if the particles touching, e.g. their distance
            high_prec dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            high_prec dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            high_prec distance = hypot(dx, dy);
            

            //cout << "at_rest, breakdown." << endl;
            ///cout << "rel_vel_x: " << rel_vel_x << endl;
            //cout << "rel_vel_y: " << rel_vel_y << endl;

            //the below decomposes relative velocity into normal and tangential components and then outputs them

            high_prec rel_vel_n = (rel_vel_x * dx + rel_vel_y * dy) / distance; //normal component
            high_prec rel_vel_t = hypot(rel_vel_x, rel_vel_y) - rel_vel_n; //tangential component

            //increase precision of output
            cout << fixed << setprecision(10);
            cout << "rel_vel_normal (pre-grav): " << rel_vel_n << endl;
            //cout << "rel_vel_t: " << rel_vel_t << endl;
            
            //create a bool that is true if the particles have a normal relative velocity smaller than epsilon
            bool at_rest = (rel_vel_n < epsilon);
            //cout << "at_rest: " << at_rest << endl;





            bool touching = (abs(distance - (particles->particle_list[i]->rad + particles->particle_list[j]->rad)) < epsilon);
            bool should_skip = at_rest && touching;

            if (touching) {
                //cout << "Particles are touching." << endl;
            }

            //3) if should_skip is true, then skip the iteration
            if (should_skip){
                cout << " zero relative velocity and touching ->Skipping interaction 1/2"<< endl;

                cout << "breakdown of touching 1/2" << endl;
                // increase precision of output
                //cout << fixed << setprecision(10);
                cout << "dx: " << dx << endl;
                cout << "dy: " << dy << endl;
                cout << "distance: " << distance << endl;
                //cout << "sum of radii: " << particles->particle_list[i]->rad + particles->particle_list[j]->rad << endl;
                //cout << "epsilon: " << epsilon << endl;
                cout << defaultfloat << endl;

                debug = true;


                //continue;
            }

            
            // Calculate distance components between particles i and j
            //high_prec dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            //high_prec dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            //high_prec distance = hypot(dx, dy);
            //if (distance < epsilon) distance = epsilon;
            //distance = sqrt(distance*distance + epsilon*epsilon);

            // Calculate gravitational force magnitude
            high_prec force = G * particles->particle_list[i]->m * particles->particle_list[j]->m / (distance * distance);

            // Force components
            high_prec fx = force * (dx / distance);
            high_prec fy = force * (dy / distance);

            // Accumulate forces for particles i and j (Newton's third law)
            net_force[i].x += fx;
            net_force[i].y += fy;
            net_force[j].x -= fx;
            net_force[j].y -= fy;
        }
    }

    // Step 3: Update velocities using the accumulated forces (first half-step)
    for (int i = 0; i < n; i++) {
        high_prec dt_i = dt * 1;
        particles->particle_list[i]->vx += 0.5 * (net_force[i].x / particles->particle_list[i]->m) * dt_i;
        particles->particle_list[i]->vy += 0.5 * (net_force[i].y / particles->particle_list[i]->m) * dt_i;

                
    }
    //calculate dx, dy
    high_prec dx3 = particles->particle_list[1]->x - particles->particle_list[0]->x;
    high_prec dy3 = particles->particle_list[1]->y - particles->particle_list[0]->y;

    //increase precision of output
    cout << fixed << setprecision(10);

    cout << "After step 3, first half-step, dx: " << dx3 << endl;
    cout << "After step 3, first half-step, dy: " << dy3 << endl;
    

    // Step 4: Update positions using the updated velocities

    //assume two particles in universe. check if they are touching

    cout << "Before step 4, dx: " << dx3 << endl;
    cout << "Before step 4, dy: " << dy3 << endl;

    cout << "distance before step 4: " << hypot(dx3, dy3) << endl;

    bool exactly_touching_pre4 = (abs(hypot(dx3, dy3) - (particles->particle_list[0]->rad + particles->particle_list[1]->rad)) < epsilon);
    cout << "exactly touching? " << (exactly_touching_pre4) << endl;

    //decompose velocity into normal and tangential components
    high_prec rel_vel_x = particles->particle_list[1]->vx - particles->particle_list[0]->vx;
    high_prec rel_vel_y = particles->particle_list[1]->vy - particles->particle_list[0]->vy;
    high_prec rel_vel = hypot(rel_vel_x, rel_vel_y);
    // Decompose relative velocity into normal component
    high_prec dx = particles->particle_list[1]->x - particles->particle_list[0]->x;
    high_prec dy = particles->particle_list[1]->y - particles->particle_list[0]->y;
    high_prec distance = hypot(dx, dy);
    high_prec rel_vel_n = (rel_vel_x * dx + rel_vel_y * dy) / distance; // normal component
    high_prec rel_vel_t = rel_vel - rel_vel_n; // tangential component
    cout << "rel_vel_normal (pre-step 4): " << rel_vel_n << endl;
    cout << "rel_vel_tangential (pre-step 4): " << rel_vel_t << endl;

    //first, initiate a debug_particles vector to store the alternatively calculated positions.
    vector<shared_ptr<Particle>> debug_particles;
    //copy the contents of particles->particle_list to debug_particles
    for (const auto& particle : particles->particle_list) {
        debug_particles.push_back(make_shared<Particle>(*particle));
    }



    for (int i = 0; i < n; i++) {
         high_prec dt_i = dt;

    
         particles->particle_list[i]->x += particles->particle_list[i]->vx * dt_i;
         particles->particle_list[i]->y += particles->particle_list[i]->vy * dt_i;

        
    


    }

    

    //an alternative way to update positions, assuming only two particles in universe is to decompose the velocity into normal and tangential components, then update the positions based on the normal component and the tangential component.

    
    high_prec nx = dx / distance;
    high_prec ny = dy / distance;

    // tangent = normal rotated +90°
    high_prec tx = -ny;
    high_prec ty =  nx;

    // particle 0
    high_prec vx0 = debug_particles[0]->vx; // use debug_particles for alternative calculation
    high_prec vy0 = debug_particles[0]->vy;
    
    high_prec v0n = vx0*nx + vy0*ny;   // normal magnitude
    high_prec v0t = vx0*tx + vy0*ty;   // tangential magnitude

    // rebuild world‑space velocity
    high_prec v0x = v0n*nx + v0t*tx;
    high_prec v0y = v0n*ny + v0t*ty;

    // advance position
    debug_particles[0]->x += v0x * dt;
    debug_particles[0]->y += v0y * dt;

    // particle 1 (same pattern)
    high_prec vx1 = debug_particles[1]->vx; // use debug_particles for alternative calculation
    high_prec vy1 = debug_particles[1]->vy;
    high_prec v1n = vx1*nx + vy1*ny;
    high_prec v1t = vx1*tx + vy1*ty;
    high_prec v1x = v1n*nx + v1t*tx;
    high_prec v1y = v1n*ny + v1t*ty;

    // advance position
    debug_particles[1]->x += v1x * dt;
    debug_particles[1]->y += v1y * dt;

    //calculate debug_rel_vel_normal and debug_rel_vel_tangential
    high_prec debug_rel_vel_x = debug_particles[1]->vx - debug_particles[0]->vx;
    high_prec debug_rel_vel_y = debug_particles[1]->vy - debug_particles[0]->vy;
    high_prec debug_rel_vel = hypot(debug_rel_vel_x, debug_rel_vel_y);
    // Decompose relative velocity into normal component
    high_prec debug_dx = debug_particles[1]->x - debug_particles[0]->x;
    high_prec debug_dy = debug_particles[1]->y - debug_particles[0]->y;
    high_prec debug_distance = hypot(debug_dx, debug_dy);   
    high_prec debug_rel_vel_n = (debug_rel_vel_x * debug_dx + debug_rel_vel_y * debug_dy) / debug_distance; // normal component
    high_prec debug_rel_vel_t = debug_rel_vel - debug_rel_vel_n; // tangential component



  


    




    //calculate dx, dy
    high_prec dx4 = particles->particle_list[1]->x - particles->particle_list[0]->x;
    high_prec dy4 = particles->particle_list[1]->y - particles->particle_list[0]->y;
    //cout << "After step 4, dx: " << dx4 << endl;
    //cout << "After step 4, dy: " << dy4 << endl;
    //reset precision
    cout << defaultfloat << endl;



    //assume two particles in universe. check if they are touching
    cout << "distance after step 4: " << hypot(dx4, dy4) << endl;
    bool exactly_touching_post4 = (abs(hypot(dx4, dy4) - (particles->particle_list[0]->rad + particles->particle_list[1]->rad)) < epsilon);
    cout << "exactly touching? " << (exactly_touching_post4) << endl;
    //decompose velocity into normal and tangential components
    rel_vel_x = particles->particle_list[1]->vx - particles->particle_list[0]->vx;
    rel_vel_y = particles->particle_list[1]->vy - particles->particle_list[0]->vy;
    rel_vel = hypot(rel_vel_x, rel_vel_y);
    // Decompose relative velocity into normal component
    dx = particles->particle_list[1]->x - particles->particle_list[0]->x;
    dy = particles->particle_list[1]->y - particles->particle_list[0]->y;
    distance = hypot(dx, dy);
    rel_vel_n = (rel_vel_x * dx + rel_vel_y * dy) / distance; // normal component
    rel_vel_t = rel_vel - rel_vel_n; // tangential component
    

    cout << "---METHODS ARE EQUIVALENT?." << endl;
    
    cout << "After step 4, the current implementation positions:" << endl;
    cout << "normal relative velocity after step 4: " << rel_vel_n << endl;
    cout << "tangential relative velocity after step 4: " << rel_vel_t << endl;
    cout << "Particle 0: (" << particles->particle_list[0]->x << ", " << particles->particle_list[0]->y << ")" << endl;
    cout << "Particle 1: (" << particles->particle_list[1]->x << ", " << particles->particle_list[1]->y << ")" << endl;

    cout << "After step 4, the alternative implementation positions:" << endl;
    cout << "normal relative velocity after alternative step 4: " << debug_rel_vel_n << endl;
    cout << "tangential relative velocity after alternative step 4: " << debug_rel_vel_t << endl;

    cout << "Particle 0: (" << debug_particles[0]->x << ", " << debug_particles[0]->y << ")" << endl;
    cout << "Particle 1: (" << debug_particles[1]->x << ", "<< debug_particles[1]->y << ")" << endl;

    

    cout << "---" << endl;




    // Step 5: Recalculate gravitational forces for the updated positions (second half-step)
    fill(net_force.begin(), net_force.end(), Vector2D{0.0, 0.0});
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // Recalculate distance components between particles i and j
            high_prec dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            high_prec dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            high_prec distance = hypot(dx, dy);

            // Skip if particles are at rest and touching
            high_prec rel_vel_x = particles->particle_list[j]->vx - particles->particle_list[i]->vx;
            high_prec rel_vel_y = particles->particle_list[j]->vy - particles->particle_list[i]->vy;
            high_prec rel_vel = hypot(rel_vel_x, rel_vel_y);
            // Decompose relative velocity into normal component
            high_prec rel_vel_n = (rel_vel_x * dx + rel_vel_y * dy) / distance; // normal component

            

            bool at_rest = (rel_vel_n < epsilon);
            bool touching = (abs(distance - (particles->particle_list[i]->rad + particles->particle_list[j]->rad)) < epsilon);
            bool should_skip = at_rest && touching;
            if (should_skip) {
                cout << " zero relative velocity and touching ->Skipping interaction 2/2"<< endl;
                //continue;
            } else {

                

                cout << "breakdown of touching 2/2" << endl;
                // increase precision of output
                cout << fixed << setprecision(10);
                cout << "dx: " << dx << endl;
                cout << "dy: " << dy << endl;
                cout << "distance: " << distance << endl;
                cout << "sum of radii: " << particles->particle_list[i]->rad + particles->particle_list[j]->rad << endl;
                cout << "epsilon: " << epsilon << endl;
                cout << defaultfloat << endl;

                if (debug) {
                    cout << "---" << endl;
                }
                
            }
            





            if (distance < epsilon) distance = epsilon;

            // Recalculate gravitational force magnitude
            high_prec force = G * particles->particle_list[i]->m * particles->particle_list[j]->m / (distance * distance);

            // Force components
            high_prec fx = force * (dx / distance);
            high_prec fy = force * (dy / distance);

            // Accumulate forces for particles i and j
            net_force[i].x += fx;
            net_force[i].y += fy;
            net_force[j].x -= fx;
            net_force[j].y -= fy;
        }
    }

    // Step 6: Update velocities again using the new forces (second half-step)
    for (int i = 0; i < n; i++) {
        high_prec dt_i = dt * 1;
        particles->particle_list[i]->vx += 0.5 * (net_force[i].x / particles->particle_list[i]->m) * dt_i;
        particles->particle_list[i]->vy += 0.5 * (net_force[i].y / particles->particle_list[i]->m) * dt_i;


    }

    //calculate post_grav rel_vel_n. assume only two particles in universe
    high_prec dx_post = particles->particle_list[1]->x - particles->particle_list[0]->x;
    high_prec dy_post = particles->particle_list[1]->y - particles->particle_list[0]->y;
    high_prec distance_post = hypot(dx_post, dy_post); 
    high_prec nx_post = dx_post / distance_post;
    high_prec ny_post = dy_post / distance_post;

    high_prec rel_vel_x_post = particles->particle_list[1]->vx - particles->particle_list[0]->vx;
    high_prec rel_vel_y_post = particles->particle_list[1]->vy - particles->particle_list[0]->vy;
    //calculate the relative velocity in the normal direction
    high_prec rel_vel_n_post = rel_vel_x_post * nx_post + rel_vel_y_post * ny_post;
    cout << "rel_vel_normal (post-grav): " << rel_vel_n_post << endl;


    // Step 7: validate conservation of energy
    high_prec TE_post = calc_TE(particles); //calc total energy after the update

    // Calculate the difference in total energy
    high_prec TE_diff = TE_post - TE_pre;
    high_prec TE_error = TE_diff / TE_pre;
    high_prec TE_error_threshold = 0.001; // Lower threshold for higher precision

    // Inform the user if the energy error is above the threshold
    if (abs(TE_error) > TE_error_threshold) {
        cout << "Warning: Energy loss (in Verlet): " << TE_error * 100 << "%" << endl;


    }

    //print the TE before and after the update

    total_TE_error_verlet += TE_post - TE_pre;

    //cout << "Verlet TE error: " << TE_post - TE_pre << endl;
    

}



void Engine::update_locations(shared_ptr<Particles> particles, shared_ptr<backed_scaler> scaler) {
    //this function will update the locations of the particles

    //1. loop through the particles
    for (int i = 0; i < particles->particle_list.size(); i++) {
        
        particles->particle_list[i]->x += particles->particle_list[i]->vx * scaler->scaler[i] * dt;
        particles->particle_list[i]->y += particles->particle_list[i]->vy * scaler->scaler[i] * dt;
        //particles->particle_list[i]->z += particles->particle_list[i]->vz * scaler->scaler[i] * dt;
    }
}

void Engine::run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {
    //this function will save the snapshots to a csv file

    // if (true){
    //     return;
    // }
    //1. define the file name, located on Inputs\rendered_scenarios

    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    //2. open the file

    ofstream file(file_name);

    //check if the file is open

    if (!file.is_open()) {
        cout << "Failed to write snapshots to cache!" << endl;
        return;
    }

    //3. write headers

    file << "step_id, particle_id,r,g,b,x,y,z,vx,vy,vz,m,rad,temp,rest" << endl;

    //4. loop through the snapshots and write the data to the file, adding the step_id

    for (int i = 0; i < particle_states->snaps.size(); i++) {
        for (int j = 0; j < particle_states->snaps[i]->particle_list.size(); j++) {
            file << std::fixed << std::setprecision(15) // Set precision to 15 decimal places
                << i << ","
                << particle_states->snaps[i]->particle_list[j]->particle_id << ","
                << (particle_states->snaps[i]->particle_list[j]->r) << ","
                << (particle_states->snaps[i]->particle_list[j]->g) << ","
                << (particle_states->snaps[i]->particle_list[j]->b) << ","
                << (particle_states->snaps[i]->particle_list[j]->x) << ","
                << (particle_states->snaps[i]->particle_list[j]->y) << ","
                << (particle_states->snaps[i]->particle_list[j]->z) << ","
                << (particle_states->snaps[i]->particle_list[j]->vx) << ","
                << (particle_states->snaps[i]->particle_list[j]->vy) << ","
                << (particle_states->snaps[i]->particle_list[j]->vz) << ","
                << (particle_states->snaps[i]->particle_list[j]->m) << ","
                << (particle_states->snaps[i]->particle_list[j]->rad) << ","
                << (particle_states->snaps[i]->particle_list[j]->temp) << ","
                << (particle_states->snaps[i]->particle_list[j]->rest) << endl;
        }
    }

    //5. save and close the file

    cout << "Saving simulation..." << endl;

    file.close();

    cout << "Snapshots saved to " << file_name << endl << endl;
    

}

bool Engine::cache_exists(shared_ptr<scenario> scenario) {
    //this function will check if a cache file exists for the scenario

    //1. define the file name, located on Inputs\rendered_scenarios

    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    //2. open the file

    ifstream file(file_name);

    //check if the file is open

    if (!file.is_open()) {
        return false;
    }

    cout << "Cache file found for " << scenario->name << endl;

    return true;
}

shared_ptr<snapshots> Engine::run_from_cache(shared_ptr<scenario> scenario) {

    //this function will read the snapshots from a csv file

    //1. define the file name, located on Inputs\rendered_scenarios

    string file_name = "Inputs/rendered_scenarios/" + scenario->name + ".csv";

    //2. open the file

    ifstream file(file_name);

    //3. read the file, using the first column to determine the step_id, store every step in a new snapshot

    shared_ptr<snapshots> particle_states = make_shared<snapshots>();

    // Create a CSV reader object
    typedef io::trim_chars<' ', '\t'> TrimPolicy;
    typedef io::double_quote_escape<',', '\"'> QuotePolicy;
    const int column_count = 15; // Adjust the column count based on your CSV structure
    io::CSVReader<column_count, TrimPolicy, QuotePolicy> in(file_name);

    // 4. Read the contents of the CSV file and store each row in a particle object. Store all particle objects in a vector.
    string col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15;
    int step_id = 0;
    shared_ptr<Particles> particles = make_shared<Particles>();

    // Discard the header row. This list is a guide to the columns in the CSV file.
    in.read_header(io::ignore_extra_column, "step_id", "particle_id", "r", "g", "b", "x", "y", "z", "vx", "vy", "vz", "m", "rad", "rest", "temp");


    while (in.read_row(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15)) {
        // Check if the step_id has changed
        if (stoi(col1) != step_id) {
            // Store the current snapshot
            particle_states->snaps.push_back(particles);
            // Create a new snapshot
            particles = make_shared<Particles>();
            step_id = stoi(col1);
        }

        // Create a new particle
        shared_ptr<Particle> particle = make_shared<Particle>();
        particle->particle_id = stoi(col2);
        particle->r = stod(col3);
        particle->g = stod(col4);
        particle->b = stod(col5);
        particle->x = stod(col6);
        particle->y = stod(col7);
        particle->z = stod(col8);
        particle->vx = stod(col9);
        particle->vy = stod(col10);
        particle->vz = stod(col11);
        particle->m = stod(col12);
        particle->rad = stod(col13);
        particle->rest = stod(col14);
        particle->temp = stod(col15);

        // Add the particle to the snapshot
        particles->particle_list.push_back(particle);

    }

    // Store the last snapshot
    particle_states->snaps.push_back(particles);

    cout << "Snapshots successfully loaded from cache: " << file_name << endl;

    return particle_states;
}


// Function to calculate the total kinetic energy (KE) and potential energy (PE) for a system of particles
high_prec Engine::calc_TE(shared_ptr<Particles> particles) {
    int n = particles->particle_list.size();
    high_prec KE = 0.0;
    high_prec PE = 0.0;
    high_prec HE = 0.0;  // Heating energy component
    high_prec epsilon = 0.001; // To avoid division by zero

    // Calculate kinetic energy and heating energy
    for (int i = 0; i < n; i++) {
        KE += 0.5 * particles->particle_list[i]->m * 
              (pow(particles->particle_list[i]->vx, 2) + pow(particles->particle_list[i]->vy, 2));
        HE += particles->particle_list[i]->temp; // Add heating energy of each particle
    }

    // Calculate potential energy
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            high_prec dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            high_prec dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            high_prec distance = hypot(dx, dy);
            distance = sqrt(distance * distance + epsilon * epsilon);
            PE += -G * particles->particle_list[i]->m * particles->particle_list[j]->m / distance;
        }
    }

    // Return the total energy (KE + PE + HE)
    return KE + PE + HE;
}

high_prec Engine::calculate_kinetic_energy(shared_ptr<Particle> p1, shared_ptr<Particle> p2) {
    high_prec KE = 0.0;
    KE += 0.5 * p1->m * (pow(p1->vx, 2) + pow(p1->vy, 2));
    KE += 0.5 * p2->m * (pow(p2->vx, 2) + pow(p2->vy, 2));
    return KE;
}

high_prec Engine::calculate_potential_energy(shared_ptr<Particle> p1, shared_ptr<Particle> p2) {
    high_prec epsilon = 0.001;
    high_prec dx = p2->x - p1->x;
    high_prec dy = p2->y - p1->y;
    high_prec distance = hypot(dx, dy);
    distance = sqrt(distance * distance + epsilon * epsilon);
    return -G * p1->m * p2->m / distance;
}

high_prec Engine::calculate_heating_energy(shared_ptr<Particle> p1, shared_ptr<Particle> p2) {
    return static_cast<high_prec>(p1->temp + p2->temp);
}

high_prec Engine::calc_TE_ij(shared_ptr<Particle> p1, shared_ptr<Particle> p2, bool verbose) {
    high_prec KE = calculate_kinetic_energy(p1, p2);
    high_prec PE = calculate_potential_energy(p1, p2);
    high_prec HE = calculate_heating_energy(p1, p2);
    
    if (verbose) {
        cout << "TE=" << KE + PE + HE
                  << ", composed of KE=" << KE
                  << ", PE=" << PE
                  << ", HE=" << HE << endl;
    }

    return KE + PE + HE;
}


momentum Engine::calc_mom(shared_ptr<Particles> particles) {
    high_prec total_momentum_x = 0.0;
    high_prec total_momentum_y = 0.0;

    // Sum up the x and y momentum components for each particle
    for (const auto& particle : particles->particle_list) {
        total_momentum_x += particle->m * particle->vx;
        total_momentum_y += particle->m * particle->vy;
    }

    return {total_momentum_x, total_momentum_y};
}

// Corrected `calc_mom_ij` for two particles
momentum Engine::calc_mom_ij(shared_ptr<Particle> p1, shared_ptr<Particle> p2) {
    high_prec total_momentum_x = p1->m * p1->vx + p2->m * p2->vx;
    high_prec total_momentum_y = p1->m * p1->vy + p2->m * p2->vy;

    return {total_momentum_x, total_momentum_y};
}

