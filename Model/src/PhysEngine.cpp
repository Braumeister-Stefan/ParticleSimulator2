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
#include <limits>
#include <cfloat>
#include <stdexcept>

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

static high_prec collissions = 0; //counter for collisions

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

    // In your run() method, after loading particles:
    cout << "=== INITIAL SYSTEM CHECK ===" << endl;
    momentum initial_system_mom = calc_mom(particles);
    cout << "Initial system momentum: x=" << initial_system_mom.x << ", y=" << initial_system_mom.y << endl;

    // Check for NaN or infinite values
    for (int i = 0; i < particles->particle_list.size(); i++) {
        auto p = particles->particle_list[i];
        if (!isfinite(p->vx) || !isfinite(p->vy) || 
            !isfinite(p->x) || !isfinite(p->y) ||
            !isfinite(p->m)) {
            cout << "ERROR: Particle " << i << " has invalid values!" << endl;
            cout << "  vx: " << p->vx << ", vy: " << p->vy << endl;
            cout << "  x: " << p->x << ", y: " << p->y << endl;
            cout << "  m: " << p->m << endl;
        }
    }
    cout << "=== END INITIAL CHECK ===" << endl << endl;

    cout << "Engine is initialized." << endl;

    // 1. Initialize the snapshots object
    shared_ptr<snapshots> particle_states = make_shared<snapshots>();

    cout << scenario->name << "'s initial states loaded" << endl << endl;

    // === ADD THIS DIAGNOSTIC CALL ===
    cout << "=== FINAL MOMENTUM ANALYSIS ===" << endl;
    check_momentum_metrics(particles, initial_system_mom);
    cout << "=== END FINAL ANALYSIS ===" << endl << endl;
    // ================================

    cout << "Total TE error (sum of parts): "
     << total_TE_error_overlap + total_TE_error_collision + total_TE_error_verlet << endl;


    cout << "Press enter to start the simulation." << endl;
    cin.ignore();
    cin.get();

    // 3. Initialize time step
    Engine::dt = scenario->dt;


    // Store initial momentum for diagnostics
    momentum initial_system_momentum = calc_mom(particles);
    cout << "Stored initial momentum: x=" << initial_system_momentum.x 
     << ", y=" << initial_system_momentum.y << endl;

    // Calculate adaptive time step based on maximum acceleration
    // high_prec max_accel = calculate_max_acceleration(particles);

    // ADD THIS:
    // Calculate adaptive time step based on maximum acceleration
    high_prec max_accel = calculate_max_acceleration(particles);
    high_prec adaptive_dt = sqrt(1e-6L / max_accel); // epsilon = 1e-6
    Engine::dt = min(scenario->dt, adaptive_dt); // Use more restrictive dt

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

        auto start_time = high_resolution_clock::now();
        update_iter = i;

        high_prec te_pre_update = calc_TE(particles);
        momentum mom_pre_update = calc_mom(particles);

        // Call update_particles
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

    cout << "Total collisions: " << collissions << endl;

    run_to_cache(scenario, particle_states);
    return particle_states;
}

void Engine::update_particles(shared_ptr<Particles> particles)
{
    // Track momentum at each substep
    momentum mom_start = calc_mom(particles);
    cout << "DEBUG: Initial momentum - x: " << mom_start.x << ", y: " << mom_start.y << endl;

    bool no_overlap = false;
    high_prec te_pre_overlap = calc_TE(particles);

    // Overlap resolution
    for (int i = 0; i < 8; i++) {
        if (no_overlap) break;
        no_overlap = true;
        overlap_iter = i;
        no_overlap = resolve_overlaps(particles);
        
        // APPLY GLOBAL MOMENTUM CORRECTION AFTER EACH OVERLAP ITERATION
        correct_global_momentum(particles, mom_start);
    }

    // Collisions
    resolve_collisions(particles);
    
    // APPLY GLOBAL MOMENTUM CORRECTION AFTER COLLISIONS
    correct_global_momentum(particles, mom_start);

    // Gravity (Velocity Verlet)
    resolve_gravity_verlet(particles);
    
    // FINAL GLOBAL MOMENTUM CORRECTION
    correct_global_momentum(particles, mom_start);

    momentum mom_final = calc_mom(particles);
    cout << "DEBUG: Final momentum change - x: " << mom_final.x - mom_start.x 
         << ", y: " << mom_final.y - mom_start.y << endl;
}


bool Engine::resolve_overlaps(shared_ptr<Particles> particles) {
    bool no_overlap = true;
    high_prec tTE_pre = calc_TE(particles);


    // STORE INITIAL MOMENTUM FOR THIS OVERLAP PASS
    momentum initial_overlap_mom = calc_mom(particles);

    high_prec error_threshold = 0.001;

    int pair_count = 0;

    //reset marginal error metrics
    s_margin_TE_error_overlap_ij_transl = 0.0;
    s_margin_TE_error_overlap_ij_corrected = 0.0;

    // IMPROVED: Use adaptive overlap resolution
    high_prec system_scale = calculate_system_scale(particles);
    high_prec adaptive_tolerance = system_scale * 1e-10;
    
    int max_iterations = 20;
    int iterations = 0;
    high_prec max_overlap = 0;
    
    do {
        no_overlap = true;
        max_overlap = 0;
        
        // Loop through particle pairs
        for (int i = 0; i < particles->particle_list.size(); i++) {
            for (int j = i + 1; j < particles->particle_list.size(); j++) {
                high_prec overlap_amount = calculate_overlap_amount(particles->particle_list[i], particles->particle_list[j]);
                if (overlap_amount > adaptive_tolerance) {
                    no_overlap = false;
                    max_overlap = max(max_overlap, overlap_amount);
                    pair_count++;

                    high_prec TE_pre_overlap = calc_TE(particles);
                    
                    resolve_overlap_ij(particles->particle_list[i], particles->particle_list[j]);

                    high_prec TE_post_overlap = calc_TE(particles);
                    high_prec TE_error_overlap = (TE_post_overlap - TE_pre_overlap) / TE_pre_overlap;
                }
            }
        }
        iterations++;
    } while (!no_overlap && iterations < max_iterations && max_overlap > adaptive_tolerance);

    high_prec tTE_post = calc_TE(particles);
    high_prec tTE_error = (tTE_post - tTE_pre) / tTE_pre;


    // APPLY MOMENTUM CORRECTION FOR THIS OVERLAP PASS
    momentum final_overlap_mom = calc_mom(particles);
    high_prec mom_error_x = final_overlap_mom.x - initial_overlap_mom.x;
    high_prec mom_error_y = final_overlap_mom.y - initial_overlap_mom.y;
    
    if (fabs(mom_error_x) > 1e-12 || fabs(mom_error_y) > 1e-12) {
        cout << "Overlap pass momentum error - dx: " << mom_error_x << ", dy: " << mom_error_y << endl;
        correct_global_momentum(particles, initial_overlap_mom);
    }

    return no_overlap;

}


void Engine::resolve_overlap_ij(shared_ptr<Particle> particle_i, shared_ptr<Particle> particle_j) 
{
    
    // STORE INITIAL MOMENTUM FOR THIS PAIR
    momentum initial_mom_ij = calc_mom_ij(particle_i, particle_j);
    
    overlap_ij_iter++;

    // Compute displacement and distance
    high_prec dx = particle_j->x - particle_i->x;
    high_prec dy = particle_j->y - particle_i->y;
    high_prec dist = hypot(dx, dy);

    // Avoid division by zero
    if (dist == 0.0) {
        dx = 0.001 * (rand() / (high_prec)RAND_MAX - 0.5);
        dy = 0.001 * (rand() / (high_prec)RAND_MAX - 0.5);
        dist = hypot(dx, dy);
    }

    high_prec overlap_amount = (particle_i->rad + particle_j->rad) - dist;
    if (overlap_amount <= 0) return;

    // Collision normal
    high_prec nx = dx / dist;
    high_prec ny = dy / dist;

    // Separate particles proportionally to their masses
    high_prec m1 = particle_i->m;
    high_prec m2 = particle_j->m;
    high_prec total_mass = m1 + m2;
    high_prec separation_buffer = 1e-5L;
    high_prec total_separation = overlap_amount + separation_buffer;
    high_prec d1 = (total_separation * m2) / total_mass;
    high_prec d2 = (total_separation * m1) / total_mass;

    particle_i->x -= nx * d1;
    particle_i->y -= ny * d1;
    particle_j->x += nx * d2;
    particle_j->y += ny * d2;

    // --- Velocity update (elastic collision) ---
    high_prec v1n = particle_i->vx * nx + particle_i->vy * ny;
    high_prec v2n = particle_j->vx * nx + particle_j->vy * ny;

    // Only process if moving toward each other
    if (v2n - v1n < 0) {
        high_prec v1t = -particle_i->vx * ny + particle_i->vy * nx;
        high_prec v2t = -particle_j->vx * ny + particle_j->vy * nx;

        high_prec v1n_new = ((m1 - m2) * v1n + 2 * m2 * v2n) / (m1 + m2);
        high_prec v2n_new = ((m2 - m1) * v2n + 2 * m1 * v1n) / (m1 + m2);

        particle_i->vx = v1n_new * nx - v1t * ny;
        particle_i->vy = v1n_new * ny + v1t * nx;
        particle_j->vx = v2n_new * nx - v2t * ny;
        particle_j->vy = v2n_new * ny + v2t * nx;
    }

    // --- Energy correction (momentum-conserving) ---
    high_prec v1x_pre = particle_i->vx;
    high_prec v1y_pre = particle_i->vy;
    high_prec v2x_pre = particle_j->vx;
    high_prec v2y_pre = particle_j->vy;

    high_prec TE_pre_ij = 0.5L * m1 * (v1x_pre*v1x_pre + v1y_pre*v1y_pre) 
                        + 0.5L * m2 * (v2x_pre*v2x_pre + v2y_pre*v2y_pre);
    high_prec TE_post_ij = 0.5L * m1 * (particle_i->vx*particle_i->vx + particle_i->vy*particle_i->vy)
                         + 0.5L * m2 * (particle_j->vx*particle_j->vx + particle_j->vy*particle_j->vy);

    high_prec delta_E = TE_post_ij - TE_pre_ij;

    // Apply correction only if energy deviation is significant
    if (fabs(delta_E) > 1e-12L)
        apply_energy_correction(particle_i, particle_j, delta_E);

    // Done — no global scaling applied

    // VERIFY MOMENTUM CONSERVATION FOR THIS PAIR
    momentum final_mom_ij = calc_mom_ij(particle_i, particle_j);
    high_prec mom_error_x = final_mom_ij.x - initial_mom_ij.x;
    high_prec mom_error_y = final_mom_ij.y - initial_mom_ij.y;
    
    if (fabs(mom_error_x) > 1e-12 || fabs(mom_error_y) > 1e-12) {
        cout << "MOMENTUM ERROR in resolve_overlap_ij: dx=" << mom_error_x 
             << ", dy=" << mom_error_y << " at iter " << overlap_ij_iter << endl;
        
        // CORRECT THE MOMENTUM ERROR
        high_prec total_mass = particle_i->m + particle_j->m;
        particle_i->vx -= mom_error_x / total_mass;
        particle_j->vx += mom_error_x / total_mass;
        particle_i->vy -= mom_error_y / total_mass;
        particle_j->vy += mom_error_y / total_mass;
    }

}

// Momentum-conserving energy correction
void Engine::apply_energy_correction(
    shared_ptr<Particle> p1,
    shared_ptr<Particle> p2,
    high_prec delta_E)
{
    // Skip if negligible
    if (fabs(delta_E) < 1e-18L) return;

    // Use existing impulse-based routine that preserves momentum
    resolve_energy_gap(p1, p2, delta_E);
}


void Engine::resolve_energy_gap(std::shared_ptr<Particle> p1,
    std::shared_ptr<Particle> p2,
    boost::multiprecision::cpp_dec_float_50 delta_E)
{
    
    // STORE INITIAL MOMENTUM
    momentum initial_mom = calc_mom_ij(p1, p2);
    
    using high_prec = boost::multiprecision::cpp_dec_float_50;

    //add to counter of energy gap corrections
    energy_gap_corrections++;

    // Extract masses
    high_prec m1 = p1->m;
    high_prec m2 = p2->m;

    // Compute collision normal
    high_prec dx = p2->x - p1->x;
    high_prec dy = p2->y - p1->y;
    high_prec distance = hypot(dx, dy);

    if (distance == 0.0) {
        return;  // Cannot compute normal, skip adjustment
    }

    high_prec nx = dx / distance;
    high_prec ny = dy / distance;

    // Compute normal and tangential directions
    high_prec tx = -ny;
    high_prec ty = nx;

    // Compute normal velocities
    high_prec v1n = p1->vx * nx + p1->vy * ny;
    high_prec v2n = p2->vx * nx + p2->vy * ny;

    // Compute tangential velocities
    high_prec v1t = p1->vx * tx + p1->vy * ty;
    high_prec v2t = p2->vx * tx + p2->vy * ty;

    // Compute reduced mass
    high_prec reduced_mass = (m1 * m2) / (m1 + m2);

    // Compute initial kinetic energy in normal and tangential directions
    high_prec KE_rel_n_initial = 0.5 * reduced_mass * (v2n - v1n) * (v2n - v1n);
    high_prec KE_rel_t_initial = 0.5 * reduced_mass * (v2t - v1t) * (v2t - v1t);

    high_prec delta_E_used = 0; // Track total energy removed

    // First, remove delta_E from the normal component
    high_prec KE_rel_n = KE_rel_n_initial;
    if (KE_rel_n < delta_E) {
        delta_E_used += KE_rel_n;
        delta_E -= KE_rel_n;
        KE_rel_n = 0;
    } else {
        KE_rel_n -= delta_E;
        delta_E_used += delta_E;
        delta_E = 0;
    }

    // Adjust normal velocity
    high_prec new_rel_vel_n = sqrt(std::max(high_prec(0.0), high_prec(2 * KE_rel_n / reduced_mass)));
    if (v2n - v1n < 0) new_rel_vel_n = -new_rel_vel_n;
    high_prec delta_rel_vel_n = new_rel_vel_n - (v2n - v1n);

    // Compute impulse and apply adjustment
    high_prec impulse_n = delta_rel_vel_n * reduced_mass;
    high_prec dvx1_n = (impulse_n / m1) * nx;
    high_prec dvy1_n = (impulse_n / m1) * ny;
    high_prec dvx2_n = (impulse_n / m2) * nx;
    high_prec dvy2_n = (impulse_n / m2) * ny;

    p1->vx -= dvx1_n;
    p1->vy -= dvy1_n;
    p2->vx += dvx2_n;
    p2->vy += dvy2_n;

    // If there's remaining energy to remove, take it from the tangential component
    high_prec KE_rel_t = KE_rel_t_initial;

    if (delta_E > 0) {
        if (KE_rel_t < delta_E) {
            delta_E_used += KE_rel_t;
            delta_E -= KE_rel_t;
            KE_rel_t = 0;
            energy_gap_corrections_incomplete++;
        } else {
            KE_rel_t -= delta_E;
            delta_E_used += delta_E;
            delta_E = 0;
        }

        // Adjust tangential velocity
        high_prec new_rel_vel_t = sqrt(max(high_prec(0.0), high_prec(2 * KE_rel_t / reduced_mass)));
        if (v2t - v1t < 0) new_rel_vel_t = -new_rel_vel_t;
        high_prec delta_rel_vel_t = new_rel_vel_t - (v2t - v1t);

        // Compute impulse for tangential adjustment
        high_prec impulse_t = delta_rel_vel_t * reduced_mass;
        high_prec dvx1_t = (impulse_t / m1) * tx;
        high_prec dvy1_t = (impulse_t / m1) * ty;
        high_prec dvx2_t = (impulse_t / m2) * tx;
        high_prec dvy2_t = (impulse_t / m2) * ty;

        p1->vx -= dvx1_t;
        p1->vy -= dvy1_t;
        p2->vx += dvx2_t;
        p2->vy += dvy2_t;
    }

    // VERIFY AND CORRECT MOMENTUM AFTER ENERGY ADJUSTMENT
    momentum final_mom = calc_mom_ij(p1, p2);
    high_prec mom_error_x = final_mom.x - initial_mom.x;
    high_prec mom_error_y = final_mom.y - initial_mom.y;
    
    if (fabs(mom_error_x) > 1e-15 || fabs(mom_error_y) > 1e-15) {
        cout << "MOMENTUM ERROR in resolve_energy_gap: dx=" << mom_error_x 
             << ", dy=" << mom_error_y << endl;
        
        // CORRECT MOMENTUM
        high_prec total_mass = p1->m + p2->m;
        p1->vx -= mom_error_x / total_mass;
        p2->vx += mom_error_x / total_mass;
        p1->vy -= mom_error_y / total_mass;
        p2->vy += mom_error_y / total_mass;
    }

}


// Add this function in PhysEngine.cpp, around line 950 (near other helper functions)
high_prec safe_momentum_error(high_prec current, high_prec initial, const string& direction) {
    if (fabs(initial) < 1e-12) {
        // If initial momentum is near zero, use absolute error
        high_prec absolute_error = current - initial;
        if (fabs(absolute_error) > 1e-12) {
            cout << "WARNING: " << direction << "-momentum error with near-zero initial: " 
                 << absolute_error << endl;
        }
        return absolute_error;
    } else {
        // Normal relative error calculation
        return (current - initial) / initial;
    }
}

// Add this function in PhysEngine.cpp, around line 965
void Engine::check_momentum_metrics(shared_ptr<Particles> particles, momentum initial_mom) {
    momentum final_mom = calc_mom(particles);
    
    cout << "=== MOMENTUM METRICS DEBUG ===" << endl;
    cout << "Initial: x=" << initial_mom.x << ", y=" << initial_mom.y << endl;
    cout << "Final:   x=" << final_mom.x << ", y=" << final_mom.y << endl;
    
    high_prec abs_error_x = final_mom.x - initial_mom.x;
    high_prec abs_error_y = final_mom.y - initial_mom.y;
    
    cout << "Absolute Error: x=" << abs_error_x << ", y=" << abs_error_y << endl;
    
    // Safe percentage calculation
    if (fabs(initial_mom.x) > 1e-12) {
        cout << "X-Momentum Error: " << (abs_error_x / initial_mom.x * 100.0) << "%" << endl;
    } else {
        cout << "X-Momentum Error: N/A (initial x-momentum is zero)" << endl;
    }
    
    if (fabs(initial_mom.y) > 1e-12) {
        cout << "Y-Momentum Error: " << (abs_error_y / initial_mom.y * 100.0) << "%" << endl;
    } else {
        cout << "Y-Momentum Error: N/A (initial y-momentum is zero)" << endl;
    }
    cout << "=== END MOMENTUM DEBUG ===" << endl;
}


high_prec Engine::heat_ij(high_prec E, shared_ptr<Particle> particle_i, shared_ptr<Particle> particle_j){
    //convert the given energy to heat
    high_prec heat = E * (particle_i->m + particle_j->m) / 2;

    //convert the heat to temperature
    high_prec temp = heat / (particle_i->m + particle_j->m);

    particle_i->temp += temp.convert_to<double>();
    particle_j->temp += temp.convert_to<double>();

    return E;
}

void Engine::resolve_collisions(shared_ptr<Particles> particles) {
    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            bool collision = check_collision(particles->particle_list[i], particles->particle_list[j]);
            
            if (collision) {
                //add collision to the counter
                collissions++;
                resolve_collision(particles->particle_list[i], particles->particle_list[j]);
            }
        }
    }
}

void Engine::resolve_collision(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    resolve_collission(particle1, particle2);
}

bool Engine::check_collision(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    //this function will check if two particles are colliding

    //1. calculate the distance between the particles
    high_prec distance = hypot(particle1->x - particle2->x, particle1->y - particle2->y);
    
    //calculate relative velocity 
    high_prec relative_velocity = hypot(particle1->vx - particle2->vx, particle1->vy - particle2->vy);

    //2. check if the distance+threshold is smaller than the sum of the radii of the particles
    high_prec sum_radii = particle1->rad + particle2->rad;
    high_prec tolerance = 1e-8L;
    if ((abs(distance - sum_radii) < tolerance) && (relative_velocity > 0)) {
        return true;
    } else {
        return false;
    }
}

bool Engine::check_overlap(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    high_prec distance = hypot(particle1->x - particle2->x, particle1->y - particle2->y);
    high_prec sum_radii = particle1->rad + particle2->rad;
    high_prec tolerance = 1e-6L; // Small tolerance for floating point
    
    // Return true if particles are overlapping (distance < sum_radii)
    // Add small buffer to prevent sticking
    return (distance < sum_radii - tolerance);
}


void Engine::resolve_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2) {
    // Store pre-collision momentum for verification
    momentum mom_pre = calc_mom_ij(particle1, particle2);

    // Calculate the normal vector
    high_prec dx = particle2->x - particle1->x;
    high_prec dy = particle2->y - particle1->y;
    high_prec distance = sqrt(dx*dx + dy*dy);
    
    if (distance == 0.0) return;
    
    high_prec nx = dx / distance;
    high_prec ny = dy / distance;

    // Calculate the tangential vector
    high_prec tx = -ny;
    high_prec ty = nx;

    // Decompose velocities into normal and tangential components
    high_prec v1n = particle1->vx * nx + particle1->vy * ny;
    high_prec v1t = particle1->vx * tx + particle1->vy * ty;
    high_prec v2n = particle2->vx * nx + particle2->vy * ny;
    high_prec v2t = particle2->vx * tx + particle2->vy * ty;

    // Use combined restitution coefficient
    high_prec combined_rest = min(particle1->rest, particle2->rest);
    
    // Masses
    high_prec m1 = particle1->m;
    high_prec m2 = particle2->m;
    high_prec total_mass = m1 + m2;

    // Calculate new normal velocities
    high_prec v1n_new, v2n_new;
    
    if (combined_rest == 0) {
        // Perfectly inelastic - common velocity
        high_prec v_final = (m1 * v1n + m2 * v2n) / total_mass;
        v1n_new = v_final;
        v2n_new = v_final;
    } else {
        // Elastic or partially elastic
        v1n_new = ((m1 - combined_rest * m2) * v1n + (1 + combined_rest) * m2 * v2n) / total_mass;
        v2n_new = ((m2 - combined_rest * m1) * v2n + (1 + combined_rest) * m1 * v1n) / total_mass;
    }

    // Tangential components remain unchanged (frictionless)
    high_prec v1t_new = v1t;
    high_prec v2t_new = v2t;

    // Reconstruct velocity vectors
    particle1->vx = v1n_new * nx + v1t_new * tx;
    particle1->vy = v1n_new * ny + v1t_new * ty;
    particle2->vx = v2n_new * nx + v2t_new * tx;
    particle2->vy = v2n_new * ny + v2t_new * ty;

    // Verify momentum conservation
    momentum mom_post = calc_mom_ij(particle1, particle2);
    high_prec mom_error_x = (mom_post.x - mom_pre.x) / (abs(mom_pre.x) + 1e-12);
    high_prec mom_error_y = (mom_post.y - mom_pre.y) / (abs(mom_pre.y) + 1e-12);
    
    if (abs(mom_error_x) > 0.01 || abs(mom_error_y) > 0.01) {
        cout << "WARNING: Momentum not conserved in collision! Error: " 
             << mom_error_x * 100 << "%, " << mom_error_y * 100 << "%" << endl;
    }


}


void Engine::resolve_gravity_verlet(std::shared_ptr<Particles> particles) {
    high_prec epsilon = 1e-6L; // small softening to avoid divide-by-zero
    int n = (int)particles->particle_list.size();

    // Total energy before integration (for diagnostics)
    high_prec TE_pre = calc_TE(particles);

    // Prepare force accumulator
    std::vector<Vector2D> net_force(n, {0.0, 0.0});

    // --- Step 1: compute pairwise forces (first evaluation) ---
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            auto pi = particles->particle_list[i];
            auto pj = particles->particle_list[j];

            high_prec dx = pj->x - pi->x;
            high_prec dy = pj->y - pi->y;

            high_prec dist2 = dx*dx + dy*dy + epsilon*epsilon;
            high_prec distance = sqrt(dist2);
            if (distance == 0) distance = epsilon;

            high_prec force = G * pi->m * pj->m / dist2;

            double fcheck = static_cast<double>(force);
            if (!std::isfinite(fcheck) || fabs(fcheck) > 1e100) {
                std::cout << "Warning: extreme gravitational force detected (first pass). Clipping to zero." << std::endl;
                force = high_prec(0);
            }

            high_prec fx = force * (dx / distance);
            high_prec fy = force * (dy / distance);

            net_force[i].x += fx;
            net_force[i].y += fy;
            net_force[j].x -= fx;
            net_force[j].y -= fy;
        }
    }

    // --- Step 2: first half-step velocity update ---
    for (int i = 0; i < n; ++i) {
        auto p = particles->particle_list[i];
        p->vx += 0.5L * (net_force[i].x / p->m) * dt;
        p->vy += 0.5L * (net_force[i].y / p->m) * dt;
    }

    // --- Step 3: position update (full step) ---
    for (int i = 0; i < n; ++i) {
        auto p = particles->particle_list[i];
        p->x += p->vx * dt;
        p->y += p->vy * dt;
    }

    // --- Step 4: recompute forces at new positions (second evaluation) ---
    std::fill(net_force.begin(), net_force.end(), Vector2D{0.0, 0.0});
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            auto pi = particles->particle_list[i];
            auto pj = particles->particle_list[j];

            high_prec dx = pj->x - pi->x;
            high_prec dy = pj->y - pi->y;

            high_prec dist2 = dx*dx + dy*dy + epsilon*epsilon;
            high_prec distance = sqrt(dist2);
            if (distance == 0) distance = epsilon;

            high_prec force = G * pi->m * pj->m / dist2;

            double fcheck = static_cast<double>(force);
            if (!std::isfinite(fcheck) || fabs(fcheck) > 1e100) {
                std::cout << "Warning: extreme gravitational force detected (second pass). Clipping to zero." << std::endl;
                force = high_prec(0);
            }

            high_prec fx = force * (dx / distance);
            high_prec fy = force * (dy / distance);

            net_force[i].x += fx;
            net_force[i].y += fy;
            net_force[j].x -= fx;
            net_force[j].y -= fy;
        }
    }

    // --- Step 5: second half-step velocity update ---
    for (int i = 0; i < n; ++i) {
        auto p = particles->particle_list[i];
        p->vx += 0.5L * (net_force[i].x / p->m) * dt;
        p->vy += 0.5L * (net_force[i].y / p->m) * dt;
    }

    // --- Diagnostics ---
    high_prec TE_post = calc_TE(particles);
    high_prec TE_diff = TE_post - TE_pre;
    high_prec TE_error = TE_pre == 0 ? TE_diff : TE_diff / TE_pre;

    if (abs(TE_error) > 0.001L) {
        std::cout << "[DEBUG_VERLET] Energy change in step: "
                  << (double)(TE_error * 100.0L)
                  << "% (TE_pre=" << (double)TE_pre
                  << ", TE_post=" << (double)TE_post << ")"
                  << std::endl;
    }

    total_TE_error_verlet += (TE_post - TE_pre);
    // --- Diagnostics (limited debug for first 100 steps) ---
    static int local_step = 0;

    if (local_step < 100) {
        high_prec TE_post = calc_TE(particles);
        high_prec TE_diff = TE_post - TE_pre;
        high_prec TE_error = TE_pre == 0 ? TE_diff : TE_diff / TE_pre;

            std::cout << "[DEBUG_VERLET_STEP] "
              << "step=" << local_step
              << ", drift=" << (double)(TE_error * 100.0L)
              << "%, TE_pre=" << (double)TE_pre
              << ", TE_post=" << (double)TE_post
              << std::endl;
    }

    local_step++;

} // ✅ function ends cleanly

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
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open cache file: " + file_name);
    }

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

// Function: calculate_total_kinetic_energy
// Purpose: Calculate the total KE of all particles in the system
high_prec Engine::calculate_total_kinetic_energy(std::shared_ptr<Particles> particles) {
    high_prec KE = 0;
    for (auto& p : particles->particle_list) {
        KE += 0.5 * p->m * (p->vx * p->vx + p->vy * p->vy);
    }
    return KE;
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

void Engine::correct_global_momentum(shared_ptr<Particles> particles, momentum initial_mom) {
    momentum current_mom = calc_mom(particles);
    high_prec mom_error_x = current_mom.x - initial_mom.x;
    high_prec mom_error_y = current_mom.y - initial_mom.y;
    
    
    if (fabs(mom_error_x) > 1e-12 || fabs(mom_error_y) > 1e-12) {
        cout << "APPLYING GLOBAL MOMENTUM CORRECTION: dx=" << mom_error_x << ", dy=" << mom_error_y << endl;
        
        high_prec total_mass = 0;
        for (auto& p : particles->particle_list) {
            total_mass += p->m;
        }
        
        // Apply small correction to all particles proportionally
        for (auto& p : particles->particle_list) {
            p->vx -= mom_error_x / total_mass;
            p->vy -= mom_error_y / total_mass;
        }
        
        // Verify correction worked
        momentum corrected_mom = calc_mom(particles);
        
    }
}