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
const double G = 6.674 * pow(10, -11); //m^3 kg^-1 s^-2

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
static double s_prev_TE = -1;
static double s_prev_step = -1;
static double s_margin_TE_error = 0.0;

static double s_margin_TE_error_overlap  = 0.0;
static double s_margin_TE_error_collision = 0.0;
static double s_margin_TE_error_integrate = 0.0;





//define static member variables
double Engine::dt = 1; //time step

// Constructor
Engine::Engine() {
    cout << "Engine initialized." << endl;
}

// Destructor
Engine::~Engine() {
    cout << "Engine destroyed." << endl;
}

double Engine::get_margin_TE_error() {
    return s_margin_TE_error;
}

double Engine::get_margin_TE_error_overlap()   { return s_margin_TE_error_overlap; }
double Engine::get_margin_TE_error_collision() { return s_margin_TE_error_collision; }
double Engine::get_margin_TE_error_integrate() { return s_margin_TE_error_integrate; }

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
    double steps_db = scenario->time / scenario->dt;
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

        double te_pre_update = calc_TE(particles);
        momentum mom_pre_update = calc_mom(particles);

        // Call update_particles
        update_particles(particles);

        double te_post_update = calc_TE(particles);
        momentum mom_post_update = calc_mom(particles);

        // Copy final state for this step
        auto particles_copy = make_unique<Particles>(*particles);
        particle_states->snaps.push_back(move(particles_copy));

        // Store the overall margin error
        double margin_TE_error = get_margin_TE_error();
        particle_states->metrics[i]->margin_TE_error = margin_TE_error;

        // Also store the three substep marginal errors
        double margin_overlap   = get_margin_TE_error_overlap();
        double margin_collision = get_margin_TE_error_collision();
        double margin_integrate = get_margin_TE_error_integrate();
        particle_states->metrics[i]->margin_TE_error_overlap   = margin_overlap;
        particle_states->metrics[i]->margin_TE_error_collision = margin_collision;
        particle_states->metrics[i]->margin_TE_error_integrate = margin_integrate;

        double te_error_update = (te_post_update - te_pre_update) / te_pre_update;
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

    run_to_cache(scenario, particle_states);
    return particle_states;
}




void Engine::update_particles(shared_ptr<Particles> particles)
{
    bool no_overlap = false;
    double te_pre_overlap = calc_TE(particles);

    // Overlap resolution
    for (int i = 0; i < 8; i++) {
        if (no_overlap) break;
        no_overlap = true;
        overlap_iter = i;
        no_overlap = resolve_overlap(particles);
    }
    double te_post_overlap = calc_TE(particles);
    double te_error_overlap = (te_post_overlap - te_pre_overlap) / te_pre_overlap;
    if (te_error_overlap > 0.05) {
        cout << "TE error overlap: " << te_error_overlap << endl;
    }
    // Store substep margin for overlap
    s_margin_TE_error_overlap = te_post_overlap - te_pre_overlap;

    // Collisions
    double te_pre_collision = calc_TE(particles);
    resolve_collisions(particles);
    double te_post_collision = calc_TE(particles);
    double te_error_collision = (te_post_collision - te_pre_collision) / te_pre_collision;
    if (te_error_collision > 0.05) {
        cout << "TE error collision: " << te_error_collision << endl;
    }
    // Store substep margin for collision
    s_margin_TE_error_collision = te_post_collision - te_pre_collision;

    // Gravity (Velocity Verlet)
    double te_pre_verlet = calc_TE(particles);
    resolve_gravity_verlet(particles);
    double te_post_verlet = calc_TE(particles);
    double te_error_verlet = (te_post_verlet - te_pre_verlet) / te_pre_verlet;
    if (te_error_verlet > 0.05) {
        cout << "TE error verlet: " << te_error_verlet << endl;
    }
    // Store substep margin for integration
    s_margin_TE_error_integrate = te_post_verlet - te_pre_verlet;

    // Now final TE for overall margin
    double current_TE_post = calc_TE(particles);
    if (s_prev_step >= 0 && update_iter == s_prev_step + 1) {
        s_margin_TE_error = current_TE_post - s_prev_TE;
    } else {
        s_margin_TE_error = 0.0;
    }
    s_prev_TE   = current_TE_post;
    s_prev_step = update_iter;
}

bool Engine::resolve_overlap(shared_ptr<Particles> particles) {

    //cout << "Resolving overlap between particles triggered" << endl;
    bool no_overlap = true;
    high_prec tTE_pre = calc_TE(particles);  // Accumulates total energy error for all particle pairs
    //momentum pre_mom = calc_mom(particles);  // Accumulates total momentum error for all particle pairs

    high_prec error_threshold = 0.001;  // Energy and momentum error threshold

    int pair_count = 0;  // To track the number of overlapping particle pairs

    // Loop through particle pairs
    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            if (check_collission(particles->particle_list[i], particles->particle_list[j], 0.0001)) {

                //no_overlap = false;

                pair_count++;

                // Call the core function to resolve overlap between particles i and j

                double TE_pre_overlap = calc_TE(particles);
                
                resolve_overlap_ij(particles->particle_list[i], particles->particle_list[j]);

                double TE_post_overlap = calc_TE(particles);
                double TE_error_overlap = (TE_post_overlap - TE_pre_overlap) / TE_pre_overlap;

                //cout << "------------------------------------------" << endl;
                //cout << "For particle pair: " << particles->particle_list[i]->particle_id << " and " << particles->particle_list[j]->particle_id << endl;
                //cout << "TE error overlap: " << TE_post_overlap - TE_pre_overlap << endl;
                //cout << "TE error overlap: " << TE_error_overlap*100 << "%" << endl;
                //cout << "------------------------------------------" << endl;

                //pause if particle_id is anything but 1 and 2
                if (particles->particle_list[i]->particle_id != 1 && particles->particle_list[i]->particle_id != 2) {
                    //cout << "Pausing for particle " << particles->particle_list[i]->particle_id << endl;
                    //cin.get();
                }

                //check if the particles are still overlapping
                if (check_collission(particles->particle_list[i], particles->particle_list[j], 0.0001)) {
                    

                    no_overlap = false;


                } else {
                    //cout << "Particles are not overlapping." << endl;
                }
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

    // Calculate initial total kinetic energy using calc_TE_ij
    high_prec TE_pre_ij = calc_TE_ij(particle_i, particle_j);

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
        dx = 0.001 * (rand() / (double)RAND_MAX - 0.5);
        dy = 0.001 * (rand() / (double)RAND_MAX - 0.5);
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
    particle_i->x = (x1_before - nx * d1).convert_to<double>();
    particle_i->y = (y1_before - ny * d1).convert_to<double>();
    particle_j->x = (x2_before + nx * d2).convert_to<double>();
    particle_j->y = (y2_before + ny * d2).convert_to<double>();

    // --------------------------------------------------------------------
    // Remove the "perfectly elastic impulse" block, because we rely on
    // resolve_energy_gap(...) to do the final energy/momentum correction.
    // --------------------------------------------------------------------

    // Re-check distances
    high_prec dx_new = particle_j->x - particle_i->x;
    high_prec dy_new = particle_j->y - particle_i->y;
    high_prec distance_new = hypot(dx_new, dy_new);

    // Overlap post separation (debug only)
    high_prec overlap_post = (particle_i->rad + particle_j->rad) - distance_new;
    //cout << "overlap_post in resolve_overlap_ij(): " << overlap_post << endl;

    // Compute TE, momentum again
    high_prec TE_post_ij = calc_TE_ij(particle_i, particle_j);
    momentum mom_post_ij = calc_mom_ij(particle_i, particle_j);

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
    high_prec rel_vel_post = (particle_j->vx - particle_i->vx) * dx_new +
                             (particle_j->vy - particle_i->vy) * dy_new;
    
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
        //cout << "Energy gap detected: " << delta_E << endl;
        resolve_energy_gap(particle_i, particle_j, delta_E);

        


    }





    high_prec distance_post_corrected = hypot(dx_new, dy_new);
    high_prec overlap_post_corrected = (particle_i->rad + particle_j->rad) - distance_post_corrected;
    high_prec rel_vel_post_corrected = (particle_j->vx - particle_i->vx) * dx_new +
                                       (particle_j->vy - particle_i->vy) * dy_new;
    
    high_prec min_restitution = min(particle_i->rest, particle_j->rest);
}



void Engine::resolve_energy_gap(std::shared_ptr<Particle> p1,
    std::shared_ptr<Particle> p2,
    boost::multiprecision::cpp_dec_float_50 delta_E)
{
    using high_prec = boost::multiprecision::cpp_dec_float_50;

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

    cout << "Available normal KE: " << KE_rel_n_initial << endl;
    cout << "Available tangential KE: " << KE_rel_t_initial << endl;

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

    p1->vx -= dvx1_n.convert_to<double>();
    p1->vy -= dvy1_n.convert_to<double>();
    p2->vx += dvx2_n.convert_to<double>();
    p2->vy += dvy2_n.convert_to<double>();

    // If there's remaining energy to remove, take it from the tangential component
    high_prec KE_rel_t = KE_rel_t_initial;

    cout << "delta_E: " << delta_E << endl;
    

    if (delta_E > 0) {
        if (KE_rel_t < delta_E) {
            std::cout << "WARNING: Not enough tangential kinetic energy either!" << std::endl;
            std::cout << "KE_rel_t (before) = " << KE_rel_t << std::endl;
            std::cout << "Remaining delta_E = " << delta_E << std::endl;
            std::cout << "Velocity components: v1t = " << v1t << ", v2t = " << v2t << std::endl;
            std::cout << "Removing as much as possible..." << std::endl;

            delta_E_used += KE_rel_t;
            delta_E -= KE_rel_t;
            KE_rel_t = 0;

            std::cout << "Press Enter to continue..." << std::endl;
            std::cin.get();
        } else {
            KE_rel_t -= delta_E;
            delta_E_used += delta_E;
            delta_E = 0;
        }

        // Adjust tangential velocity
        high_prec new_rel_vel_t = sqrt(std::max(high_prec(0.0), high_prec(2 * KE_rel_t / reduced_mass)));
        if (v2t - v1t < 0) new_rel_vel_t = -new_rel_vel_t;
        high_prec delta_rel_vel_t = new_rel_vel_t - (v2t - v1t);

        

        // Compute impulse for tangential adjustment
        high_prec impulse_t = delta_rel_vel_t * reduced_mass;
        high_prec dvx1_t = (impulse_t / m1) * tx;
        high_prec dvy1_t = (impulse_t / m1) * ty;
        high_prec dvx2_t = (impulse_t / m2) * tx;
        high_prec dvy2_t = (impulse_t / m2) * ty;

        p1->vx -= dvx1_t.convert_to<double>();
        p1->vy -= dvy1_t.convert_to<double>();
        p2->vx += dvx2_t.convert_to<double>();
        p2->vy += dvy2_t.convert_to<double>();
    }

    cout << "tangential KE after it is used: " << KE_rel_t << endl;
    cout << "normal KE after it is used: " << KE_rel_n << endl;
    cout << "delta_E after it is resolved: " << delta_E << endl;

    // ===================== FINAL CHECK =====================
    v1n = p1->vx * nx + p1->vy * ny;
    v2n = p2->vx * nx + p2->vy * ny;
    high_prec KE_rel_n_final = 0.5 * reduced_mass * (v2n - v1n) * (v2n - v1n);

    v1t = p1->vx * tx + p1->vy * ty;
    v2t = p2->vx * tx + p2->vy * ty;
    high_prec KE_rel_t_final = 0.5 * reduced_mass * (v2t - v1t) * (v2t - v1t);

    

    if (delta_E > 0) {
        std::cout << "WARNING: Energy gap not fully resolved!" << std::endl;
        std::cout << "delta_E = " << delta_E << std::endl;
    }
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

            bool collission = check_collission(particles->particle_list[i], particles->particle_list[j], -0.000001);
            
            
            

            if (collission) {
                //resolve the collission

                resolve_collission(particles->particle_list[i], particles->particle_list[j]);

                
            }
            

        }
    }

    
    
    
}


bool Engine::check_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2, double threshold) {
    //this function will check if two particles are colliding

    //1. calculate the distance between the particles
    double distance = hypot(particle1->x - particle2->x, particle1->y - particle2->y);
    
    


    //2. check if the distance+threshold is smaller than the sum of the radii of the particles & relative velocity is 
    if ((distance + threshold)< particle1->rad + particle2->rad ) {

        //calculate relative velocity 
        high_prec relative_velocity = hypot(particle1->vx - particle2->vx, particle1->vy - particle2->vy);
       // high_prec velocity_threshold = compute_velocity_threshold(particle1, particle2, dt);



        //if (relative_velocity < velocity_threshold) {

            //cout << "Relative velocity is too low, assuming no collision" << endl;
            //cout << "Relative velocity: " << relative_velocity << endl;
            //cout << "Velocity threshold: " << velocity_threshold << endl;

            

            

            

            //return false;
        //} else {
            //cout << "relative velocity above threshold, assuming collision:" << velocity_threshold << endl;
            //cout << "relative velocity: " << relative_velocity << endl;

            //check if threshold is negative
            //if (threshold < 0) {
                //cout << "Threshold is negative" << endl;

            //}
        //}

        //print high_prec values for distance, radius sum and threshold

        high_prec distance_high_prec = distance;
        high_prec radius_sum_high_prec = particle1->rad + particle2->rad;
        high_prec threshold_high_prec = threshold;

        // cout << "Distance: " << distance_high_prec << endl;
        // cout << "Radius sum: " << radius_sum_high_prec << endl;
        // cout << "Threshold: " << threshold_high_prec << endl;

        // //print relative velocity, using hypot to calculate the magnitude
        
        // cout << "Relative velocity: " << relative_velocity << endl;




        return true;
    }
    else {
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

    if (combined_rest == 0) {
        // Calculate common velocity after collision (conservation of momentum)
        high_prec v_final = (particle1->m * v1n + particle2->m * v2n) / total_mass;
        v1n_new = v_final;
        v2n_new = v_final;
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
    particle1->vx = (v1n_vec.x + v1t_vec.x).convert_to<double>();
    particle1->vy = (v1n_vec.y + v1t_vec.y).convert_to<double>();
    particle2->vx = (v2n_vec.x + v2t_vec.x).convert_to<double>();
    particle2->vy = (v2n_vec.y + v2t_vec.y).convert_to<double>();

    //convert the energy loss to heat

    // Calculate the energy loss due to restitution
    high_prec energy_loss = initial_kinetic_energy - final_kinetic_energy;
    //cout << "Energy loss: " << energy_loss << endl;

    



    


    //convert the energy loss to heat
    heat_ij(energy_loss, particle1, particle2);

        

    

    //calculate TE after the collision
    high_prec TE_post = calc_TE_ij(particle1, particle2);

    //print the TE before and after the collision
    //cout << "Collision TE error: " << TE_post - TE_pre << endl;
}


void Engine::resolve_gravity_euler(shared_ptr<Particles> particles) {//abandoned technique
    // This function resolves the gravitational attraction between particles using Euler's method. Not sufficiently accurate for energy conservation.
    

    // Define a small value to avoid division by zero

    high_prec epsilon = 0.0000001;

    // Loop through the particles
    for (int i = 0; i < particles->particle_list.size(); i++) {
        for (int j = i + 1; j < particles->particle_list.size(); j++) {
            // Skip self-interaction
            if (i == j) continue;

            //calculate TE before the update
            //high_prec TE_pre = calc_TE_ij(particles->particle_list[i], particles->particle_list[j]);
            

            // Calculate the distance between particles i and j
            high_prec dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            high_prec dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            high_prec distance = hypot(dx, dy);

            if (distance < epsilon) {
                
                distance = epsilon;
            }

            // Calculate the gravitational force magnitude
            high_prec force = G * particles->particle_list[i]->m * particles->particle_list[j]->m / (distance * distance);

            // Calculate the components of the force
            high_prec fx = force * (dx / distance);
            high_prec fy = force * (dy / distance);

            // Calculate KE, PE, and TE for the system of i and j before the force update, all are magnitudes
            high_prec KE_pre = 0.5 * particles->particle_list[i]->m * (pow(particles->particle_list[i]->vx, 2) + pow(particles->particle_list[i]->vy, 2)) +
                               0.5 * particles->particle_list[j]->m * (pow(particles->particle_list[j]->vx, 2) + pow(particles->particle_list[j]->vy, 2));
            high_prec PE_pre = -G * particles->particle_list[i]->m * particles->particle_list[j]->m / distance;
            high_prec TE_pre = KE_pre + PE_pre;



            // Update the velocities of both particles to conserve momentum
            particles->particle_list[i]->vx += (fx / particles->particle_list[i]->m).convert_to<double>() * dt;
            particles->particle_list[i]->vy += (fy / particles->particle_list[i]->m).convert_to<double>() * dt;

            particles->particle_list[j]->vx -= (fx / particles->particle_list[j]->m).convert_to<double>() * dt;
            particles->particle_list[j]->vy -= (fy / particles->particle_list[j]->m).convert_to<double>() * dt;

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
    double epsilon = 0.001; // To avoid division by zero
    int n = particles->particle_list.size();

    double TE_pre = calc_TE(particles); //calc total energy before the update


    // Step 1: Initialize net forces for all particles
    vector<Vector2D_double> net_force(n, {0.0, 0.0});


    
    // Step 2: Calculate and accumulate gravitational forces for all pairs
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            
            // Calculate distance components between particles i and j
            double dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            double dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            double distance = hypot(dx, dy);
            //if (distance < epsilon) distance = epsilon;
            distance = sqrt(distance*distance + epsilon*epsilon);

            // Calculate gravitational force magnitude
            double force = G * particles->particle_list[i]->m * particles->particle_list[j]->m / (distance * distance);

            // Force components
            double fx = force * (dx / distance);
            double fy = force * (dy / distance);

            // Accumulate forces for particles i and j (Newton's third law)
            net_force[i].x += fx;
            net_force[i].y += fy;
            net_force[j].x -= fx;
            net_force[j].y -= fy;
        }
    }

    // Step 3: Update velocities using the accumulated forces (first half-step)
    for (int i = 0; i < n; i++) {
        double dt_i = dt * 1;
        particles->particle_list[i]->vx += 0.5 * (net_force[i].x / particles->particle_list[i]->m) * dt_i;
        particles->particle_list[i]->vy += 0.5 * (net_force[i].y / particles->particle_list[i]->m) * dt_i;
    }

    // Step 4: Update positions using the updated velocities
    for (int i = 0; i < n; i++) {
        double dt_i = dt * 1;
        particles->particle_list[i]->x += particles->particle_list[i]->vx * dt_i;
        particles->particle_list[i]->y += particles->particle_list[i]->vy * dt_i;
    }


    // Step 5: Recalculate gravitational forces for the updated positions (second half-step)
    fill(net_force.begin(), net_force.end(), Vector2D_double{0.0, 0.0});
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // Recalculate distance components between particles i and j
            double dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            double dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            double distance = hypot(dx, dy);
            if (distance < epsilon) distance = epsilon;

            // Recalculate gravitational force magnitude
            double force = G * particles->particle_list[i]->m * particles->particle_list[j]->m / (distance * distance);

            // Force components
            double fx = force * (dx / distance);
            double fy = force * (dy / distance);

            // Accumulate forces for particles i and j
            net_force[i].x += fx;
            net_force[i].y += fy;
            net_force[j].x -= fx;
            net_force[j].y -= fy;
        }
    }

    // Step 6: Update velocities again using the new forces (second half-step)
    for (int i = 0; i < n; i++) {
        double dt_i = dt * 1;
        particles->particle_list[i]->vx += 0.5 * (net_force[i].x / particles->particle_list[i]->m) * dt_i;
        particles->particle_list[i]->vy += 0.5 * (net_force[i].y / particles->particle_list[i]->m) * dt_i;
    }

    // Step 7: validate conservation of energy
    double TE_post = calc_TE(particles); //calc total energy after the update

    // Calculate the difference in total energy
    double TE_diff = TE_post - TE_pre;
    double TE_error = TE_diff / TE_pre;
    double TE_error_threshold = 0.001; // Lower threshold for higher precision

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

    if (true){
        return;
    }
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
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->r) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->g) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->b) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->x) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->y) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->z) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->vx) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->vy) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->vz) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->m) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->rad) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->temp) << ","
                << static_cast<double>(particle_states->snaps[i]->particle_list[j]->rest) << endl;
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

    cout << "Snapshots successfully loaded from cache: " << file_name << std::endl;

    return particle_states;
}


// Function to calculate the total kinetic energy (KE) and potential energy (PE) for a system of particles
double Engine::calc_TE(shared_ptr<Particles> particles) {
    int n = particles->particle_list.size();
    double KE = 0.0;
    double PE = 0.0;
    double HE = 0.0;  // Heating energy component
    double epsilon = 0.001; // To avoid division by zero

    // Calculate kinetic energy and heating energy
    for (int i = 0; i < n; i++) {
        KE += 0.5 * particles->particle_list[i]->m * 
              (pow(particles->particle_list[i]->vx, 2) + pow(particles->particle_list[i]->vy, 2));
        HE += particles->particle_list[i]->temp; // Add heating energy of each particle
    }

    // Calculate potential energy
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            double dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            double distance = hypot(dx, dy);
            distance = sqrt(distance * distance + epsilon * epsilon);
            PE += -G * particles->particle_list[i]->m * particles->particle_list[j]->m / distance;
        }
    }

    // Return the total energy (KE + PE + HE)
    return KE + PE + HE;
}

double Engine::calc_TE_ij(shared_ptr<Particle> p1, shared_ptr<Particle> p2) {
    double KE = 0.0;
    double PE = 0.0;
    double HE = 0.0;  // Heating energy component
    double epsilon = 0.001; // To avoid division by zero

    // Calculate kinetic energy for both particles
    KE += 0.5 * p1->m * (pow(p1->vx, 2) + pow(p1->vy, 2));
    KE += 0.5 * p2->m * (pow(p2->vx, 2) + pow(p2->vy, 2));

    // Include heating energy for each particle
    HE += p1->temp + p2->temp;

    // Calculate potential energy between the two particles
    double dx = p2->x - p1->x;
    double dy = p2->y - p1->y;
    double distance = hypot(dx, dy);
    distance = sqrt(distance * distance + epsilon * epsilon);
    PE += -G * p1->m * p2->m / distance;

    // Return the total energy (KE + PE + HE)
    return KE + PE + HE;
}

momentum Engine::calc_mom(shared_ptr<Particles> particles) {
    double total_momentum_x = 0.0;
    double total_momentum_y = 0.0;

    // Sum up the x and y momentum components for each particle
    for (const auto& particle : particles->particle_list) {
        total_momentum_x += particle->m * particle->vx;
        total_momentum_y += particle->m * particle->vy;
    }

    return {total_momentum_x, total_momentum_y};
}

// Corrected `calc_mom_ij` for two particles
momentum Engine::calc_mom_ij(shared_ptr<Particle> p1, shared_ptr<Particle> p2) {
    double total_momentum_x = p1->m * p1->vx + p2->m * p2->vx;
    double total_momentum_y = p1->m * p1->vy + p2->m * p2->vy;

    return {total_momentum_x, total_momentum_y};
}

