#ifndef PHYSENGINE_H
#define PHYSENGINE_H

// Standard libraries
#include <memory>
#include <vector>
#include <string>

// Internal libraries
#include "InitStructs.h"
#include "MathUtils.h"

// Using declarations
using namespace std;
using high_prec = boost::multiprecision::cpp_dec_float_50;

class Engine {
public:
    // Momentum struct
    struct momentum {
        high_prec x, y;
    };

    // Static member variable
    static high_prec dt;

    // Constructor and Destructor
    Engine();
    ~Engine();

    // Main simulation functions
    shared_ptr<snapshots> run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles);
    void update_particles(shared_ptr<Particles> particles);
    bool resolve_overlaps(shared_ptr<Particles> particles);
    void resolve_overlap_ij(shared_ptr<Particle> particle_i, shared_ptr<Particle> particle_j);
    void resolve_collisions(shared_ptr<Particles> particles);
    void resolve_collision(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    void resolve_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    void resolve_gravity_verlet(shared_ptr<Particles> particles);
    void update_locations(shared_ptr<Particles> particles, shared_ptr<backed_scaler> scaler);

    // Energy and momentum calculations
    high_prec calc_TE(shared_ptr<Particles> particles);
    high_prec calc_TE_ij(shared_ptr<Particle> p1, shared_ptr<Particle> p2, bool verbose = false);
    momentum calc_mom(shared_ptr<Particles> particles);
    momentum calc_mom_ij(shared_ptr<Particle> p1, shared_ptr<Particle> p2);
    high_prec calculate_total_kinetic_energy(std::shared_ptr<Particles> particles);

    // Collision detection
    bool check_overlap(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    bool check_collision(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    bool check_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);

    // Energy management
    void resolve_energy_gap(shared_ptr<Particle> p1, shared_ptr<Particle> p2, high_prec delta_E);
    void apply_energy_correction(shared_ptr<Particle> p1, shared_ptr<Particle> p2, high_prec delta_E);
    high_prec heat_ij(high_prec E, shared_ptr<Particle> particle_i, shared_ptr<Particle> particle_j);

    // Cache operations
    void run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states);
    bool cache_exists(shared_ptr<scenario> scenario);
    shared_ptr<snapshots> run_from_cache(shared_ptr<scenario> scenario);


    // In the public section of Engine class in PhysEngine.h, add:
    void check_momentum_metrics(shared_ptr<Particles> particles, momentum initial_mom);

    // Metric getters
    high_prec get_margin_TE_error();
    high_prec get_margin_TE_error_overlap();
    high_prec get_margin_TE_error_collision();
    high_prec get_margin_TE_error_integrate();
    high_prec get_margin_TE_error_overlap_ij_transl();
    high_prec get_margin_TE_error_overlap_ij_corrected();

private:
    // Helper functions for energy calculations
    high_prec calculate_kinetic_energy(shared_ptr<Particle> p1, shared_ptr<Particle> p2);
    high_prec calculate_potential_energy(shared_ptr<Particle> p1, shared_ptr<Particle> p2);
    high_prec calculate_heating_energy(shared_ptr<Particle> p1, shared_ptr<Particle> p2);

    high_prec calculate_max_acceleration(shared_ptr<Particles> particles) {
        high_prec max_accel = 0;
        for (int i = 0; i < particles->particle_list.size(); i++) {
            high_prec accel_magnitude = hypot(particles->particle_list[i]->vx, particles->particle_list[i]->vy) / dt;
            if (accel_magnitude > max_accel) {
                max_accel = accel_magnitude;
            }
        }
        return max_accel;
    }

    high_prec calculate_system_scale(shared_ptr<Particles> particles) {
        high_prec max_distance = 0;
        for (int i = 0; i < particles->particle_list.size(); i++) {
            for (int j = i + 1; j < particles->particle_list.size(); j++) {
                high_prec dist = hypot(particles->particle_list[i]->x - particles->particle_list[j]->x,
                                    particles->particle_list[i]->y - particles->particle_list[j]->y);
                max_distance = max(max_distance, dist);
            }
        }
        return max_distance > 0 ? max_distance : 1.0;
    }

    high_prec calculate_overlap_amount(shared_ptr<Particle> p1, shared_ptr<Particle> p2) {
        high_prec dist = hypot(p1->x - p2->x, p1->y - p2->y);
        high_prec sum_radii = p1->rad + p2->rad;
        return max(high_prec(0.0), sum_radii - dist);
    }

    // ADD THIS METHOD DECLARATION:
    void correct_global_momentum(shared_ptr<Particles> particles, momentum initial_mom);

};

#endif // PHYSENGINE_H