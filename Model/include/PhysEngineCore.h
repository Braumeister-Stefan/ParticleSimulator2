// ===============================
// PhysEngineCore.h   (PhysEngineCore SET)
// ===============================
#ifndef PHYSENGINECORE_H
#define PHYSENGINECORE_H

#include <memory>
#include <vector>
#include <string>
#include <chrono>

#include "InitStructs.h"
#include "MathUtils.h"

class EngineCore {
public:
    using high_prec = double;

    // Single source-of-truth dt for the entire simulation
    static high_prec dt;

    struct StepEnergyTrace {
        high_prec TE0 = 0;
        high_prec TE1 = 0;
        high_prec TE2 = 0;
        high_prec TE3 = 0;

        high_prec dE_overlap   = 0;
        high_prec dE_collision = 0;
        high_prec dE_verlet    = 0;

        high_prec rel_overlap   = 0;
        high_prec rel_collision = 0;
        high_prec rel_verlet    = 0;

        high_prec KE = 0;
        high_prec PE = 0;
        high_prec HE = 0;
        high_prec TE = 0;

        // Stage timing (seconds)
        double collision_seconds = 0;
        double verlet_seconds    = 0;
    };

    EngineCore();
    virtual ~EngineCore() = default;

    // Adaptive dt (MUST match old behavior exactly)
    void configure_dt(std::shared_ptr<scenario> scenario, std::shared_ptr<Particles> particles);

    // One full physics step: overlaps -> collisions -> verlet (+ energies if trace)
    void step(std::shared_ptr<Particles> particles, StepEnergyTrace* trace = nullptr);

    // Metric getters (kept for wrapper)
    high_prec get_margin_TE_error();
    high_prec get_margin_TE_error_overlap();
    high_prec get_margin_TE_error_collision();
    high_prec get_margin_TE_error_integrate();
    high_prec get_margin_TE_error_overlap_ij_transl();
    high_prec get_margin_TE_error_overlap_ij_corrected();

    // Energy
    high_prec calc_TE(std::shared_ptr<Particles> particles);
    high_prec calc_TE_ij(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, bool verbose = false);
    high_prec calculate_total_kinetic_energy(std::shared_ptr<Particles> particles);

    high_prec calculate_system_kinetic_energy(std::shared_ptr<Particles> particles);
    high_prec calculate_system_potential_energy(std::shared_ptr<Particles> particles);
    high_prec calculate_system_heating_energy(std::shared_ptr<Particles> particles);

    // Overlaps / Collisions / Gravity
    bool resolve_overlaps(std::shared_ptr<Particles> particles);
    void resolve_overlap_ij(std::shared_ptr<Particle> particle_i, std::shared_ptr<Particle> particle_j);

    void resolve_collisions(std::shared_ptr<Particles> particles);
    void resolve_collision(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2);
    void resolve_collission(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2);

    void resolve_gravity_verlet(std::shared_ptr<Particles> particles);

    // Misc
    void update_locations(std::shared_ptr<Particles> particles, std::shared_ptr<backed_scaler> scaler);
    void set_overlap_beta(std::shared_ptr<scenario> scenario);
    void set_collision_distance_tolerance(std::shared_ptr<scenario> scenario);

    // Collision detection
    bool check_overlap(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2);
    bool check_collision(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2);
    bool check_collission(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2);

    // Energy management
    void resolve_energy_gap(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, high_prec delta_E);
    void apply_energy_correction(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, high_prec delta_E);
    high_prec heat_ij(high_prec E, std::shared_ptr<Particle> particle_i, std::shared_ptr<Particle> particle_j);

    // Timing helper
    using clock_t = std::chrono::high_resolution_clock;
    static double seconds_between(const clock_t::time_point& start, const clock_t::time_point& end);

    high_prec calculate_overlap_amount(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);

    // Public counters (referenced by wrapper build)
    static int overlap_iter;
    static int overlap_ij_iter;
    static int update_iter;

    static int energy_gap_corrections;
    static int energy_gap_corrections_incomplete;

    static long long collisions;

    static high_prec collision_distance_tolerance_;

private:
    // MUST match old calculate_max_acceleration semantics exactly:
    // accel_magnitude = hypot(vx,vy)/dt    (dt already set to scenario->dt before calling)
    high_prec calculate_max_acceleration(std::shared_ptr<Particles> particles);

    high_prec calculate_kinetic_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);
    high_prec calculate_potential_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);
    high_prec calculate_heating_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);

    // Rank 1: cached post-drift forces reused as next step's pre-drift forces
    std::vector<Vector2D> cached_force_;
    bool forces_valid_ = false;

private:
    static high_prec margin_TE_error;
    static high_prec margin_TE_error_overlap;
    static high_prec margin_TE_error_collision;
    static high_prec margin_TE_error_integrate;
    static high_prec margin_TE_error_overlap_ij_transl;
    static high_prec margin_TE_error_overlap_ij_corrected;
};

#endif // PHYSENGINECORE_H
