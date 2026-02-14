// ===============================
// PhysEngineCore.cpp  (PhysEngineCore SET)
// ===============================
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdlib>

#include "../include/PhysEngineCore.h"

using std::cout;
using std::endl;

using high_prec = EngineCore::high_prec;

// PHYSICS CONSTANTS
static const high_prec kGravG = 6.674e-11;

// TIMESTEP CONSTANTS
static const high_prec kAdaptiveDtAccelEps = 1e-6;

//iteration
int EngineCore::update_iter = 0;

// COLLISION CONSTANTS
high_prec EngineCore::collision_distance_tolerance_ = 1e-5;
//static const high_prec kContactBiasBeta = 0.5;     // larger = stronger

//to set kContactBiasBeta based on a value in scenario object
static high_prec kContactBiasBeta = 0.1;

static const high_prec kContactBiasSlop = 0.0;       // tolerance for small overlaps

// VELOCITY VERLET CONSTANTS
static const high_prec kGravitySofteningEps = 1e-6;
static const double    kGravityForceClipAbs = 1e100;

// ENERGY CALCULATION CONSTANTS
static const high_prec kEnergySofteningEps = 1e-3;


// =======================
// LOCAL HELPERS
// =======================
static inline double safe_double_from_high_prec(const high_prec& v) {
    return v;
}

static inline high_prec abs_hp(high_prec v) { return v < 0.0 ? -v : v; }

static inline high_prec safe_rel_error(high_prec post, high_prec pre) {
    high_prec diff = post - pre;
    if (pre == 0.0) return diff;
    return diff / pre;
}

// =======================
// HOT-LOOP HELPERS (MINIMAL, NO API CHANGES)
// =======================
static inline bool check_collission_ptr(const Particle* particle1, const Particle* particle2) {
    high_prec dx = particle2->x - particle1->x;
    high_prec dy = particle2->y - particle1->y;

    high_prec dist2 = dx * dx + dy * dy;
    high_prec sum_radii = particle1->rad + particle2->rad;

    // Squared-distance early reject:
    // old: (distance - sum_radii) < tol  <=> distance < sum_radii + tol
    const high_prec thresh = sum_radii + EngineCore::collision_distance_tolerance_;
    if (!(dist2 < thresh * thresh)) return false;

    if (dist2 == 0.0) return false;

    high_prec distance = sqrt(dist2);

    high_prec nx = dx / distance;
    high_prec ny = dy / distance;

    high_prec rvx = particle2->vx - particle1->vx;
    high_prec rvy = particle2->vy - particle1->vy;
    high_prec v_rel_n = rvx * nx + rvy * ny;

    return (v_rel_n < 0.0);
}

static inline void resolve_collission_ptr(Particle* particle1, Particle* particle2) {
    high_prec KE_pre = 0.5L * particle1->m * (particle1->vx * particle1->vx + particle1->vy * particle1->vy)
                     + 0.5L * particle2->m * (particle2->vx * particle2->vx + particle2->vy * particle2->vy);

    high_prec dx = particle2->x - particle1->x;
    high_prec dy = particle2->y - particle1->y;
    high_prec distance = sqrt(dx * dx + dy * dy);
    if (distance == 0.0) return;

    high_prec nx = dx / distance;
    high_prec ny = dy / distance;

    high_prec tx = -ny;
    high_prec ty = nx;

    high_prec v1n = particle1->vx * nx + particle1->vy * ny;
    high_prec v1t = particle1->vx * tx + particle1->vy * ty;
    high_prec v2n = particle2->vx * nx + particle2->vy * ny;
    high_prec v2t = particle2->vx * tx + particle2->vy * ty;

    high_prec combined_rest = std::min(particle1->rest, particle2->rest);

    high_prec m1 = particle1->m;
    high_prec m2 = particle2->m;
    high_prec total_mass = m1 + m2;

    // compute relative normal speed (should be < 0 due to check_collision)
    high_prec v_rel_n = (v2n - v1n);
    if (v_rel_n >= 0.0) return;

    // penetration-based bias velocity (stronger), only applied because v_rel_n < 0
    high_prec bias_vel = 0.0;
    if (EngineCore::dt != 0.0) {
        high_prec penetration = (particle1->rad + particle2->rad) - distance;
        if (penetration > kContactBiasSlop) {
            bias_vel = kContactBiasBeta * (penetration - kContactBiasSlop) / EngineCore::dt;
        }
    }

    // baseline target separating relative speed from restitution
    // (v_rel_post = -e * v_rel_pre)
    high_prec base_target = -combined_rest * v_rel_n; // >= 0

    // enforce a minimum separation speed to unwind overlap
    high_prec v_rel_target = std::max(base_target, bias_vel);

    // momentum conservation along normal, with specified relative speed
    high_prec P = m1 * v1n + m2 * v2n;
    high_prec v1n_new = (P - m2 * v_rel_target) / total_mass;
    high_prec v2n_new = v1n_new + v_rel_target;

    high_prec v1t_new = v1t;
    high_prec v2t_new = v2t;

    particle1->vx = v1n_new * nx + v1t_new * tx;
    particle1->vy = v1n_new * ny + v1t_new * ty;
    particle2->vx = v2n_new * nx + v2t_new * tx;
    particle2->vy = v2n_new * ny + v2t_new * ty;

    high_prec KE_post = 0.5L * particle1->m * (particle1->vx * particle1->vx + particle1->vy * particle1->vy)
                      + 0.5L * particle2->m * (particle2->vx * particle2->vx + particle2->vy * particle2->vy);

    high_prec energy_lost = KE_pre - KE_post;
    // Only add heat if energy_lost is significant (avoid numerical flicker)
    const high_prec EPSILON = 1e-10;
    if (std::abs(energy_lost) > EPSILON) {
        high_prec temp = energy_lost / 2.0;
        particle1->temp += temp;
        particle2->temp += temp;
    }
}


// =======================
// STATIC MEMBERS
// =======================
high_prec EngineCore::dt = 1;

long long EngineCore::collisions = 0;

high_prec EngineCore::margin_TE_error = 0;
high_prec EngineCore::margin_TE_error_collision = 0;
high_prec EngineCore::margin_TE_error_integrate = 0;


// =======================
// TIMING
// =======================
double EngineCore::seconds_between(const EngineCore::clock_t::time_point& start,
                                   const EngineCore::clock_t::time_point& end) {
    return std::chrono::duration<double>(end - start).count();
}

// =======================
// LIFECYCLE
// =======================
EngineCore::EngineCore() {
    cout << "EngineCore initialized." << endl;
}

// =======================
// METRIC GETTERS
// =======================
high_prec EngineCore::get_margin_TE_error() { return margin_TE_error; }
high_prec EngineCore::get_margin_TE_error_collision() { return margin_TE_error_collision; }
high_prec EngineCore::get_margin_TE_error_integrate() { return margin_TE_error_integrate; }


// =======================
// HELPERS (DT)
// =======================
high_prec EngineCore::calculate_max_acceleration(std::shared_ptr<Particles> particles) {
    // EXACT old semantics:
    // accel_magnitude = hypot(vx,vy) / dt
    // dt MUST already be set to scenario->dt before calling (configure_dt does that).
    if (dt == 0.0) return 0.0;

    high_prec max_accel = 0;
    for (int i = 0; i < (int)particles->particle_list.size(); i++) {
        high_prec vx = particles->particle_list[i]->vx;
        high_prec vy = particles->particle_list[i]->vy;
        high_prec speed = sqrt(vx * vx + vy * vy);
        high_prec accel_magnitude = speed / dt;
        if (accel_magnitude > max_accel) max_accel = accel_magnitude;
    }
    return max_accel;
}

// =======================
// ADAPTIVE DT
// =======================
void EngineCore::configure_dt(std::shared_ptr<scenario> scenario, std::shared_ptr<Particles> particles) {
    // 1) EXACT old first assignment
    dt = scenario->dt;

    // 2) EXACT old max_accel call (uses dt internally)
    high_prec max_accel = calculate_max_acceleration(particles);

    // 3) EXACT old adaptive dt computation and clamp
    //    adaptive_dt = sqrt(1e-6L / max_accel)
    //    dt = min(scenario->dt, adaptive_dt)
    if (max_accel != 0.0) {
        high_prec adaptive_dt = sqrt(kAdaptiveDtAccelEps / max_accel);

        const high_prec base_dt = scenario->dt;

        // IMPORTANT: avoid std::min on mixed types (this is what typically introduces functional differences)
        dt = (adaptive_dt < base_dt) ? adaptive_dt : base_dt;
    } else {
        // If max_accel == 0, old behavior effectively keeps dt = scenario->dt
        dt = scenario->dt;
    }
}

void EngineCore::set_overlap_beta(std::shared_ptr<scenario> scenario) {
    //to set kContactBiasBeta based on a value in scenario object
    kContactBiasBeta = scenario->beta;
}

void EngineCore::set_collision_distance_tolerance(std::shared_ptr<scenario> scenario) {
    collision_distance_tolerance_ = scenario->collision_distance_tolerance;
}

// =======================
// 1 CORE STEP
// =======================
void EngineCore::step(std::shared_ptr<Particles> particles, StepEnergyTrace* trace) {

    if (trace) {
        trace->TE1 = calc_TE(particles);
    }

    if (trace) {
        auto t0 = clock_t::now();
        resolve_collisions(particles);
        auto t1 = clock_t::now();
        trace->collision_seconds = seconds_between(t0, t1);
    } else {
        resolve_collisions(particles);
    }

    if (trace) {
        trace->TE2 = calc_TE(particles);
        trace->dE_collision  = trace->TE2 - trace->TE1;
        trace->rel_collision = safe_rel_error(trace->TE2, trace->TE1);
    }

    if (trace) {
        auto t0 = clock_t::now();
        resolve_gravity_verlet(particles);
        auto t1 = clock_t::now();
        trace->verlet_seconds = seconds_between(t0, t1);
    } else {
        resolve_gravity_verlet(particles);
    }

    if (trace) {
        trace->TE3 = calc_TE(particles);
        trace->dE_verlet  = trace->TE3 - trace->TE2;
        trace->rel_verlet = safe_rel_error(trace->TE3, trace->TE2);

        trace->KE = calculate_system_kinetic_energy(particles);
        trace->PE = calculate_system_potential_energy(particles);
        trace->HE = calculate_system_heating_energy(particles);
        trace->TE = trace->KE + trace->PE + trace->HE;
    }
}

// =======================
// ENERGY
// =======================
high_prec EngineCore::calculate_system_kinetic_energy(std::shared_ptr<Particles> particles) {
    high_prec KE = 0.0;
    for (const auto& p : particles->particle_list) {
        KE += 0.5 * p->m * (p->vx * p->vx + p->vy * p->vy);
    }
    return KE;
}

high_prec EngineCore::calculate_system_heating_energy(std::shared_ptr<Particles> particles) {
    high_prec HE = 0.0;
    for (const auto& p : particles->particle_list) {
        HE += static_cast<high_prec>(p->temp);
    }
    return HE;
}

high_prec EngineCore::calculate_system_potential_energy(std::shared_ptr<Particles> particles) {
    int n = (int)particles->particle_list.size();
    high_prec PE = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            high_prec dx = particles->particle_list[j]->x - particles->particle_list[i]->x;
            high_prec dy = particles->particle_list[j]->y - particles->particle_list[i]->y;
            high_prec distance = sqrt(dx * dx + dy * dy + kEnergySofteningEps * kEnergySofteningEps);
            PE += -kGravG * particles->particle_list[i]->m * particles->particle_list[j]->m / distance;
        }
    }
    return PE;
}

high_prec EngineCore::calc_TE(std::shared_ptr<Particles> particles) {
    high_prec KE = calculate_system_kinetic_energy(particles);
    high_prec PE = calculate_system_potential_energy(particles);
    high_prec HE = calculate_system_heating_energy(particles);
    return KE + PE + HE;
}

high_prec EngineCore::calculate_total_kinetic_energy(std::shared_ptr<Particles> particles) {
    return calculate_system_kinetic_energy(particles);
}

high_prec EngineCore::calculate_kinetic_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2) {
    high_prec KE = 0.0;
    KE += 0.5 * p1->m * (p1->vx * p1->vx + p1->vy * p1->vy);
    KE += 0.5 * p2->m * (p2->vx * p2->vx + p2->vy * p2->vy);
    return KE;
}

high_prec EngineCore::calculate_potential_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2) {
    high_prec dx = p2->x - p1->x;
    high_prec dy = p2->y - p1->y;
    high_prec distance = sqrt(dx * dx + dy * dy + kEnergySofteningEps * kEnergySofteningEps);
    return -kGravG * p1->m * p2->m / distance;
}

high_prec EngineCore::calculate_heating_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2) {
    return static_cast<high_prec>(p1->temp + p2->temp);
}

high_prec EngineCore::calc_TE_ij(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2, bool verbose) {
    (void)verbose;
    high_prec KE = calculate_kinetic_energy(p1, p2);
    high_prec PE = calculate_potential_energy(p1, p2);
    high_prec HE = calculate_heating_energy(p1, p2);
    return KE + PE + HE;
}

high_prec EngineCore::calculate_overlap_amount(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2) {
    high_prec dx = p1->x - p2->x;
    high_prec dy = p1->y - p2->y;
    high_prec dist = sqrt(dx * dx + dy * dy);
    high_prec sum_radii = p1->rad + p2->rad;
    return sum_radii - dist;
}

high_prec EngineCore::heat_ij(high_prec E, std::shared_ptr<Particle> particle_i, std::shared_ptr<Particle> particle_j) {
    high_prec heat = E * (particle_i->m + particle_j->m) / 2;
    high_prec temp = heat / (particle_i->m + particle_j->m);

    particle_i->temp += temp;
    particle_j->temp += temp;

    return E;
}

bool EngineCore::check_overlap(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2) {
    return calculate_overlap_amount(particle1, particle2) > 0.0;
}

bool EngineCore::check_collission(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2) {
    // keep public API behavior; use ptr helper internally
    return check_collission_ptr(particle1.get(), particle2.get());
}

void EngineCore::resolve_collisions(std::shared_ptr<Particles> particles) {
    // MINIMAL: avoid shared_ptr refcount churn in inner loops by using raw pointers
    const int n = (int)particles->particle_list.size();
    auto& list = particles->particle_list;

    for (int i = 0; i < n; i++) {
        Particle* pi = list[i].get();
        for (int j = i + 1; j < n; j++) {
            Particle* pj = list[j].get();
            if (check_collission_ptr(pi, pj)) {
                collisions++;
                // MINIMAL: resolve via ptr helper to avoid shared_ptr copies
                resolve_collission_ptr(pi, pj);
            }
        }
    }
}

void EngineCore::resolve_collission(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2) {
    // keep public API behavior; use ptr helper internally
    resolve_collission_ptr(particle1.get(), particle2.get());
}

// Local helper: compute all-pairs gravity forces into an output buffer
static void compute_gravity_forces(
    std::vector<Vector2D>& forces,
    std::vector<std::shared_ptr<Particle>>& list,
    int n)
{
    std::fill(forces.begin(), forces.end(), Vector2D{0.0, 0.0});
    for (int i = 0; i < n; ++i) {
        Particle* pi = list[i].get();
        for (int j = i + 1; j < n; ++j) {
            Particle* pj = list[j].get();

            high_prec dx = pj->x - pi->x;
            high_prec dy = pj->y - pi->y;

            high_prec dist2 = dx * dx + dy * dy + kGravitySofteningEps * kGravitySofteningEps;
            high_prec distance = sqrt(dist2);
            if (distance == 0) distance = kGravitySofteningEps;

            high_prec force = kGravG * pi->m * pj->m / dist2;

            high_prec fx = force * (dx / distance);
            high_prec fy = force * (dy / distance);

            forces[i].x += fx;
            forces[i].y += fy;
            forces[j].x -= fx;
            forces[j].y -= fy;
        }
    }
}

void EngineCore::resolve_gravity_verlet(std::shared_ptr<Particles> particles) {
    const int n = (int)particles->particle_list.size();
    auto& list = particles->particle_list;

    // Invalidate cache if particle count changed
    if ((int)cached_force_.size() != n) {
        cached_force_.assign(n, Vector2D{0.0, 0.0});
        forces_valid_ = false;
    }

    // --- First half-kick ---
    // On step 0, compute forces from scratch.  On subsequent steps, reuse
    // the post-drift forces cached at the end of the previous step.
    if (!forces_valid_) {
        compute_gravity_forces(cached_force_, list, n);
        forces_valid_ = true;
    }

    for (int i = 0; i < n; ++i) {
        Particle* p = list[i].get();
        p->vx += 0.5 * (cached_force_[i].x / p->m) * dt;
        p->vy += 0.5 * (cached_force_[i].y / p->m) * dt;
    }

    // --- Drift ---
    for (int i = 0; i < n; ++i) {
        Particle* p = list[i].get();
        p->x += p->vx * dt;
        p->y += p->vy * dt;
    }

    // --- Compute forces at new positions and cache for next step ---
    compute_gravity_forces(cached_force_, list, n);

    // --- Second half-kick ---
    for (int i = 0; i < n; ++i) {
        Particle* p = list[i].get();
        p->vx += 0.5 * (cached_force_[i].x / p->m) * dt;
        p->vy += 0.5 * (cached_force_[i].y / p->m) * dt;
    }
}
