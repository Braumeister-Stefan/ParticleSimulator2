// ===============================
// PhysEngineCore.cpp  (PhysEngineCore SET)
// ===============================
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <cstdlib>
#include <chrono>
#include <memory>
#include <array>

#include "../include/PhysEngineCore.h"

using std::cout;
using std::endl;

using high_prec = EngineCore::high_prec;

// =======================
// NON-STRUCT PARAMETERS (TUNING KNOBS)
// =======================
// Barnes–Hut opening angle (smaller => more accurate, slower)
static const high_prec kBHTheta = 0.9;

// Smallest cell half-size allowed before we stop subdividing and use a bucket
static const high_prec kBHMinHalf = 1e-8;

// Small padding to avoid zero-size root bounds
static const high_prec kBHBoundsPad = 1e-9;

// =======================
// PHYSICS CONSTANTS
// =======================
static const high_prec kGravG = 6.674e-11;

// TIMESTEP CONSTANTS
static const high_prec kAdaptiveDtAccelEps = 1e-6;

// iteration
int EngineCore::update_iter = 0;

// COLLISION CONSTANTS
high_prec EngineCore::collision_distance_tolerance_ = 1e-5;

// to set kContactBiasBeta based on a value in scenario object
static high_prec kContactBiasBeta = 0.1;

static const high_prec kContactBiasSlop = 0.0;       // tolerance for small overlaps

// VELOCITY VERLET CONSTANTS
static const high_prec kGravitySofteningEps = 1e-6;
static const double    kGravityForceClipAbs = 1e100; // unused (kept)

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
// HOT-LOOP HELPERS
// =======================
static inline bool check_collission_ptr(const Particle* particle1, const Particle* particle2) {
    high_prec dx = particle2->x - particle1->x;
    high_prec dy = particle2->y - particle1->y;
    high_prec dz = particle2->z - particle1->z;

    high_prec dist2 = dx * dx + dy * dy + dz * dz;
    high_prec sum_radii = particle1->rad + particle2->rad;

    const high_prec thresh = sum_radii + EngineCore::collision_distance_tolerance_;
    if (!(dist2 < thresh * thresh)) return false;
    if (dist2 == 0.0) return false;

    high_prec distance = sqrt(dist2);

    high_prec nx = dx / distance;
    high_prec ny = dy / distance;
    high_prec nz = dz / distance;

    high_prec rvx = particle2->vx - particle1->vx;
    high_prec rvy = particle2->vy - particle1->vy;
    high_prec rvz = particle2->vz - particle1->vz;
    high_prec v_rel_n = rvx * nx + rvy * ny + rvz * nz;

    return (v_rel_n < 0.0);
}

static inline void resolve_collission_ptr(Particle* particle1, Particle* particle2) {
    high_prec KE_pre = 0.5L * particle1->m * (particle1->vx * particle1->vx + particle1->vy * particle1->vy + particle1->vz * particle1->vz)
                     + 0.5L * particle2->m * (particle2->vx * particle2->vx + particle2->vy * particle2->vy + particle2->vz * particle2->vz);

    high_prec dx = particle2->x - particle1->x;
    high_prec dy = particle2->y - particle1->y;
    high_prec dz = particle2->z - particle1->z;
    high_prec distance = sqrt(dx * dx + dy * dy + dz * dz);
    if (distance == 0.0) return;

    high_prec nx = dx / distance;
    high_prec ny = dy / distance;
    high_prec nz = dz / distance;

    // Normal-direction velocity components (scalar projections)
    high_prec v1n = particle1->vx * nx + particle1->vy * ny + particle1->vz * nz;
    high_prec v2n = particle2->vx * nx + particle2->vy * ny + particle2->vz * nz;

    // Tangential velocity components (full 3D vectors)
    high_prec v1tx = particle1->vx - v1n * nx;
    high_prec v1ty = particle1->vy - v1n * ny;
    high_prec v1tz = particle1->vz - v1n * nz;
    high_prec v2tx = particle2->vx - v2n * nx;
    high_prec v2ty = particle2->vy - v2n * ny;
    high_prec v2tz = particle2->vz - v2n * nz;

    high_prec combined_rest = std::min(particle1->rest, particle2->rest);

    high_prec m1 = particle1->m;
    high_prec m2 = particle2->m;
    high_prec total_mass = m1 + m2;

    high_prec v_rel_n = (v2n - v1n);
    if (v_rel_n >= 0.0) return;

    high_prec bias_vel = 0.0;
    if (EngineCore::dt != 0.0) {
        high_prec penetration = (particle1->rad + particle2->rad) - distance;
        if (penetration > kContactBiasSlop) {
            bias_vel = kContactBiasBeta * (penetration - kContactBiasSlop) / EngineCore::dt;
        }
    }

    high_prec base_target = -combined_rest * v_rel_n; // >= 0
    high_prec v_rel_target = std::max(base_target, bias_vel);

    high_prec P = m1 * v1n + m2 * v2n;
    high_prec v1n_new = (P - m2 * v_rel_target) / total_mass;
    high_prec v2n_new = v1n_new + v_rel_target;

    // Reconstruct full 3D velocities from new normal + unchanged tangential
    particle1->vx = v1n_new * nx + v1tx;
    particle1->vy = v1n_new * ny + v1ty;
    particle1->vz = v1n_new * nz + v1tz;
    particle2->vx = v2n_new * nx + v2tx;
    particle2->vy = v2n_new * ny + v2ty;
    particle2->vz = v2n_new * nz + v2tz;

    high_prec KE_post = 0.5L * particle1->m * (particle1->vx * particle1->vx + particle1->vy * particle1->vy + particle1->vz * particle1->vz)
                      + 0.5L * particle2->m * (particle2->vx * particle2->vx + particle2->vy * particle2->vy + particle2->vz * particle2->vz);

    high_prec energy_lost = KE_pre - KE_post;
    const high_prec EPSILON = 1e-10;
    if (std::abs(energy_lost) > EPSILON) {
        high_prec temp = energy_lost / 2.0;
        particle1->temp += temp;
        particle2->temp += temp;
    }
}

// =======================
// DIRECT ALL-PAIRS GRAVITY - DEPRECATED
// =======================
static void compute_gravity_forces(
    std::vector<Vector2D>& forces,
    std::vector<std::shared_ptr<Particle>>& list,
    int n,
    long long* out_pair_calcs)
{
    std::fill(forces.begin(), forces.end(), Vector2D{0.0, 0.0, 0.0});

    if (!out_pair_calcs) {
        for (int i = 0; i < n; ++i) {
            Particle* pi = list[i].get();
            for (int j = i + 1; j < n; ++j) {
                Particle* pj = list[j].get();

                high_prec dx = pj->x - pi->x;
                high_prec dy = pj->y - pi->y;
                high_prec dz = pj->z - pi->z;

                high_prec dist2 = dx * dx + dy * dy + dz * dz + kGravitySofteningEps * kGravitySofteningEps;
                high_prec distance = sqrt(dist2);
                if (distance == 0) distance = kGravitySofteningEps;

                high_prec force = kGravG * pi->m * pj->m / dist2;

                high_prec fx = force * (dx / distance);
                high_prec fy = force * (dy / distance);
                high_prec fz = force * (dz / distance);

                forces[i].x += fx;
                forces[i].y += fy;
                forces[i].z += fz;
                forces[j].x -= fx;
                forces[j].y -= fy;
                forces[j].z -= fz;
            }
        }
        return;
    }

    for (int i = 0; i < n; ++i) {
        Particle* pi = list[i].get();
        for (int j = i + 1; j < n; ++j) {
            Particle* pj = list[j].get();

            ++(*out_pair_calcs);

            high_prec dx = pj->x - pi->x;
            high_prec dy = pj->y - pi->y;
            high_prec dz = pj->z - pi->z;

            high_prec dist2 = dx * dx + dy * dy + dz * dz + kGravitySofteningEps * kGravitySofteningEps;
            high_prec distance = sqrt(dist2);
            if (distance == 0) distance = kGravitySofteningEps;

            high_prec force = kGravG * pi->m * pj->m / dist2;

            high_prec fx = force * (dx / distance);
            high_prec fy = force * (dy / distance);
            high_prec fz = force * (dz / distance);

            forces[i].x += fx;
            forces[i].y += fy;
            forces[i].z += fz;
            forces[j].x -= fx;
            forces[j].y -= fy;
            forces[j].z -= fz;
        }
    }
}

// =======================
// BARNES–HUT TREE HELPERS
// =======================
// Octant bit-mask constants: bit0=east(x>=cx), bit1=north(y>=cy), bit2=up(z>=cz)
static const int kOctantEast  = 1; // bit 0: x >= node.cx
static const int kOctantNorth = 2; // bit 1: y >= node.cy
static const int kOctantUp    = 4; // bit 2: z >= node.cz

// Octant index: combined bit-mask of East, North, Up flags
static inline int bh_octant(const BHNode& node, high_prec x, high_prec y, high_prec z) {
    return (x >= node.cx ? kOctantEast : 0) | (y >= node.cy ? kOctantNorth : 0) | (z >= node.cz ? kOctantUp : 0);
}

static inline void bh_child_center(const BHNode& node, int q, high_prec& out_cx, high_prec& out_cy, high_prec& out_cz) {
    const high_prec h2 = node.half * 0.5;
    out_cx = node.cx + (q & kOctantEast  ? h2 : -h2);
    out_cy = node.cy + (q & kOctantNorth ? h2 : -h2);
    out_cz = node.cz + (q & kOctantUp    ? h2 : -h2);
}

static inline void bh_add_mass_com(BHNode& node, high_prec m, high_prec x, high_prec y, high_prec z) {
    const high_prec M0 = node.mass;
    const high_prec M1 = M0 + m;
    if (M1 == 0.0) return;

    node.comx = (node.comx * M0 + x * m) / M1;
    node.comy = (node.comy * M0 + y * m) / M1;
    node.comz = (node.comz * M0 + z * m) / M1;
    node.mass = M1;
}

// collision prune: min distance from point to node cube vs (ri + node.max_rad + tol)
static inline bool bh_prune_collision(const BHNode& node, const Particle* pi) {
    high_prec dx = abs_hp(pi->x - node.cx) - node.half;
    high_prec dy = abs_hp(pi->y - node.cy) - node.half;
    high_prec dz = abs_hp(pi->z - node.cz) - node.half;
    if (dx < 0.0) dx = 0.0;
    if (dy < 0.0) dy = 0.0;
    if (dz < 0.0) dz = 0.0;

    const high_prec minDist2 = dx * dx + dy * dy + dz * dz;
    const high_prec R = pi->rad + node.max_rad + EngineCore::collision_distance_tolerance_;
    return (minDist2 > R * R);
}

static void bh_insert(BHNode& node, int idx, std::vector<std::shared_ptr<Particle>>& list) {
    Particle* p = list[idx].get();

    // keep max radius for collision pruning
    node.max_rad = std::max(node.max_rad, p->rad);

    // If too small, bucket everything
    if (node.half <= kBHMinHalf) {
        if (node.particle_index >= 0) {
            node.bucket.push_back(node.particle_index);
            node.particle_index = -1;
        }
        node.bucket.push_back(idx);
        bh_add_mass_com(node, p->m, p->x, p->y, p->z);
        return;
    }

    // If node has bucket, keep bucketing
    if (!node.bucket.empty()) {
        node.bucket.push_back(idx);
        bh_add_mass_com(node, p->m, p->x, p->y, p->z);
        return;
    }

    // Empty leaf
    if (!node.has_children() && node.particle_index < 0) {
        node.particle_index = idx;
        node.mass = p->m;
        node.comx = p->x;
        node.comy = p->y;
        node.comz = p->z;
        return;
    }

    // Leaf with one particle -> subdivide and reinsert both
    if (!node.has_children() && node.particle_index >= 0) {
        const int old_idx = node.particle_index;
        node.particle_index = -1;

        // reset aggregate, then re-accumulate via inserts into children
        node.mass = 0.0;
        node.comx = 0.0;
        node.comy = 0.0;
        node.comz = 0.0;

        auto insert_into_child = [&](int id) {
            Particle* pp = list[id].get();
            const int q = bh_octant(node, pp->x, pp->y, pp->z);
            if (!node.child[q]) {
                high_prec ccx, ccy, ccz;
                bh_child_center(node, q, ccx, ccy, ccz);
                node.child[q] = std::make_unique<BHNode>(ccx, ccy, ccz, node.half * 0.5);
            }
            bh_insert(*node.child[q], id, list);
            bh_add_mass_com(node, pp->m, pp->x, pp->y, pp->z);
        };

        insert_into_child(old_idx);
        insert_into_child(idx);
        return;
    }

    // Internal node: insert into one child
    const int q = bh_octant(node, p->x, p->y, p->z);
    if (!node.child[q]) {
        high_prec ccx, ccy, ccz;
        bh_child_center(node, q, ccx, ccy, ccz);
        node.child[q] = std::make_unique<BHNode>(ccx, ccy, ccz, node.half * 0.5);
    }

    bh_insert(*node.child[q], idx, list);

    //
    bh_add_mass_com(node, p->m, p->x, p->y, p->z);
}

static std::unique_ptr<BHNode> build_bh_tree(std::vector<std::shared_ptr<Particle>>& list, int n) {
    if (n <= 0) return nullptr;

    high_prec minx = list[0]->x, maxx = list[0]->x;
    high_prec miny = list[0]->y, maxy = list[0]->y;
    high_prec minz = list[0]->z, maxz = list[0]->z;

    for (int i = 1; i < n; ++i) {
        const Particle* p = list[i].get();
        if (p->x < minx) minx = p->x;
        if (p->x > maxx) maxx = p->x;
        if (p->y < miny) miny = p->y;
        if (p->y > maxy) maxy = p->y;
        if (p->z < minz) minz = p->z;
        if (p->z > maxz) maxz = p->z;
    }

    const high_prec cx = (minx + maxx) * 0.5;
    const high_prec cy = (miny + maxy) * 0.5;
    const high_prec cz = (minz + maxz) * 0.5;

    high_prec half = std::max({maxx - minx, maxy - miny, maxz - minz}) * 0.5;
    if (half <= 0.0) half = kBHBoundsPad;
    half += kBHBoundsPad;

    auto root = std::make_unique<BHNode>(cx, cy, cz, half);

    for (int i = 0; i < n; ++i) {
        bh_insert(*root, i, list);
    }
    return root;
}

static inline void add_force_pair(
    const Particle* pi,
    const Particle* pj,
    Vector2D& outFi,
    long long* out_interactions)
{
    high_prec dx = pj->x - pi->x;
    high_prec dy = pj->y - pi->y;
    high_prec dz = pj->z - pi->z;

    high_prec dist2 = dx * dx + dy * dy + dz * dz + kGravitySofteningEps * kGravitySofteningEps;
    high_prec dist  = sqrt(dist2);
    if (dist == 0.0) dist = kGravitySofteningEps;

    high_prec F = kGravG * pi->m * pj->m / dist2;

    outFi.x += F * (dx / dist);
    outFi.y += F * (dy / dist);
    outFi.z += F * (dz / dist);

    if (out_interactions) ++(*out_interactions);
}

static inline void add_force_approx(
    const Particle* pi,
    high_prec mass,
    high_prec comx,
    high_prec comy,
    high_prec comz,
    Vector2D& outFi,
    long long* out_interactions)
{
    high_prec dx = comx - pi->x;
    high_prec dy = comy - pi->y;
    high_prec dz = comz - pi->z;

    high_prec dist2 = dx * dx + dy * dy + dz * dz + kGravitySofteningEps * kGravitySofteningEps;
    high_prec dist  = sqrt(dist2);
    if (dist == 0.0) dist = kGravitySofteningEps;

    high_prec F = kGravG * pi->m * mass / dist2;

    outFi.x += F * (dx / dist);
    outFi.y += F * (dy / dist);
    outFi.z += F * (dz / dist);

    if (out_interactions) ++(*out_interactions);
}

static void bh_accumulate_force(
    const BHNode* node,
    int i,
    std::vector<std::shared_ptr<Particle>>& list,
    Vector2D& outFi,
    long long* out_interactions)
{
    if (!node) return;
    if (node->mass == 0.0) return;

    const Particle* pi = list[i].get();

    // Degenerate bucket => direct sum
    if (!node->bucket.empty()) {
        for (int idx : node->bucket) {
            if (idx == i) continue;
            add_force_pair(pi, list[idx].get(), outFi, out_interactions);
        }
        return;
    }

    // Leaf with single particle
    if (!node->has_children()) {
        const int j = node->particle_index;
        if (j < 0 || j == i) return;
        add_force_pair(pi, list[j].get(), outFi, out_interactions);
        return;
    }

    // Barnes–Hut acceptance
    high_prec dx = node->comx - pi->x;
    high_prec dy = node->comy - pi->y;
    high_prec dz = node->comz - pi->z;

    high_prec dist2 = dx * dx + dy * dy + dz * dz + kGravitySofteningEps * kGravitySofteningEps;
    high_prec dist  = sqrt(dist2);
    if (dist == 0.0) dist = kGravitySofteningEps;

    const high_prec s = node->half * 2.0; // side length
    if ((s / dist) < kBHTheta) {
        // Approximate this node as one body at COM
        add_force_approx(pi, node->mass, node->comx, node->comy, node->comz, outFi, out_interactions);
        return;
    }

    // Open node
    for (int q = 0; q < 8; ++q) {
        if (node->child[q]) {
            bh_accumulate_force(node->child[q].get(), i, list, outFi, out_interactions);
        }
    }
}

static void compute_gravity_forcesBH(
    std::vector<Vector2D>& forces,
    std::vector<std::shared_ptr<Particle>>& list,
    int n,
    long long* out_pair_calcs)
{
    std::fill(forces.begin(), forces.end(), Vector2D{0.0, 0.0, 0.0});
    if (n <= 0) return;

    auto root = build_bh_tree(list, n);
    if (!root) return;

    for (int i = 0; i < n; ++i) {
        Vector2D Fi{0.0, 0.0, 0.0};
        bh_accumulate_force(root.get(), i, list, Fi, out_pair_calcs);
        forces[i] = Fi;
    }
}

// =======================
// COLLISIONS VIA BH TREE 
// =======================
static void bh_collide_query(
    const BHNode* node,
    int i,
    Particle* pi,
    std::vector<std::shared_ptr<Particle>>& list,
    long long* out_checks,
    long long* out_hits)
{
    if (!node) return;
    if (bh_prune_collision(*node, pi)) return;

    if (!node->bucket.empty()) {
        for (int j : node->bucket) {
            if (j <= i) continue;
            Particle* pj = list[j].get();
            if (out_checks) ++(*out_checks);
            if (check_collission_ptr(pi, pj)) {
                if (out_hits) ++(*out_hits);
                EngineCore::collisions++;
                resolve_collission_ptr(pi, pj);
            }
        }
        return;
    }

    if (!node->has_children()) {
        const int j = node->particle_index;
        if (j <= i) return;
        Particle* pj = list[j].get();
        if (out_checks) ++(*out_checks);
        if (check_collission_ptr(pi, pj)) {
            if (out_hits) ++(*out_hits);
            EngineCore::collisions++;
            resolve_collission_ptr(pi, pj);
        }
        return;
    }

    for (int q = 0; q < 4; ++q) {
        if (node->child[q]) {
            bh_collide_query(node->child[q].get(), i, pi, list, out_checks, out_hits);
        }
    }
}

// =======================
// UNIFIED COLLISIONS
// =======================
static void resolve_collisions_unified(
    std::vector<std::shared_ptr<Particle>>& list,
    int n,
    long long* out_checks,
    long long* out_hits)
{
    auto root = build_bh_tree(list, n);
    if (!root) return;

    if (!out_checks || !out_hits) {
        for (int i = 0; i < n; ++i) {
            Particle* pi = list[i].get();
            bh_collide_query(root.get(), i, pi, list, nullptr, nullptr);
        }
        return;
    }

    for (int i = 0; i < n; ++i) {
        Particle* pi = list[i].get();
        bh_collide_query(root.get(), i, pi, list, out_checks, out_hits);
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
    //cout << "EngineCore initialized." << endl;
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
    if (dt == 0.0) return 0.0;

    high_prec max_accel = 0;
    for (int i = 0; i < (int)particles->particle_list.size(); i++) {
        high_prec vx = particles->particle_list[i]->vx;
        high_prec vy = particles->particle_list[i]->vy;
        high_prec vz = particles->particle_list[i]->vz;
        high_prec speed = sqrt(vx * vx + vy * vy + vz * vz);
        high_prec accel_magnitude = speed / dt;
        if (accel_magnitude > max_accel) max_accel = accel_magnitude;
    }
    return max_accel;
}

// =======================
// ADAPTIVE DT
// =======================
void EngineCore::configure_dt(std::shared_ptr<scenario> scenario, std::shared_ptr<Particles> particles) {
    dt = scenario->dt;

    high_prec max_accel = calculate_max_acceleration(particles);

    if (max_accel != 0.0) {
        high_prec adaptive_dt = sqrt(kAdaptiveDtAccelEps / max_accel);

        const high_prec base_dt = scenario->dt;
        dt = (adaptive_dt < base_dt) ? adaptive_dt : base_dt;
    } else {
        dt = scenario->dt;
    }
}

void EngineCore::set_overlap_beta(std::shared_ptr<scenario> scenario) {
    kContactBiasBeta = scenario->beta;
}

void EngineCore::set_collision_distance_tolerance(std::shared_ptr<scenario> scenario) {
    collision_distance_tolerance_ = scenario->collision_distance_tolerance;
}

// =======================
// 1 CORE STEP
// Uses cached_force_ as "a(t)" from previous step.
// On step 0: compute initial forces once.
// Each step: drift -> compute new forces (cached_force_) -> second half-kick.
// =======================
void EngineCore::step(std::shared_ptr<Particles> particles, StepEnergyTrace* trace) {

    const int n = (int)particles->particle_list.size();
    auto& list = particles->particle_list;

    if ((int)cached_force_.size() != n) cached_force_.resize(n);

    const bool need_initial_force = (EngineCore::update_iter == 0);

    if (trace) {
        trace->collision_seconds = 0.0;
        trace->verlet_seconds = 0.0;
        trace->gravity_force_seconds = 0.0;
        trace->integration_seconds = 0.0;
        trace->overhead_seconds = 0.0;
        trace->total_seconds = 0.0;
        trace->timing_residual_seconds = 0.0;

        trace->collision_checks   = 0;
        trace->collision_hits     = 0;
        trace->collision_hit_pct  = 0.0;
        trace->gravity_pair_calcs = 0;

        long long collision_checks = 0;
        long long collision_hits   = 0;
        long long gravity_pairs    = 0;

        auto t_total0 = clock_t::now();

        // overhead: TE1
        {
            auto t0 = clock_t::now();
            trace->TE1 = calc_TE(particles);
            auto t1 = clock_t::now();
            trace->overhead_seconds += seconds_between(t0, t1);
        }

        // collisions (counted)
        {
            auto t0 = clock_t::now();
            resolve_collisions_unified(list, n, &collision_checks, &collision_hits);
            auto t1 = clock_t::now();
            trace->collision_seconds = seconds_between(t0, t1);
        }

        // overhead: TE2 + bookkeeping
        {
            auto t0 = clock_t::now();
            trace->TE2 = calc_TE(particles);
            trace->dE_collision  = trace->TE2 - trace->TE1;
            trace->rel_collision = safe_rel_error(trace->TE2, trace->TE1);
            auto t1 = clock_t::now();
            trace->overhead_seconds += seconds_between(t0, t1);
        }

        // gravity + verlet (cached pre-force)
        {
            // initial force only once
            if (need_initial_force) {
                auto t0 = clock_t::now();
                compute_gravity_forcesBH(cached_force_, list, n, &gravity_pairs);
                auto t1 = clock_t::now();
                trace->gravity_force_seconds += seconds_between(t0, t1);
            }

            // integration: first half-kick + drift
            {
                auto t0 = clock_t::now();
                verlet_half_kick(list, cached_force_, n, dt);
                for (int i = 0; i < n; ++i) {
                    Particle* p = list[i].get();
                    p->x += p->vx * dt;
                    p->y += p->vy * dt;
                    p->z += p->vz * dt;
                }
                auto t1 = clock_t::now();
                trace->integration_seconds += seconds_between(t0, t1);
            }

            // force (post) ALWAYS: becomes cached_force_ for next step
            {
                auto t0 = clock_t::now();
                compute_gravity_forcesBH(cached_force_, list, n, &gravity_pairs);
                auto t1 = clock_t::now();
                trace->gravity_force_seconds += seconds_between(t0, t1);
            }

            // integration: second half-kick
            {
                auto t0 = clock_t::now();
                verlet_half_kick(list, cached_force_, n, dt);
                auto t1 = clock_t::now();
                trace->integration_seconds += seconds_between(t0, t1);
            }

            trace->verlet_seconds = trace->gravity_force_seconds + trace->integration_seconds;
        }

        // overhead: TE3 + energy breakdown
        {
            auto t0 = clock_t::now();

            trace->TE3 = calc_TE(particles);
            trace->dE_verlet  = trace->TE3 - trace->TE2;
            trace->rel_verlet = safe_rel_error(trace->TE3, trace->TE2);

            trace->KE = calculate_system_kinetic_energy(particles);
            trace->PE = calculate_system_potential_energy(particles);
            trace->HE = calculate_system_heating_energy(particles);
            trace->TE = trace->KE + trace->PE + trace->HE;

            auto t1 = clock_t::now();
            trace->overhead_seconds += seconds_between(t0, t1);
        }

        auto t_total1 = clock_t::now();
        trace->total_seconds = seconds_between(t_total0, t_total1);

        trace->timing_residual_seconds =
            trace->total_seconds
            - (trace->collision_seconds
               + trace->gravity_force_seconds
               + trace->integration_seconds
               + trace->overhead_seconds);

        trace->collision_checks  = collision_checks;
        trace->collision_hits    = collision_hits;
        trace->collision_hit_pct = (collision_checks > 0)
            ? (100.0 * (double)collision_hits / (double)collision_checks)
            : 0.0;

        trace->gravity_pair_calcs = gravity_pairs;
        return;
    }

    // ---- No-trace path ----
    resolve_collisions_unified(list, n, nullptr, nullptr);

    if (need_initial_force) {
        compute_gravity_forcesBH(cached_force_, list, n, nullptr);
    }

    verlet_half_kick(list, cached_force_, n, dt);

    for (int i = 0; i < n; ++i) {
        Particle* p = list[i].get();
        p->x += p->vx * dt;
        p->y += p->vy * dt;
        p->z += p->vz * dt;
    }

    compute_gravity_forcesBH(cached_force_, list, n, nullptr);
    verlet_half_kick(list, cached_force_, n, dt);
}

// =======================
// ENERGY
// =======================
high_prec EngineCore::calculate_system_kinetic_energy(std::shared_ptr<Particles> particles) {
    high_prec KE = 0.0;
    for (const auto& p : particles->particle_list) {
        KE += 0.5 * p->m * (p->vx * p->vx + p->vy * p->vy + p->vz * p->vz);
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
            high_prec dz = particles->particle_list[j]->z - particles->particle_list[i]->z;
            high_prec distance = sqrt(dx * dx + dy * dy + dz * dz + kEnergySofteningEps * kEnergySofteningEps);
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
    KE += 0.5 * p1->m * (p1->vx * p1->vx + p1->vy * p1->vy + p1->vz * p1->vz);
    KE += 0.5 * p2->m * (p2->vx * p2->vx + p2->vy * p2->vy + p2->vz * p2->vz);
    return KE;
}

high_prec EngineCore::calculate_potential_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2) {
    high_prec dx = p2->x - p1->x;
    high_prec dy = p2->y - p1->y;
    high_prec dz = p2->z - p1->z;
    high_prec distance = sqrt(dx * dx + dy * dy + dz * dz + kEnergySofteningEps * kEnergySofteningEps);
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
    high_prec dz = p1->z - p2->z;
    high_prec dist = sqrt(dx * dx + dy * dy + dz * dz);
    high_prec sum_radii = p1->rad + p2->rad;
    return sum_radii - dist;
}

high_prec EngineCore::heat_ij(high_prec E, std::shared_ptr<Particle> particle_i, std::shared_ptr<Particle> particle_j) {
    high_prec heat = E * (particle_i->m + particle_j->m) / 2;
    high_prec temp = heat / (particle_i->m + particle_i->m);

    particle_i->temp += temp;
    particle_j->temp += temp;

    return E;
}

bool EngineCore::check_overlap(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2) {
    return calculate_overlap_amount(particle1, particle2) > 0.0;
}

bool EngineCore::check_collision(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2) {
    return check_collission_ptr(particle1.get(), particle2.get());
}

bool EngineCore::check_collission(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2) {
    return check_collission_ptr(particle1.get(), particle2.get());
}

void EngineCore::resolve_collisions(std::shared_ptr<Particles> particles) {
    auto& list = particles->particle_list;
    const int n = (int)list.size();
    resolve_collisions_unified(list, n, nullptr, nullptr);
}

void EngineCore::resolve_collision(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2) {
    resolve_collission_ptr(particle1.get(), particle2.get());
}

void EngineCore::resolve_collission(std::shared_ptr<Particle> particle1, std::shared_ptr<Particle> particle2) {
    resolve_collission_ptr(particle1.get(), particle2.get());
}

// Helper for velocity Verlet half-kick
void EngineCore::verlet_half_kick(std::vector<std::shared_ptr<Particle>>& list,
                                  const std::vector<Vector2D>& forces,
                                  int n,
                                  high_prec dt)
{
    for (int i = 0; i < n; ++i) {
        Particle* p = list[i].get();
        p->vx += 0.5 * (forces[i].x / p->m) * dt;
        p->vy += 0.5 * (forces[i].y / p->m) * dt;
        p->vz += 0.5 * (forces[i].z / p->m) * dt;
    }
}
