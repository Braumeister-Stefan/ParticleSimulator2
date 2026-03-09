#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <memory>
#include <array>   // <-- ADD THIS

using namespace std;
using high_prec = double;

struct Particle {
    int particle_id;
    high_prec r;
    high_prec g;
    high_prec b;
    int rgb;
    high_prec x;
    high_prec y;
    high_prec z;
    high_prec vx;
    high_prec vy;
    high_prec vz;
    high_prec m;
    high_prec rad;
    high_prec rest;
    high_prec temp = 0;
};

struct Particles {
    vector<shared_ptr<Particle>> particle_list;

    Particles() = default;

    Particles(const Particles& other) {
        particle_list.reserve(other.particle_list.size());
        for (const auto& particle : other.particle_list) {
            particle_list.push_back(make_shared<Particle>(*particle));
        }
    }
};

// time scaler for backtracking
struct backed_scaler {
    vector<double> scaler;
};

// =======================================================
// Barnes–Hut node (needed by compute_gravity_forcesBH)
// =======================================================
struct BHNode {
    high_prec cx = 0.0;   // node center x
    high_prec cy = 0.0;   // node center y
    high_prec cz = 0.0;   // node center z
    high_prec half = 0.0; // half side-length

    high_prec mass = 0.0; // total mass in node
    high_prec comx = 0.0; // center-of-mass x
    high_prec comy = 0.0; // center-of-mass y
    high_prec comz = 0.0; // center-of-mass z

    //collision specs
    high_prec max_rad = 0.0; // max radius of any particle in this node (for collision detection)

    // Leaf storage options:
    int particle_index = -1;     // single particle leaf (optional)
    std::vector<int> bucket;     // multi-particle leaf (optional, if you used it)

    // Children: octree (8 nodes) indexed by bit-mask: bit0=east(x>=cx), bit1=north(y>=cy), bit2=up(z>=cz)
    std::array<std::unique_ptr<BHNode>, 8> child;

    BHNode(high_prec in_cx, high_prec in_cy, high_prec in_cz, high_prec in_half)
        : cx(in_cx), cy(in_cy), cz(in_cz), half(in_half) {}

    inline bool has_children() const {
        return (child[0] || child[1] || child[2] || child[3] ||
                child[4] || child[5] || child[6] || child[7]);
    }
};

#endif // PARTICLES_H
