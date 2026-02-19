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
    high_prec half = 0.0; // half side-length

    high_prec mass = 0.0; // total mass in node
    high_prec comx = 0.0; // center-of-mass x
    high_prec comy = 0.0; // center-of-mass y

    //collision specs
    high_prec max_rad = 0.0; // max radius of any particle in this node (for collision detection)

    // Leaf storage options:
    int particle_index = -1;     // single particle leaf (optional)
    std::vector<int> bucket;     // multi-particle leaf (optional, if you used it)

    // Children: 0=NW,1=NE,2=SW,3=SE (as your helpers imply)
    std::array<std::unique_ptr<BHNode>, 4> child;

    BHNode(high_prec in_cx, high_prec in_cy, high_prec in_half)
        : cx(in_cx), cy(in_cy), half(in_half) {}

    inline bool has_children() const {
        return (child[0] || child[1] || child[2] || child[3]);
    }
};

#endif // PARTICLES_H
