

#ifndef PARTICLES_H
#define PARTICLES_H


#include <vector>
#include <memory>

#include <boost/multiprecision/cpp_dec_float.hpp>



using namespace std;
using namespace boost::multiprecision;
using high_prec = cpp_dec_float_50;

struct Particle {
    int particle_id;
    high_prec r; //rgb values
    high_prec g;
    high_prec b;
    int rgb; //rgb code
    high_prec x; //position values
    high_prec y;
    high_prec z;
    high_prec vx; //velocity values
    high_prec vy;
    high_prec vz;
    high_prec m; //mass
    high_prec rad; //radius of sphere
    high_prec rest; //restitution parameter of sphere
    high_prec temp=0; //temperature of the particle


};


struct Particles {
    vector<shared_ptr<Particle>> particle_list;

    //default constructor
    Particles() = default;

    //copy constructor
    Particles(const Particles& other) {
        particle_list.reserve(other.particle_list.size());
        for (const auto& particle : other.particle_list) {
            particle_list.push_back(make_shared<Particle>(*particle));
        }
    }
};


//time scaler for backtracking
struct backed_scaler {
    vector<double> scaler;
};

#endif // PARTICLES_H