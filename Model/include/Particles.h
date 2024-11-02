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
    double r; //rgb values
    double g;
    double b;
    int rgb; //rgb code
    double x; //position values
    double y;
    double z;
    double vx; //velocity values
    double vy;
    double vz;
    double m; //mass
    double rad; //radius of sphere
    double rest; //restitution parameter of sphere
    double temp; //temperature of the particle


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