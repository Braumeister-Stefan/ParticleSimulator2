#ifndef PARTICLES_H
#define PARTICLES_H


#include <vector>
#include <memory>

using namespace std;

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


struct backed_scaler {
    vector<double> scaler;
};

#endif // PARTICLES_H