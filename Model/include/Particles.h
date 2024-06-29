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
    vector<unique_ptr<Particle>> particle_list;
};

#endif // PARTICLES_H