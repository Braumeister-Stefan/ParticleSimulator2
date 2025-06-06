#ifndef PHYS_ENGINE_H
#define PHYS_ENGINE_H

//Standard libraries
#include <iostream>
#include <memory>

//Internal libraries
#include "InitStructs.h"

//external libraries
#include <boost/multiprecision/cpp_dec_float.hpp>



using namespace std;
using namespace boost::multiprecision;

using high_prec = cpp_dec_float_50;


class Engine {
public:

    //declare static member variables
    static high_prec dt;

    //declare particle_states as a static member variable
    shared_ptr<snapshots> particle_states;


    // Constructor
    Engine();

    // Destructor
    ~Engine();

    //structs

    //structure with momentum x,y vector

    struct momentum {
        high_prec x;
        high_prec y;

        // Define arithmetic operators
        momentum operator+(const momentum& other) const {
            return { x + other.x, y + other.y };
        }

        momentum operator-(const momentum& other) const {
            return { x - other.x, y - other.y };
        }

        momentum operator*(const high_prec& scalar) const {
            return { x * scalar, y * scalar };
        }

        momentum operator/(const high_prec& scalar) const {
            return { x / scalar, y / scalar };
        }

        momentum operator/(const momentum& other) const {
            return { x / other.x, y / other.y };
        }
    };

    //Function to run the simulation and return snapshots of each time step
    shared_ptr<snapshots> run(shared_ptr<scenario> scenario, shared_ptr<Particles> particles);

    //Function to update the particles
    void update_particles(shared_ptr<Particles> particles);

    //function to access margin TE errors
    high_prec get_margin_TE_error();
    high_prec get_margin_TE_error_overlap();
    high_prec get_margin_TE_error_collision();
    high_prec get_margin_TE_error_integrate();
    high_prec get_overlap_iters_in_step();
    high_prec get_margin_TE_error_overlap_ij_transl();
    high_prec get_margin_TE_error_overlap_ij_corrected();



    //Functions to resolve overlap between particles
    bool resolve_overlaps(shared_ptr<Particles> particles);
    void resolve_overlap_ij(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    

    //heating functions 
    high_prec heat_ij(high_prec E, shared_ptr<Particle> particle_i, shared_ptr<Particle> particle_j);

    //functions to resolve energy gap after inelastic collission (currently also used as stop gap for energy conservation in resolve_ovverlap_ij)
    void resolve_energy_gap(shared_ptr<Particle> particle_i, shared_ptr<Particle> particle_j, high_prec delta_E);
    

    

    //Functions to resolve collissions between particles
    void resolve_collisions(shared_ptr<Particles> particles);


    void resolve_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);

    bool check_collission(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    bool check_overlap(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);



    //Functions to resolve gravitational attraction between particles
    void resolve_gravity_euler(shared_ptr<Particles> particles);
    void resolve_gravity_verlet(shared_ptr<Particles> particles);

    //Function to update the locations of the particles
    void update_locations(shared_ptr<Particles> particles, shared_ptr<backed_scaler> scaler);

    //Function to interact with cache
    void run_to_cache(shared_ptr<scenario> scenario, shared_ptr<snapshots> snapshots);
    bool cache_exists(shared_ptr<scenario> scenario);
    shared_ptr<snapshots> run_from_cache(shared_ptr<scenario> scenario);

    //Function to validate TE
    high_prec calc_TE(shared_ptr<Particles> particles);
    high_prec calc_TE_ij(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2, bool verbose = false);
    high_prec calculate_kinetic_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);
    high_prec calculate_potential_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);
    high_prec calculate_heating_energy(std::shared_ptr<Particle> p1, std::shared_ptr<Particle> p2);
    

    momentum calc_mom(shared_ptr<Particles> particles); 
    momentum calc_mom_ij(shared_ptr<Particle> particle1, shared_ptr<Particle> particle2);
    
};

#endif
