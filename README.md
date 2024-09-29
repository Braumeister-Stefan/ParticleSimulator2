# ParticleSimulator2

## Project description

This project aims to improve upon the lessons learned from ParticleSimulator1. The objective is to simulate in 2 dimensions, the kinetic and gravitational interactions of n spheroids in an optimized way.

## Objectives

### Primary objectives

*	Separate code in header and source files
*	Use classes instead of structs - 29/9: Classes seem to differ mostly in access rights and are only usefull when I learn about inheritance and polymorphism. I will stick to structs for now.
*	Use smart pointers. Each particle should only be defined once. -29/9: Smart pointers used but confusion on when to use unique_ptr and when to use shared_ptr.
*	Implement a quadtree (store the O(n^n) design for comparison)  - 29/9: NOT IMPLEMENTED

### Secondary objectives

* Loading of objects/run, user interface to navigate options, using .csv instead of .txt - 29/9: csv implemented, changing parameters is much user friendlier now.(see devnotes)
*	Implement a simple helper function that can be wrapped around any function and prints run time to measure efficiency
*	Allow for particle radius and particle mass to be individualized - 29/9: physics seem reasonable but not validated yet
*	Separate computation and play-back of simulations - DONE

## Recycling from ParticleSimulator1

The elements below are to be recycled from ParticleSimulator1, in accordance with the primary and secondary objectives:
*	Kinetic Collission function (including detection, backtracking and resolution)
*	Complex object generation (sphere only is OK) (including storage of complex object)
*	GNUplot plotting engine
*	Kinetic energy, momentum, fps trackers

## Project structure

*	Main.cpp file

Header files:
*	ObjHandler.h : this will contain all functionality in generating, saving and unpacking complex objects
*	PhysEngine.h: this will initialize the simulation, run through the timesteps, update particle locations and save Particle information & meta information (fps etc) at the end of the simulation
*	ParticlePlotter.h: this will set up the plotting engine and plot particles for every timestep
*	Interfacer.h: this will read the input parameters, define scenarios and allow the user to choose which scenario to run.
* Particles.h: this will define the individual, grouped and time-variant grouped particles and meta information per timestep.
* Utils.h: this will contain all math heavy functions  called upon by PhysEngine 

Source files:
* ObjHandler.cpp
* PhysEngine.cpp
* ParticlePlotter.cpp
* Interfacer.cpp
* Particles.cpp
* Utils.cpp


## DevNotes

Next steps (29/9/24):
*	Complex object rendering, Saving and retrieving from cache (conditional on cache param in scenario specified). Support for circle only is fine
        should include manual inputs into complex name groups "with group defined particles are rendered from the same complex object" to allow for easier user input. Satisfies the GUI criteria for the project.
* Validation of momentum and kinetic energy conservation
* number of particles, fps tracker on the plot handled through the metrics and plotter classes


Issues:
*	When to use unique_ptr and when to use shared_ptr? I see the logic of pointers giving the benefit of clear memory allocation and deallocation. However, they require complex syntax. I will try to use them in computationally demanding parts of the code, but not simple ones such as reading and writing inputs.


