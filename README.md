# ParticleSimulator2

## Project description

This project aims to improve upon the lessons learned from ParticleSimulator1. The objective is to simulate in 2 dimensions, the kinetic and gravitational interactions of n spheroids in an optimized way.

## Objectives

### Primary objectives

*	Separate code in header and source files [Done]
*	Use classes instead of structs - 29/9: Classes seem to differ mostly in access rights and are only usefull when I learn about inheritance and polymorphism. I will stick to structs for now. [Abandoned]
*	Use smart pointers. Each particle should only be defined once.  [DONE]
*	Implement a quadtree (store the O(n^n) design for comparison)  - [DONE] 
### Secondary objectives

* Loading of objects/run, user interface to navigate options, using .csv instead of .txt - 29/9: csv implemented, changing parameters is much user friendlier now.(see devnotes) [DONE]
*	Implement a simple helper function that can be wrapped around any function and prints run time to measure efficiency [Done]
*	Allow for particle radius and particle mass to be individualized - 29/9: physics seem reasonable but not validated yet [DONE]
*	Separate computation and play-back of simulations - [DONE]

## Recycling from ParticleSimulator1

The elements below are to be recycled from ParticleSimulator1, in accordance with the primary and secondary objectives:
*	Kinetic Collission function (including detection, backtracking and resolution) [DONE]
*	Complex object generation (sphere only is OK) (including storage of complex object) [DONE]
*	GNUplot plotting engine [DONE]
*	Kinetic energy, momentum, fps trackers [DONE]

## Project structure

*	Main.cpp file

Header files:
* ObjHandler.h : this will contain all functionality in generating, saving and unpacking complex objects
* PhysEngineWrapper.h: this will initialize the simulation, run through the timesteps, measure energy and time differences with benchmark
* PhysEngineCore.h: for a given timestep, this file will handle parameter choices, collissions, integration, update particle locations and save Particle information & meta information at the end of the simulation. 
* ParticlePlotter.h: this will set up the plotting engine and plot particles for every timestep, it will also assign RGB values conditional on temperature.
* Interfacer.h: this will read the input parameters, define scenarios and allow the user to choose which scenario to run.
* Particles.h: this will define the individual, grouped and time-variant grouped particles and meta information per timestep.
* PhysMetrics.h: this will define all functionality related to the validation (Both physics and efficiency related) metrics.
* MathUtils.h: this will contain all math heavy functions  called upon by PhysEngine 

Source files:
* ObjHandler.cpp
* PhysEngineWrapper.cpp
* PhysEngineCore.cpp
* ParticlePlotter.cpp
* Interfacer.cpp
* PhysMetrics.cpp
* Particles.cpp
* MathUtils.cpp










