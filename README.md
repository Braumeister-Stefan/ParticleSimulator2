# ParticleSimulator2 (2026)
_By Stephane Mertens de Wilmars (smdw1997@gmail.com)_

## Project description
ParticleSimulator2 is a 2D physics engine that supports Newtonian gravity and collissions for n-body systems. It's complexity is **O(nlogn)** through implementation of the Barnes-Hut approximation algorithm.

Below are some examples of simulations that are possible with this engine. These and more are stored in [Media](Media).
### Planet Core Formation
WIP
### Simulation of Hydrogen Cloud Collapse
WIP
### Impact of a moon and a planet
WIP
### Stable rotation of a moon around a star 
WIP


## Objectives
I worked on this project to improve upon the lessons learned from ParticleSimulator1 and build a solid C++ foundation. 

Below are the main objectives I set myself to reach before considering the project a success. After 1 year and 8 months, all objectives were reached. I tried to log my learning journey and the project evolution towards these objectives in [DevNotes](DevNotes.md).

### Primary objectives

*	Separate code in header and source files [DONE]
*	Use smart pointers. Each particle should only be defined once.  [DONE]
*	Implement a quadtree (store the O(n^n) design for comparison)   [DONE]

### Secondary objectives

* Loading of objects/run, user interface to navigate options, using .csv instead of .txt [DONE]
*	Implement a simple helper function that can be wrapped around any function and prints run time to measure efficiency [DONE]
*	Allow for particle radius and particle mass to be individualized [DONE]
*	Separate computation and play-back of simulations [DONE]

## Project structure

*	Main.cpp: This is the main file from which the code runs. This is what you should build locally.

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

## User Guide
This section provides a rough guide for those willing to play around with the tool. I believe it could be of use for anyone developing space-based games in c++ or for undergraduate physics/compsci students. Happy to respond to any inquiries over email if there is interest.
### First use

### Logging system

### Main Inputs
All below inputs are accessible through scenario_inputs.csv and object_inputs.csv. To generate a new scenario/object, simply add a line with the desired inputs for your run.

WIP
### Further Inputs
The below is a (not exhaustive) list of parameters that can be changed only by editing the code itself. This is advised only for those wishing to extend/tweak the project.

### Limitations & Warnings
WIP











