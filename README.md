# ParticleSimulator2

## Project description

This project aims to improve upon the lessons learned from ParticleSimulator1. The objective is to simulate in 2 dimensions, the kinetic and gravitational interactions of n spheroids in an optimized way.

## Objectives

### Primary objectives

*	Separate code in header and source files [Done]
*	Use classes instead of structs - 29/9: Classes seem to differ mostly in access rights and are only usefull when I learn about inheritance and polymorphism. I will stick to structs for now. [Abandoned]
*	Use smart pointers. Each particle should only be defined once. -29/9: Smart pointers used but confusion on when to use unique_ptr and when to use shared_ptr. [DONE]
*	Implement a quadtree (store the O(n^n) design for comparison)  - 29/9: NOT IMPLEMENTED. 24/10: Needs some consideration whether all key components (overlap resolution, collission detection, velocity verlet integration are compatible)

### Secondary objectives

* Loading of objects/run, user interface to navigate options, using .csv instead of .txt - 29/9: csv implemented, changing parameters is much user friendlier now.(see devnotes) [DONE]
*	Implement a simple helper function that can be wrapped around any function and prints run time to measure efficiency
*	Allow for particle radius and particle mass to be individualized - 29/9: physics seem reasonable but not validated yet [DONE]
*	Separate computation and play-back of simulations - [DONE]

## Recycling from ParticleSimulator1

The elements below are to be recycled from ParticleSimulator1, in accordance with the primary and secondary objectives:
*	Kinetic Collission function (including detection, backtracking and resolution) [DONE]
*	Complex object generation (sphere only is OK) (including storage of complex object) [DONE]
*	GNUplot plotting engine [DONE]
*	Kinetic energy, momentum, fps trackers 

## Project structure

*	Main.cpp file

Header files:
* ObjHandler.h : this will contain all functionality in generating, saving and unpacking complex objects
* PhysEngine.h: this will initialize the simulation, run through the timesteps, update particle locations and save Particle information & meta information (fps etc) at the end of the simulation
* ParticlePlotter.h: this will set up the plotting engine and plot particles for every timestep
* Interfacer.h: this will read the input parameters, define scenarios and allow the user to choose which scenario to run.
* Particles.h: this will define the individual, grouped and time-variant grouped particles and meta information per timestep.
* PhysMetrics.h: this will define all functionality related to the validation (Both physics and efficiency related) metrics.
* Utils.h: this will contain all math heavy functions  called upon by PhysEngine 

Source files:
* ObjHandler.cpp
* PhysEngine.cpp
* ParticlePlotter.cpp
* Interfacer.cpp
* PhysMetrics.cpp
* Particles.cpp
* Utils.cpp


## DevNotes

Update (25/1):

After working on and off for the last 4 months on this project, I managed to get inelastic collissions working, so that particles behave in a physically realistic way (clumping together while retaining linear and angular momentum) and converting KE to heat, which shows visually as particles brightening. An issue I was not able to solve is that seemingly, energy is created from nothing. I will take a break (though not abandon just yet) this project to focus on my mathematics learning (MST125 at Open University) and perhaps start another python/c++ sideproject. 

Next steps (25/1/25)
My best guesses is that the issue of increasing energy is that this has something to do with overlap resolution  An issue I was not able to solve is that seemingly, energy is created from nothing once particles are collided (leading to disintegration of the system). My best guesses is that this has something to do with overlap resolution 
Possible next steps to look into a solution is to investigate why it seems that after some time a 2-body inelastic colliding system seems to converge and not build up additional energy. Another clue is that this convergence seems to happen much later when the timestep is smaller. I would guess that the way I correct for energy difference after geometrically separating particles has issues in the math. Two alternative things to explore is if the issue is numeral inaccuracy (particle structure is still double datatype even though calculations are almost all high_prec datatype) or a wrong order of overlap resolution, collission resolution and verlet integration. (see .PNG for illustration, note the initial drop is expected due to conversion of KE to HE).

Next steps (29/9/24):
*	Complex object rendering, Saving and retrieving from cache (conditional on cache param in scenario specified). Support for circle only is fine [DONE]
        should include manual inputs into complex name groups "with group defined particles are rendered from the same complex object" to allow for easier user input. Satisfies the GUI criteria for the project. [DONE]
        4/10: Would be nice to have GUI on 1/3 of screen and plot on 2/3 of screen.
* Validation of momentum and kinetic energy conservation
* number of particles, fps tracker on the plot handled through the metrics and plotter classes

Notes (24/10/24): 
* Lost 2 days of work due to a deleting some files. Managed to recover but teaches the valuable lesson to push small incremental changes to git instead of everything at once.
* Currently able to simulate moderately complex (n=50, 5000 steps) simulations, limited to perfectly elastic conditions though.
         * Proof of Concept gif is added to repo
* Math and physics needed to understand the current PhysEngine is not very high, but to come up with the right solution required multiple evenings of reading, youtube and use of LLMs.
* Chatgpt-4o is a great tool, but will only provide the elementary building blocks once the problem & solution is clearly described. For future problems, it would be a good idea to take some time drawing out a problem and solution on paper first, before starting to generate code and incrementally improve it.
* Some doubt on whether implementing quadtree method is 1) within scope of my skillset, 2) has sufficient pay off versus simpler efficiency increasing techniques (e.g. only checking for overlaps on nearby particles).
* Some (Arbitrary) threshold on artificial energy needs to be set. I will put it at 1% meaning that simulations that exceed this threshold should not be accepted.
* There will always be a trade-off between precision and computational speed. In case of the high_prec data structure, this is a worthwile tradeoff.

Next steps (24/10/24):
* Wrap up some minor points set out above in 29/9.
        * number of particles, fps tracker (technically correct but shifts too fast to see).
        * validation of momentum
        * GUI on 1/3 screen and plot other 2/3 screen. Darkmode (if not too hard)
* Saving of rendered simulations should only save down the 1/dt-th frames to make saving simulations less of a wait.
* Gravitational attraction has a dampening factor (ɛ̝) that prevents division by (near zero). This needs to be made relative to the size of particles so as to allow smaller particles to function properly.
* Inelastic (and partially inelastic) collissions need to be programmed. To my understanding this involves a transfer of kinetic energy to heat which needs to be programmed. Can be strictly visual for now, e.g. adding heat as a Particle property and linearly increasing r,g,b values on the plot based on its value.
* Complex object "SPINNING CLOUD (name is WIP)" to be made, which is same as CIRCLE but with angular momentum. Needs to be tested so that it remains stable (for a certain mass) over time if not affected by external forces. This will allow me to make meaningfull scenarios which should push against the computational bounderies of an O(n2) method.




