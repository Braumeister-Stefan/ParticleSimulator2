# ParticleSimulator2

## Project description

This project aims to improve upon the lessons learned from ParticleSimulator1. The objective is to simulate in 2 dimensions, the kinetic and gravitational interactions of n spheroids in an optimized way.

## Objectives

### Primary objectives

*	Separate code in header and source files [Done]
*	Use classes instead of structs - 29/9: Classes seem to differ mostly in access rights and are only usefull when I learn about inheritance and polymorphism. I will stick to structs for now. [Abandoned]
*	Use smart pointers. Each particle should only be defined once. -29/9: Smart pointers used but confusion on when to use unique_ptr and when to use shared_ptr. [DONE]
*	Implement a quadtree (store the O(n^n) design for comparison)  - 29/9: NOT IMPLEMENTED. 24/10: Needs some consideration whether all key components (overlap resolution, collission detection, velocity verlet integration are compatible) -6/2/25: collission detection gets a paralel grid method to speed up overlap detection?

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


## DevNotes

### Notes (14/02/26)
This branch will focus on the implementation of Barnes-Hut approximation algorithm to reduce the complexity from O(n^2) to O(nlogn). Any results are to be carefully validated against BH_validationset scenario to evaluate speedup and error between the baseline and a BH quadtree. Proper logging of time and energy is set up to make this process more rigorous

### Notes (08/02/26)
Performance of the simulation is sped up 5000 times, which well exceeds the 5 times goal I set. Further simulation optimization would need to implement a quadtree and/or spatial grid for collission detection, which is the last goal I set at the beginning, so that's exciting. For evidence, see proofofconcept01-02-26.gif.

Some secondary changes
* RGB assignment is reworked. Now it should give more apparant heating visuals as is shown in the POC. Choosing the right parameters for most appealing/informative visual seems to be more of an art.
* Due to n>1k, GNUPlotter was struggling with visualisation. Plotter is now made more efficient and can skip frames in a more customizable way. Still slightly jittery plotting, but this has no bearing on the actual simulated events so is not a priority as of now.

Lessons learned:
* Some time (2.5h) was lost on git versioning problems. Lesson learned is to always keep sufficient additional backups and when a problem arrises on git, to not deprioritize in favour of more interesting work.
* Each run has a different execution time, and even with 6 runs, there is still a considerable variation due to other processes running/clocking times. A robust approach would need to test each change multiple times and have a null-hypothesis of there not being a significant improvement. Only when this hypothesis can be rejected should the change be kept.
For measurement of computational gains for which 10x speed ups are expected, a robust approach is not needed and can be eye-balled. For expected gains of 1.5x speedups, it is sufficient to take the average time of a handful of runs (say 6). Where optimizations of 1.2x or less are under consideration, a robust and standardized approach would be necessary. Using the above rules of thumb, I only kept optimizations that gave >1.5x improvements.

### Notes (01/02/26)
Following much headaches the physics engine now is able to run simulations based on inelastic collissions with a moderate complexity (n=60). Looking back on my C++ journey, I have started the initial work on ParticleSimulator1 (which was limited to kinetic interactions and ignored gravity) 1 year and 10 months ago. Reflecting back on this, I would say that I did learn basic C++ concept more towards the beginning of the project, but learning shifted almost completely away from this in the last year, towards grasping of geometry, basic physics principles and a bit of visual/game development. As an example, initially I had to engage on broad goals, plotting of particles, reading of input files, separation in header/source files etc, while more recently I focussed exclusively on different ways to deal with overlap resolution. Both are good learning, but in the next stage of the project (see below in next steps) I will work on computational efficiency, which should bring me back to the computer science/C++ side of learning. Below I will discuss in order 1) recap of the inelastic issue and eventual fix, 2) lessons learned, 3) Next steps.

Issue & Fix (01/02/26)

Gravity in each time step will attract particles. For perfectly elastic collissions, each overlap corresponds with an immediate rebounce, hence naive overlap corrections were fine. For imperfectly elastic particles however, this rebounce doesn't happen, hence overlaps will occur in each timesteps and hence any errors introduced by overlap corrections explode over time. I struggled months with how to properly handling overlap correction, with the constraint that any correction needs to respect the conservation of energy. The ultimate solution is to not separate overlap resolution and collission resolution, and simply handle overlap resolution as part of collissions. The core change to do so is when two particles overlap, to introduce an extra outwards velocity component which is larger depending on the size of the overlap. The drawback is that a parameter is to be selected on how strong this bias is. A too large bias seems to give visual appearance of elasticity (which shouldnt be for inelastic particles), while a too small one leads to overlaps growing larger. 0.5 seems to be a good value balancing both for a range of mass values (more mass is more gravity, hence requiring a higher parameter). As a result, benchmark simulations for perfectly inelastic particles now have energy errors well under 1%, which means this topic is closed off. See for a proof of the results the file in this branche: proofofconcept01-02-26.gif.

Looking back at it, there are 3 components leading to this happy (though dragged out) resolution: 1) introduction of a clear benchmark and for each change testing again against this benchmark. This avoids solving one problem and creating 2 new ones, which especially in AI assisted programming is an endemic issue 2) I first encountered this issue exactly one year ago, hence AI has improved significantly. 3) I worked with 2 different freelancers through Fiverr (more on that below), which while it provided mixed results allowed some mistakes to float up and be fixed. 


Lessons Learned (01/02/26).

As I realised the time commitment combined with the little progress over the last year was starting to drag on me, I decided to look for a freelancer on Fiverr platform. This gave mixed results, but overall I think it was a valuable experience (hired help twice). Success is fully defined by how good instructions are and its not always easy to come up with a clear statement of work. For an overview of my approach in managing Fiverr Freelancers, See the readme on https://github.com/Braumeister-Stefan/ParticleSimulator2/tree/Fiverr_branch2 and the statements of work (SOW) on that branch for phase 1 & 2.

Looking back, I spent hundreds of hours on this project, which gave me a lot of learning but I should be aware of the opportunity costs as well. A weekend spent on (a small fraction of) this project means a weekend I can't develop something else, work on my health or follow a course. I want to at least work out the most obvious computational inefficiencies, and create a cool capstone scenario to demonstrate the result of my work. After that, I need to pause and think of the pros/cons whether I continue with implementing a quadcore or call it a wrap and move on to something else.

Next Steps (01/02/26)

Major:
* 3 benchmarks "Planet + Moon System" has been created for different simulation lengths, ranging from 5 hours to 2 minutes and stored for each energy error and runtime in a benchmark. As I am comfortable with the Physics being correct, next I want to get rid of any efficiencies in the code that slow down the simulation, measuring for each change improvements in runtime (and ensure no new energy errors are introduced!). As I was not at all focussed on computational efficiency in my own work, prompts to AI and Fiverr collaborations I do expect there should be some low hanging fruits. Definition of success is at least a 80% reduction in computation time. To be scientifically rigourous I should run simulations multiple times to get a stable average simulation time, but at this stage this seems overkill. [DONE, 500,000% improvement, well over target]

Minor:
* Better handling of RGB values. Currently very slow and visually not appealing. In essence, I want to have no change in colour with mild heating and only large heating should give a visual cue. (Think a brown orb, impacted and the particles on the side of impact should become bright/reddish). [DONE]
* Simulation outcomes can now be loaded as objects, which allows iterative construction of more complex objects. This is necessary to create multi-particle objects such as planets, as it is (in my understanding) an unsolved problem in mathematics how many smaller circles fit into a larger one, hence allowing gravity to do its work acts as a numerical solver to this packing problem. A change that would make this more useful is if I can edit the (system) velocities & positions of these objects in the same interface as for new objects. [DONE]


### Next steps (25/1/25)
My best guesses is that the issue of increasing energy is that this has something to do with overlap resolution  An issue I was not able to solve is that seemingly, energy is created from nothing once particles are collided (leading to disintegration of the system). My best guesses is that this has something to do with overlap resolution 
Possible next steps to look into a solution is to investigate why it seems that after some time a 2-body inelastic colliding system seems to converge and not build up additional energy. Another clue is that this convergence seems to happen much later when the timestep is smaller. I would guess that the way I correct for energy difference after geometrically separating particles has issues in the math. Two alternative things to explore is if the issue is numeral inaccuracy (particle structure is still double datatype even though calculations are almost all high_prec datatype) or a wrong order of overlap resolution, collission resolution and verlet integration. (see .PNG for illustration, note the initial drop is expected due to conversion of KE to HE).

After working on and off for the last 4 months on this project, I managed to get inelastic collissions working, so that particles behave in a physically realistic way (clumping together while retaining linear and angular momentum) and converting KE to heat, which shows visually as particles brightening. An issue I was not able to solve is that seemingly, energy is created from nothing. I will take a break (though not abandon just yet) from this project to focus on my mathematics learning (MST125 at Open University) and perhaps start another python/c++ sideproject. 

### Next steps (29/9/24):
*	Complex object rendering, Saving and retrieving from cache (conditional on cache param in scenario specified). Support for circle only is fine [DONE]
        should include manual inputs into complex name groups "with group defined particles are rendered from the same complex object" to allow for easier user input. Satisfies the GUI criteria for the project. [DONE]
        4/10: Would be nice to have GUI on 1/3 of screen and plot on 2/3 of screen.
* Validation of momentum and kinetic energy conservation [DONE]
* number of particles, fps tracker on the plot handled through the metrics and plotter classes [DONE]

### Notes (24/10/24): 
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




