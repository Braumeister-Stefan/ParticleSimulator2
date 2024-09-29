//This class can be called by the model and requires a scenario & snapshots object as input. It will write the snapshots to a csv file, read to snapshots object again, and handle related functions.

#ifndef PARTICLE_DUMPER_H
#define PARTICLE_DUMPER_H

#include <iostream>
#include <vector>
#include <memory>
#include "../include/Particles.h"

using namespace std;

class ParticleDumper {
public:
    // Constructor
    ParticleDumper();

    // Destructor
    ~ParticleDumper();


    // Functions

    //write the snapshots to a csv file

    //read the snapshots from a csv file

    //check for the existence of a csv file names snapshots.csv, if it does not exist, create it

    //TO implement

};

#endif // OBJ_HANDLER_H