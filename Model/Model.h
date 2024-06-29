#ifndef PARTICLEMODEL_H
#define PARTICLEMODEL_H

#include <iostream>
#include <memory>
#include "include/Interfacer.h"
#include "include/PhysEngine.h"
#include "include/ParticlePlotter.h"
#include "include/ObjHandler.h"

// Use the standard namespace
using namespace std;

class ParticleModel {
public:
    // Constructor
    ParticleModel();

    // Destructor
    ~ParticleModel();


    // Submodels
    unique_ptr<Interfacer> interfacer;
    unique_ptr<Engine> engine;
    unique_ptr<Plotter> plotter;
    unique_ptr<ObjHandler> obj_handler;

};

#endif // PARTICLEMODEL_H

