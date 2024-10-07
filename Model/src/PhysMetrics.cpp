#include "../include/PhysMetrics.h"
#include "../include/InitStructs.h"

#include <memory>

//namespaces
using namespace std;

// Constructor
Metrics::Metrics() {
    cout << "Metrics computor initialized." << endl;
}

// Destructor
Metrics::~Metrics() {
    cout << "Metrics computor destroyed." << endl;
}

//define the compute_metrics function

shared_ptr<snapshots> Metrics::compute_metrics(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {
    //this function computes the metrics for the simulation. Used for validation purposes.
    

    cout << "Computing metrics for " << particle_states->snaps.size() << " snapshots." << endl;

    //loop through the snapshots and compute the metrics
    for (int i = 0; i < particle_states->snaps.size(); i++) {
        
        
    }

    cout << "Metrics computed." << endl << endl;

    //average the fps metrics every 100 snapshots to smooth the curve
    for (int i = 0; i < particle_states->metrics.size(); i++) {
        if (i % 10 == 0) {
            double avg_fps = 0;
            for (int j = 0; j < 10; j++) {
                avg_fps += particle_states->metrics[i + j]->fps;
            }
            avg_fps /= 10;
            particle_states->metrics[i]->fps = avg_fps;
        }
    }

    //PLACEHOLDER, if metrics is a null pointer, create a new metrics object fill with zeros
    if (particle_states->metrics.size() == 0) {
        for (int i = 0; i < particle_states->snaps.size(); i++) {
            particle_states->metrics.push_back(make_shared<test_metrics_t>());
            particle_states->metrics[i]->fps = 0;
            particle_states->metrics[i]->memory = 0;
            particle_states->metrics[i]->cpu = 0;
            particle_states->metrics[i]->gpu = 0;
        }
    }


    return particle_states;
}