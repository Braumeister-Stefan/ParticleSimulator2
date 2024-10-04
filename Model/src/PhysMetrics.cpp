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

shared_ptr<test_metrics> Metrics::compute_metrics(shared_ptr<scenario> scenario, shared_ptr<snapshots> particle_states) {
    
    //initialize the metrics object
    shared_ptr<test_metrics> metrics = make_shared<test_metrics>();

    cout << "Computing metrics for " << particle_states->snaps.size() << " snapshots." << endl;

    //loop through the snapshots and compute the metrics
    for (int i = 0; i < particle_states->snaps.size(); i++) {
        

        //initialize the metrics object
        unique_ptr<test_metrics_t> metric = make_unique<test_metrics_t>();

        //to implement
        //compute the metrics for the snapshot
        metric->time = 0.1;
        metric->memory = 0.1;
        metric->cpu = 0.1;
        metric->gpu = 0.1;

        //add the metrics to the metrics object
        metrics->metrics.push_back(move(metric));
    }

    cout << "Metrics computed." << endl << endl;

    return metrics;
}