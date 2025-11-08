
#ifndef SNAPSHOTS_H
#define SNAPSHOTS_H

#include <vector>
#include <memory>
#include "Particles.h"  // Needed because snapshots contains Particles

using namespace std;

class snapshots {
public:
    vector<shared_ptr<Particles>> snaps;

    snapshots() = default;
    ~snapshots() = default;
};

#endif
