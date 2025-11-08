
#ifndef BACKED_SCALER_H
#define BACKED_SCALER_H

#include <vector>
#include <memory>

using namespace std;

class backed_scaler {
public:
    vector<double> scaler; // matches how you used scaler[i]
    
    backed_scaler() = default;
    ~backed_scaler() = default;
};

#endif
