//
//  Scaler.hpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/11/24.
//

#ifndef Scaler_hpp
#define Scaler_hpp

#include <stdio.h>
#include "PSMC.hpp"

class Scaler {
    
public:
    
    int num_temporal_bins;
    vector<double> old_grids;
    vector<double> new_grids;
    vector<double> mutation_counts;
    vector<double> branch_length;
    
    Scaler();
    
    void rescale(PSMC &psmc);
    
    void rewrite_temporal_grids(PSMC &psmc);
    
    void add_mutations(int count, double theta, double t);
    
    void add_branch_length(double t, double l);
    
    void compute_new_grids();
};

#endif /* Scaler_hpp */
