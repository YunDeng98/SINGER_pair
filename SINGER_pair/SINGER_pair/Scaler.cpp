//
//  Scaler.cpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/11/24.
//

#include "Scaler.hpp"

void Scaler::rescale(PSMC &psmc) {
    num_temporal_bins = psmc.num_temporal_bins;
    temporal_grids = psmc.temporal_grids;
}

void Scaler::add_mutations(int count, double theta, double t) {
    double overlap;
    for (int k = 0; k < num_temporal_bins; k++) {
        if (temporal_grids[k] > t) {
            break;
        }
        overlap = (min(t, temporal_grids[k+1]) - temporal_grids[k])/t;
        mutation_counts[k] += overlap*count/theta;
    }
}

void Scaler::add_branch_length(double t, double l) {
    double overlap;
    for (int k = 0; k < num_temporal_bins; k++) {
        if (temporal_grids[k] > t) {
            break;
        }
        overlap = min(t, temporal_grids[k+1]) - temporal_grids[k];
        branch_length[k] += overlap*l;
    }
}
