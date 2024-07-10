//
//  Scaler.cpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/11/24.
//

#include "Scaler.hpp"

Scaler::Scaler() {}

void Scaler::rescale(PSMC &psmc) {
    num_temporal_bins = psmc.num_temporal_bins;
    mutation_counts.resize(num_temporal_bins);
    branch_length.resize(num_temporal_bins);
    old_grids.push_back(0.0);
    for (double x : psmc.mid_points) {
        old_grids.push_back(x);
    }
    // double m = 0;
    for (int i = 0; i < psmc.num_het_sites.size(); i++) {
        if (psmc.num_het_sites[i] > 0) {
            // m = psmc.thetas[i]/(psmc.spatial_grids[i+1] - psmc.spatial_grids[i]);
            add_mutations(psmc.num_het_sites[i], 2*psmc.thetas[i], psmc.posterior_averages[i]);
        }
    }
    for (int i = 0; i < psmc.posterior_averages.size() - 1; i++) {
        add_branch_length(psmc.posterior_averages[i], 1);
    }
    cout << accumulate(psmc.num_het_sites.begin(), psmc.num_het_sites.end(), 0.0) << endl;
    cout << accumulate(mutation_counts.begin(), mutation_counts.end(), 0.0) << endl;
    cout << accumulate(branch_length.begin(), branch_length.end(), 0.0) << endl;
    // cout << accumulate(psmc.posterior_averages.begin(), psmc.posterior_averages.end(), 0.0) << endl;
    compute_new_grids();
    rewrite_temporal_grids(psmc);
}

void Scaler::add_mutations(int count, double theta, double t) {
    assert(t < old_grids.back());
    double overlap;
    for (int k = 0; k < num_temporal_bins; k++) {
        if (old_grids[k] > t) {
            break;
        }
        overlap = (min(t, old_grids[k+1]) - old_grids[k])/t;
        mutation_counts[k] += overlap*count/theta;
    }
    // cout << accumulate(mutation_counts.begin(), mutation_counts.end(), 0.0) << endl;
}

void Scaler::add_branch_length(double t, double l) {
    double overlap;
    for (int k = 0; k < num_temporal_bins; k++) {
        if (old_grids[k] > t) {
            break;
        }
        overlap = min(t, old_grids[k+1]) - old_grids[k];
        branch_length[k] += overlap*l;
    }
}

void Scaler::compute_new_grids() {
    new_grids.resize(old_grids.size());
    new_grids[0] = 0;
    for (int i = 1; i < new_grids.size(); i++) {
        new_grids[i] = new_grids[i-1] + mutation_counts[i-1]/ branch_length[i-1]*(old_grids[i] - old_grids[i-1]);
    }
}

void Scaler::rewrite_temporal_grids(PSMC &psmc) {
    for (int i = 0; i < num_temporal_bins; i++) {
        psmc.mid_points[i] = new_grids[i+1];
    }
}
