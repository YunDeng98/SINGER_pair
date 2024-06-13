//
//  PSMC.cpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/10/24.
//

#include "PSMC.hpp"

PSMC::PSMC() {}

void PSMC::initialize(VCF_reader &vr, Rate_map &rm) {
    double m = rm.mean_rate();
    Ne = vr.Ne;
    start = vr.start;
    end = vr.end;
    bin_size = ceil(unit_theta/m/Ne/2);
    int n = ceil((end - start + 2)/bin_size);
    forward_probs.resize(n, vector<double>(num_temporal_bins, 0.0));
    backward_probs.resize(n, vector<double>(num_temporal_bins, 0.0));
    emission_probs.resize(n, vector<double>(num_temporal_bins, 0.0));
    spatial_grids.resize(n + 1);
    num_het_sites.resize(n);
    posterior_averages.resize(n);
    thetas.resize(n);
    rhos.resize(n);
    temporal_grids.resize(num_temporal_bins + 1);
    mid_points.resize(num_temporal_bins);
    forward_factors.resize(num_temporal_bins);
    backward_factors.resize(num_temporal_bins);
    diagonals.resize(num_temporal_bins);
    lower_diagonals.resize(num_temporal_bins);
    upper_diagonals.resize(num_temporal_bins);
    lower_sums.resize(num_temporal_bins);
    upper_sums.resize(num_temporal_bins);
    create_spatial_bins();
    create_temporal_bins();
    compute_thetas(rm);
    compute_rhos(m);
}

void PSMC::create_spatial_bins() {
    for (int i = 0; i < spatial_grids.size(); i++) {
        spatial_grids[i] = min(start - 1 + bin_size*i, end + 1);
    }
}

void PSMC::create_temporal_bins() {
    temporal_grids[0] = 0;
    double q = 0;
    for (int i = 1; i < num_temporal_bins; i++) {
        q = 1 - (i+0.0)/num_temporal_bins;
        temporal_grids[i] = -log(q);
    }
    temporal_grids[num_temporal_bins] = numeric_limits<double>::infinity();
    for (int i = 0; i < num_temporal_bins; i++) {
        q = 1 - (i+0.5)/num_temporal_bins;
        mid_points[i] = -log(q);
    }
}

void PSMC::compute_num_het_sites(VCF_reader &vr) {
    num_het_sites = vr.compute_num_het_sites(spatial_grids);
}

void PSMC::compute_thetas(Rate_map &rm) {
    for (int i = 0; i < thetas.size(); i++) {
        thetas[i] = rm.segment_distance(spatial_grids[i], spatial_grids[i+1])*Ne;
    }
}

void PSMC::compute_rhos(double m) {
    for (int i = 0; i < rhos.size(); i++) {
        rhos[i] = m*Ne*(spatial_grids[i+1] - spatial_grids[i]);
    }
}

void PSMC::compute_emission_probs() {
    double l = 0;
    double e = 0;
    for (int i = 0; i < thetas.size(); i++) {
        for (int j = 0; j < temporal_grids.size() - 1; j++) {
            l = 2*mid_points[j];
            e = exp(-thetas[i]*l)*pow(thetas[i]*l, num_het_sites[i]);
            emission_probs[i][j] = e;
        }
    }
}

void PSMC::forward_algorithm() {
    compute_forward_factors();
    for (int k = 0; k < num_temporal_bins; k++) {
        forward_probs[0][k] = 1.0/num_temporal_bins;
    }
    forward_emit(0);
    for (int i = 0; i < thetas.size() - 1; i++) {
        forward_transition(i);
        forward_emit(i+1);
    }
}

void PSMC::backward_algorithm() {
    compute_backward_factors();
    int n = (int) thetas.size();
    for (int k = 0; k < num_temporal_bins; k++) {
        backward_probs[n - 1][k] = 1;
    }
    backward_emit(n - 1);
    for (int j = n - 1; j > 0; j--) {
        backward_transition(j);
        backward_emit(j-1);
    }
}

void PSMC::posterior_average() {
    double ws;
    double wt;
    for (int i = 0; i < spatial_grids.size() - 1; i++) {
        ws = 0;
        wt = 0;
        for (int k = 0; k < num_temporal_bins; k++) {
            ws += forward_probs[i][k]*backward_probs[i][k];
            wt += forward_probs[i][k]*backward_probs[i][k]*mid_points[k];
        }
        posterior_averages[i] = wt/ws;
    }
}

void PSMC::forward_emit(int i) {
    double ws = 0;
    for (int k = 0; k < num_temporal_bins; k++) {
        forward_probs[i][k] *= emission_probs[i][k];
        ws += forward_probs[i][k];
    }
    for (int k = 0; k < num_temporal_bins; k++) {
        forward_probs[i][k] /= ws;
    }
}

void PSMC::backward_emit(int j) {
    double ws = 0;
    for (int k = 0; k < num_temporal_bins; k++) {
        backward_probs[j][k] *= emission_probs[j][k];
        ws += backward_probs[j][k];
    }
    for (int k = 0; k < num_temporal_bins; k++) {
        backward_probs[j][k] /= ws;
    }
}

double PSMC::psmc_cdf(double rho, double s, double t) {
    double l;
    double pre_factor;
    if (t <= s) {
        l = 2*t;
    } else {
        l = 2*s;
    }
    if (l == 0) {
        pre_factor = rho;
    } else {
        pre_factor = (1 - exp(-rho*l))/l;
    }
    double integral;
    double cdf;
    if (t <= s) {
        integral = 2*t + 2*exp(-t) - 2;
    } else {
        integral = 2*s + exp(-t) + exp(-t) - 2*exp(s-t);
    }
    cdf = pre_factor*integral;
    return cdf;
}

double PSMC::psmc_prob(double rho, double s, double t1, double t2) {
    double prob;
    double base;
    double l = 2*s;
    if (t1 < s and t2 > s) {
        base = exp(-rho*l);
    } else {
        base = 0;
    }
    double gap;
    double uq = 0;
    double lq = 0;
    uq = psmc_cdf(rho, s, t2);
    lq = psmc_cdf(rho, s, t1);
    gap = uq - lq;
    prob = base + gap;
    assert(!isnan(prob) and prob >= 0 and prob <= 1);
    return prob;
}

void PSMC::compute_forward_diagonals(int i) {
    double diag;
    for (int k = 0; k < num_temporal_bins; k++) {
        diag = psmc_prob(rhos[i], mid_points[k], temporal_grids[k], temporal_grids[k+1]);
        diagonals[k] = diag;
    }
}

void PSMC::compute_forward_lower_diagonals(int i) {
    upper_diagonals[0] = 0;
    double lower_diag;
    for (int k = 0; k < num_temporal_bins - 1; k++) {
        lower_diag = psmc_prob(rhos[i], mid_points[k+1], temporal_grids[k], temporal_grids[k+1]);
        lower_diagonals[k] = lower_diag;
    }
}

void PSMC::compute_forward_upper_diagonals(int i) {
    upper_diagonals[0] = 0;
    double upper_diag;
    for (int k = 1; k < num_temporal_bins; k++) {
        upper_diag = psmc_prob(rhos[i], mid_points[k-1], temporal_grids[k], temporal_grids[k+1]);
        upper_diagonals[k] = upper_diag;
    }
}

void PSMC::compute_forward_lower_sums(int i) {
    lower_sums[0] = 0;
    for (int k = 1; k < num_temporal_bins; k++) {
        lower_sums[k] = upper_diagonals[k]*forward_probs[i][k-1] + forward_factors[k]*lower_sums[k-1];
    }
}

void PSMC::compute_forward_upper_sums(int i) {
    partial_sum(forward_probs[i].rbegin(), forward_probs[i].rend()-1, upper_sums.rbegin()+1);
}

void PSMC::compute_backward_diagonals(int j) {
    double diag;
    for (int k = 0; k < num_temporal_bins; k++) {
        diag = psmc_prob(rhos[j], mid_points[k], temporal_grids[k], temporal_grids[k+1]);
        diagonals[k] = diag;
    }
}

void PSMC::compute_backward_lower_diagonals(int j) {
    upper_diagonals[0] = 0;
    double upper_diag;
    for (int k = 0; k < num_temporal_bins - 1; k++) {
        upper_diag = psmc_prob(rhos[j], mid_points[k+1], temporal_grids[k], temporal_grids[k+1]);
        upper_diagonals[k] = upper_diag;
    }
}

void PSMC::compute_backward_upper_diagonals(int j) {
    upper_diagonals[0] = 0;
    double upper_diag;
    for (int k = 1; k < num_temporal_bins; k++) {
        upper_diag = psmc_prob(rhos[j], mid_points[k-1], temporal_grids[k], temporal_grids[k+1]);
        upper_diagonals[k] = upper_diag;
    }
}

void PSMC::compute_backward_lower_sums(int j) {
    lower_sums[0] = 0;
    for (int k = 1; k < num_temporal_bins; k++) {
        lower_sums[k] = lower_sums[k-1] + backward_probs[j][k-1]*lower_diagonals[k-1];
    }
}

void PSMC::compute_backward_upper_sums(int j) {
    upper_sums[num_temporal_bins - 1] = 0;
    for (int k = num_temporal_bins - 2; k >= 0; k--) {
        upper_sums[k] = upper_sums[k+1]*backward_factors[k] + backward_probs[j][k+1]*upper_diagonals[k+1];
    }
}

void PSMC::compute_forward_factors() {
    forward_factors[0] = 0;
    for (int j = 1; j < num_temporal_bins; j++) {
        forward_factors[j] = (exp(-temporal_grids[j]) - exp(-temporal_grids[j+1]))/(exp(-temporal_grids[j-1]) - exp(-temporal_grids[j]));
    }
}

void PSMC::compute_backward_factors() {
    backward_factors[num_temporal_bins - 1] = 0;
    for (int j = 0; j < num_temporal_bins - 1; j++) {
        backward_factors[j] = psmc_prob(rhos[0], mid_points[j], temporal_grids[num_temporal_bins - 1], temporal_grids[num_temporal_bins])/psmc_prob(rhos[0], mid_points[j+1], temporal_grids[num_temporal_bins - 1], temporal_grids[num_temporal_bins]);
    }
}

void PSMC::forward_transition(int i) {
    compute_forward_diagonals(rhos[i]);
    compute_forward_lower_diagonals(rhos[i]);
    compute_forward_upper_diagonals(rhos[i]);
    compute_forward_lower_sums(i);
    compute_forward_upper_sums(i);
    for (int k = 0; k < num_temporal_bins; k++) {
        forward_probs[i+1][k] = diagonals[k]*forward_probs[i][k] + lower_diagonals[k]*upper_sums[k] + lower_sums[k];
    }
    // cout << accumulate(forward_probs[i+1].begin(), forward_probs[i+1].end(), 0.0) << endl;
}

void PSMC::backward_transition(int j) {
    compute_backward_diagonals(rhos[j]);
    compute_backward_lower_diagonals(rhos[j]);
    compute_backward_upper_diagonals(rhos[j]);
    compute_backward_lower_sums(j);
    compute_backward_upper_sums(j);
    for (int k = 0; k < num_temporal_bins; k++) {
        backward_probs[j-1][k] = diagonals[k]*backward_probs[j][k] + upper_sums[k] + lower_sums[k];
    }
    cout << accumulate(backward_probs[j-1].begin(), backward_probs[j-1].end(), 0.0) << endl;
}

void PSMC::write_posterior_average(string filename) {
    ofstream file;
    file.open(filename);
    for (int i = 0; i < spatial_grids.size() - 1; i++) {
        file << spatial_grids[i] <<  " " << spatial_grids[i+1] << " " << posterior_averages[i] << endl;
    }
    return;
}
