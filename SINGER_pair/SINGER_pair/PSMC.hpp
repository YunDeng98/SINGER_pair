//
//  PSMC.hpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/10/24.
//

#ifndef PSMC_hpp
#define PSMC_hpp

#include <stdio.h>
#include <cmath>
#include <cassert>
#include "Rate_map.hpp"
#include "VCF_reader.hpp"

class PSMC {
    
public:
    
    double Ne;
    double unit_theta = 2.5e-2;
    int num_temporal_bins = 30;
    int bin_size;
    double start;
    double end;
    vector<double> temporal_grids;
    vector<double> mid_points;
    
    vector<double> diagonals;
    vector<double> upper_diagonals;
    vector<double> lower_diagonals;
    vector<double> lower_sums;
    vector<double> upper_sums;
    vector<double> forward_factors;
    vector<double> backward_factors;
    
    vector<vector<double>> forward_probs;
    vector<vector<double>> backward_probs;
    
    vector<double> spatial_grids;
    vector<double> num_het_sites;
    vector<double> thetas;
    vector<double> rhos;
    vector<vector<double>> emission_probs;
    
    vector<double> posterior_averages;
    
    PSMC();
    
    void initialize(VCF_reader &vr, Rate_map &rm);
    
    void create_spatial_bins();
    
    void create_temporal_bins();
    
    void compute_num_het_sites(VCF_reader &vr);
    
    void compute_thetas(Rate_map &rm);
    
    void compute_rhos(double m);
    
    void compute_emission_probs();
    
    void forward_algorithm();
    
    void backward_algorithm();
    
    void posterior_average();

// private:
    
    double psmc_cdf(double rho, double s, double t);
    
    double psmc_prob(double rho, double s, double t1, double t2);
    
    void compute_forward_diagonals(int i);
    
    void compute_forward_lower_diagonals(int i);
    
    void compute_forward_upper_diagonals(int i);
    
    void compute_forward_lower_sums(int i);
    
    void compute_forward_upper_sums(int i);
    
    void compute_backward_diagonals(int i);
    
    void compute_backward_lower_diagonals(int i);
    
    void compute_backward_upper_diagonals(int i);
    
    void compute_backward_lower_sums(int i);
    
    void compute_backward_upper_sums(int i);
    
    void compute_forward_factors();
    
    void compute_backward_factors();
    
    void forward_transition(int i);
    
    void backward_transition(int j);
    
    void forward_emit(int i);
    
    void backward_emit(int j);
    
    void write_posterior_average(string filename);
};

#endif /* PSMC_hpp */
