//
//  VCF_reader.hpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/10/24.
//

#ifndef VCF_reader_hpp
#define VCF_reader_hpp

#include <stdio.h>
#include "Rate_map.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

class VCF_reader {
    
public:
    
    Rate_map rm;
    double Ne = 0;
    double start = INT_MAX;
    double end = 0;
    double sequence_length = 0;
    double num_het_sites = 0;
    vector<double> het_sites = {};
    
    void read_vcf(string filename);
    
    void compute_Ne(double m);
    
    vector<double> compute_num_het_sites(vector<double> coordinates);
};

#endif /* VCF_reader_hpp */
