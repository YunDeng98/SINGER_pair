//
//  VCF_reader.cpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/10/24.
//

#include "VCF_reader.hpp"

void VCF_reader::read_vcf(string filename) {
    string vcf_file = filename;
    ifstream file(vcf_file);
    string line;
    int num_individuals = 0;
    double prev_pos = -1;
    vector<double> genotypes = {0, 0};
    while (getline(file, line)) {
        if (line.substr(0, 6) == "#CHROM") {
            istringstream iss(line);
            vector<string> fields;
            string field;
            while (iss >> field) {
                fields.push_back(field);
            }
            num_individuals = (int) fields.size() - 9;
            assert(num_individuals == 1);
            continue;
        } else if (line[0] == '#') {
            continue; // skip these header lines
        }
        istringstream iss(line);
        string chrom, id, ref, alt, qual, filter, info, format, genotype;
        double pos;
        iss >> chrom >> pos >> id >> ref >> alt >> qual >> filter >> info >> format;
        start = min(start, pos);
        end = max(end, pos);
        if (pos == prev_pos) {continue;} // skip multi-allelic sites
        if (ref.size() > 1 or alt.size() > 1) {
            continue;
        } // skip multi-allelic sites or structural variant
        
        streampos old_pos = file.tellg();
        string next_line;
        if (getline(file, next_line)) {
            istringstream next_iss(next_line);
            string next_chrom;
            int next_pos;
            next_iss >> next_chrom >> next_pos;
            if (next_pos == pos) {
                prev_pos = pos;
                continue;
            }
            file.seekg(old_pos);
        }
        while (iss >> genotype) {
            if (genotype[0] == '1') {
                genotypes[0] = 1;
            } else {
                genotypes[0] = 0;
            }
            if (genotype[2] == '1') {
                genotypes[1] = 1;
            } else {
                genotypes[1] = 0;
            }
        }
        int genotype_sum = accumulate(genotypes.begin(), genotypes.end(), 0.0);
        if (genotype_sum == 1) {
            het_sites.push_back(pos);
            num_het_sites += 1;
        }
    }
}

void VCF_reader::compute_Ne(double m) {
    Ne = num_het_sites/(end - start)/m/2;
}

vector<double> VCF_reader::compute_num_het_sites(vector<double> coordinates) {
    int n = (int) coordinates.size() - 1;
    vector<double> counts = vector<double>(n);
    int i = 0;
    int j = 0;
    double lb, ub, pos;
    while (i < n) {
        pos = het_sites[j];
        lb = coordinates[i];
        ub = coordinates[i+1];
        if (pos > lb and pos <= ub) {
            counts[i] += 1;
            j++;
        } else {
            i++;
        }
    }
    return counts;
}
