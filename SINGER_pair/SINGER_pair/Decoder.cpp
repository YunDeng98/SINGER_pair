//
//  Decoder.cpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/10/24.
//

#include "Decoder.hpp"

Decoder::Decoder() {}

void Decoder::posterior_average_decode(string vcf_filename, string output_filename, string mut_map_filename) {
    VCF_reader vr = VCF_reader();
    Rate_map rm = Rate_map();
    rm.load_map(mut_map_filename);
    double m = rm.mean_rate();
    vr.read_vcf(vcf_filename);
    vr.compute_Ne(m);
    PSMC psmc = PSMC();
    psmc.initialize(vr, rm);
    psmc.compute_num_het_sites(vr);
    psmc.compute_emission_probs();
    psmc.forward_algorithm();
    psmc.backward_algorithm();
    psmc.posterior_average();
    /*
    for (int i = 0; i < 5; i++) {
        Scaler scaler = Scaler();
        scaler.rescale(psmc);
        psmc.posterior_average();
    }
     */
    psmc.write_posterior_average(output_filename);
}
