//
//  Decoder.hpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/10/24.
//

#ifndef Decoder_hpp
#define Decoder_hpp

#include <stdio.h>
#include "PSMC.hpp"
#include "Scaler.hpp"

class Decoder {
    
public:
    
    Decoder();
    
    void full_decode(string vcf_filename, string output_filename, string mut_map_filename);
    
    void posterior_average_decode(string vcf_filename, string output_filename, string mut_map_filename);
    
};

#endif /* Decoder_hpp */
