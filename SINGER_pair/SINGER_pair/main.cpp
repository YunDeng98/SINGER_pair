//
//  main.cpp
//  SINGER_pair
//
//  Created by Yun Deng on 6/10/24.
//

#include <iostream>
#include "Decoder.hpp"

int main(int argc, const char * argv[]) {
    string vcf_filename, output_filename, mut_map_filename;
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        if (arg == "-vcf_filename") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -vcf_filename flag cannot be empty. " << endl;
                exit(1);
            }
            vcf_filename = argv[++i];
        }
        else if (arg == "-output_filename") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -output_filename flag cannot be empty. " << endl;
                exit(1);
            }
            output_filename = argv[++i];
        }
        else if (arg == "-mut_map_filename") {
            if (i + 1 > argc || argv[i+1][0] == '-') {
                cerr << "Error: -mut_map_filename flag cannot be empty. " << endl;
                exit(1);
            }
            mut_map_filename = argv[++i];
        }
        else {
            cerr << "Error: Unknown flag. " << arg << endl;
            exit(1);
        }
    }
    Decoder decoder = Decoder();
    decoder.posterior_average_decode(vcf_filename, output_filename, mut_map_filename);
    return 0;
}
