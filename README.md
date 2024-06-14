# SINGER_pair
SINGER when specifically optimized to 2 sequences 

To use this special version of SINGER on a pair of sequences (including unphased individual), simply do:

```
singer_pair -vcf_filename vcf_file -output_filename output -mut_map_filename mut_map
```
The mutation map file should be of the following format, that each row contains the start of the bin, end of the bin, and the rate in the bin:

```
0 1000 2e-8
1000 2000 1e-8
2000 5000 2e-8
```
the map must **fully cover** the VCF file from the first to the last variant
