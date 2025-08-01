## Test datasets

This directory contains small test datasets derived from real samples and [minimal metadata](./test_metadata.xlsx) for the samples.  These datasets include:

| dataset (R1 file) | property |
| ------- | -------- | 
| empty_fastq_R1 | An empty fastq.  Pipeline should not fail on empty input |
| D_andersoni_R1_001| Dermacentor andersoni - no pathogens |
| D_variabilis_R1_001| Dermacentor variabilis - no pathogens |
| I_pacificus_R1_001| Ixodes pacificus - no pathogens |
| I_scapularis_R1_001| Ixodes scapularis - no pathogens |
| Borrelia_burgdoreri_B31_1000_S13_L001_R1_001| A dataset from a tick positive for B. burgdorferi B31, 1000 total reads |
| Borrelia_burgdoreri_B31_100_S13_L001_R1_001| A dataset from a tick positive for B. burgdorferi B31, 100 total reads |
| Borrelia_miyamotoi_R1_001.fastq | A dataset from a tick positive for B. miyamotoi |
| pathogen_pos_S96_L001_R1_001| A dataset from a positive control sample positive for multiple pathogens |

These datasets were collected using the [./create_test_fastq script](./create_test_fastq) in this directory.  Note that these test dataset files are named using various typical Illumina naming conventions (e.g. including L001 or not) and correspond to compressed and non-compressed fastq.

These datasets can be used as input by using the test profile, for example:

```
nextflow run stenglein-lab/tick_surveillance -profile singularity,test
```
