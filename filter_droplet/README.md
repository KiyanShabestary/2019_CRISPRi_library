## filter_droplet

This script aims to filter sequencing data obtained after microfluidic-based sorting of lactate producers.

Data is treated as follow:

1. Low number of reads (<32) are set to 0.
2. sgRNAs with at least one sample with 10-fold difference within replicates are removed.
3. Relative data is computed per sample.
4. Duplicate sgRNAs are removed.
5. Ratios sorted/unsorted (obtained from gDNA) are computed.
6. sgRNAs in each top fraction are counted to make a final list.

To run the script, go to 2019_CRISPRi_library/filter_droplet and execute the following command in a terminal:

`python filter_droplet.py`

The script uses python3 and requires installation of pandas.
