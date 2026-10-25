# Scripts in ToA
Scripts used in the paper "Tree of Archaeplastida as a showcase for systematic artifacts in nuclear and plastid phylogenomics".

## [`Rhodelphis_LowerCov_RANDOM.py`](./Rhodelphis_LowerCov_RANDOM.py)
This script was used to generate the `nuc276ΔP_cov50%` alignments in the paper.
To lower the coverage of *Rhodelphis*, the script works as follows:

1. Gets a list of the indices of valid amino acid positions in the sequence of the target OTU
1. Randomly selects positions from this list using `random.sample()` from the Python standard library
1. Outputs a new sequence in which the non-selected positions are replaced by gaps

The script has the following options:

```shellsession
$ python3 Rhodelphis_LowerCov_RANDOM.py -h
usage: Rhodelphis_LowerCov_RANDOM.py [-h] -f FASTA_FILE -r SELECT_RATE -n TRIAL_NUM [-s SEED] -o OUTPUT_PREFIX

Lower the coverage of the Rhodelphis sequence in a phylogenomic alignment by random selection.

options:
  -h, --help            show this help message and exit
  -f FASTA_FILE, --fasta_file FASTA_FILE
                        input FASTA file
  -r SELECT_RATE, --select_rate SELECT_RATE
                        percentage of valid positions to retain (0-100)
  -n TRIAL_NUM, --trial_num TRIAL_NUM
                        number of replicates to generate (>= 1)
  -s SEED, --seed SEED  seed value for the random sampling
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        prefix of the output files
```

The script was executed as follows:

```shellsession
$ python3 Rhodelphis_LowerCov_RANDOM.py \
     -f nuc276.fasta -r 50 -n 10 -s 100 -o Rhodelphis_RandomSample
```

The actual files used in the analyses are available as supplementary data.
