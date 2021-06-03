# SwarmTCR
## Purpose
SwarmTCR predicts T-cell receptor (TCR) specificity using the 'distance' between TCRs. SwarmTCR is an optimized adaption of TCRdist (Dash et al. 2017) that uses both protein sequence identity of the complementary determining region (CDR) loops and particle swarm optimization. Distance is determined in the same manner as the original TCRdist methodology, where the alignment (BLOSUM62 matrix) values at each CDR loop are multiplied by a weight, the sum of all CDR loop values is the distance. Unlike TCRdist, SwarmTCR optimizes the weight for each CDR loop, making it specific to the repertoire tested. SwarmTCR is functional on both single-cell and deep-sequencing data, input parameters for each will be detailed below.

## Installation
```console
cd ~/SwarmTCR
g++  -Wall -O3 -fopenmp -o swarmTCR util.C pso.C swarmTCR.C
```
The ```-fopenmp``` option allows the computation of the sequence alignments and optimization to run in parallel on a multicore machine using the OpenMP library. It is optional but highly recommended as it speed things up.

## Single-cell usage
```console
cd ~/SwarmTCR
./swarmTCR -r <reference set> -i <training sample set> -t <test sample set> -1 <TCRdist output file> -2 <SwarmTCR output file> -n <number of iterations> -s <swarm size>
```
## Single-cell example usage
```console
cd ~/SwarmTCR
./swarmTCR -r C_ref_Jun19_gilgfvftl_0.txt -i C_tra_Jun19_gilgfvftl_0.txt-t C_val_Jun19_gilgfvftl_0.txt -1 TCRdist_out.txt -2 SwarmTCR_out.txt -n 20 -s 25
```

## Bulk sequencing usage
```console
cd ~/SwarmTCR
./swarmTCR -r <training reference set> -i <training sample set> -t <test sample set> -x <test reference set> -1 <TCRdist output file> -2 <SwarmTCR output file> -n <number of iterations> -s <swarm size>
```
## Bulk sequencing example usage
```console
cd ~/SwarmTCR
./swarmTCR -r C_refTra_Jun19_gilgfvftl_0.txt -i C_tra_Jun19_gilgfvftl_0.txt -t C_tes_Jun19_gilgfvftl_0.txt -x C_refTes_Jun19_gilgfvftl_0.txt -1 TCRdist_out.txt -2 SwarmTCR_out.txt -n 20 -s 25
```

The ```-x``` is only used if a test reference set is availabe. In our implementation this was only used for deep-sequencing (see methods).

**Input files**

Each row in this file contains a TCR ID, boolean flag to infer a positve or negative label, and complete CDR loop protein data. The flag is either a ```1``` or ```0```, ```1``` indicates a positive label, ```0``` a negative label.

**Bulk sequencing input**

Bulk sequencing input will use either the alpha or beta chain only, there is no need to specify chain type.
```
<TCR id> <flag> <CDR1> <CDR2> <CDR2.5> <CDR3>
...
```

**Single-cell input**

Single-cell input uses both alpha and beta chains, alpha CDRs are listed first, then beta.
```
<TCR id> <flag> <CDR1_alpha> <CDR2_alpha> <CDR2.5_alpha> <CDR3_alpha> <CDR1_beta> <CDR2_beta> <CDR2.5_beta> <CDR3_beta>
...
```

**Output files**

There are a total of 3 ouput files:
1) TCRdist precision recall values
2) SwarmTCR precision recall values
3) SwarmTCR scores (```SwarmTCR output file```.scores). This file contains the nearest-neighbor TCR ID, label, and SwarmTCR distance. 

Precision recall files are formatted as follows:

Column 1: Precision values

Column 2: Recall values
  
  
Score file is formatted as follows:

Column 1: Nearest-neighbor match (TCR ID) 

Column 2: Binary label (0 = negative, 1 = positive)

Column 3: Distance value (actual distance = distance * -1)

