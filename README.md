# SwarmTCR
## Purpose
SwarmTCR predicts T-cell receptor (TCR) specificity using the 'distance' between TCRs. SwarmTCR is an optimized adaption of TCRdist (Dash et al. 2017) that uses both protein sequence identity of the complementary determining region (CDR) loops and particle swarm optimization. Distance is determined in the same manner as the original TCRdist methodology, where the alignment (BLOSUM62 matrix) values at each CDR loop are multiplied by a weight, the sum of all CDR loop values is the distance. Unlike TCRdist, SwarmTCR optimizes the weight for each CDR loop, making it specific to the repertoire tested. SwarmTCR is functional on both single-cell and deep-sequencing data, input parameters for each will be detailed below.

## Installation
```console
cd ~/SwarmTCR
g++  -Wall -O3 -fopenmp -o swarmTCR util.C pso.C swarmTCR.C
```
The ```-fopenmp``` option allows the computation of the sequence alignments and optimization to run in parallel on a multicore machine using the OpenMP library. It is optional but highly recommended as it speed things up. For most Mac OS X users, the OpenMP library needs to be installed separately. To do that, install Homebrew (https://brew.sh/) if not already installed, and do the following:

```console
brew install llvm
brew install libomp
/usr/local/Cellar/llvm/12.0.0_1/bin/clang++  -Wall -O3 -fopenmp -o swarmTCR util.C pso.C swarmTCR.C -L/usr/local/Cellar/libomp/12.0.0/lib
```
12.0.0 in the last command needs to be replaced with the correct version in your system.

## Single-cell usage
```console
cd ~/SwarmTCR
cd Example_data
./swarmTCR -r <reference set> -i <training sample set> -t <test sample set> -1 <TCRdist output file> -2 <SwarmTCR output file> -n <number of iterations> -s <swarm size>
```
## Single-cell example usage
```console
cd ~/SwarmTCR
cd Example_data
./swarmTCR -r SC_GILGFVFTL_Reference.txt -i SC_GILGFVFTL_Train_Sample.txt -t SC_GILGFVFTL_Test_Sample.txt -1 TCRdist_out.txt -2 SwarmTCR_out.txt -n 20 -s 25
```

## Bulk sequencing usage
```console
cd ~/SwarmTCR
cd Example_data
./swarmTCR -r <training reference set> -i <training sample set> -t <test sample set> -x <test reference set> -1 <TCRdist output file> -2 <SwarmTCR output file> -n <number of iterations> -s <swarm size>
```
## Bulk sequencing example usage
```console
cd ~/SwarmTCR
cd Example_data
./swarmTCR -r BS_YVLDHLIVV_Train_Reference.txt -i BS_YVLDHLIVV_Train_Sample.txt -t BS_YVLDHLIVV_Test_Sample.txt -x BS_YVLDHLIVV_Test_Reference.txt -1 TCRdist_out.txt -2 SwarmTCR_out.txt -n 20 -s 25
```

The ```-x``` is only used if a test reference set is availabe. In our implementation this was only used for deep sequencing (see methods).

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

## Using pre-trained weights to perform classification
In order to use SwarmTCR on previously uncharacterized TCR sequences,
you can use the program SwarmTCRClassify.

** Compiling SwarmTCRClassify **

To compile the program, simply type:
```
g++  -Wall -O3 -fopenmp -o swarmTCRClassify util.C pso.C swarmTCRClassify.C
```

** Running SwarmTCRClassify **

The program can be run like this:
```
./swarmTCRClassify -w REFERENCE_WEIGHTS -i INPUT_TCRs
```

These two files are as follows:
- REFERENCE_WEIGHTS. This file should contain a header on the first line, and in the following lines it should have the full filename (i.e., path + filename) of each reference file, epitope name, and pre-trained weights. For example:
```
FILE_NAME EPITOPE CDR1A_W CDR2A_W CDR2.5A_W CDR3A_W CDR1B_W CDR2B_W CDR2.5B_W CDR3B_W
./GILGFVFTL/GILGFVFTL_Classifier-Reference.txt GILGFVFTL 0.140705 0 0 0 1 0 0.959139 1
./GLCTLVAML/GLCTLVAML_Classifier-Reference.txt GLCTLVAML 1 0.074943 0.611607 0.989865 1 1 0.314176 1
```
- INPUT_TCRs. This file lists the TCRs we want to classify and is structured as:
```
TCR_ID CDR1A CDR2A CDR2.5A CDR3A CDR1B CDR2B CDR2.5B CDR3B
iedb_ylepgpvtv_111 DSAIYN IQSSQRE DKSSGR AVLSSGGSNYKLT SGHTA FQGTGA PEGSV ASSFIGGTDTQY
iedb_ylepgpvtv_112 DSAIYN IQSSQRE DKSSGR AVLSSGGSNYKLTF SGHTA FQGTGA PEGSV ASSFIGGTDTQYF
```
Please notice that both files begin with a header line.

The program writes to standard output the TCR_ID, the best (least negative) score, the corresponding epitope, and best matching TCR in the reference. You can use the second column (score) to threshold the results. High-confidence results will have scores close to zero.
