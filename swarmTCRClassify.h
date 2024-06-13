// This file is part of the swarmTCRClassify program
// Copyright (c) 2024 Dario Ghersi
// Version 20240609


#ifndef _swarmTCRClassify_
#define _swarmTCRClassify_

//////////////////////////////////////////////////////////////////////
// CONSTANTS                                                        //
//////////////////////////////////////////////////////////////////////

#define USAGE "\nUsage: swarmTCRClassify -w REFERENCE_WEIGHTS -i INPUT_TCRs\n"
#define NUMCDR 8 // maximum number of CDR regions (single-cell seq.)
#define GAPPENCDR3 -8
#define GAPPENCDR12 -4
#define BLOSUM_SIZE 24

/////////////////////////////////////////////////////////////////////
// STRUCTURES                                                      //
/////////////////////////////////////////////////////////////////////

struct TCRData {
  char **sampleID;
  char **cdrSeqs[NUMCDR];
  unsigned int **cdrSeqsNum[NUMCDR];
  unsigned int *lengthCDR[NUMCDR];
  unsigned int numTCRs;
};

struct Results {
  double *scores;
  char **TCRID;
  char *epitope;
  char **refTCRID;
  unsigned int numTCRs;
};

struct RefWeightsData {
  char **fileNames;
  char **epitopes;
  double *weights[NUMCDR];
  unsigned int numRefs;
};

//////////////////////////////////////////////////////////////////////
// CLASSES                                                          //
//////////////////////////////////////////////////////////////////////

class Parameters {

 public:
  char *weightsFileName; // weights for each epitope
  char *inputTCRsFileName; // input TCR set

  Parameters(char **, int);
};

#endif
