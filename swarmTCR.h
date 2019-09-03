// This file is part of the swarmTCR program
// Copyright (c) 2018 Dario Ghersi and Ryan Ehrlich
// Version 20181223

using namespace std;

#ifndef _swarmTCR_
#define _swarmTCR_

//////////////////////////////////////////////////////////////////////
// CONSTANTS                                                        //
//////////////////////////////////////////////////////////////////////

#define USAGE "\nUsage: swarmTCR -r REFERENCE -i TRAINING -t TEST -n ITERATIONS -s SWARM_SIZE -1 OUTPUT_STD -2 OUTPUT_OPT [-x TEST_REFERENCE]\n"
#define NUMCDR 8 // maximum number of CDR regions (single-cell seq.)
#define GAPPENCDR3 -8
#define GAPPENCDR12 -4
#define PRECISION_BINS 350
#define BLOSUM_SIZE 24

//////////////////////////////////////////////////////////////////////
// STRUCTURES                                                       //
//////////////////////////////////////////////////////////////////////

struct PrecRec {
  double precision;
  double recall;
};

//////////////////////////////////////////////////////////////////////

struct Results {
  double *scores;
  int *labels;
  unsigned int numCells;
};
  
//////////////////////////////////////////////////////////////////////

struct SingleCell {
  char **sampleID;
  unsigned int *epitope; // 1 for the epitope of interest, 0 otherwise
  char **cdrSeqs[NUMCDR];
  unsigned int **cdrSeqsNum[NUMCDR];
  unsigned int *lengthCDR[NUMCDR];
  unsigned int numCells;
};

//////////////////////////////////////////////////////////////////////
// CLASSES                                                          //
//////////////////////////////////////////////////////////////////////

class Parameters {

 public:
  char *refSetFileName; // reference set
  char *trainingSetFileName; // training set
  char *testRefFileName; // optional reference file for testing
  char *testSetFileName; // test set
  char *outputStdFileName; // output PR file for standard weights
  char *outputOptFileName; // output PR file for optimized weights
  unsigned int numIterations; // PSO parameter
  unsigned int swarmSize; // PSO parameter
  bool useTestRef; // boolean for the optional reference test file

  Parameters(char **, int);
};

#endif
