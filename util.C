// This file is part of the swarmTCR program
// Copyright (c) 2020 Dario Ghersi and Ryan Ehrlich
// Version 20200227

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string.h>
#include "util.h"
#include "BLOSUM62.h"
#include "swarmTCR.h"

using namespace std;

extern unsigned int CURRCDR;
double epsilon = numeric_limits<double>::epsilon();

//////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                        //
//////////////////////////////////////////////////////////////////////

double averagePrecision(struct Results *res, char *outFileName)
{
  // calculate the average precision
 
  double ap = 0.0;

  // open the output file
  fstream outFile;
  if (outFileName != 0) {
    outFile.open(outFileName, fstream::out);
    outFile << "Precision\tRecall\n";
  }
  
  // get the range of scores
  double min = res->scores[0], max = res->scores[0];
  for (unsigned int i = 1; i < res->numCells; i++) {
    
    if (min > res->scores[i]) {
      min = res->scores[i];
    }
    
    if (max < res->scores[i]) {
      max = res->scores[i];
    }
  }

  // calculate the increment
  double threshold = max;
  double increment = (max - min) / PRECISION_BINS;
  struct PrecRec precRec;
  if (increment < 0) {
    increment = - increment;
  }

  // calculate the average precision
  double oldRec = 0.0;
  double oldPrec = 0.0;
  for (unsigned int i = 0; i <= PRECISION_BINS; i++) {
    calculatePrecRecall(res, threshold, &precRec);
    ap += precRec.precision * (precRec.recall - oldRec);

    // (optional) print the precision/recall to file
    if (outFileName != 0 && ((fabs(oldRec - precRec.recall) > epsilon) ||
			     (fabs(oldPrec - precRec.precision) > epsilon))) {
      outFile << precRec.precision << "\t" << precRec.recall << endl;
    }
    oldRec = precRec.recall;
    oldPrec = precRec.precision;
    threshold -= increment;
  }

  // close the output file
  if (outFileName != 0) {
    outFile.close();
  }
  
  return ap;
}

//////////////////////////////////////////////////////////////////////

void calculatePrecRecall(struct Results *res, double threshold,
                          struct PrecRec *precRec)
{
  // calculate precision and recall

  unsigned int tp = 0, numRetrieved = 0, numRelevant = 0;

  for (unsigned int i = 0; i < res->numCells; i++) {
    if (res->scores[i] >= threshold) {
      numRetrieved += 1;
      if (res->labels[i] == 1) {
        tp++;
      }
    }
    if (res->labels[i] == 1) {
      numRelevant++;
    }
  }

  precRec->precision = (double) tp / numRetrieved;
  precRec->recall = (double) tp / numRelevant;
}

//////////////////////////////////////////////////////////////////////

void checkCommandLineArgs(char **argv, int argc)
{
  // check all the parameters have been provided

  bool err = false;

  if (!cmdOptionExists(argv, argv+argc, "-r")) {
    cerr << "Reference set file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-i")) {
    cerr << "Training set file missing\n";
  }
  if (!cmdOptionExists(argv, argv+argc, "-t")) {
    cerr << "Test set file missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-n")) {
    cerr << "Max. number of iterations missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-s")) {
    cerr << "Swarm size missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-1")) {
    cerr << "Output file for standard weights missing\n";
    err = true;
  }
  if (!cmdOptionExists(argv, argv+argc, "-2")) {
    cerr << "Output file for optimized weights missing\n";
    err = true;
  }
 
  if (err) {
    cout << USAGE;
    exit(1);
  }
}

//////////////////////////////////////////////////////////////////////

bool cmdOptionExists(char **begin, char **end,
                     const string & option)
{
  return std::find(begin, end, option) != end;
}

//////////////////////////////////////////////////////////////////////

void convertAA2Num(char *string, unsigned int *numbers)
{
  // Convert a string of amino acids to numbers, using
  // the indexes in BLOSUM62.h

  bool found = false;

  for (unsigned int i = 0; i < strlen(string); i++) {
    found = false;
    for (unsigned int j = 0; j < NUM_AA; j++) { 
      if (string[i] == BLOSUM_AA[j]) {
        numbers[i] = j;
        found = true;
      }
    }
    if (!found) {
      cerr << "Unrecognized amino acid: " << string[i] << endl;
      exit(1);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void freeSingleCellData(struct SingleCell *sc)
{
  // free the memory allocated for single cell structures

  free(sc->epitope);
  for (unsigned int i = 0; i < sc->numCells; i++) {
    free(sc->sampleID[i]);

    for (unsigned int j = 0; j < CURRCDR; j++) {
      free(sc->cdrSeqs[j][i]);
      free(sc->cdrSeqsNum[j][i]);
    }
  }

  for (unsigned int j = 0; j < CURRCDR; j++) {
    free(sc->lengthCDR[j]);
  }

  free(sc->sampleID);
}

//////////////////////////////////////////////////////////////////////

char *getCmdOption(char **begin, char **end,
                   const string & option)
{
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////

int globalAlign(unsigned int *seqA, unsigned int lenA,
		unsigned int *seqB, unsigned int lenB, int gapPen,
		int **MAT)
{
  // Apply the Needleman-Wunsch algorithm to align two sequences
  // and return the optimal alignment score

  int **dpMat;
  int globalScore;
  
  // allocate memory for the score matrix
  dpMat = (int **) malloc(sizeof(int *) * (lenB + 1));
  for (unsigned int i = 0; i <= lenB; i++) {
    dpMat[i] = (int *) malloc(sizeof(int) * (lenA + 1));
  }

  // initialize the matrix
  for (unsigned int i = 0; i <= lenA; i++) {
    dpMat[0][i] = i * gapPen;
  }

  for (unsigned int j = 0; j <= lenB; j++) {
    dpMat[j][0] = j * gapPen;
  }

  // dynamic programming step
  for (unsigned int i = 1; i <= lenB; i++) {
    for (unsigned int j = 1; j <= lenA; j++) {
      dpMat[i][j] = maximum(dpMat[i - 1][j - 1] +\
			    MAT[seqB[i - 1]][seqA[j - 1]],
                            dpMat[i - 1][j] + gapPen,\
			    dpMat[i][j - 1] + gapPen);
    }
  }

  globalScore = dpMat[lenB][lenA];

  // free memory for the alignment matrix
  for (unsigned int i = 0; i <= lenB; i++) {
    free(dpMat[i]);
  }
  free(dpMat);

  return globalScore;
}

//////////////////////////////////////////////////////////////////////

int maximum(int a, int b, int c)
{
  // return the maximum of three numbers

  return max(max(a, b), c);
}

//////////////////////////////////////////////////////////////////////

void modifyBLOSUM(int **newBLOSUM) {
  // build a modified blosum matrix as described in Paul Thomas 2017

  for (unsigned int i = 0; i < BLOSUM_SIZE; i++) {
    for (unsigned int j = 0; j < BLOSUM_SIZE; j++) {
      if (i == j) {
        newBLOSUM[i][j] = 0;
      }
      else if (BLOSUM62[i][j] < 0) {
	      newBLOSUM[i][j] = -4;
      }
      else {
	      newBLOSUM[i][j] = -(4 - BLOSUM62[i][j]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

struct Results *nearestNeighbor(struct SingleCell *ref, struct SingleCell *sc, 
                                double *weights, int **MAT)
{
  // return the nearest neighbor for each receptor

  double score = 0.0;

  // memory allocation for the results
  struct Results *res;
  res = (struct Results *) malloc(sizeof(struct Results));
  res->scores = (double *) malloc(sizeof(double) * sc->numCells);
  res->labels = (int *) malloc(sizeof(int) * sc->numCells);
  
  res->numCells = sc->numCells;

  // process each cell
  for (unsigned int i = 0; i < sc->numCells; i++) {

    double maxScore = -10000.0;
    double cdrScore;
    double gapPenalty = GAPPENCDR12;

    for (unsigned int j = 0; j < ref->numCells; j++) {
      score = 0.0;
      for (unsigned int region = 0; region < CURRCDR; region++) {
	  
	      if (region == 3 || region == 7) { // region-specific gap pen.
	        gapPenalty = GAPPENCDR3;
	      }
	      else {
	        gapPenalty = GAPPENCDR12;
	      }

        // compute the alignment score for the region
	      cdrScore = globalAlign(sc->cdrSeqsNum[region][i],\
                               sc->lengthCDR[region][i],\
                               ref->cdrSeqsNum[region][j],\
                               ref->lengthCDR[region][j], gapPenalty, MAT);

        score += cdrScore * weights[region];
      }
    
      if (score > maxScore) {
        maxScore = score;
      }
    }

    // update the results
    res->labels[i] = sc->epitope[i];
    res->scores[i] = maxScore;
  }

  return res;
}

//////////////////////////////////////////////////////////////////////

struct SingleCell *readSingleCellData(string infileName)
{
  // store single cell data for each receptor

  string line;
  unsigned int cellCount = 0, cdrCount = 0;
  struct SingleCell *sc;
  sc = (struct SingleCell *) malloc(sizeof(struct SingleCell));

  // open the input file
  fstream infile;
  infile.open(infileName.c_str(), fstream::in);

  // complain if the file doesn't exist
  if (! infile.good()) {
    cerr << "Can't open " << infileName << endl;
    exit(1);
  }

  // count the cdr regions (single cell vs. bulk sequencing)
  string header, token;
  CURRCDR = -2; // to avoid counting the SAMPLE_ID and EPITOPE as CDRs
  getline(infile, header);
  stringstream ss(header);
  while (getline(ss, token, ' ')) {
    CURRCDR++;
  }

  // count how many receptors there are
  unsigned int numRec = 0;
  while (getline(infile, line)) {
    numRec++;
  }

  sc->numCells = numRec;

  // allocate the memory needed to store the receptors
  sc->sampleID = (char **) malloc(sizeof(char *) * numRec);
  sc->epitope = (unsigned int *) malloc(sizeof(unsigned int) * numRec);
  for (unsigned int i = 0; i < CURRCDR; i++) {
    sc->cdrSeqs[i] = (char **) malloc(sizeof(char *) * numRec);
    sc->cdrSeqsNum[i] = (unsigned int **) malloc(sizeof(unsigned int *) * numRec);
    sc->lengthCDR[i] = (unsigned int *) malloc(sizeof(unsigned int *) * numRec);
  }

  // rewind the file
  infile.clear();
  infile.seekg(0);

  // load the single cell data into the structure
  getline(infile, line); // skip the header

  while (getline(infile, line)) {
    cdrCount = 0;
    stringstream ss(line);

    // sample ID
    getline(ss, token, ' ');
    sc->sampleID[cellCount] = (char *) malloc(sizeof(char) * \
					      token.length() + 1);
    strcpy(sc->sampleID[cellCount], token.c_str());
    
    // epitope
    getline(ss, token, ' ');
    sc->epitope[cellCount] = atoi(token.c_str());
    
    while (getline(ss, token, ' ')) {
      sc->cdrSeqs[cdrCount][cellCount] = (char *) malloc(sizeof(char) * \
							 token.length() + 1);
      sc->cdrSeqsNum[cdrCount][cellCount] = (unsigned int *) \
	malloc(sizeof(unsigned int) * token.length());

      strcpy(sc->cdrSeqs[cdrCount][cellCount], token.c_str()); // add the string of the CDR
      sc->lengthCDR[cdrCount][cellCount] = token.length();

      convertAA2Num(sc->cdrSeqs[cdrCount][cellCount],\
		    sc->cdrSeqsNum[cdrCount][cellCount]); // add the number version
      
      cdrCount++;
    }

    cellCount++;
  }

  infile.close();

  return sc;
}
