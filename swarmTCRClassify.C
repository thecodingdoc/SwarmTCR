// This file is part of the swarmTCR program
// Copyright (c) 2024 Dario Ghersi
// Version: 20240610

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream> 
#include "util.h"
#include "swarmTCRClassify.h"

//////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES                                                 //
//////////////////////////////////////////////////////////////////////

unsigned int CURRCDR = 8; // actual number of cdr regions (8 for
                          // single cell, 4 for bulk sequencing)

//////////////////////////////////////////////////////////////////////
// CONSTRUCTORS                                                     //
//////////////////////////////////////////////////////////////////////

Parameters::Parameters(char **argv, int argc)
{
  // parse the command-line arguments

  // reference set file name
  weightsFileName = getCmdOption(argv, argv + argc, "-w");

  // training set file name
  inputTCRsFileName = getCmdOption(argv, argv + argc, "-i");
}

//////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                        //
//////////////////////////////////////////////////////////////////////

void checkCommandLineArgsClassify(char **argv, int argc)
{
  // check all the parameters have been provided

  bool err = false;

  if (!cmdOptionExists(argv, argv+argc, "-w")) {
    cerr << "Weight file missing\n";
    err = true;
  }

  if (!cmdOptionExists(argv, argv+argc, "-i")) {
    cerr << "Input TCRs file missing\n";
    err = true;
  }
    
  if (err) {
    cout << USAGE;
    exit(1);
  }
}

//////////////////////////////////////////////////////////////////////

void freeTCRData(struct TCRData *input)
{
  // free the memory allocated for input structures

  for (unsigned int i = 0; i < input->numTCRs; i++) {
    free(input->sampleID[i]);

    for (unsigned int j = 0; j < CURRCDR; j++) {
      free(input->cdrSeqs[j][i]);
      free(input->cdrSeqsNum[j][i]);
    }
  }

  for (unsigned int j = 0; j < CURRCDR; j++) {
    free(input->lengthCDR[j]);
  }

  free(input->sampleID);
}

//////////////////////////////////////////////////////////////////////

struct Results *nearestNeighbor(struct TCRData *ref, struct TCRData *input, 
                                double *weights, int **MAT)
{
  // return the nearest neighbor for each receptor

  double score = 0.0;

  // memory allocation for the results
  struct Results *res;
  res = (struct Results *) malloc(sizeof(struct Results));
  res->scores = (double *) malloc(sizeof(double) * input->numTCRs);
  res->refTCRID = (char **) malloc(sizeof(char **) * input->numTCRs);
  res->TCRID = (char **) malloc(sizeof(char **) * input->numTCRs);
  
  res->numTCRs = input->numTCRs;

  // process each cell
  for (unsigned int i = 0; i < input->numTCRs; i++) {

    double maxScore = -10000.0;
    double cdrScore;
    double gapPenalty = GAPPENCDR12;

    // memory allocation for the TCR ids
    res->TCRID[i] =  (char *) malloc(sizeof (char *));
    res->refTCRID[i] =  (char *) malloc(sizeof (char *));

    char *refID = NULL;
    res->TCRID[i] = input->sampleID[i];
    for (unsigned int j = 0; j < ref->numTCRs; j++) {
      score = 0.0;
      for (unsigned int region = 0; region < CURRCDR; region++) {
	  
	      if (region == 3 || region == 7) { // region-specific gap pen.
	        gapPenalty = GAPPENCDR3;
	      }
	      else {
	        gapPenalty = GAPPENCDR12;
	      }

	      //  compute the alignment score for the region
	      cdrScore = globalAlign(input->cdrSeqsNum[region][i],\
                               input->lengthCDR[region][i],\
                               ref->cdrSeqsNum[region][j],\
                               ref->lengthCDR[region][j], gapPenalty, MAT);

        score += cdrScore * weights[region];
      }
    
      if (score > maxScore) {
        maxScore = score;
	refID = ref->sampleID[j];
      }
    }

    // update the results
    res->scores[i] = maxScore;
    res->refTCRID[i] = refID;
  }

  return res;
}

//////////////////////////////////////////////////////////////////////

struct TCRData *readTCRData(string infileName)
{
  // store TCR data that needs to be classified

  string line;
  unsigned int tcrCount = 0, cdrCount = 0;
  struct TCRData *sc;
  sc = (struct TCRData *) malloc(sizeof(struct TCRData));

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
  CURRCDR = -1; // to avoid counting the SAMPLE_ID as CDRs
  getline(infile, header);
  stringstream ss(header);
  while (getline(ss, token, ' ')) {
    CURRCDR++;
  }

  // count how many receptors there are
  unsigned int numTCRs = 0;
  while (getline(infile, line)) {
    numTCRs++;
  }

  sc->numTCRs = numTCRs;
  // allocate the memory needed to store the receptors
  sc->sampleID = (char **) malloc(sizeof(char *) * numTCRs);
  for (unsigned int i = 0; i < CURRCDR; i++) {
    sc->cdrSeqs[i] = (char **) malloc(sizeof(char *) * numTCRs);
    sc->cdrSeqsNum[i] = (unsigned int **) malloc(sizeof(unsigned int *) * numTCRs);
    sc->lengthCDR[i] = (unsigned int *) malloc(sizeof(unsigned int *) * numTCRs);
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
    sc->sampleID[tcrCount] = (char *) malloc(sizeof(char) * \
                                              token.length() + 1);
    strcpy(sc->sampleID[tcrCount], token.c_str());
    
    while (getline(ss, token, ' ')) {
      sc->cdrSeqs[cdrCount][tcrCount] = (char *) malloc(sizeof(char) * \
                                                         token.length() + 1);
      sc->cdrSeqsNum[cdrCount][tcrCount] = (unsigned int *) \
        malloc(sizeof(unsigned int) * token.length());

      strcpy(sc->cdrSeqs[cdrCount][tcrCount], token.c_str()); // add the string of the CDR
      sc->lengthCDR[cdrCount][tcrCount] = token.length();

      convertAA2Num(sc->cdrSeqs[cdrCount][tcrCount],\
                    sc->cdrSeqsNum[cdrCount][tcrCount]); // add the number version
      
      cdrCount++;
    }

    tcrCount++;
  }

  infile.close();

  return sc;
}

//////////////////////////////////////////////////////////////////////

struct RefWeightsData *readRefWeightsData(string infileName) {

  string line;
  unsigned int refCount = 0, cdrCount = 0;
  struct RefWeightsData *refWeights;
  refWeights = (struct RefWeightsData *) malloc(sizeof(struct RefWeightsData));

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
  CURRCDR = -2; // to avoid counting the filename and epitope as CDRs
  getline(infile, header);
  stringstream ss(header);
  while (getline(ss, token, ' ')) {
    CURRCDR++;
  }

  // count how many receptors there are
  unsigned int numRefs = 0;
  while (getline(infile, line)) {
    numRefs++;
  }

  refWeights->numRefs = numRefs;
  // allocate the memory needed to store the reference files and weights
  refWeights->fileNames = (char **) malloc(sizeof(char *) * numRefs);
  refWeights->epitopes = (char **) malloc(sizeof(char *) * numRefs);
  for (unsigned int i = 0; i < CURRCDR; i++) {
    refWeights->weights[i] = (double *) malloc(sizeof(double) * numRefs);
  }

  // rewind the file
  infile.clear();
  infile.seekg(0);

  // load the single cell data into the structure
  getline(infile, line); // skip the header

  refCount = 0;
  while (getline(infile, line)) {
    cdrCount = 0;
    stringstream ss(line);

    // filename
    getline(ss, token, ' ');
    refWeights->fileNames[refCount] = (char *) malloc(sizeof(char) * \
						      token.length() + 1);
    strcpy(refWeights->fileNames[refCount], token.c_str());

    // epitope
    getline(ss, token, ' ');
    refWeights->epitopes[refCount] = (char *) malloc(sizeof(char) * \
						     token.length() + 1);
    strcpy(refWeights->epitopes[refCount], token.c_str());

    // weights
    while (getline(ss, token, ' ')) {
      sscanf(token.c_str(), "%lf", &refWeights->weights[cdrCount][refCount]);

      cdrCount++;
    }

    refCount++;
  }

  infile.close();

  return refWeights;
}

//////////////////////////////////////////////////////////////////////

void printBestMatch(struct Results **res, unsigned int numRefs) {

  // print the header
  cout << "TCR_ID EPITOPE SCORE REF_TCR_ID" << endl;

  // cycle through each TCR
  for (unsigned int i = 0; i < res[0]->numTCRs; i++) {
    unsigned int index = -1;
    double bestScore = -1E6;

    // find the highest score
    for (unsigned int j = 0; j < numRefs; j++) {
      if (res[j]->scores[i] > bestScore) {
	  index = j;
	  bestScore = res[j]->scores[i];
	}
    }

    // print the results
    cout << res[0]->TCRID[i] << " " << res[index]->epitope << " " <<
      bestScore << " " << res[index]->refTCRID[i] << endl;
  }
}

//////////////////////////////////////////////////////////////////////
// MAIN PROGRAM                                                     //
//////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

  // make sure all parameters are there
  checkCommandLineArgsClassify(argv, argc);

  // get the parameters
  Parameters p(argv, argc);

  // create a modified BLOSUM matrix
  int **newBLOSUM;
  newBLOSUM = (int **) malloc(sizeof(int *) * BLOSUM_SIZE);
  for (unsigned int i = 0; i < BLOSUM_SIZE; i++) {
    newBLOSUM[i] = (int *) malloc(sizeof(int) * BLOSUM_SIZE);
  }
  modifyBLOSUM(newBLOSUM);
  
  // read the input data
  struct TCRData *inputData;
  inputData = readTCRData(p.inputTCRsFileName);

  // read the reference and weights data
  struct RefWeightsData *refWeightsData;
  refWeightsData = readRefWeightsData(p.weightsFileName);

  // read the references
  struct TCRData **refData = (struct TCRData **) malloc(sizeof(struct TCRData *) * refWeightsData->numRefs);
  for (unsigned int i = 0; i < refWeightsData->numRefs; i++) {
    refData[i] = readTCRData(refWeightsData->fileNames[i]);
  }

  // perform nearest neighbor for each reference set
  struct Results **res = (struct Results **) malloc(sizeof(struct Results *) * refWeightsData->numRefs);
  for (unsigned int i = 0; i < refWeightsData->numRefs; i++) {
    double weights[NUMCDR];
    for (unsigned j = 0; j < NUMCDR; j++) {
      weights[j] = refWeightsData->weights[j][i];
    }
    res[i] = nearestNeighbor(refData[i], inputData, weights, newBLOSUM);
    res[i]->epitope = refWeightsData->epitopes[i];
  }

  // print the results
  printBestMatch(res, refWeightsData->numRefs);
  
  // free the memory allocated for the modified BLOSUM
  for (unsigned int i = 0; i < BLOSUM_SIZE; i++) {
    free(newBLOSUM[i]);
  }
  free(newBLOSUM);

  // free the memory for input data and reference data
  freeTCRData(inputData);
  for (unsigned int i = 0; i < refWeightsData->numRefs; i++) {
    free(refData[i]);
  }
  free(refData);

  // free the memory for the results array
  for (unsigned int i = 0; i < refWeightsData->numRefs; i++) {
    free(res[i]);
  }
  free(res);
  
  return 0;
}
