// This file is part of the swarmTCR program
// Copyright (c) 2018 Dario Ghersi and Ryan Ehrlich
// Version: 20181222

#include <cstdlib>
#include <iostream>
#include "swarmTCR.h"
#include "pso.h"
#include "util.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES                                                 //
//////////////////////////////////////////////////////////////////////

unsigned int CURRCDR = 8; // actual number of cdr regions (8 for
                          // single cell, 4 for deep sequencing)


//////////////////////////////////////////////////////////////////////
// CONSTRUCTORS                                                     //
//////////////////////////////////////////////////////////////////////

Parameters::Parameters(char **argv, int argc)
{
  // parse the command-line arguments

  // reference set file name
  refSetFileName = getCmdOption(argv, argv + argc, "-r");

  // training set file name
  trainingSetFileName = getCmdOption(argv, argv + argc, "-i");

  // test set file name
  testSetFileName = getCmdOption(argv, argv + argc, "-t");

  // test reference file name
  testRefFileName = getCmdOption(argv, argv + argc, "-x");
  if (testRefFileName) {
    useTestRef = true;
  }
  else {
    useTestRef = false;
  }

  // output PR file standard weights
  outputStdFileName = getCmdOption(argv, argv + argc, "-1");

  // output PR file optimized weights
  outputOptFileName = getCmdOption(argv, argv + argc, "-2");

  // number of iterations
  numIterations = atoi(getCmdOption(argv, argv + argc, "-n"));

  // swarm size
  swarmSize = atoi(getCmdOption(argv, argv + argc, "-s"));
}

//////////////////////////////////////////////////////////////////////
// MAIN PROGRAM                                                     //
//////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {

  // make sure all parameters are there
  checkCommandLineArgs(argv, argc);

  // get the parameters
  Parameters p(argv, argc);

  // create a modified BLOSUM matrix
  int **newBLOSUM;
  newBLOSUM = (int **) malloc(sizeof(int *) * BLOSUM_SIZE);
  for (unsigned int i = 0; i < BLOSUM_SIZE; i++) {
    newBLOSUM[i] = (int *) malloc(sizeof(int) * BLOSUM_SIZE);
  }
  modifyBLOSUM(newBLOSUM);
  
  // read the single cell data
  struct SingleCell *refSetData, *trainingSetData, *testSetData,
                    *testRefData;
  refSetData = readSingleCellData(p.refSetFileName);
  trainingSetData = readSingleCellData(p.trainingSetFileName);
  testSetData = readSingleCellData(p.testSetFileName);
  if (p.useTestRef) {
    testRefData = readSingleCellData(p.testRefFileName);
  }
  else {
    testRefData = refSetData;
  }

  // run PSO
  double *optimalWeights;
  optimalWeights = runPSO(refSetData, trainingSetData, newBLOSUM,\
                          p.swarmSize, p.numIterations);

  // test with nearest-neighbor on the training set
  double weights[] = {1.0, 1.0, 1.0, 3.0, 1.0, 1.0, 1.0, 3.0};
  struct Results *res;

  // optional test reference data
  if (p.useTestRef) {
    res = nearestNeighbor(testRefData, testSetData, weights, newBLOSUM);
  }
  else {
    res = nearestNeighbor(refSetData, testSetData, weights, newBLOSUM);
  }

  // calculate average precision with the original weights
  double ap = averagePrecision(res, p.outputStdFileName);
  cout << "\n---------------------------------------\n";
  cout << "Test results with original weights: " << ap << endl;

  // calculate average precision with the optimized weights
  res =  nearestNeighbor(refSetData, testSetData, optimalWeights, newBLOSUM);
  ap = averagePrecision(res, p.outputOptFileName);
  cout << "\nTest results with optimized weights: " << ap << endl;

  // free memory for optimal weights 
  free(optimalWeights);

  // free the memory allocated for the modified BLOSUM
  for (unsigned int i = 0; i < BLOSUM_SIZE; i++) {
    free(newBLOSUM[i]);
  }
  free(newBLOSUM);

  // free the memory for single cell data
  freeSingleCellData(trainingSetData);
  freeSingleCellData(testSetData);
  free(trainingSetData);
  free(testSetData);

  // free the memory for the results
  free(res->scores);
  free(res->labels);
  free(res);

  return 0;
}
