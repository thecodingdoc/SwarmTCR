// This file is part of the swarmTCR program
// Copyright (c) 2018 Dario Ghersi and Ryan Ehrlich

#include <cstdlib>
#include <iostream>
#include <vector>
#include "pso.h"
#include "util.h"

extern unsigned int CURRCDR;

//////////////////////////////////////////////////////////////////////
// CONSTRUCTORS                                                     //
//////////////////////////////////////////////////////////////////////

Particle::Particle(struct SingleCell *refSingleCell,\
                   struct SingleCell *trainingSingleCell,\
                   int **newBLOSUM)
{

  // assign the pointers to the reference and training set
  ref = refSingleCell;
  training = trainingSingleCell;

  // initialize the position and velocity using random numbers [0, 1]
  for (unsigned int i = 0; i < CURRCDR; i++) {
    position[i] = (double) rand() / RAND_MAX;
    velocity[i] = (double) rand() / RAND_MAX;
  }

  // assign a very low initial best score
  bestScore = -1000000.0;

  mat = newBLOSUM;
}

//////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                        //
//////////////////////////////////////////////////////////////////////

double Particle::runNearestNeighbor()
{

  Results *res;

  // compute the predictions
  res = nearestNeighbor(ref, training, position, mat);
  
  // calculate the average precision
  double ap = averagePrecision(res, 0);
  if (ap > bestScore) {
    bestScore = ap;
    for (unsigned int i = 0; i < CURRCDR; i++) {
      bestPos[i] = position[i];
    }
  }
  free(res);

  return ap;
}

//////////////////////////////////////////////////////////////////////

double * runPSO(struct SingleCell *refSingleCell,\
         struct SingleCell *trainingSingleCell, int **newBLOSUM,
         unsigned int swarmSize, unsigned int numIterations)
{
  
  srand(time(NULL)); // initialize the random number generator

  // allocate memory for the final optimal position
  double *globalBestPos;
  globalBestPos = (double *) malloc(sizeof(double) * CURRCDR);

  // create and initialize the vector of particles
  vector<Particle> particles;
  for (unsigned int i = 0; i < swarmSize; i++) {
    particles.push_back(Particle(refSingleCell, trainingSingleCell,\
                                 newBLOSUM));
   particles[i].runNearestNeighbor();
  }

  // iterate the process
  double globalBest = -10000.0;
  for (unsigned int i = 0; i < numIterations; i++) {

    // find the global best score
    for (unsigned int j = 0; j < swarmSize; j++) {
      if (particles[j].bestScore > globalBest) {
	      globalBest = particles[j].bestScore;
	
        for (unsigned int k = 0; k < CURRCDR; k++) {
	        globalBestPos[k] = particles[j].bestPos[k];
	      }
      }
    }

    // compute the velocity for each particle,
    // update the position, and compute the objective function
    #pragma omp parallel for
    for (unsigned int j = 0; j < swarmSize; j++) {
      for (unsigned int k = 0; k < CURRCDR; k++) {
        // update the velocity
        double r1 = (double) rand() / RAND_MAX;
        double r2 = (double) rand() / RAND_MAX;
        particles[j].velocity[k] = INERTIA * particles[j].velocity[k] +\
                                   C1 * r1 * (particles[j].bestPos[k] - particles[j].position[k]) +\
                                   C2 * r2 * (globalBestPos[k] - particles[j].position[k]);
        // update the position
        particles[j].position[k] = particles[j].position[k] + particles[j].velocity[k];

        // make sure the position is bounded between 0 and 1 inclusive
        if (particles[j].position[k] < 0.0) {
          particles[j].position[k] = 0.0;
        }
        else if (particles[j].position[k] > 1.0) {
          particles[j].position[k] = 1.0;
        }
      }
      particles[j].runNearestNeighbor();
    }

    // print the current state
    cout << "#" << i << " " << globalBest << " -- ";
    for (unsigned int k = 0; k < CURRCDR; k++) {
      cout << globalBestPos[k] << " ";
    }
    cout << endl;
  }

  return globalBestPos;
}
