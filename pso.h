// This file is part of the swarmTCR program
// Copyright (c) 2018 Dario Ghersi and Ryan Ehrlich

#include <string>
#include "swarmTCR.h"

#ifndef _pso_
#define _pso_

#define INERTIA 0.5
#define C1 2.0
#define C2 2.0

class Particle {

 public:
  double position[NUMCDR];
  double velocity[NUMCDR];
  double bestPos[NUMCDR];
  double bestScore;
  int **mat;

  struct SingleCell *ref;
  struct SingleCell *training;

  Particle(struct SingleCell *, struct SingleCell *, int **);

  double runNearestNeighbor();
};

//////////////////////////////////////////////////////////////////////

double * runPSO(struct SingleCell *, struct SingleCell *, int **,\
                unsigned int, unsigned int);

#endif
