// This file is part of the swarmTCR program
// Copyright (c) 2018 Dario Ghersi and Ryan Ehrlich

using namespace std;

#ifndef _util_
#define _util_

double averagePrecision(struct Results *, char *);
void calculatePrecRecall(struct Results *, double, struct PrecRec *);
void checkCommandLineArgs(char **, int);
bool cmdOptionExists(char **, char **, const string &);
void convertAA2Num(char *, unsigned int *);
void freeSingleCellData(struct SingleCell *);
char *getCmdOption(char **, char **, const string & );
int globalAlign(unsigned int *, unsigned int, unsigned int *,
                unsigned int, int, int **);
int maximum(int, int, int);
void modifyBLOSUM(int **);
struct Results * nearestNeighbor(struct SingleCell *, struct SingleCell *, double *, int **);
struct SingleCell * readSingleCellData(string);

#endif
