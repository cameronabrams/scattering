/*
 *   Structure factor for a periodic system
 *
 *   (c) 2014 cameron f abrams
 *
 *
 */

#ifndef _SKDATASPACE_H_
#define _SKDATASPACE_H_
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

typedef struct SKBIN {
  int i;
  double s;
  int c;
} SkBin;

typedef struct SKDATASPACESTRUCT {
  int n; // number of atoms
  int nframes; // number of trajectory frames
  double c[3]; // orthorectangular cell size
  double * b; // atom scattering lengths (A)
  double ** x; // coordinates; x[0][i] is x-coordinate of atom i
  
  int nk[3]; // number of k-vectors in each dimension

  double ** co; // cosine-terms:  co[i][j] = cos(2*j*M_PI/c[i] * x[i][this_atom]) 
  double ** si; // sine-terms:    si[i][j] = sin(2*j*M_PI/c[i] * x[i][this_atom])
  double *** C; // cosine-terms aggregator
  double *** S; // sine-terms aggregator
  
  FILE * dump;

  SkBin * Sk;
  double dq;
  double qmax;
  int nq;

  int frame; // current frame number
  
} SkDataSpace;

SkDataSpace * NewSkDataSpace ( int n, int nframes, int nkx, int nky, int nkz );
void FreeSkDataSpace ( SkDataSpace * sds );

double SkDataSpace_setdq ( SkDataSpace * sds, double dq );
double * SkDataSpace_getCellSizeAddr ( SkDataSpace * sds );
double * SkDataSpace_getScatLengthAddr ( SkDataSpace * sds );
double * SkDataSpace_getCoordAddr ( SkDataSpace * sds, int d );

void SkDataSpace_scaleCoords ( SkDataSpace * sds, int d );

int SkDataSpace_setFrame ( SkDataSpace * sds, int f );

void SkDataSpace_report ( SkDataSpace * sds );
void SkDataSpace_updateSk ( SkDataSpace * sds );
void SkDataSpace_outputSk ( SkDataSpace * sds );
#endif
