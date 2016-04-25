
#ifndef BICGSTAB_H_DEFINED
#define BICGSTAB_H_DEFINED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "tools.h"

#define DEFAULTACC pow(10,-7)	//	for the Biconjugate gradient method		//
#define INTERACTIVE_ACCURACY 0

//	Solves complex Linear system Kx = b, K is Kernel, x is density, size is the dimension of the square matrix Kernel	//
int CompBiCG(int size, double *density, double *densityimg, double **Kernel, double **Kernelimg, double *b, double *bimg);
//	Solves complex linear system in parallel using opnemp	//
int omp_CompBiCG(int size, double *density, double *densityimg, double **Kernel, double **Kernelimg, \
		double *b, double *bimg, int chunk);


#endif
