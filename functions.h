
#ifndef FUNCTIONS_H_DEFINED
#define FUNCTIoNS_H_DEFINED

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "PDE.h"
#include "frame.h"
#include "constants.h"
#include "tools.h"

#define PSI_ETA 0.012	//	For nonhomogeneous mode, this d etermines how smooth the point source is, the smaller the more concentrated 	//
#define DIST_MODE 0 	//	0 for circle(s), 1 for level set function	//
#define NONHOM_MODE 0	//	1 for nonhomogeneous mode, 0 for homogeneous mode.	//
#define MOMENT_ORDER 1
#define AVG_KERNEL "cosine"	//	hat or cosine, this is for the discretized delta function	//

extern int N;
extern double H;
extern int SIZE;
extern double WAVE;

//	Old way of computing the kernel, new way computes the gradient analytically	//
double cal_partial(int mode, double x1, double y1, double x2, double y2);
double cal_regpartial(double x1, double y1, double x2, double y2, double tau);
double cal_delta(double epsilon, double x);
double cal_nonhom(double zstarx, double zstary, double A[]);

//	Old ways of computing the derivatives of distance functions	//
double gradientdx(double x, double y);
double gradientdy(double x, double y);
double laplace(double x, double y);
//	Fundamental solution	//
double phi(double xi, double yi, double xj, double yj);
double phiimg(double xi, double yi, double xj, double yj);
//	Point source for nonhomogeneous equation	//
double psi(double xj, double yj);
double psiimg(double xj, double yj);

double distance(double x, double y);

//	Extrapolation weight to evaluate boundary values	//
double PolyWeight(double t, double tau);
double SineWeight(double t, double tau);

#endif
