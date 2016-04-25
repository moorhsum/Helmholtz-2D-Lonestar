
#include "functions.h"

double cal_partial(int mode, double xi, double yi, double xj, double yj)
{
	//	Center gradient of the fundamental solution on the second variable	//
	switch (mode)
	{
		case 1:
			return (phi(xi, yi, xj+H, yj)-phi(xi, yi, xj-H, yj))/(2.0*H);
			break;
		case 2:
			return (phi(xi, yi, xj, yj+H)-phi(xi, yi, xj, yj-H))/(2.0*H);
			break;
		case 3:
			return (phiimg(xi, yi, xj+H, yj) - phiimg(xi, yi, xj-H, yj))/(2.0*H);
			break;
		case 4:
			return (phiimg(xi, yi, xj, yj+H) - phiimg(xi, yi, xj, yj-H))/(2.0*H);
			break;
		default:
			printf("Error in partial mode.\n");
			exit(0);
	}
}

double gradientdx(double x, double y)
{	return (distance( x+H,y) - distance(x-H,y))/(2.0*H);	}
double gradientdy(double x, double y)
{	return (distance(x, y+H)- distance(x, y-H))/(2.0*H);	}
double laplace(double x, double y)
{
	double Dplusx, Dminusx, Dplusy, Dminusy;
	Dplusx = (distance(x+H, y) - distance(x,y))/H;
	Dminusx = (distance(x,y) - distance(x-H,y))/H;
	Dplusy = (distance(x, y+H) - distance(x,y))/H;
	Dminusy = (distance(x,y) - distance(x,y-H))/H;
	return ( Dplusx-Dminusx + Dplusy-Dminusy )/H;
}

double cal_delta(double epsilon, double x)
{
	if (strcmp(AVG_KERNEL, "hat") == 0)
	{
		if (fabs(x) <= epsilon)
			return 1./epsilon - fabs(x)/(epsilon*epsilon);
		else
			return 0.0;
	}
	else if (strcmp(AVG_KERNEL, "cosine") == 0)
	{
	
		if (fabs(x) <= epsilon)
			return ( 1/(2.0*epsilon)) * (1.0 + cos(PI*x/epsilon));
		else
			return 0.0;
	}
	printf("The kernel is neither cosine nor hat. Weird\n");
	exit(0);
}

double phi(double xi, double yi, double xj, double yj)
{
	double norm;

	norm = sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) );
	if (norm == 0.0)
	{
		printf("At (%lf,%lf) with (%lf, %lf), norm = 0.\n", xi, yi, xj, yj);
		return 0.0;
	}

	if (LAPLACE==1)
		return log(norm)/(2.0*PI);
	else
		return 0.25*BesselY0(WAVE*norm);
}
double phiimg(double xi, double yi, double xj, double yj)
{
	double norm;
	norm = sqrt( (xi-xj)*(xi-xj) + (yi-yj)*(yi-yj) );

	if (norm == 0.0)
	{
		printf("At (%lf,%lf) with (%lf, %lf), norm = 0.\n", xi, yi, xj, yj);
		return 0.0;
	}

	if (LAPLACE==1)
		return log(norm)/(2.0*PI);
	else
		return -0.25*BesselJ0(WAVE*norm);
}

double psi(double xj, double yj)
{
	double norm;
	static double eta = PSI_ETA;
	static double neta = (-1.0) * PSI_ETA * PSI_ETA;
	
	norm = fabs( pow(xj-Wavex, 2.0) + pow(yj-Wavey, 2.0) );

	if (norm > 0.01)
		return 0.0;
	else
		return exp(norm/neta)/(eta*eta);	//	e^(-(|x|/eta)^2) / eta		//
}

double psiimg(double xj, double yj)
{
	static int flag = 0;
	if (flag ==0)
	{
		printf("psiimg has not been written yet.\n");
		flag++;
	}
	return 0.0;
}


double cal_nonhom(double zstarx, double zstary, double A[])
{
	static double Waveh = (WAVEXINT_U - WAVEXINT_L)/WaveN;
	double wx, wy;
	double realphi, imgphi, realpsi, imgpsi;
	double sum = 0.0;
	double sumimg = 0.0;
	int i, j;

	for (i = 0; i < WaveN; i++)
	{
		wy = WAVEYINT_L + Waveh * (double) (i);
		for (j = 0; j < WaveN; j++)
		{
			wx = WAVEXINT_L + Waveh * (double)(j);
			if ( distance(wx, wy) >= 0)
			{
				sum += 0.0;
			}
			else
			{
				realphi = phi(zstarx, zstary, wx, wy);
				imgphi = phiimg(zstarx, zstary, wx, wy);
				realpsi = psi(wx, wy);
				imgpsi = psiimg(wx, wy);
				sum += ( realphi*realpsi - imgphi*imgpsi);
				sumimg += (imgphi*realpsi + realphi*imgpsi);
				if ( (sum != sum) || (sum+1.0==sum) )
				{
					printf("Nonhomogeneous part at (%lf,%lf) with (%lf,%lf) is %lf.\n", wx, wy, zstarx, zstary, sum);
					exit(0);
				}
			}
		}
	}
	sum = sum * Waveh * Waveh;
	sumimg = sumimg * Waveh * Waveh;
	A[0] = sum;
	A[1] = sumimg;
	return sum;
}



double distance(double x, double y)
{
	static int ini = 0;
	static double **distance;
	static double **real;
	static double *coordinate;
	static double distance_h;
	int i, j;
	int indx, indy;
	double templl, templh, temphl, temphh, templ, temph;
	double zx, zy;
	double temp;
	double d1, d2, d3;
	static int matrixsize;
	char distfilename[50];


	FILE *fpdist;
	
	if (DIST_MODE==0)	//	Test on special case of a circle, with exact distance values	//
	{
	
		d1 = radius - sqrt(fabs( pow(x-Cx,2.0) + pow(y-Cy,2.0) ));	
		return d1;
		d2 = radius2 - sqrt(fabs( pow(x-Cx2,2.0) + pow(y-Cy2,2.0) ));
		return min(d1, -1.0*d2);
		d3 = radius3 - sqrt(fabs( pow(x-Cx3,2.0) + pow(y-Cy3,2.0) ));
		return max(max(d1, d2), d3);
	}
	else		//	Use distance file. The distance file is from redistancing a level set somewhere else	//
	{
		if (ini==0)
		{
			sprintf(distfilename, "distance%d.txt", N);
			if ((fpdist = fopen(distfilename, "r"))==NULL)
			{
				printf("Open %s error.\n", distfilename);
				exit(0);
			}
			matrixsize = SIZE;
			distance_h = (INT_U-INT_L)/matrixsize;
			real = (double **) malloc((int)(SIZE) * sizeof(double *));
			distance = (double **) malloc(matrixsize * sizeof(double *));
			coordinate = (double *) malloc(matrixsize * sizeof(double));
			
			for (i = 0; i < SIZE; i++)
			{
				real[i] = (double *) malloc((int)(SIZE) * sizeof(double));
				if (real[i]==NULL)
				{
					printf("NULL pointer at real[%d].\n", i);
					exit(0);
				}
			}

			for (i = 0 ; i < matrixsize; i++)
			{
				distance[i] = (double *) malloc(matrixsize * sizeof(double));
				if (distance[i]==NULL)
				{
					printf("NULL pointer at distance[%d].\n", i);
					exit(0);
				}
				coordinate[i] = INT_L + i * distance_h;
			}
			for (i = 0; i < matrixsize; i++)
			{
				for (j = 0; j < matrixsize; j++)
				{
					fscanf(fpdist, "%lf", &distance[i][j]);
				}
			}
			if (matrixsize < SIZE)	//	THIS WAS FOR INTRODUCING A FLEXIBILITY THAT WE HAVE FINER GRID THAN THE DISTANCE FILE CAN OFFER
			{
				for ( i = 0; i < SIZE; i++)
				{
					zy = INT_L + H * i;
					for (j = 0; j < SIZE; j++)
					{
						zx = INT_L + H * j;
						indx = (int) ( (zx-INT_L)*matrixsize/(INT_U-INT_L) );
						indy = (int) ( (zy-INT_L)*matrixsize/(INT_U-INT_L) );

						templl = distance[indy][indx];
						if ((indy == matrixsize-1)||(indx==matrixsize-1))
						{
							temphl = distance[indy][indx];
							templh = distance[indy][indx];
							temphh = distance[indy][indx];
						}
						if (indy < matrixsize-1)
						{
							templh = distance[indy+1][indx];
						}
						if (indx < matrixsize-1)
						{
							temphl = distance[indy][indx+1];
						}
						if ((indy < matrixsize-1)&&(indx < matrixsize-1))
						{
							temphh = distance[indy+1][indx+1];
						}

					   templ = templl + (temphl - templl) * (zx - coordinate[indx]) / distance_h;
					   temph = templh + (temphh - templh) * (zx - coordinate[indx]) / distance_h;
					   real[i][j] = templ + (temph - templ) * (zy - coordinate[indy])/ distance_h;
					}
				}
			}
			ini++;
		}
		if (matrixsize >= SIZE)
		{
			if (ini==1)
			{
				ini++;
			}
			indx = (int) ( (x-INT_L)*matrixsize/(INT_U-INT_L) );
			indy = (int) ( (y-INT_L)*matrixsize/(INT_U-INT_L) );
			if (indy < 0)
				printf("indx = %d, indy = %d.\n", indx, indy);
			return -1.0 * distance[indy][indx];
		}
		else
		{
			indx = (int) ( (x-INT_L)*N/(INT_U-INT_L) );
			indy = (int) ( (y-INT_L)*N/(INT_U-INT_L) );
			return -1.0 * real[indy][indx];
		}
	}
}

//	Extrapolation weight for evaluation on boundary (polynomial)
double PolyWeight(double t, double tau)
{
	double term1, term2, term3, tau2, tau3, tau4;
	if ( (t <= tau) || (t >= 1.0) )
		return 0.0;
	else
	{
		tau2 = tau * tau;
		tau3 = tau2 * tau;
		tau4 = tau2 * tau2;
//		return 2./(1.-tau) - 4.*fabs(t-0.5*(1.+tau))/((1.-tau)*(1.-tau));
//		return 6. * (t-1.)*(t-tau)/((tau-1.)*(tau-1.)*(tau-1.));
//		return -30.*(t-1.)*(t-1.)*(t-tau)*(t-tau)/pow(tau-1,5.0);

		tau2 = tau * tau;
		if (1)
		{
			tau3 = tau2 * tau;
			tau4 = tau2 * tau2;

			term1 = -6. * (3. + 8.*tau + 3.* tau2)*t*t;
			term2 = 4. * (5. + 16.*tau + 16. * tau2 + 5. * tau3)*t;
			term3 = -1.* (5. + 20.*tau + 34.*tau2 + 20.*tau3 + 5.*tau4);

			return 210. * (t-1.)*(t-1.)*(t-tau)*(t-tau)*(term1 + term2 + term3)/pow(tau-1.,9.0);
		}

		term1 = 7. * (1+tau)*t;
		term2 = -2. * (2. + 3. * tau + 2. * tau2);

		return 60. * (t-1.)*(t-1.)*(t-tau)*(t-tau)*(term1 + term2)/pow(tau-1., 7.0);

//		term1 = 420. * ( 1. + 3.*tau + tau2 ) * t * t;
//		term2 = -60. * ( 8. + 27.*tau + 27.*tau2 + 8.*tau3 ) * t;
//		term3 = 60. * ( 2. + 8.*tau + 15.*tau2 + 8.*tau3 + 2.*tau4 );
//		return ( t - 1. ) * ( t - tau ) * ( term1 + term2 + term3 ) / pow(tau-1., 7.0);
	}
}

//	Extrapolation weight for evaluation on boundary (sine weight)
double SineWeight(double t, double tau)
{
	double term1, term2, term3, tau2, PI2, PI3;
	double denom;

	if ( (t <= tau) || (t >= 1.0) )
		return 0.0;
	if (MOMENT_ORDER == 2)
	{
		tau2 = tau * tau;
		PI2 = PI * PI;
		PI3 = PI2 * PI;

		term1 = 2.*(20. - 40.*tau + 9.*PI2*tau + 20.*tau2) * (1. - cos(2.*PI*(t-tau)/(1.-tau))) / (3.*PI2 - 31.);
		term2 = 9.*(PI2 + PI2 * tau) * (tau-1.)/32. * (cos(PI*(t-tau)/(1.-tau)) - cos(3.*PI*(t-tau)/(1.-tau)));
		term3 = -9. * (3.*PI + PI3 - 6.*PI*tau + 4.*PI3*tau + 3.*PI*tau2 + PI3*tau2 ) * \
			(3.*sin(PI*(t-tau)/(1.-tau)) - sin(3.*PI*(t-tau)/(1.-tau))) / ( 16. * (3. * PI2 - 31.) );

		return (term1 + term2 + term3) / ( (tau-1.)*(tau-1.)*(tau-1.) );
	}
	else	//	FIRST MOMENT	//
	{
		denom = tau - 1.;
		term1 = PI * pow( cos( PI * (t-1.) /(2.*denom) ), 3.0);
		term2 = 7. + tau + 3. * (5. + 3. * tau) * cos(PI * (tau-t)/denom);
		term3 = sin( PI * (t-1.)/(2. * denom))/(denom*denom);
		return term1 * term2 * term3;
	}
}

