#ifndef INDEX_V2_H
#define INDEX_V2_H

#include "Matrix.h"

int qg2index(int q,int g,int Q,int G);


int sq2index(int s,int q,const int *S,int Q);


int sqg2index(int s,int q,int g,const int *S,int Q,int G);


int qq2index(int q1,int q2,int Q);


void makeSigma(int g,int G,std::vector<std::vector<double> > &Sigma,const int Q,
	       const double gamma2,const double *tau2,
	       const double *a,const double *sigma2g,
	       const double *r);


void makeSigma(int g,int G,std::vector<std::vector<double> > &Sigma,
	       const std::vector<int> &on,const int Q,
	       const double gamma2,const double *tau2,
	       const double *a,const double *sigma2g,
	       const double *r);


double nuGibbs(double *nu,int Q,int G,const int *S,double gamma2,
	       const double *tau2Rho,const double *a,const double *rho,
	       const double *sigma2,const double *phi,
	       const int *psi,const double *x,
	       const int *delta,const double *Delta,Random &ran,int draw);  


double DeltaGibbs(int g,double *Delta,int Q,int G,const int *S,double c2,
		  const double *tau2R,const double *b,const double *r,
		  const double *sigma2,const double *phi,
		  const int *psi,const double *x,
		  const int *delta,const double *nu,Random &ran,int draw);  


double DeltaGibbs(double *Delta,int Q,int G,const int *S,double c2,
		  const double *tau2R,const double *b,const double *r,
		  const double *sigma2,const double *phi,
		  const int *psi,const double *x,
		  const int *delta,const double *nu,Random &ran,int draw);

#endif
