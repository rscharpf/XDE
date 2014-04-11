#ifndef INDEX_V2_H
#define INDEX_V2_H

#include "Matrix_v2.h"

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



void updateMRF1perfect_onedelta(int g,vector<int> &valueLower,
				vector<int> &valueUpper,
				const vector<double> &potOn,
				const vector<double> &potOff,
				const vector<vector<int> > &neighbour,
				double eta0,double omega0,double kappa,
				Random &ran);


double perfectMRF1_onedelta(int *delta,int G,
			    const vector<vector<int> > &neighbour,
			    const vector<double> &potOn,
			    const vector<double> &potOff,
			    double eta0,double omega0,double kappa,
			    unsigned int *seed,int draw);


void updateMRF2perfect_onedelta(int g,vector<int> &valueLower,
				vector<int> &valueUpper,
				const vector<double> &potOn,
				const vector<double> &potOff,
				const vector<vector<int> > &neighbour,
				double alpha,double beta,Random &ran);



double perfectMRF2_onedelta(int *delta,int G,
			    const vector<vector<int> > &neighbour,
			    const vector<double> &potOn,
			    const vector<double> &potOff,
			    double alpha,double beta,
			    unsigned int *seed,int draw);


void updateMRF2perfect(int q,int g,int Q,int G,vector<int> &valueLower,
		       vector<int> &valueUpper,const vector<double> &potOn,
		       const vector<double> &potOff,
		       const vector<vector<int> > &neighbour,
		       double alpha,double beta,double betag,
		       Random &ran);


double perfectMRF2(int *delta,int Q,int G,
		   const vector<vector<int> > &neighbour,
		   const vector<double> &potOn,
		   const vector<double> &potOff,
		   double alpha,double beta,
		   double betag,unsigned int *seed,int draw);


double OmegaGibbs(double df,const vector<vector<vector<double> > > &D,
		  const vector<int> &oldClique,
		  const vector<vector<int> > &oldComponents,
		  int Q,int G,const double *Delta,
		  const double *r,const double *sigma2,
		  const double *tau2R,const double *b,
		  vector<vector<vector<double> > > &Omega,
		  Random &ran,int draw);

double DeltaStarGibbs(const vector<int> &oldClique,
		      const vector<vector<int> > &oldComponents,
		      int Q, int G, const int *S,
		      double *Delta,const double *r,
		      const double *sigma2,const double *phi,
		      const double *tau2R,const double *b,
		      const double *nu,const int *delta,
		      const int *psi,const double *x,
		      const vector<vector<vector<double> > > &Omega,
		      Random &ran,int draw);

void transformGraph(const int *nClique,const int *oldClique,const int *nOldComponents,
		    const int *oldComponents,vector<int> &oldCliqueTransformed,
		    vector<vector<int> > &oldComponentsTransformed);

void transformOmega(const int *nClique,const int *nOldComponents,
		    const int *nNewComponents,const double *Omega,
		    vector<vector<vector<double> > > &OmegaTransformed);

void inverseTransformOmega(const vector<vector<vector<double> > > &OmegaTransformed,
			  double *Omega);


#endif
