#ifndef POTENTIAL_V2_H
#define POTENTIAL_V2_H

#include <vector>

using namespace std;

#include "Random.h"
#include "Matrix_v2.h"


double potentialXqg(int q,
		    int g,
		    int Q,
		    int G,
		    const int *S,
		    const double *x,
		    const int *psi,
		    const double *nu,
		    const int *delta,
		    const double *Delta,
		    const double *sigma2,
		    const double *phi);


double potentialXg(int g, 
		   int Q,
		   int G,
		   const int *S,
		   const double *x,
		   const int *psi,
		   const double *nu,
		   const int *delta,
		   const double *Delta,
		   const double *sigma2,
		   const double *phi);



double potentialX(int Q,
		  int G,
		  const int *S,
		  const double *x,
		  const int *psi,
		  const double *nu,
		  const int *delta,
		  const double *Delta,
		  const double *sigma2,
		  const double *phi);


double potentialNug(int g,
		    int Q,
		    int G,
		    const double *nu,
		    double gamma2,
		    const double *a,
		    const double *rho,
		    const double *tau2Rho,
		    const double *sigma2);


double potentialNu(int Q,
		   int G,
		   const double *nu,
		   double gamma2,
		   const double *a,
		   const double *rho,
		   const double *tau2Rho,
		   const double *sigma2);


double potentialDDeltag(int g,
			int Q,
			int G,
			const std::vector<int> &on,
			const double *Deltag,
			double c2,
			const double *b,
			const double *r,
			const double *tau2R,
			const double *sigma2g);


double potentialDDeltag(int g,
			int Q,
			int G,
			const int *delta,
			const double *Delta,
			double c2,
			const double *b,
			const double *r,
			const double *tau2R,
			const double *sigma2);


double potentialDDelta(int Q,
		       int G,
		       const int *delta,
		       const double *Delta,
		       double c2,
		       const double *b,
		       const double *r,
		       const double *tau2R,
		       const double *sigma2);


double potentialDDeltaStar_HyperInverseWishart(const double *Delta,
					       const double *b,
					       const double *sigma2,
					       const double *tau2R,
					       const double *r,
					       int Q,int G,
					       const vector<vector<vector<double> > > &Omega,
					       const vector<int> &oldClique,
					       const vector<vector<int> > &oldComponents);

double potentialDDeltaStar_HyperInverseWishart(int gene,
					       const double *Delta,
					       const double *b,
					       const double *sigma2,
					       const double *tau2R,
					       const double *r,
					       int Q,int G,
					       const vector<vector<vector<double> > > &Omega,
					       const vector<int> &oldClique,
					       const vector<vector<int> > &oldComponents);



double potentialOmega_HyperInverseWishart(const vector<vector<vector<double> > > &Omega,
					  const vector<vector<vector<double> > > &D,
					  double df,
					  const vector<int> &oldClique,
					  const vector<vector<int> > &oldComponents);


double potentialA(int Q,
		  const double *a,
		  double pA0,
		  double pA1,
		  double alphaA,
		  double betaA);


double potentialB(int Q,
		  const double *b,
		  double pB0,
		  double pB1,
		  double alphaB,
		  double betaB);


double potentialR(int Q,
		  const double *r,
		  double nuR);


double potentialRho(int Q,
		    const double *rho,
		    double nuRho);


double potentialSigma2qg(int q,
			 int g,
			 int Q,
			 int G,
			 const double *sigma2,
			 const double *l,
			 const double *t);


double potentialSigma2(int Q,
		       int G,
		       const double *sigma2,
		       const double *l,
		       const double *t);


double potentialPhiqg(int q,
		      int g,
		      int Q,
		      int G,
		      const double *phi,
		      const double *lambda,
		      const double *theta);


double potentialPhi(int Q,
		    int G,
		    const double *phi,
		    const double *lambda,
		    const double *theta);



double potentialDeltag(int g,
		       int Q,
		       int G,
		       const int *deltag,
		       const double *xi);


double potentialDeltag_onedelta(int g,
				int Q,
				int G,
				const int *delta,
				const double *xi);


double potentialDelta(int Q,
		      int G,
		      const int *delta,
		      const double *xi);


double potentialDelta_onedelta(int Q,
			       int G,
			       const int *delta,
			       const double *xi);


double potentialXi(int Q,
		   const double *xi,
		   double alphaXi,
		   double betaXi);


double potentialXi_onedelta(const double *xi,
			    double alphaXi,
			    double betaXi);


double potentialC2(void);


double potentialGamma2(void);


double potentialTau2Rho(void);


double potentialTau2R(void);


double potentialT(void);


double potentialL(void);


double potentialTheta(void);


double potentialLambda(void);


double potentialDelta_MRF1_onedelta(int Q,
				    int G,
				    const int *delta,
				    const vector<vector<int> > &neighbour,
				    double eta0,
				    double omega0,
				    double kappa);
 

double potentialDelta_MRF2_onedelta(int Q,
				    int G,
				    const int *delta,
				    const vector<vector<int> > &neighbour,
				    double alpha,
				    double beta);


double potentialDelta_MRF2(int Q,
			   int G,
			   const int *delta,
			   const vector<vector<int> > &neighbour,
			   double alpha,
			   double beta,
			   double betag);


double potentialEta0(double eta0,double alphaEta,double betaEta);


double potentialOmega0(double omega0,double pOmega0,double lambdaOmega);


double potentialKappa(double kappa,double lambdaKappa);


double potentialAlpha(void);


double potentialBeta(void);


double potentialBetag(void);

#endif
