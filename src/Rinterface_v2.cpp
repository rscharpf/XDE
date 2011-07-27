#include <R.h>

#include "Update_v2.h"

extern "C" {
  
  void updateA(int *seed,
	       int *nTry,
	       int *nAccept,
	       double *epsilon,
	       double *a,
	       int *Q,
	       int *G,
	       const double *nu,
	       double *gamma2,
	       const double *rho,
	       const double *sigma2,
	       const double *tau2Rho,
	       double *pA0,
	       double *pA1,
	       double *alphaA,
	       double *betaA) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateA(&seedU,*nTry,nAccept,*epsilon,a,*Q,*G,nu,*gamma2,rho,sigma2,
	    tau2Rho,*pA0,*pA1,*alphaA,*betaA);

    *seed = (int) seedU;

    return;
  }


  void updateB(int *seed,
	       int *nTry,
	       int *nAccept,
	       double *epsilon,
	       double *b,
	       int *Q,
	       int *G,
	       const int *delta,
	       const double *Delta,
	       double *c2,
	       const double *r,
	       const double *sigma2,
	       const double *tau2R,
	       double *pB0,
	       double *pB1,
	       double *alphaB,
	       double *betaB) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateB(&seedU,*nTry,nAccept,*epsilon,b,*Q,*G,delta,Delta,*c2,r,sigma2,
	    tau2R,*pB0,*pB1,*alphaB,*betaB);

    *seed = (int) seedU;

    return;
  }


  void updateTau2RhoNu(int *seed,
		       int *nTry,
		       int *nAccept,
		       double *epsilon,
		       double *tau2Rho,
		       double *nu,
		       int *Q,
		       int *G,
		       const int *S,
		       const double *x,
		       const int *psi,
		       const int *delta,
		       const double *Delta,
		       double *gamma2,
		       const double *rho,
		       const double *sigma2,
		       const double *phi,
		       const double *a) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateTau2RhoNu(&seedU,*nTry,nAccept,*epsilon,tau2Rho,nu,*Q,*G,S,x,
		    psi,delta,Delta,*gamma2,rho,sigma2,phi,a);

    *seed = (int) seedU;

    return;
  }



  void updateTau2RDDelta(int *seed,
			 int *nTry,
			 int *nAccept,
			 double *epsilon,
			 double *tau2R,
			 double *Delta,
			 int *Q,
			 int *G,
			 const int *S,
			 const double *x,
			 const int *psi,
			 const double *nu,
			 const int *delta,
			 double *c2,
			 const double *r,
			 const double *sigma2,
			 const double *phi,
			 const double *b) {
    unsigned int seedU = (unsigned int) *seed;
    
    updateTau2RDDelta(&seedU,*nTry,nAccept,*epsilon,tau2R,Delta,*Q,*G,S,x,
		      psi,nu,delta,*c2,r,sigma2,phi,b);

    *seed = (int) seedU;

    return;
  }



} // extern "C"
    
