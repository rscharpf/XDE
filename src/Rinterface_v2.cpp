#include <R.h>

#include "Update_v2.h"

extern "C" {
  
  void RupdateA(int *seed,
		double *epsilon,
		int *nTry,
		int *nAccept,
		double *a,
		int *Q,
		int *G,
		double *pA0,
		double *pA1,
		double *alphaA,
		double *betaA,
		double *nu,
		double *gamma2,
		double *rho,
		double *tau2Rho,
		double *sigma2) {
    unsigned int seedU = (unsigned int) *seed;

    updateA(&seedU,*epsilon,*nTry,nAccept,a,*Q,*G,*pA0,*pA1,*alphaA,
	    *betaA,nu,*gamma2,rho,tau2Rho,sigma2);

    *seed = (int) seedU;

    return;
  }


} // extern "C"
    
