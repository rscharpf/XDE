#ifndef UPDATE_V2_H
#define Update_V2_H

#include <vector>

#include "Random.h"
#include "Potential_v2.h"


inline void updateA(unsigned int *seed,
		    double epsilon,
		    int nTry,
		    int *nAccept,
		    double *a,
		    int Q,
		    int G,
		    double pA0,
		    double pA1,
		    double alphaA,
		    double betaA,
		    const double *nu,
		    double gamma2,
		    const double *rho,
		    const double *tau2Rho,
		    const double *sigma2) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    int q = (int) (ran.Unif01() * Q);
    
    double oldValue = a[q];
    double p0 = 0.0;
    double p1 = 0.0;
    if (oldValue > 0.0 && oldValue < 1.0)
      {
	if (pA0 > 0.0 && oldValue - epsilon < 0.0) 
	  p0 = (epsilon - oldValue) / (2.0 * epsilon);
	if (pA1 > 0.0 && oldValue + epsilon > 1.0) 
	  p1 = (oldValue + epsilon - 1.0) / (2.0 * epsilon);
      }
    
    double newValue;
    double lower = 0.0;
    double upper = 0.0;
    double u = ran.Unif01();
    if (u < p0)
      newValue = 0.0;
    else if (u < p0 + p1)
      newValue = 1.0;
    else
      {
	lower = oldValue - epsilon;
	upper = oldValue + epsilon;
	if (lower < 0.0) lower = 0.0;
	if (upper > 1.0) upper = 1.0;
	newValue = lower + (upper - lower) * ran.Unif01();
      }
    
    
    double p0Back = 0.0;
    double p1Back = 0.0;
    if (newValue > 0.0 && newValue < 1.0)
      {
	if (pA0 > 0.0 && newValue - epsilon < 0.0) 
	  p0Back = (epsilon - newValue) / (2.0 * epsilon);
	if (pA1 > 0.0 && newValue + epsilon > 1.0) 
	  p1Back = (newValue + epsilon - 1.0) / (2.0 * epsilon);
      }
    
    double lowerBack = 0.0;
    double upperBack = 1.0;
    if (oldValue > 0.0 && oldValue < 1.0)
      {
	lowerBack = newValue - epsilon;
	upperBack = newValue + epsilon;
	if (lowerBack < 0.0) lowerBack = 0.0;
	if (upperBack > 1.0) upperBack = 1.0;
      }
    
    
    double pot = 0.0;
    if (oldValue == 0.0)  // then (newValue > 0.0 && newValue < 1.0)
      {
	pot -= - log(1.0);
	pot -= - log(1.0 / (upper - lower));
	pot += - log(p0Back);
      }
    else if (oldValue == 1.0) // then (newValue > 0.0 && newValue < 1.0)
      {
	pot -= - log(1.0);
	pot -= - log(1.0 / (upper - lower));
	pot += - log(p1Back);
      }
    else  // then (oldValue > 0.0 && oldValue < 1.0)
      {
	if (newValue == 0.0)
	  {
	    pot -= - log(p0);
	    pot += - log(1.0);
	    pot += - log(1.0 / (upperBack - lowerBack));
	  }
	else if (newValue == 1.0)
	  {
	    pot -= - log(p1);
	    pot += - log(1.0);
	    pot += - log(1.0 / (upperBack - lowerBack));
	  }
	else  // then (newValue > 0.0 && newValue < 1.0)
	  {
	    pot -= - log(1.0 - p0 - p1);
	    pot -= - log(1.0 / (upper - lower));
	    pot += - log(1.0 - p0Back - p1Back);
	    pot += - log(1.0 / (upperBack - lowerBack));
	  }
      }
    
    pot -= potentialA(Q,a,pA0,pA1,alphaA,betaA);
    pot -= potentialNu(Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);

    a[q] = newValue;
    
    pot += potentialA(Q,a,pA0,pA1,alphaA,betaA);
    pot += potentialNu(Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    
    a[q] = oldValue;
    
    if (ran.Unif01() <= exp(- pot))
      {
	a[q] = newValue;
	(*nAccept)++;
      }
  }
  
  *seed = ran.ChangeSeed(*seed);

  return;
}
		       



#endif
