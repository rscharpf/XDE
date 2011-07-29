#ifndef UPDATE_V2_H
#define Update_V2_H

#include <vector>

#include "Random.h"
#include "Potential_v2.h"
#include "Utility_v2.h"


/*
// order of variables:
   unsigned int *seed,
   int nTry,
   int *nAccept,
   double epsilon,
   the variables to be changed,
   int Q,
   int G,
   const int *S
   const double *x,
   const int *psi,
   const double *nu,
   const int *delta,
   const double *Delta,
   double c2,
   double gamma2,
   const double *r,
   const double *rho,
   const double *sigma2,
   const double *phi,
   const double *t,
   const double *l,
   const double *theta,
   const double *lambda,
   const double *tau2R,
   const double *tau2Rho,
   const double *xi,
   const double *a,
   const double *b,
   double pA0,
   double pA1,
   double alphaA,
   double betaA,
   double pB0,
   double pB1,
   double alphaB,
   double betaB
*/

inline void updateA(unsigned int *seed,
		    int nTry,
		    int *nAccept,
		    double epsilon,
		    double *a,
		    int Q,
		    int G,
		    const double *nu,
		    double gamma2,
		    const double *rho,
		    const double *sigma2,
		    const double *tau2Rho,
		    double pA0,
		    double pA1,
		    double alphaA,
		    double betaA) {
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
    
    if (ran.Unif01() <= exp(- pot)) {
      a[q] = newValue;
      (*nAccept)++;
    }
  }
  
  *seed = ran.ChangeSeed(*seed);

  return;
}
		       





inline void updateB(unsigned int *seed,
		    int nTry,
		    int *nAccept,
		    double epsilon,
		    double *b,
		    int Q,
		    int G,
		    const int *delta,
		    const double *Delta,
		    double c2,
		    const double *r,
		    const double *sigma2,
		    const double *tau2R,
		    double pB0,
		    double pB1,
		    double alphaB,
		    double betaB) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    int q = (int) (ran.Unif01() * Q);
    
    double oldValue = b[q];
    double p0 = 0.0;
    double p1 = 0.0;
    if (oldValue > 0.0 && oldValue < 1.0)
      {
	if (pB0 > 0.0 && oldValue - epsilon < 0.0) 
	  p0 = (epsilon - oldValue) / (2.0 * epsilon);
	if (pB1 > 0.0 && oldValue + epsilon > 1.0) 
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
	if (pB0 > 0.0 && newValue - epsilon < 0.0) 
	  p0Back = (epsilon - newValue) / (2.0 * epsilon);
	if (pB1 > 0.0 && newValue + epsilon > 1.0) 
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
    
    pot -= potentialB(Q,b,pB0,pB1,alphaB,betaB);
    pot -= potentialDDelta(Q,G,delta,Delta,c2,b,r,tau2R,sigma2);

    b[q] = newValue;
    
    pot += potentialB(Q,b,pB0,pB1,alphaB,betaB);
    pot += potentialDDelta(Q,G,delta,Delta,c2,b,r,tau2R,sigma2);
    
    b[q] = oldValue;
    
    if (ran.Unif01() <= exp(- pot)) {
      b[q] = newValue;
      (*nAccept)++;
    }
  }
  
  *seed = ran.ChangeSeed(*seed);

  return;
}
		       







inline void updateTau2RhoNu(unsigned int *seed,
			    int nTry,
			    int *nAccept,
			    double epsilon,
			    double *tau2Rho,
			    double *nu,
			    int Q,
			    int G,
			    const int *S,
			    const double *x,
			    const int *psi,
			    const int *delta,
			    const double *Delta,
			    double gamma2,
			    const double *rho,
			    const double *sigma2,
			    const double *phi,
			    const double *a) {
  Random ran(*seed);
  
  if (Q > 1) {
    int k;
    for (k = 0; k < nTry; k++) {
      int q = (int) (Q * ran.Unif01());
      int p = (int) ((Q - 1) * ran.Unif01());
      if (p >= q) p++;
      
      double upper = 1.0 + epsilon;
      double lower = 1.0 / upper;
      
      double u = lower + (upper - lower) * ran.Unif01();
      double *oldValues = (double *) calloc(Q,sizeof(double));
      double *newValues = (double *) calloc(Q,sizeof(double));
      
      int i;
      for (i = 0; i < Q; i++)
	{
	  oldValues[i] = tau2Rho[i];
	  newValues[i] = tau2Rho[i];
	}
      
      newValues[q] *= u;
      newValues[p] /= u;
      
      double prod = 1.0;
      for (i = 0; i < Q; i++)
	prod *= newValues[i];
      
      prod = exp(log(prod) / Q);
      for (i = 0; i < Q; i++)
	newValues[i] /= prod;
      
      double pot = - log(1.0 / (u * u));

      //
      // propose new values for nu from full conditionals
      //
      
      double *newNu = (double *) calloc(Q * G,sizeof(double));
      double *newNuOldTau2 = (double *) calloc(Q * G,sizeof(double));

      pot += nuGibbs(newNu,Q,G,S,gamma2,newValues,a,rho,sigma2,
		     phi,psi,x,delta,Delta,ran);
      pot -= nuGibbs(newNuOldTau2,Q,G,S,gamma2,oldValues,a,rho,sigma2,
		     phi,psi,x,delta,Delta,ran);
      
      pot -= potentialTau2Rho();
      pot += potentialTau2Rho();


      if (ran.Unif01() <= exp(- pot)) {
	int q;
	for (q = 0; q < Q; q++)
	  tau2Rho[q] = newValues[q];
	int k;
	for (k = 0; k < Q * G; k++)
	  nu[k] = newNu[k];
	
	(*nAccept)++;
      }
      else {
	int k;
	for (k = 0; k < Q * G; k++)
	  nu[k] = newNuOldTau2[k];
      }

      free(newNu);
      free(newNuOldTau2);
      free(oldValues);
      free(newValues);
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}



  
inline void updateTau2RDDelta(unsigned int *seed,
			      int nTry,
			      int *nAccept,
			      double epsilon,
			      double *tau2R,
			      double *Delta,
			      int Q,
			      int G,
			      const int *S,
			      const double *x,
			      const int *psi,
			      const double *nu,
			      const int *delta,
			      double c2,
			      const double *r,
			      const double *sigma2,
			      const double *phi,
			      const double *b) {
  Random ran(*seed);

  if (Q > 1) {
    int k;
    for (k = 0; k < nTry; k++) {
      int q = (int) (Q * ran.Unif01());
      int p = (int) ((Q - 1) * ran.Unif01());
      if (p >= q) p++;
      
      double upper = 1.0 + epsilon;
      double lower = 1.0 / upper;
      
      double u = lower + (upper - lower) * ran.Unif01();
      double *oldValues = (double *) calloc(Q,sizeof(double));
      double *newValues = (double *) calloc(Q,sizeof(double));
      
      int i;
      for (i = 0; i < Q; i++)
	{
	  oldValues[i] = tau2R[i];
	  newValues[i] = tau2R[i];
	}
      
      newValues[q] *= u;
      newValues[p] /= u;
      
      double prod = 1.0;
      for (i = 0; i < Q; i++)
	prod *= newValues[i];
      
      prod = exp(log(prod) / Q);
      for (i = 0; i < Q; i++)
	newValues[i] /= prod;
      
      double pot = - log(1.0 / (u * u));

      //
      // propose new values for Delta from full conditionals
      //
      
      double *newDDelta = (double *) calloc(Q * G,sizeof(double));
      double *newDDeltaOldTau2 = (double *) calloc(Q * G,sizeof(double));

      pot += DeltaGibbs(newDDelta,Q,G,S,c2,newValues,b,r,sigma2,
			phi,psi,x,delta,nu,ran);
      pot -= DeltaGibbs(newDDeltaOldTau2,Q,G,S,c2,oldValues,b,r,sigma2,
			phi,psi,x,delta,nu,ran);
      
      pot -= potentialTau2R();
      pot += potentialTau2R();

      if (ran.Unif01() <= exp(- pot)) {
	int q;
	for (q = 0; q < Q; q++)
	  tau2R[q] = newValues[q];

	int k;
	for (k = 0; k < Q * G; k++)
	  Delta[k] = newDDelta[k];
	
	(*nAccept)++;
      }
      else { 
	int k;
	for (k = 0; k < Q * G; k++)
	  Delta[k] = newDDeltaOldTau2[k];
      }

      free(newDDelta);
      free(newDDeltaOldTau2);
      free(oldValues);
      free(newValues);
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}



inline void updateNu(unsigned int *seed,
		     int *nAccept,
		     double *nu,
		     int Q,
		     int G,
		     const int *S,
		     const double *x,
		     const int *psi,
		     const int *delta,
		     const double *Delta,
		     double gamma2,
		     const double *rho,
		     const double *sigma2,
		     const double *phi,
		     const double *tau2Rho,
		     const double *a) {
  Random ran(*seed);
  
  //
  // propose new values for nu from full conditionals
  //
  
  nuGibbs(nu,Q,G,S,gamma2,tau2Rho,a,rho,sigma2,
	  phi,psi,x,delta,Delta,ran);
  (*nAccept)++;

  *seed = ran.ChangeSeed(*seed);

  return;
}




  
inline void updateDelta(unsigned int *seed,
			int *nAccept,
			double *Delta,
			int Q,
			int G,
			const int *S,
			const double *x,
			const int *psi,
			const double *nu,
			const int *delta,
			double c2,
			const double *r,
			const double *sigma2,
			const double *phi,
			const double *tau2R,
			const double *b) {
  Random ran(*seed);

  DeltaGibbs(Delta,Q,G,S,c2,tau2R,b,r,sigma2,phi,psi,x,delta,nu,ran);
  (*nAccept)++;

  *seed = ran.ChangeSeed(*seed);

  return;
}




#endif
