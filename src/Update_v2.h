#ifndef UPDATE_V2_H
#define Update_V2_H

#include <vector>

#include "Random.h"
#include "Potential_v2.h"
#include "Utility_v2.h"
#include "Cholesky.h"

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
   double betaB,
   double nuR,
   double nuRho,
   double c2Max
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




inline void updateC2(unsigned int *seed,
		     int nTry,
		     int *nAccept,
		     double *c2,
		     int Q,
		     int G,
		     const int *delta,
		     const double *Delta,
		     const double *r,
		     const double *sigma2,
		     const double *tau2R,
		     const double *b,
		     double c2Max) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    (*nAccept)++;

    //
    // set prior parameters
    //
    
    double s = -1.0;
    double lambda = 0.0;
    
    //
    // update parameters based on available observations
    //
    
    int g;
    for (g = 0; g < G; g++)
      {
	int nOn = 0;
	vector<int> on(Q,0);
	int q;
	for (q = 0; q < Q; q++) {
	  int kqg = qg2index(q,g,Q,G);
	  on[q] = delta[kqg];
	  nOn += delta[kqg];
	}
	
	if (nOn >= 1)
	  {
	    vector<vector<double> > varInv;
	    makeSigma(varInv,on,Q,1.0,tau2R,b,sigma2 + g * Q,r);
	    
	    //
	    // pick out the elements of Delta that are active
	    //
	    
	    vector<double> DeltaActive(nOn,0.0);
	    int qq = 0;
	    for (q = 0; q < Q; q++)
	      {
		if (on[q] == 1)
		  {
		    int kqg = qg2index(q,g,Q,G);
		    DeltaActive[qq] = Delta[kqg];
		    qq++;
		  }
	      }
	    
	    vector<vector<double> > var;
	    inverse(varInv,var);
	    
	    double value = quadratic(var,DeltaActive);
	    lambda += 0.5 * value;
	    s += 0.5 * nOn;
	  }
      }
    
    //
    // Draw new value
    //
    
    double newValue;
    if (s > 0.0) // there exist at least one delta == 1
      {
	int nTry = 0;
	do
	  {
	    nTry++;
	    newValue = ran.InverseGamma(s,lambda);
	  }
	while (newValue > c2Max && nTry < 100);
	if (nTry == 100) {
	  newValue = *c2;   // propose an unchanged value!
	  (*nAccept)--;
	}
      }
    else if (s == 0.0)
      {
	double fmax;
	if (lambda < c2Max)
	  fmax = exp(-1.0)/lambda;
	else
	  fmax = exp(-lambda/c2Max) / c2Max;
	int nTry = 0;
	int accept = 0;
	do
	  {
	    nTry++;
	    newValue = c2Max * ran.Unif01();
	    double alpha = (exp(-lambda/newValue)/newValue) / fmax;
	    accept = (ran.Unif01() <= alpha);
	  }
	while (accept == 0 && nTry < 100);
	if (nTry == 100) {
	  newValue = *c2;   // propose an unchanged value!
	  (*nAccept)--;
	}
      }
    else if (s == -0.5)
      {
	double fmax;
	if (lambda < 0.5 * c2Max)
	  fmax = exp(-0.5)/sqrt(2*lambda);
	else
	  fmax = exp(-lambda/c2Max) / sqrt(c2Max);
	int nTry = 0;
	int accept = 0;
	do 
	  {
	    nTry++;
	    newValue = c2Max * ran.Unif01();
	    double alpha = (exp(-lambda/newValue)/sqrt(newValue)) / fmax;
	    accept = (ran.Unif01() < alpha);
	  }
	while (accept == 0 && newValue < 100);
	if (nTry == 100) {
	  newValue = *c2;   // propose an unchanged value!
	  (*nAccept)--;
	}
      }
    else
      newValue = c2Max * ran.Unif01();
    
    //
    // Set new value
    //
    
    *c2 = newValue;
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}




inline void updateGamma2(unsigned int *seed,
			 int *nAccept,
			 double *gamma2,
			 int Q,
			 int G,
			 const double *nu,
			 const double *rho,
			 const double *sigma2,
			 const double *tau2Rho,
			 const double *a) {
  Random ran(*seed);

  //
  // set prior parameters
  //

  double s = -1.0;
  double lambda = 0.0;

  //
  // update parameters based on available observations
  //

  int g;
  for (g = 0; g < G; g++)
    {
      vector<vector<double> > varInv;
      makeSigma(varInv,Q,1.0,tau2Rho,a,sigma2 + g * Q,rho);
      
      vector<vector<double> > var;
      inverse(varInv,var);

      vector<double> nuActive(Q,0.0);
      int q;
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	nuActive[q] = nu[kqg];
      }
      
      double value = quadratic(var,nuActive);
      lambda += 0.5 * value;
      s += 0.5 * Q;
    }
  
  //
  // Draw new value
  //
  
  double newValue = ran.InverseGamma(s,lambda);
  
  //
  // Set new value
  //
  
  *gamma2 = newValue;

  (*nAccept)++;

  *seed = ran.ChangeSeed(*seed);

  return;
}





inline void updateRC2(unsigned int *seed,
		      int nTry,
		      int *nAccept,
		      double epsilon,
		      double *r,
		      double *c2,
		      int Q,
		      int G,
		      const int *delta,
		      const double *Delta,
		      const double *sigma2,
		      const double *tau2R,
		      const double *b,
		      double nuR,
		      double c2Max) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    vector<vector<double> > oldR;
    vector<vector<double> > newR;
    oldR.resize(Q);
    newR.resize(Q);
    int p,q;
    for (p = 0; p < Q; p++)
      {
	oldR[p].resize(Q);
	newR[p].resize(Q);
	for (q = 0; q < Q; q++)
	  {
	    int kqq = qq2index(p,q,Q);
	    oldR[p][q] = r[kqq];
	    newR[p][q] = r[kqq];
	  }
      }
    
    //
    // propose new values for r
    //  
    
    double pot = 0.0;
    double u = epsilon * ran.Norm01();
    
    //
    // draw what element to change
    //
    
    vector<double> prob(Q);
    for (p = 0; p < Q; p++)
      prob[p] = 1.0 / ((double) Q);
    int pp = ran.Discrete(prob);
    prob.resize(Q - 1);
    for (p = 0; p < Q - 1; p++)
      prob[p] = 1.0 / ((double) (Q - 1));
    
    int qq = ran.Discrete(prob);
    qq += (qq >= pp);
    
    //
    // compute potential new correlation value
    //
    
    newR[pp][qq] = oldR[pp][qq] * exp(u) / 
      (1 - oldR[pp][qq] + oldR[pp][qq] * exp(u));
    newR[qq][pp] = newR[pp][qq];
    double *newRVector = (double *) calloc(Q * Q,sizeof(double));
    qq = 0;
    for (p = 0; p < Q; p++)
      for (q = 0; q < Q; q++) {
	newRVector[qq] = newR[p][q];
	qq++;
      }    

    //
    // check that potential new correlation matrix is positive definite
    //
    
    int err = 0;
    Cholesky chol(newR,err);
    if (err == 0) {
      
      //
      // compute Jacobian determinant
      //
      
      double y = oldR[pp][qq];
      double xtilde = log(y) - log(1.0 - y) + u;
      double pot1;
      if (xtilde <= 0.0)
	pot1 = - log(1.0 + exp(xtilde));
      else
	pot1 = - xtilde - log(1.0 + exp(- xtilde));
      double potdytildedxtilde = - xtilde - 2.0 * pot1;
      
      double potdxtildedy = log(1.0 - y) + log(y);
      
      pot += potdytildedxtilde + potdxtildedy;
      
      //
      // if any of the proposed new values are negative, reject the proposal
      //
      
      int isNeg = 0;
      for (p = 0; p < Q; p++)
	for (q = 0; q < Q; q++)
	  isNeg += (newR[p][q] < 0.0);
      
      if (isNeg == 0) {

	//
	// propose new value for c2
	//
	
	double oldC2 = *c2;
	double newC2;
	
	double s = -1.0;
	double lambda = 0.0;
	
	int g;
	for (g = 0; g < G; g++)
	  {
	    int nOn = 0;
	    vector<int> on(Q,0);
	    for (q = 0; q < Q; q++) {
	      int kqg = qg2index(q,g,Q,G);
	      on[q] = delta[kqg];
	      nOn += delta[kqg];
	    }

	    vector<vector<double> > varInv;
	    makeSigma(varInv,on,Q,1.0,tau2R,b,sigma2 + g * Q,newRVector);
	    
	    vector<vector<double> > var;
	    inverse(varInv,var);
      
	    //
	    // pick out elements of Delta that are active
	    //

	    vector<double> DeltaActive(nOn,0.0);
	    int qq = 0;
	    for (q = 0; q < Q; q++)
	      {
		if (on[q] == 1)
		  {
		    int kqg = qg2index(q,g,Q,G);
		    DeltaActive[qq] = Delta[kqg];
		    qq++;
		  }
	      }
	    
	    double value = quadratic(var,DeltaActive);
	    lambda += 0.5 * value;
	    s += 0.5 * nOn;
	  }
	
	if (s > 0.0)
	  {
	    newC2 = ran.InverseGamma(s,lambda);
	    pot -= ran.PotentialInverseGamma(s,lambda,newC2);
	  }
	else
	  {
	    newC2 = c2Max * ran.Unif01();
	    pot -= - log(1.0 / c2Max);
	  }
  
	//
	// compute potential for inverse proposal for c2
	//
	
	s = -1.0;
	lambda = 0.0;
	
	for (g = 0; g < G; g++)
	  {
	    int nOn = 0;
	    vector<int> on(Q,0);
	    for (q = 0; q < Q; q++) {
	      int kqg = qg2index(q,g,Q,G);
	      on[q] = delta[kqg];
	      nOn += delta[kqg];
	    }

	    vector<vector<double> > varInv;
	    makeSigma(varInv,on,Q,1.0,tau2R,b,sigma2 + g * Q,r);

	    vector<vector<double> > var;
	    inverse(varInv,var);
      
	    //
	    // pick out elements of Delta that are active
	    //

	    vector<double> DeltaActive(nOn,0.0);
	    int qq = 0;
	    for (q = 0; q < Q; q++)
	      {
		if (on[q] == 1)
		  {
		    int kqg = qg2index(q,g,Q,G);
		    DeltaActive[qq] = Delta[kqg];
		    qq++;
		  }
	      }

	    double value = quadratic(var,DeltaActive);
	    lambda += 0.5 * value;
	    s += 0.5 * nOn;
	  }
	
	if (s > 0.0)
	  pot += ran.PotentialInverseGamma(s,lambda,oldC2);
	else
	  pot += - log(1.0 / c2Max);
  
	//
	// compute potentials for new and old states
	//
	
	if (newC2 < c2Max) { // otherwise do not accept
	  pot -= potentialR(Q,r,nuR);
	  pot -= potentialC2();
	  pot -= potentialDDelta(Q,G,delta,Delta,*c2,b,r,tau2R,sigma2);

	  pot += potentialR(Q,newRVector,nuR);
	  pot += potentialC2();
	  pot += potentialDDelta(Q,G,delta,Delta,newC2,b,newRVector,tau2R,sigma2);
	  
	  if (ran.Unif01() <= exp(- pot)) {
	    (*nAccept)++;
	    
	    *c2 = newC2;
	    int kk;
	    for (kk = 0; kk < Q * Q; kk++)
	      r[kk] = newRVector[kk];
	  }
	}
      }
    }
    
    free(newRVector);
  }
  
  *seed = ran.ChangeSeed(*seed);
  
  return;
}




inline void updateRhoGamma2(unsigned int *seed,
			    int nTry,
			    int *nAccept,
			    double epsilon,
			    double *rho,
			    double *gamma2,
			    int Q,
			    int G,
			    const double *nu,
			    const double *sigma2,
			    const double *tau2Rho,
			    const double *a,
			    double nuRho) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    vector<vector<double> > oldRho;
    vector<vector<double> > newRho;
    oldRho.resize(Q);
    newRho.resize(Q);
    int p,q;
    for (p = 0; p < Q; p++)
      {
	oldRho[p].resize(Q);
	newRho[p].resize(Q);
	for (q = 0; q < Q; q++)
	  {
	    int kqq = qq2index(p,q,Q);
	    oldRho[p][q] = rho[kqq];
	    newRho[p][q] = rho[kqq];
	  }
      }
    
    double pot = 0.0;
    double u = epsilon * ran.Norm01();
    
    //
    // draw what element to change
    //
    
    vector<double> prob(Q);
    for (p = 0; p < Q; p++)
      prob[p] = 1.0 / ((double) Q);
    int pp = ran.Discrete(prob);
    prob.resize(Q - 1);
    for (p = 0; p < Q - 1; p++)
      prob[p] = 1.0 / ((double) (Q - 1));
    
    int qq = ran.Discrete(prob);
    qq += (qq >= pp);
    
    //
    // compute potential new correlation value
    //
    
    newRho[pp][qq] = oldRho[pp][qq] * exp(u) / 
      (1 - oldRho[pp][qq] + oldRho[pp][qq] * exp(u));
    newRho[qq][pp] = newRho[pp][qq];
    double *newRhoVector = (double *) calloc(Q * Q,sizeof(double));
    qq = 0;
    for (p = 0; p < Q; p++)
      for (q = 0; q < Q; q++) {
	newRhoVector[qq] = newRho[p][q];
	qq++;
      }    
    
    // check that potential new correlation matrix is positive definite
    //
    
    int err = 0;
    Cholesky chol(newRho,err);
    if (err == 0) {

      //
      // compute Jacobian determinant
      //
      
      double y = oldRho[pp][qq];
      double xtilde = log(y) - log(1.0 - y) + u;
      double pot1;
      if (xtilde <= 0.0)
	pot1 = - log(1.0 + exp(xtilde));
      else
	pot1 = - xtilde - log(1.0 + exp(- xtilde));
      double potdytildedxtilde = - xtilde - 2.0 * pot1;
      
      double potdxtildedy = log(1.0 - y) + log(y);
      
      pot += potdytildedxtilde + potdxtildedy;
      
      //
      // if any of the proposed new values are negative, reject the proposal
      //
      
      int isNeg = 0;
      for (p = 0; p < Q; p++)
	for (q = 0; q < Q; q++)
	  isNeg += (newRho[p][q] < 0.0);

      if (isNeg == 0) {
	
	//
	// propose new values for gamma2
	//
	
	double oldGamma2 = *gamma2;
	double newGamma2;
	
	double s = -1.0;
	double lambda = 0.0;
	int g;
	for (g = 0; g < G; g++) {
	  vector<vector<double> > varInv;
	  makeSigma(varInv,Q,1.0,tau2Rho,a,sigma2 + g * Q,newRhoVector);

	  vector<vector<double> > var;
	  inverse(varInv,var);
      
	  vector<double> nuActive(Q,0.0);
	  for (q = 0; q < Q; q++) {
	    int kqg = qg2index(q,g,Q,G);
	    nuActive[q] = nu[kqg];
	  }
	  
	  double value = quadratic(var,nuActive);
	  lambda += 0.5 * value;
	  s += 0.5 * Q;
	}
	
	newGamma2 = ran.InverseGamma(s,lambda);
	pot -= ran.PotentialInverseGamma(s,lambda,newGamma2);
	
	//
	// compute potential for inverse proposal for gamma2
	//
	
	s = -1.0;
	lambda = 0.0;
	for (g = 0; g < G; g++) {
	  vector<vector<double> > varInv;
	  makeSigma(varInv,Q,1.0,tau2Rho,a,sigma2 + g * Q,rho);
	  
	  vector<vector<double> > var;
	  inverse(varInv,var);
      
	  vector<double> nuActive(Q);
	  for (q = 0; q < Q; q++) {
	    int kqg = qg2index(q,g,Q,G);
	    nuActive[q] = nu[kqg];
	  }
	  
	  double value = quadratic(var,nuActive);
	  lambda += 0.5 * value;
	  s += 0.5 * Q;
	}
	
	pot += ran.PotentialInverseGamma(s,lambda,oldGamma2);
	
	//
	// compute potential for new and old states
	//
	
	pot -= potentialRho(Q,rho,nuRho);
	pot -= potentialGamma2();
	pot -= potentialNu(Q,G,nu,*gamma2,a,rho,tau2Rho,sigma2);

	pot += potentialRho(Q,newRhoVector,nuRho);
	pot += potentialGamma2();
	pot += potentialNu(Q,G,nu,newGamma2,a,newRhoVector,tau2Rho,sigma2);

	if (ran.Unif01() <= exp(- pot))
	  {
	    (*nAccept)++;
	    
	    *gamma2 = newGamma2;
	    int kk;
	    for (kk = 0; kk < Q * Q; kk++)
	      rho[kk] = newRhoVector[kk];
	  }
      }
    }
    
    free(newRhoVector);
  }
    
  *seed = ran.ChangeSeed(*seed);
  
  return;
}






inline void updateSigma2(unsigned int *seed,
			 int nTry,
			 int *nAccept,
			 double epsilon,
			 double *sigma2,
			 int Q,
			 int G,
			 const int *S,
			 const double *x,
			 const int *psi,
			 const double *nu,
			 const int *delta,
			 const double *Delta,
			 double c2,
			 double gamma2,
			 const double *r,
			 const double *rho,
			 const double *phi,
			 const double *t,
			 const double *l,
			 const double *tau2R,
			 const double *tau2Rho,
			 const double *a,
			 const double *b) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {

    //
    // propose new value
    //

    int q = (int) (ran.Unif01() * Q);
    int g = (int) (ran.Unif01() * G);

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    int kqg = qg2index(q,g,Q,G);
    double oldValue = sigma2[kqg];
    double newValue = oldValue * u;

    double pot = - log(1.0 / u);

    //
    // subtract potential for old state
    //

    std::vector<int> on(Q,0);
    int qq;
    for (qq = 0; qq < Q; qq++) {
      int index = qg2index(qq,g,Q,G);
      on[qq] = delta[index];
    }

    pot -= potentialSigma2qg(q,g,Q,G,sigma2,l,t);
    pot -= potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    pot -= potentialNug(Q,nu + g * Q,gamma2,a,rho,tau2Rho,sigma2 + g * Q);
    pot -= potentialDDeltag(Q,on,Delta + g * Q,c2,b,r,tau2R,sigma2 + g * Q);

    //
    // add potential for new state
    //

    sigma2[kqg] = newValue;
    pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);
    pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    pot += potentialNug(Q,nu + g * Q,gamma2,a,rho,tau2Rho,sigma2 + g * Q);
    pot += potentialDDeltag(Q,on,Delta + g * Q,c2,b,r,tau2R,sigma2 + g * Q);
    sigma2[kqg] = oldValue;

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      sigma2[kqg] = newValue;
      (*nAccept)++;
    }
  }
    
  *seed = ran.ChangeSeed(*seed);
  
  return;
}















#endif
