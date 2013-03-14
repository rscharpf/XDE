//C++ code that is wrapped by Rinterface.cpp

#include <vector>
#include <iostream>
#include <fstream>

#include "Random_v2.h"
#include "Cholesky.h"
#include "Potential_v2.h"
#include "Utility_v2.h"

void updateANu(unsigned int *seed,
	       int nTry,
	       int *nAccept,
	       double epsilon,
	       double *a,
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

    //
    // propose new values for nu from full conditionals
    //

    double *newNu = (double *) calloc(Q * G,sizeof(double));

    a[q] = newValue;
    pot -= nuGibbs(newNu,Q,G,S,gamma2,tau2Rho,a,rho,sigma2,
		   phi,psi,x,delta,Delta,ran,1);
    a[q] = oldValue;
    pot += nuGibbs(nu,Q,G,S,gamma2,tau2Rho,a,rho,sigma2,
		   phi,psi,x,delta,Delta,ran,1);


    pot -= potentialA(Q,a,pA0,pA1,alphaA,betaA);
    pot -= potentialNu(Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    a[q] = newValue;
    pot += potentialA(Q,a,pA0,pA1,alphaA,betaA);
    pot += potentialNu(Q,G,newNu,gamma2,a,rho,tau2Rho,sigma2);
    pot += potentialX(Q,G,S,x,psi,newNu,delta,Delta,sigma2,phi);
    a[q] = oldValue;

    cout << "value of pot in updateANu : " << pot << endl;

    if (ran.Unif01() <= exp(- pot)) {
      a[q] = newValue;
      int k;
      for (k = 0; k < Q * G; k++)
	nu[k] = newNu[k];

      (*nAccept)++;
    }

    free(newNu);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}






void updateBDDelta(unsigned int *seed,
		   int nTry,
		   int *nAccept,
		   double epsilon,
		   double *b,
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

    //
    // propose new values for Delta from full conditionals
    //

    double *newDDelta = (double *) calloc(Q * G,sizeof(double));

    b[q] = newValue;
    pot -= DeltaGibbs(newDDelta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);
    b[q] = oldValue;
    pot += DeltaGibbs(Delta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);


    pot -= potentialB(Q,b,pB0,pB1,alphaB,betaB);
    pot -= potentialDDelta(Q,G,delta,Delta,c2,b,r,tau2R,sigma2);
    pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    b[q] = newValue;
    pot += potentialB(Q,b,pB0,pB1,alphaB,betaB);
    pot += potentialDDelta(Q,G,delta,newDDelta,c2,b,r,tau2R,sigma2);
    pot += potentialX(Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);
    b[q] = oldValue;

    if (ran.Unif01() <= exp(- pot)) {
      b[q] = newValue;
      int k;
      for (k = 0; k < Q * G; k++)
	Delta[k] = newDDelta[k];

      (*nAccept)++;
    }

    free(newDDelta);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}








void updateTau2RhoNu(unsigned int *seed,
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

      pot -= nuGibbs(newNu,Q,G,S,gamma2,newValues,a,rho,sigma2,
		     phi,psi,x,delta,Delta,ran,1);
      pot += nuGibbs(nu,Q,G,S,gamma2,oldValues,a,rho,sigma2,
		     phi,psi,x,delta,Delta,ran,1);

      pot -= potentialTau2Rho();
      pot -= potentialNu(Q,G,nu,gamma2,a,rho,oldValues,sigma2);
      pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

      pot += potentialTau2Rho();
      pot += potentialNu(Q,G,newNu,gamma2,a,rho,newValues,sigma2);
      pot += potentialX(Q,G,S,x,psi,newNu,delta,Delta,sigma2,phi);

      if (ran.Unif01() <= exp(- pot)) {
	int q;
	for (q = 0; q < Q; q++)
	  tau2Rho[q] = newValues[q];
	int k;
	for (k = 0; k < Q * G; k++)
	  nu[k] = newNu[k];

	(*nAccept)++;
      }

      free(newNu);
      free(oldValues);
      free(newValues);
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}




void updateTau2RDDelta(unsigned int *seed,
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

      pot -= DeltaGibbs(newDDelta,Q,G,S,c2,newValues,b,r,sigma2,
			phi,psi,x,delta,nu,ran,1);
      pot += DeltaGibbs(Delta,Q,G,S,c2,oldValues,b,r,sigma2,
			phi,psi,x,delta,nu,ran,1);


      pot -= potentialTau2R();
      pot -= potentialDDelta(Q,G,delta,Delta,c2,b,r,tau2R,sigma2);
      pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

      pot += potentialTau2R();
      pot += potentialDDelta(Q,G,delta,newDDelta,c2,b,r,newValues,sigma2);
      pot += potentialX(Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);

      if (ran.Unif01() <= exp(- pot)) {
	int q;
	for (q = 0; q < Q; q++)
	  tau2R[q] = newValues[q];

	int k;
	for (k = 0; k < Q * G; k++)
	  Delta[k] = newDDelta[k];

	(*nAccept)++;
      }

      free(newDDelta);
      free(oldValues);
      free(newValues);
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}



void updateNu(unsigned int *seed,
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
	  phi,psi,x,delta,Delta,ran,1);
  (*nAccept) += Q * G;

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateDDelta(unsigned int *seed,
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

  DeltaGibbs(Delta,Q,G,S,c2,tau2R,b,r,sigma2,phi,psi,x,delta,nu,ran,1);
  (*nAccept) += Q * G;

  *seed = ran.ChangeSeed(*seed);

  return;
}




void updateC2(unsigned int *seed,
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
	    makeSigma(g,G,varInv,on,Q,1.0,tau2R,b,sigma2,r);

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




void updateC2DDelta(unsigned int *seed,
		    int nTry,
		    int *nAccept,
		    double epsilon,
		    double *c2,
		    double *Delta,
		    int Q,
		    int G,
		    const int *S,
		    const double *x,
		    const int *psi,
		    const double *nu,
		    const int *delta,
		    const double *r,
		    const double *sigma2,
		    const double *phi,
		    const double *tau2R,
		    const double *b,
		    double c2Max) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;

    double u = lower + (upper - lower) * ran.Unif01();
    double oldValue = *c2;
    double newValue = oldValue * u;

    if (newValue > c2Max) return;

    double pot = - log(1.0 / u);

    double *newDDelta = (double *) calloc(Q * G,sizeof(double));

    pot -= DeltaGibbs(newDDelta,Q,G,S,newValue,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);
    pot += DeltaGibbs(Delta,Q,G,S,oldValue,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);

    pot -= potentialC2();
    pot -= potentialDDelta(Q,G,delta,Delta,oldValue,b,r,tau2R,sigma2);
    pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    pot += potentialC2();
    pot += potentialDDelta(Q,G,delta,newDDelta,newValue,b,r,tau2R,sigma2);
    pot += potentialX(Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);

    if (ran.Unif01() <= exp(- pot)) {
      *c2 = newValue;
      int kk;
      for (kk = 0; kk < Q * G; kk++)
	Delta[kk] = newDDelta[kk];

      (*nAccept)++;
    }

    free(newDDelta);
  }



  return;
}




void updateGamma2(unsigned int *seed,
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
      makeSigma(g,G,varInv,Q,1.0,tau2Rho,a,sigma2,rho);

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





void updateGamma2Nu(unsigned int *seed,
		    int nTry,
		    int *nAccept,
		    double epsilon,
		    double *gamma2,
		    double *nu,
		    int Q,
		    int G,
		    const int *S,
		    const double *x,
		    const int *psi,
		    const int *delta,
		    const double *Delta,
		    const double *rho,
		    const double *sigma2,
		    const double *phi,
		    const double *tau2Rho,
		    const double *a) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;

    double u = lower + (upper - lower) * ran.Unif01();
    double oldValue = *gamma2;
    double newValue = oldValue * u;

    double pot = - log(1.0 / u);

    double *newNu = (double *) calloc(Q * G,sizeof(double));

    pot -= nuGibbs(newNu,Q,G,S,newValue,tau2Rho,a,rho,sigma2,
		   phi,psi,x,delta,Delta,ran,1);
    pot += nuGibbs(nu,Q,G,S,oldValue,tau2Rho,a,rho,sigma2,
		   phi,psi,x,delta,Delta,ran,1);


    pot -= potentialGamma2();
    pot -= potentialNu(Q,G,nu,oldValue,a,rho,tau2Rho,sigma2);
    pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    pot += potentialGamma2();
    pot += potentialNu(Q,G,newNu,newValue,a,rho,tau2Rho,sigma2);
    pot += potentialX(Q,G,S,x,psi,newNu,delta,Delta,sigma2,phi);

    if (ran.Unif01() <= exp(- pot)) {
      *gamma2 = newValue;
      int k;
      for (k = 0; k < Q * G; k++)
	nu[k] = newNu[k];

      (*nAccept)++;
    }

    free(newNu);
  }

  return;
}





void updateRDDelta(unsigned int *seed,
		   int nTry,
		   int *nAccept,
		   double epsilon,
		   double *r,
		   double *Delta,
		   int Q,
		   int G,
		   const int *S,
		   const double *x,
		   const int *psi,
		   const double *nu,
		   const int *delta,
		   double c2,
		   const double *sigma2,
		   const double *phi,
		   const double *tau2R,
		   const double *b,
		   double nuR) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    vector<vector<double> > oldR;
    vector<vector<double> > newR;
    oldR.resize(Q);
    newR.resize(Q);
    int p,q;
    for (p = 0; p < Q; p++) {
      oldR[p].resize(Q);
      newR[p].resize(Q);
    }
    for (p = 0; p < Q; p++) {
      oldR[p][p] = 1.0;
      newR[p][p] = 1.0;
      for (q = p + 1; q < Q; q++) {
	int kqq = qq2index(p,q,Q);
	oldR[p][q] = r[kqq];
	newR[p][q] = r[kqq];
	oldR[q][p] = r[kqq];
	newR[q][p] = r[kqq];
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

    newR[pp][qq] = oldR[pp][qq] * exp(u) /
      (1 - oldR[pp][qq] + oldR[pp][qq] * exp(u));
    newR[qq][pp] = newR[pp][qq];
    double *newRVector = (double *) calloc(Q * (Q - 1) / 2,sizeof(double));
    for (p = 0; p < Q; p++)
      for (q = p + 1; q < Q; q++) {
	int kqq = qq2index(p,q,Q);
	newRVector[kqq] = newR[p][q];
      }

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
	// propose new values for Delta from full conditional
	//

	double *newDDelta = (double *) calloc(Q * G,sizeof(double));

	pot -= DeltaGibbs(newDDelta,Q,G,S,c2,tau2R,b,newRVector,sigma2,
			  phi,psi,x,delta,nu,ran,1);
	pot += DeltaGibbs(Delta,Q,G,S,c2,tau2R,b,r,sigma2,
			  phi,psi,x,delta,nu,ran,1);


	pot -= potentialR(Q,r,nuR);
	pot -= potentialDDelta(Q,G,delta,Delta,c2,b,r,tau2R,sigma2);
	pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

	pot += potentialR(Q,newRVector,nuR);
	pot += potentialDDelta(Q,G,delta,newDDelta,c2,b,newRVector,
			       tau2R,sigma2);
	pot += potentialX(Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);

	if (ran.Unif01() <= exp(- pot)) {
	  int kk;
	  for (kk = 0; kk < Q * (Q - 1) / 2; kk++)
	    r[kk] = newRVector[kk];

	  for (kk = 0; kk < Q * G; kk++)
	    Delta[kk] = newDDelta[kk];

	  (*nAccept)++;
	}

	free(newDDelta);
      }
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateRC2(unsigned int *seed,
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
    for (p = 0; p < Q; p++) {
      oldR[p].resize(Q);
      newR[p].resize(Q);
    }
    for (p = 0; p < Q; p++) {
      oldR[p][p] = 1.0;
      newR[p][p] = 1.0;
      for (q = p + 1; q < Q; q++) {
	int kqq = qq2index(p,q,Q);
	oldR[p][q] = r[kqq];
	newR[p][q] = r[kqq];
	oldR[q][p] = r[kqq];
	newR[q][p] = r[kqq];
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
    double *newRVector = (double *) calloc(Q * (Q - 1) / 2,sizeof(double));
    for (p = 0; p < Q; p++)
      for (q = p + 1; q < Q; q++) {
	int kqq = qq2index(p,q,Q);
	newRVector[kqq] = newR[p][q];
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
	    makeSigma(g,G,varInv,on,Q,1.0,tau2R,b,sigma2,newRVector);

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
	    makeSigma(g,G,varInv,on,Q,1.0,tau2R,b,sigma2,r);

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
	    for (kk = 0; kk < Q * (Q - 1) / 2; kk++)
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



void updateRhoNu(unsigned int *seed,
		 int nTry,
		 int *nAccept,
		 double epsilon,
		 double *rho,
		 double *nu,
		 int Q,
		 int G,
		 const int *S,
		 const double *x,
		 const int *psi,
		 const int *delta,
		  const double *Delta,
		  double gamma2,
		  const double *sigma2,
		  const double *phi,
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
     for (p = 0; p < Q; p++) {
       oldRho[p].resize(Q);
       newRho[p].resize(Q);
     }
     for (p = 0; p < Q; p++) {
       oldRho[p][p] = 1.0;
       newRho[p][p] = 1.0;
       for (q = p + 1; q < Q; q++)
	 {
	   int kqq = qq2index(p,q,Q);
	   oldRho[p][q] = rho[kqq];
	   newRho[p][q] = rho[kqq];
	   oldRho[q][p] = rho[kqq];
	   newRho[q][p] = rho[kqq];
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
     double *newRhoVector = (double *) calloc(Q * (Q - 1) / 2,sizeof(double));
     for (p = 0; p < Q; p++)
       for (q = p + 1; q < Q; q++) {
	 int kqq = qq2index(p,q,Q);
	 newRhoVector[kqq] = newRho[p][q];
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
	 // propose new values for nu from full conditional
	 //

	 double *newNu = (double *) calloc(Q * G,sizeof(double));

	 pot -= nuGibbs(newNu,Q,G,S,gamma2,tau2Rho,a,newRhoVector,sigma2,
			phi,psi,x,delta,Delta,ran,1);
	 pot += nuGibbs(nu,Q,G,S,gamma2,tau2Rho,a,rho,sigma2,
			phi,psi,x,delta,Delta,ran,1);


	 pot -= potentialRho(Q,rho,nuRho);
	 pot -= potentialNu(Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
	 pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

	 pot += potentialRho(Q,newRhoVector,nuRho);
	 pot += potentialNu(Q,G,newNu,gamma2,a,newRhoVector,tau2Rho,sigma2);
	 pot += potentialX(Q,G,S,x,psi,newNu,delta,Delta,sigma2,phi);

	 if (ran.Unif01() <= exp(- pot)) {
	   int kk;
	   for (kk = 0; kk < Q * (Q - 1) / 2; kk++)
	     rho[kk] = newRhoVector[kk];

	   for (kk = 0; kk < Q * G; kk++)
	     nu[kk] = newNu[kk];

	   (*nAccept)++;
	 }

	 free(newNu);
       }
     }
   }

   *seed = ran.ChangeSeed(*seed);

   return;
 }





 void updateRhoGamma2(unsigned int *seed,
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
     for (p = 0; p < Q; p++) {
       oldRho[p].resize(Q);
       newRho[p].resize(Q);
     }
     for (p = 0; p < Q; p++) {
       oldRho[p][p] = 1.0;
       newRho[p][p] = 1.0;
       for (q = p + 1; q < Q; q++) {
	 int kqq = qq2index(p,q,Q);
	 oldRho[p][q] = rho[kqq];
	 newRho[p][q] = rho[kqq];
	 oldRho[q][p] = rho[kqq];
	 newRho[q][p] = rho[kqq];
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
     double *newRhoVector = (double *) calloc(Q * (Q - 1) / 2,sizeof(double));
     for (p = 0; p < Q; p++)
       for (q = p + 1; q < Q; q++) {
	 int kqq = qq2index(p,q,Q);
	 newRhoVector[kqq] = newRho[p][q];
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
	  makeSigma(g,G,varInv,Q,1.0,tau2Rho,a,sigma2,newRhoVector);

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
	  makeSigma(g,G,varInv,Q,1.0,tau2Rho,a,sigma2,rho);

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
	    for (kk = 0; kk < Q * (Q - 1) / 2; kk++)
	      rho[kk] = newRhoVector[kk];
	  }
      }
    }

    free(newRhoVector);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}






void updateSigma2(unsigned int *seed,
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
    pot -= potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    pot -= potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);

    //
    // add potential for new state
    //

    sigma2[kqg] = newValue;
    pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);
    pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    pot += potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    pot += potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);
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







void updateSigma2_HyperInverseWishart(unsigned int *seed,
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
				      double gamma2,
				      const double *r,
				      const double *rho,
				      const double *phi,
				      const double *t,
				      const double *l,
				      const double *tau2R,
				      const double *tau2Rho,
				      const double *a,
				      const double *b,
				      const vector<vector<vector<double> > > &Omega,
				      const vector<int> &oldClique,
				      const vector<vector<int> > &oldComponents) {
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
    pot -= potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    pot -= potentialDDeltaStar_HyperInverseWishart(g,Delta,b,sigma2,tau2R,r,Q,G,Omega,oldClique,oldComponents);

    //
    // add potential for new state
    //

    sigma2[kqg] = newValue;
    pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);
    pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    pot += potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    pot += potentialDDeltaStar_HyperInverseWishart(g,Delta,b,sigma2,tau2R,r,Q,G,Omega,oldClique,oldComponents);
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







void updatePhi(unsigned int *seed,
	       int nTry,
	       int *nAccept,
	       double epsilon,
	       double *phi,
	       int Q,
	       int G,
	       const int *S,
	       const double *x,
	       const int *psi,
	       const double *nu,
	       const int *delta,
	       const double *Delta,
	       const double *sigma2,
	       const double *theta,
	       const double *lambda) {
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
    double oldValue = phi[kqg];
    double newValue = oldValue * u;

    double pot = - log(1.0 / u);

    //
    // subtract potential for old state
    //

    pot -= potentialPhiqg(q,g,Q,G,phi,lambda,theta);
    pot -= potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    //
    // add potential for new state
    //

    phi[kqg] = newValue;
    pot += potentialPhiqg(q,g,Q,G,phi,lambda,theta);
    pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    phi[kqg] = oldValue;

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      phi[kqg] = newValue;
      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}






void updateTheta(unsigned int *seed,
		 int nTry,
		 int *nAccept,
		 double epsilon,
		 double *theta,
		 int Q,
		 int G,
		 const double *phi,
		 const double *lambda) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {

    //
    // propose new value
    //

    int q = (int) (ran.Unif01() * Q);

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = theta[q];
    double newValue = oldValue * u;

    double pot = - log(1.0 / u);

    //
    // subtract potential for old state
    //

    pot -= potentialTheta();
    int g;
    for (g = 0; g < G; g++)
      pot -= potentialPhiqg(q,g,Q,G,phi,lambda,theta);

    //
    // add potential for new state
    //

    theta[q] = newValue;
    pot += potentialTheta();
    for (g = 0; g < G; g++)
      pot += potentialPhiqg(q,g,Q,G,phi,lambda,theta);
    theta[q] = oldValue;

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      theta[q] = newValue;
      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}






void updateLambda(unsigned int *seed,
		  int nTry,
		  int *nAccept,
		  double epsilon,
		  double *lambda,
		  int Q,
		  int G,
		  const double *phi,
		  const double *theta) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {

    //
    // propose new value
    //

    int q = (int) (ran.Unif01() * Q);

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = lambda[q];
    double newValue = oldValue * u;

    double pot = - log(1.0 / u);

    //
    // subtract potential for old state
    //

    pot -= potentialLambda();
    int g;
    for (g = 0; g < G; g++)
      pot -= potentialPhiqg(q,g,Q,G,phi,lambda,theta);

    //
    // add potential for new state
    //

    lambda[q] = newValue;
    pot += potentialLambda();
    for (g = 0; g < G; g++)
      pot += potentialPhiqg(q,g,Q,G,phi,lambda,theta);
    lambda[q] = oldValue;

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      lambda[q] = newValue;
      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}






void updateT(unsigned int *seed,
	     int nTry,
	     int *nAccept,
	     double epsilon,
	     double *t,
	     int Q,
	     int G,
	     const double *sigma2,
	     const double *l) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {

    //
    // propose new value
    //

    int q = (int) (ran.Unif01() * Q);

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = t[q];
    double newValue = oldValue * u;

    double pot = - log(1.0 / u);

    //
    // subtract potential for old state
    //

    pot -= potentialT();
    int g;
    for (g = 0; g < G; g++)
      pot -= potentialSigma2qg(q,g,Q,G,sigma2,l,t);

    //
    // add potential for new state
    //

    t[q] = newValue;
    pot += potentialT();
    for (g = 0; g < G; g++)
      pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);
    t[q] = oldValue;

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      t[q] = newValue;
      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateL(unsigned int *seed,
	     int nTry,
	     int *nAccept,
	     double epsilon,
	     double *l,
	     int Q,
	     int G,
	     const double *sigma2,
	     const double *t) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {

    //
    // propose new value
    //

    int q = (int) (ran.Unif01() * Q);

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = l[q];
    double newValue = oldValue * u;

    double pot = - log(1.0 / u);

    //
    // subtract potential for old state
    //

    pot -= potentialL();
    int g;
    for (g = 0; g < G; g++)
      pot -= potentialSigma2qg(q,g,Q,G,sigma2,l,t);

    //
    // add potential for new state
    //

    l[q] = newValue;
    pot += potentialL();
    for (g = 0; g < G; g++)
      pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);
    l[q] = oldValue;

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      l[q] = newValue;
      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateXi(unsigned int *seed,
	      int *nAccept,
	      double *xi,
	      int Q,
	      int G,
	      const int *delta,
	      double alphaXi,
	      double betaXi) {
  Random ran(*seed);

  int q;
  for (q = 0; q < Q; q++) {
    //
    // set prior parameters
    //

    double alpha = alphaXi;
    double beta = betaXi;

    //
    // update parameters based on available observations
    //

    int g;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      if (delta[kqg] == 1)
	alpha += 1.0;
      else
	beta += 1.0;
    }

    //
    // Draw new value
    //

    double newValue = ran.Beta(alpha,beta);
    xi[q] = newValue;
    (*nAccept)++;
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}






void updateXi_onedelta(unsigned int *seed,
		       int *nAccept,
		       double *xi,
		       int Q,
		       int G,
		       const int *delta,
		       double alphaXi,
		       double betaXi) {
  Random ran(*seed);

  int q,g;
  for (g = 0; g < G; g++) {
    int nOn = 0;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      nOn += delta[kqg];
    }
    if (nOn != 0 && nOn != Q) {
      cout << "Error found in function \"updateXi_onedelta\":" << endl;
      cout << "All delta's for any gene must be equal." << endl;
      cout << "For gene \"" << g << "\" this requirement is not fulfilled." <<
	endl;
      cout << "Aborting." << endl;
      exit(-1);
    }
  }

  double alpha = alphaXi;
  double beta = betaXi;

  //
  // update parameters based on available observations
  //

  for (g = 0; g < G; g++) {
    int kqg = qg2index(0,g,Q,G);
    if (delta[kqg] == 1)
      alpha += 1.0;
    else
      beta += 1.0;
  }

  //
  // Draw new value
  //

  double newValue = ran.Beta(alpha,beta);
  for (q = 0; q < Q; q++)
    xi[q] = newValue;
  (*nAccept)++;

  *seed = ran.ChangeSeed(*seed);

  return;
}




void updateDeltaDDelta(unsigned int *seed,
		       int nTry,
		       int *nAccept,
		       int *delta,
		       double *Delta,
		       int Q,
		       int G,
		       const int *S,
		       const double *x,
		       const int *psi,
		       const double *nu,
		       double c2,
		       const double *r,
		       const double *sigma2,
		       const double *phi,
		       const double *tau2R,
		       const double *xi,
		       const double *b) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    int q = (int) (ran.Unif01() * Q);
    int g = (int) (ran.Unif01() * G);

    int kqg = qg2index(q,g,Q,G);
    int oldDelta = delta[kqg];
    int newDelta = 1 - oldDelta;

    //
    // propose new values for Delta from full conditionals
    //

    double *newDDelta = (double *) calloc(Q * G,sizeof(double));

    delta[kqg] = newDelta;
    pot -= DeltaGibbs(g,newDDelta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);
    delta[kqg] = oldDelta;
    pot += DeltaGibbs(g,Delta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);

    delta[kqg] = oldDelta;
    pot -= potentialDeltag(g,Q,G,delta,xi);
    pot -= potentialDDeltag(g,Q,G,delta,Delta,c2,b,r,tau2R,sigma2);
    pot -= potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    delta[kqg] = newDelta;
    pot += potentialDeltag(g,Q,G,delta,xi);
    pot += potentialDDeltag(g,Q,G,delta,newDDelta,c2,b,r,tau2R,sigma2);
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);
    delta[kqg] = oldDelta;

    if (ran.Unif01() <= exp(- pot)) {
      delta[kqg] = newDelta;
      int qq;
      for (qq = 0; qq < Q; qq++) {
	int index = qg2index(qq,g,Q,G);
	if (delta[index] == 1)
	  Delta[index] = newDDelta[index];
      }

      (*nAccept)++;
    }

    free(newDDelta);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}



void updateDeltaDDelta_onedelta(unsigned int *seed,
				int nTry,
				int *nAccept,
				int *delta,
				double *Delta,
				int Q,
				int G,
				const int *S,
				const double *x,
				const int *psi,
				const double *nu,
				double c2,
				const double *r,
				const double *sigma2,
				const double *phi,
				const double *tau2R,
				const double *xi,
				const double *b) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    int g = (int) (ran.Unif01() * G);

    int nOn = 0;
    int q;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      nOn += delta[kqg];
    }
    if (nOn != 0 && nOn != Q) {
      cout << "Error found in function \"updateDeltaDDelta_onedelta\":" << endl;
      cout << "All delta's for any gene must be equal." << endl;
      cout << "For gene \"" << g << "\" this requirement is not fulfilled." <<
	endl;
      cout << "Aborting." << endl;
      exit(-1);
    }

    int kqg = qg2index(0,g,Q,G);
    int oldDelta = delta[kqg];
    int newDelta = 1 - oldDelta;

    //
    // propose new values for Delta from full conditionals
    //

    double *newDDelta = (double *) calloc(Q * G,sizeof(double));

    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = newDelta;
    }
    pot -= DeltaGibbs(g,newDDelta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }
    pot += DeltaGibbs(g,Delta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);


    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }
    pot -= potentialDeltag_onedelta(g,Q,G,delta,xi);
    pot -= potentialDDeltag(g,Q,G,delta,Delta,c2,b,r,tau2R,sigma2);
    pot -= potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = newDelta;
    }
    pot += potentialDeltag_onedelta(g,Q,G,delta,xi);
    pot += potentialDDeltag(g,Q,G,delta,newDDelta,c2,b,r,tau2R,sigma2);
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }

    if (ran.Unif01() <= exp(- pot)) {
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	delta[kqg] = newDelta;
	if (newDelta == 1)
	  Delta[kqg] = newDDelta[kqg];
      }

      (*nAccept)++;
    }

    free(newDDelta);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}




void updateLSigma2(unsigned int *seed,
		   int nTry,
		   int *nAccept,
		   double epsilon,
		   double *l,
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

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = l[q];
    double newValue = oldValue * u;
    double *oldSigma2 = (double *) calloc(G,sizeof(double));
    double *newSigma2 = (double *) calloc(G,sizeof(double));
    int g;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      oldSigma2[g] = sigma2[kqg];
      newSigma2[g] = oldSigma2[g] + oldValue * (u - 1.0);
    }

    double pot = - log(1.0 / u);

    //
    // subtract potential for old state
    //

    pot -= potentialL();
    for (g = 0; g < G; g++) {
      std::vector<int> on(Q,0);
      int qq;
      for (qq = 0; qq < Q; qq++) {
	int index = qg2index(qq,g,Q,G);
	on[qq] = delta[index];
      }

      pot -= potentialSigma2qg(q,g,Q,G,sigma2,l,t);
      pot -= potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
      pot -= potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
      pot -= potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);
    }

    //
    // add potential for new state
    //

    l[q] = newValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      sigma2[kqg] = newSigma2[g];
    }
    pot += potentialL();
    for (g = 0; g < G; g++) {
      std::vector<int> on(Q,0);
      int qq;
      for (qq = 0; qq < Q; qq++) {
	int index = qg2index(qq,g,Q,G);
	on[qq] = delta[index];
      }

      pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);
      pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
      pot += potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
      pot += potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);
    }
    l[q] = oldValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      sigma2[kqg] = oldSigma2[g];
    }

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      l[q] = newValue;
      for (g = 0; g < G; g++) {
	int kqg = qg2index(q,g,Q,G);
	sigma2[kqg] = newSigma2[g];
      }

      (*nAccept)++;
    }

    free(oldSigma2);
    free(newSigma2);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateLSigma2_HyperInverseWishart(unsigned int *seed,
				       int nTry,
				       int *nAccept,
				       double epsilon,
				       double *l,
				       double *sigma2,
				       int Q,
				       int G,
				       const int *S,
				       const double *x,
				       const int *psi,
				       const double *nu,
				       const int *delta,
				       const double *Delta,
				       double gamma2,
				       const double *r,
				       const double *rho,
				       const double *phi,
				       const double *t,
				       const double *tau2R,
				       const double *tau2Rho,
				       const double *a,
				       const double *b,
				       const vector<vector<vector<double> > > &Omega,
				       const vector<int> &oldClique,
				       const vector<vector<int> > &oldComponents) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {

    //
    // propose new value
    //

    int q = (int) (ran.Unif01() * Q);

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = l[q];
    double newValue = oldValue * u;
    double *oldSigma2 = (double *) calloc(G,sizeof(double));
    double *newSigma2 = (double *) calloc(G,sizeof(double));
    int g;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      oldSigma2[g] = sigma2[kqg];
      newSigma2[g] = oldSigma2[g] + oldValue * (u - 1.0);
    }

    double pot = - log(1.0 / u);

    //
    // subtract potential for old state
    //

    pot -= potentialL();
    for (g = 0; g < G; g++) {
      std::vector<int> on(Q,0);
      int qq;
      for (qq = 0; qq < Q; qq++) {
	int index = qg2index(qq,g,Q,G);
	on[qq] = delta[index];
      }

      pot -= potentialSigma2qg(q,g,Q,G,sigma2,l,t);
      pot -= potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
      pot -= potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    }
    pot -= potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,Q,G,Omega,oldClique,oldComponents);

    //
    // add potential for new state
    //

    l[q] = newValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      sigma2[kqg] = newSigma2[g];
    }
    pot += potentialL();
    for (g = 0; g < G; g++) {
      std::vector<int> on(Q,0);
      int qq;
      for (qq = 0; qq < Q; qq++) {
	int index = qg2index(qq,g,Q,G);
	on[qq] = delta[index];
      }

      pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);
      pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
      pot += potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    }
    pot += potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,Q,G,Omega,oldClique,oldComponents);

    l[q] = oldValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      sigma2[kqg] = oldSigma2[g];
    }

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      l[q] = newValue;
      for (g = 0; g < G; g++) {
	int kqg = qg2index(q,g,Q,G);
	sigma2[kqg] = newSigma2[g];
      }

      (*nAccept)++;
    }

    free(oldSigma2);
    free(newSigma2);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateTSigma2(unsigned int *seed,
		   int nTry,
		   int *nAccept,
		   double epsilon,
		   double *t,
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

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = t[q];
    double newValue = oldValue * u;
    double *oldSigma2 = (double *) calloc(G,sizeof(double));
    double *newSigma2 = (double *) calloc(G,sizeof(double));
    int g;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      oldSigma2[g] = sigma2[kqg];
      newSigma2[g] = sqrt(u) * (oldSigma2[g] - l[q]) + l[q];
    }

    double pot = - (((double) G) / 2.0 - 1.0) * log(u);

    //
    // subtract potential for old state
    //

    pot -= potentialT();
    for (g = 0; g < G; g++) {
      std::vector<int> on(Q,0);
      int qq;
      for (qq = 0; qq < Q; qq++) {
	int index = qg2index(qq,g,Q,G);
	on[qq] = delta[index];
      }

      pot -= potentialSigma2qg(q,g,Q,G,sigma2,l,t);
      pot -= potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
      pot -= potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
      pot -= potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);
    }

    //
    // add potential for new state
    //

    t[q] = newValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      sigma2[kqg] = newSigma2[g];
    }
    pot += potentialT();
    for (g = 0; g < G; g++) {
      std::vector<int> on(Q,0);
      int qq;
      for (qq = 0; qq < Q; qq++) {
	int index = qg2index(qq,g,Q,G);
	on[qq] = delta[index];
      }

      pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);
      pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
      pot += potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
      pot += potentialDDeltag(g,Q,G,on,Delta,c2,b,r,tau2R,sigma2);
    }
    t[q] = oldValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      sigma2[kqg] = oldSigma2[g];
    }

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      t[q] = newValue;
      for (g = 0; g < G; g++) {
	int kqg = qg2index(q,g,Q,G);
	sigma2[kqg] = newSigma2[g];
      }

      (*nAccept)++;
    }

    free(oldSigma2);
    free(newSigma2);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateTSigma2_HyperInverseWishart(unsigned int *seed,
				       int nTry,
				       int *nAccept,
				       double epsilon,
				       double *t,
				       double *sigma2,
				       int Q,
				       int G,
				       const int *S,
				       const double *x,
				       const int *psi,
				       const double *nu,
				       const int *delta,
				       const double *Delta,
				       double gamma2,
				       const double *r,
				       const double *rho,
				       const double *phi,
				       const double *l,
				       const double *tau2R,
				       const double *tau2Rho,
				       const double *a,
				       const double *b,
				       const vector<vector<vector<double> > > &Omega,
				       const vector<int> &oldClique,
				       const vector<vector<int> > &oldComponents) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {

    //
    // propose new value
    //

    int q = (int) (ran.Unif01() * Q);

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = t[q];
    double newValue = oldValue * u;
    double *oldSigma2 = (double *) calloc(G,sizeof(double));
    double *newSigma2 = (double *) calloc(G,sizeof(double));
    int g;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      oldSigma2[g] = sigma2[kqg];
      newSigma2[g] = sqrt(u) * (oldSigma2[g] - l[q]) + l[q];
    }

    double pot = - (((double) G) / 2.0 - 1.0) * log(u);

    //
    // subtract potential for old state
    //

    pot -= potentialT();
    for (g = 0; g < G; g++) {
      std::vector<int> on(Q,0);
      int qq;
      for (qq = 0; qq < Q; qq++) {
	int index = qg2index(qq,g,Q,G);
	on[qq] = delta[index];
      }

      pot -= potentialSigma2qg(q,g,Q,G,sigma2,l,t);
      pot -= potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
      pot -= potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    }
    pot -= potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,Q,G,Omega,oldClique,oldComponents);


    //
    // add potential for new state
    //

    t[q] = newValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      sigma2[kqg] = newSigma2[g];
    }
    pot += potentialT();
    for (g = 0; g < G; g++) {
      std::vector<int> on(Q,0);
      int qq;
      for (qq = 0; qq < Q; qq++) {
	int index = qg2index(qq,g,Q,G);
	on[qq] = delta[index];
      }

      pot += potentialSigma2qg(q,g,Q,G,sigma2,l,t);
      pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
      pot += potentialNug(g,Q,G,nu,gamma2,a,rho,tau2Rho,sigma2);
    }
    pot += potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,Q,G,Omega,oldClique,oldComponents);

    t[q] = oldValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      sigma2[kqg] = oldSigma2[g];
    }

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      t[q] = newValue;
      for (g = 0; g < G; g++) {
	int kqg = qg2index(q,g,Q,G);
	sigma2[kqg] = newSigma2[g];
      }

      (*nAccept)++;
    }

    free(oldSigma2);
    free(newSigma2);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateLambdaPhi(unsigned int *seed,
		     int nTry,
		     int *nAccept,
		     double epsilon,
		     double *lambda,
		     double *phi,
		     int Q,
		     int G,
		     const int *S,
		     const double *x,
		     const int *psi,
		     const double *nu,
		     const int *delta,
		     const double *Delta,
		     const double *sigma2,
		     const double *theta) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {

    //
    // propose new value
    //

    int q = (int) (ran.Unif01() * Q);

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = lambda[q];
    double newValue = oldValue * u;
    double *oldPhi = (double *) calloc(G,sizeof(double));
    double *newPhi = (double *) calloc(G,sizeof(double));
    int g;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      oldPhi[g] = phi[kqg];
      newPhi[g] = oldPhi[g] + oldValue * (u - 1.0);
    }

    double pot = - log(1.0 / u);

    //
    // subtract potential for old state
    //

    pot -= potentialLambda();
    for (g = 0; g < G; g++) {
      pot -= potentialPhiqg(q,g,Q,G,phi,lambda,theta);
      pot -= potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    }

    //
    // add potential for new state
    //

    lambda[q] = newValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      phi[kqg] = newPhi[g];
    }
    pot += potentialLambda();
    for (g = 0; g < G; g++) {
      pot += potentialPhiqg(q,g,Q,G,phi,lambda,theta);
      pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    }
    lambda[q] = oldValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      phi[kqg] = oldPhi[g];
    }

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      lambda[q] = newValue;
      for (g = 0; g < G; g++) {
	int kqg = qg2index(q,g,Q,G);
	phi[kqg] = newPhi[g];
      }

      (*nAccept)++;
    }

    free(oldPhi);
    free(newPhi);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}






void updateThetaPhi(unsigned int *seed,
		    int nTry,
		    int *nAccept,
		    double epsilon,
		    double *theta,
		    double *phi,
		    int Q,
		    int G,
		    const int *S,
		    const double *x,
		    const int *psi,
		    const double *nu,
		    const int *delta,
		    const double *Delta,
		    const double *sigma2,
		    const double *lambda) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {

    //
    // propose new value
    //

    int q = (int) (ran.Unif01() * Q);

    double upper = 1.0 + epsilon;
    double lower = 1.0 / upper;
    double u = lower + (upper - lower) * ran.Unif01();

    double oldValue = theta[q];
    double newValue = oldValue * u;
    double *oldPhi = (double *) calloc(G,sizeof(double));
    double *newPhi = (double *) calloc(G,sizeof(double));
    int g;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      oldPhi[g] = phi[kqg];
      newPhi[g] = sqrt(u) * (oldPhi[g] - lambda[q]) + lambda[q];
    }

    double pot = - (((double) G) / 2.0 - 1.0) * log(u);

    //
    // subtract potential for old state
    //

    pot -= potentialTheta();
    for (g = 0; g < G; g++) {
      pot -= potentialPhiqg(q,g,Q,G,phi,lambda,theta);
      pot -= potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    }

    //
    // add potential for new state
    //

    theta[q] = newValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      phi[kqg] = newPhi[g];
    }
    pot += potentialTheta();
    for (g = 0; g < G; g++) {
      pot += potentialPhiqg(q,g,Q,G,phi,lambda,theta);
      pot += potentialXqg(q,g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    }
    theta[q] = oldValue;
    for (g = 0; g < G; g++) {
      int kqg = qg2index(q,g,Q,G);
      phi[kqg] = oldPhi[g];
    }

    //
    // accept or reject proposal
    //

    if (ran.Unif01() <= exp(- pot)) {
      theta[q] = newValue;
      for (g = 0; g < G; g++) {
	int kqg = qg2index(q,g,Q,G);
	phi[kqg] = newPhi[g];
      }

      (*nAccept)++;
    }

    free(oldPhi);
    free(newPhi);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}






void updateDeltaDDelta_MRF1_onedelta(unsigned int *seed,
				     int nTry,
				     int *nAccept,
				     int *delta,
				     double *Delta,
				     int Q,
				     int G,
				     const int *S,
				     const double *x,
				     const int *psi,
				     const double *nu,
				     double c2,
				     const double *r,
				     const double *sigma2,
				     const double *phi,
				     const double *tau2R,
				     const double *b,
				     const vector<vector<int> > &neighbour,
				     double eta0,
				     double omega0,
				     double kappa) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    int g = (int) (ran.Unif01() * G);

    int nOn = 0;
    int q;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      nOn += delta[kqg];
    }
    if (nOn != 0 && nOn != Q) {
      cout << "Error found in function \"updateDeltaDDelta_MRF1_onedelta\":" <<
	endl;
      cout << "All delta's for any gene must be equal." << endl;
      cout << "For gene \"" << g << "\" this requirement is not fulfilled." <<
	endl;
      cout << "Aborting." << endl;
      exit(-1);
    }

    int kqg = qg2index(0,g,Q,G);
    int oldDelta = delta[kqg];
    int newDelta = 1 - oldDelta;

    //
    // propose new values for Delta from full conditionals
    //

    double *newDDelta = (double *) calloc(Q * G,sizeof(double));

    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = newDelta;
    }
    pot -= DeltaGibbs(g,newDDelta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }
    pot += DeltaGibbs(g,Delta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);


    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }
    pot -= potentialDelta_MRF1_onedelta(Q,G,delta,neighbour,eta0,omega0,kappa);
    pot -= potentialDDeltag(g,Q,G,delta,Delta,c2,b,r,tau2R,sigma2);
    pot -= potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = newDelta;
    }
    pot += potentialDelta_MRF1_onedelta(Q,G,delta,neighbour,eta0,omega0,kappa);
    pot += potentialDDeltag(g,Q,G,delta,newDDelta,c2,b,r,tau2R,sigma2);
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }

    if (ran.Unif01() <= exp(- pot)) {
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	delta[kqg] = newDelta;
	if (newDelta == 1)
	  Delta[kqg] = newDDelta[kqg];
      }

      (*nAccept)++;
    }

    free(newDDelta);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateDeltaDDelta_MRF2_onedelta(unsigned int *seed,
				     int nTry,
				     int *nAccept,
				     int *delta,
				     double *Delta,
				     int Q,
				     int G,
				     const int *S,
				     const double *x,
				     const int *psi,
				     const double *nu,
				     double c2,
				     const double *r,
				     const double *sigma2,
				     const double *phi,
				     const double *tau2R,
				     const double *b,
				     const vector<vector<int> > &neighbour,
				     double alpha,
				     double beta) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    int g = (int) (ran.Unif01() * G);

    int nOn = 0;
    int q;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      nOn += delta[kqg];
    }
    if (nOn != 0 && nOn != Q) {
      cout << "Error found in function \"updateDeltaDDelta_MRF2_onedelta\":" <<
	endl;
      cout << "All delta's for any gene must be equal." << endl;
      cout << "For gene \"" << g << "\" this requirement is not fulfilled." <<
	endl;
      cout << "Aborting." << endl;
      exit(-1);
    }

    int kqg = qg2index(0,g,Q,G);
    int oldDelta = delta[kqg];
    int newDelta = 1 - oldDelta;

    //
    // propose new values for Delta from full conditionals
    //

    double *newDDelta = (double *) calloc(Q * G,sizeof(double));

    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = newDelta;
    }
    pot -= DeltaGibbs(g,newDDelta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }
    pot += DeltaGibbs(g,Delta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);


    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }
    pot -= potentialDelta_MRF2_onedelta(Q,G,delta,neighbour,alpha,beta);
    pot -= potentialDDeltag(g,Q,G,delta,Delta,c2,b,r,tau2R,sigma2);
    pot -= potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = newDelta;
    }
    pot += potentialDelta_MRF2_onedelta(Q,G,delta,neighbour,alpha,beta);
    pot += potentialDDeltag(g,Q,G,delta,newDDelta,c2,b,r,tau2R,sigma2);
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }

    if (ran.Unif01() <= exp(- pot)) {
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	delta[kqg] = newDelta;
	if (newDelta == 1)
	  Delta[kqg] = newDDelta[kqg];
      }

      (*nAccept)++;
    }

    free(newDDelta);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateDeltaDDelta_MRF2(unsigned int *seed,
			    int nTry,
			    int *nAccept,
			    int *delta,
			    double *Delta,
			    int Q,
			    int G,
			    const int *S,
			    const double *x,
			    const int *psi,
			    const double *nu,
			    double c2,
			    const double *r,
			    const double *sigma2,
			    const double *phi,
			    const double *tau2R,
			    const double *b,
			    const vector<vector<int> > &neighbour,
			    double alpha,
			    double beta,
			    double betag) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    int q = (int) (ran.Unif01() * Q);
    int g = (int) (ran.Unif01() * G);

    int kqg = qg2index(q,g,Q,G);
    int oldDelta = delta[kqg];
    int newDelta = 1 - oldDelta;

    //
    // propose new values for Delta from full conditionals
    //

    double *newDDelta = (double *) calloc(Q * G,sizeof(double));

    delta[kqg] = newDelta;
    pot -= DeltaGibbs(g,newDDelta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);

    delta[kqg] = oldDelta;
    pot += DeltaGibbs(g,Delta,Q,G,S,c2,tau2R,b,r,sigma2,
		      phi,psi,x,delta,nu,ran,1);


    delta[kqg] = oldDelta;
    pot -= potentialDelta_MRF2(Q,G,delta,neighbour,alpha,beta,betag);
    pot -= potentialDDeltag(g,Q,G,delta,Delta,c2,b,r,tau2R,sigma2);
    pot -= potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    delta[kqg] = newDelta;
    pot += potentialDelta_MRF2(Q,G,delta,neighbour,alpha,beta,betag);
    pot += potentialDDeltag(g,Q,G,delta,newDDelta,c2,b,r,tau2R,sigma2);
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);
    delta[kqg] = oldDelta;


    if (ran.Unif01() <= exp(- pot)) {
      delta[kqg] = newDelta;
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	if (delta[kqg] == 1)
	  Delta[kqg] = newDDelta[kqg];
      }

      (*nAccept)++;
    }

    free(newDDelta);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateEta0Omega0Kappa_MRF1_onedelta(unsigned int *seed,
					 int nTry,
					 int *nAccept,
					 double epsilonEta0,
					 double epsilonOmega0,
					 double epsilonKappa,
					 double *eta0,
					 double *omega0,
					 double *kappa,
					 int Q,
					 int G,
					 const int *delta,
					 const vector<vector<int> > &neighbour,
					 double alphaEta,
					 double betaEta,
					 double pOmega0,
					 double lambdaOmega,
					 double lambdaKappa) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    //
    // propose new value(s)
    //

    double oldEta0 = *eta0;
    double newEta0 = *eta0;
    if (epsilonEta0 > 0.0) {
      double upper = oldEta0 + epsilonEta0;
      double lower = oldEta0 - epsilonEta0;
      if (upper > 1.0) upper = 1.0;
      if (lower < 0.0) lower = 0.0;
      newEta0 = lower + (upper - lower) * ran.Unif01();
      pot -= - log(1.0 / (upper - lower));

      upper = newEta0 + epsilonEta0;
      lower = newEta0 - epsilonEta0;
      if (upper > 1.0) upper = 1.0;
      if (lower < 0.0) lower = 0.0;
      pot += - log(1.0 / (upper - lower));
    }

    double oldOmega0 = *omega0;
    double newOmega0 = *omega0;
    if (epsilonOmega0 > 0.0) {
      if (oldOmega0 == 0.0) {
	double upper = epsilonOmega0;
	double lower = 0.0;
	newOmega0 = lower + (upper - lower) * ran.Unif01();
	pot -= - log(1.0 / (upper - lower));

	double prob0 = - (newOmega0 - epsilonOmega0);
	pot += - log(prob0);
      }
      else {
	double prob0 = - (oldOmega0 - epsilonOmega0);
	if (prob0 < 0.0) prob0 = 0.0;
	double upper = oldOmega0 + epsilonOmega0;
	double lower = oldOmega0 - epsilonOmega0;
	if (lower < 0.0) lower = 0.0;
	double u = ran.Unif01();
	if (u < prob0) {
	  newOmega0 = 0.0;
	  pot -= - log(prob0);
	  double upper = epsilonOmega0;
	  double lower = 0.0;
	  pot += - log(1.0 / (upper - lower));
	}
	else {
	  newOmega0 = lower + (upper - lower) * ran.Unif01();
	  pot -= - log(1.0 - prob0);
	  pot -= - log(1.0 / (upper - lower));

	  double prob0 = - (newOmega0 - epsilonOmega0);
	  if (prob0 < 0.0) prob0 = 0.0;
	  double upper = newOmega0 + epsilonOmega0;
	  double lower = newOmega0 - epsilonOmega0;
	  if (lower < 0.0) lower = 0.0;
	  pot += - log(1.0 - prob0);
	  pot += - log(1.0 / (upper - lower));
	}
      }
    }

    double oldKappa = *kappa;
    double newKappa = *kappa;
    if (epsilonKappa > 0.0) {
      double upper = oldKappa + epsilonKappa;
      double lower = oldKappa - epsilonKappa;
      if (lower < 0.0) lower = 0.0;
      newKappa = lower + (upper - lower) * ran.Unif01();
      pot -= - log(1.0 / (upper - lower));

      upper = newKappa + epsilonKappa;
      lower = newKappa - epsilonKappa;
      if (lower < 0.0) lower = 0.0;
      pot += - log(1.0 / (upper - lower));
    }


    //    cout << "eta0: " << newEta0 << ", omega0: " << newOmega0 << ", kappa: " << newKappa << endl;
    int *dd = (int *) calloc(G,sizeof(int));
    vector<double> potZero(G,0.0);
    unsigned int dummy = 1;
    unsigned int seedPerfect = ran.ChangeSeed(dummy);
    perfectMRF1_onedelta(dd,G,neighbour,potZero,potZero,newEta0,newOmega0,
			 newKappa,&seedPerfect,1);
    ran.ChangeSeed(seedPerfect);
    int *deltaTemp = (int *) calloc(Q * G,sizeof(int));
    int q,g;
    for (g = 0; g < G; g++)
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	deltaTemp[kqg] = dd[g];
      }


    pot -= potentialEta0(oldEta0,alphaEta,betaEta);
    pot -= potentialOmega0(oldOmega0,pOmega0,lambdaOmega);
    pot -= potentialKappa(oldKappa,lambdaKappa);
    pot -= potentialDelta_MRF1_onedelta(Q,G,delta,neighbour,oldEta0,
					oldOmega0,oldKappa);
    pot -= potentialDelta_MRF1_onedelta(Q,G,deltaTemp,neighbour,newEta0,
					newOmega0,newKappa);

    pot += potentialEta0(newEta0,alphaEta,betaEta);
    pot += potentialOmega0(newOmega0,pOmega0,lambdaOmega);
    pot += potentialKappa(newKappa,lambdaKappa);
    pot += potentialDelta_MRF1_onedelta(Q,G,delta,neighbour,newEta0,
					newOmega0,newKappa);
    pot += potentialDelta_MRF1_onedelta(Q,G,deltaTemp,neighbour,oldEta0,
					oldOmega0,oldKappa);

    free(dd);
    free(deltaTemp);

    if (ran.Unif01() < exp(- pot)) {
      *eta0 = newEta0;
      *omega0 = newOmega0;
      *kappa = newKappa;

      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateAlphaBeta_MRF2_onedelta(unsigned int *seed,
				   int nTry,
				   int *nAccept,
				   double epsilonAlpha,
				   double epsilonBeta,
				   double *alpha,
				   double *beta,
				   int Q,
				   int G,
				   const int *delta,
				   const vector<vector<int> > &neighbour) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    //
    // propose new value(s)
    //

    double oldAlpha = *alpha;
    double newAlpha = *alpha;
    if (epsilonAlpha > 0.0) {
      newAlpha += epsilonAlpha * ran.Norm01();
    }

    double oldBeta = *beta;
    double newBeta = *beta;
    if (epsilonBeta > 0.0) {
      newBeta += epsilonBeta * ran.Norm01();
      if (newBeta < 0.0) return;
    }



    //    cout << "alpha: " << newAlpha << ", beta: " << newBeta << endl;
    int *dd = (int *) calloc(G,sizeof(int));
    vector<double> potZero(G,0.0);
    unsigned int dummy = 1;
    unsigned int seedPerfect = ran.ChangeSeed(dummy);
    perfectMRF2_onedelta(dd,G,neighbour,potZero,potZero,newAlpha,newBeta,
			 &seedPerfect,1);
    ran.ChangeSeed(seedPerfect);
    int *deltaTemp = (int *) calloc(Q * G,sizeof(int));
    int q,g;
    for (g = 0; g < G; g++)
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	deltaTemp[kqg] = dd[g];
      }


    pot -= potentialAlpha();
    pot -= potentialBeta();
    pot -= potentialDelta_MRF2_onedelta(Q,G,delta,neighbour,oldAlpha,oldBeta);
    pot -= potentialDelta_MRF2_onedelta(Q,G,deltaTemp,neighbour,newAlpha,newBeta);

    pot += potentialAlpha();
    pot += potentialBeta();
    pot += potentialDelta_MRF2_onedelta(Q,G,delta,neighbour,newAlpha,newBeta);
    pot += potentialDelta_MRF2_onedelta(Q,G,deltaTemp,neighbour,oldAlpha,oldBeta);

    free(dd);
    free(deltaTemp);

    if (ran.Unif01() < exp(- pot)) {
      *alpha = newAlpha;
      *beta = newBeta;

      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateAlphaBetaBetag_MRF2(unsigned int *seed,
			       int nTry,
			       int *nAccept,
			       double epsilonAlpha,
			       double epsilonBeta,
			       double epsilonBetag,
			       double *alpha,
			       double *beta,
			       double *betag,
			       int Q,
			       int G,
			       const int *delta,
			       const vector<vector<int> > &neighbour) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    //
    // propose new value(s)
    //

    double oldAlpha = *alpha;
    double newAlpha = *alpha;
    if (epsilonAlpha > 0.0) {
      newAlpha += epsilonAlpha * ran.Norm01();
    }

    double oldBeta = *beta;
    double newBeta = *beta;
    if (epsilonBeta > 0.0) {
      newBeta += epsilonBeta * ran.Norm01();
      if (newBeta < 0.0) return;
    }


    double oldBetag = *betag;
    double newBetag = *betag;
    if (epsilonBetag > 0.0) {
      newBetag += epsilonBetag * ran.Norm01();
      if (newBetag < 0.0) return;
    }


    //    cout << "alpha: " << newAlpha << ", beta: " << newBeta << ", betag: " <<
    //      newBetag << endl;
    int *dd = (int *) calloc(Q * G,sizeof(int));
    vector<double> potZero(Q * G,0.0);
    unsigned int dummy = 1;
    unsigned int seedPerfect = ran.ChangeSeed(dummy);
    perfectMRF2(dd,Q,G,neighbour,potZero,potZero,newAlpha,newBeta,
		newBetag,&seedPerfect,1);
    ran.ChangeSeed(seedPerfect);

    pot -= potentialAlpha();
    pot -= potentialBeta();
    pot -= potentialBetag();
    pot -= potentialDelta_MRF2(Q,G,delta,neighbour,oldAlpha,oldBeta,oldBetag);
    pot -= potentialDelta_MRF2(Q,G,dd,neighbour,newAlpha,newBeta,newBetag);

    pot += potentialAlpha();
    pot += potentialBeta();
    pot += potentialBetag();
    pot += potentialDelta_MRF2(Q,G,delta,neighbour,newAlpha,newBeta,newBetag);
    pot += potentialDelta_MRF2(Q,G,dd,neighbour,oldAlpha,oldBeta,oldBetag);

    free(dd);

    if (ran.Unif01() < exp(- pot)) {
      *alpha = newAlpha;
      *beta = newBeta;
      *betag = newBetag;

      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}




void updateOmega_HyperInverseWishart(unsigned int *seed,
				     int *nAccept,
				     vector<vector<vector<double> > > &Omega,
				     int Q,
				     int G,
				     const double *Delta,
				     const double *r,
				     const double *sigma2,
				     const double *tau2R,
				     const double *b,
				     double df,
				     const vector<vector<vector<double> > > &D,
				     const vector<int> &oldClique,
				     const vector<vector<int> > &oldComponents) {
  Random ran(*seed);

  //
  // propose new values from full conditional
  //

  vector<vector<vector<double> > > OmegaOld(Omega);

  double pot = - OmegaGibbs(df,D,oldClique,oldComponents,Q,G,Delta,
			  r,sigma2,tau2R,b,Omega,ran,1);

  /*
  pot += OmegaGibbs(df,D,oldClique,oldComponents,Q,G,Delta,
		    r,sigma2,tau2R,b,OmegaOld,ran,0);

  pot += potentialOmega_HyperInverseWishart(Omega,D,df,oldClique,oldComponents);
  pot += potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,Q,G,Omega,oldClique,oldComponents);

  pot -= potentialOmega_HyperInverseWishart(OmegaOld,D,df,oldClique,oldComponents);
  pot -= potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,r,Q,G,
						 OmegaOld,oldClique,oldComponents);

  cout << "UpdateOmega::Potential: " << pot << endl;
  */

  (*nAccept)++;

  *seed = ran.ChangeSeed(*seed);

  return;
}



void updateDDeltaStar_HyperInverseWishart(unsigned int *seed,
					  int *nAccept,
					  double *Delta,
					  int Q,
					  int G,
					  const int *S,
					  const double *x,
					  const int *psi,
					  const double *nu,
					  const int *delta,
					  const double *r,
					  const double *sigma2,
					  const double *phi,
					  const double *tau2R,
					  const double *b,
					  const vector<vector<vector<double> > > &Omega,
					  const vector<int> &oldClique,
					  const vector<vector<int> > &oldComponents) {

  Random ran(*seed);

  //
  // propose new values from full conditional
  //

  double *oldValues = Delta;
  double *newValues = (double *) calloc(Q * G,sizeof(double));

  double pot = DeltaStarGibbs(oldClique,oldComponents,Q,G,S,newValues,
			      r,sigma2,phi,tau2R,b,nu,delta,psi,x,Omega,ran,1);

  /*
  double pot2 = DeltaStarGibbs(oldClique,oldComponents,Q,G,S,oldValues,
			       r,sigma2,phi,tau2R,b,nu,delta,psi,x,Omega,ran,0);

  double pot3 = potentialDDeltaStar_HyperInverseWishart(newValues,b,sigma2,tau2R,r,Q,G,Omega,oldClique,oldComponents);
  double pot4 = potentialX(Q,G,S,x,psi,nu,delta,newValues,sigma2,phi);

  double pot5 = potentialDDeltaStar_HyperInverseWishart(oldValues,b,sigma2,tau2R,r,Q,G,Omega,oldClique,oldComponents);
  double pot6 = potentialX(Q,G,S,x,psi,nu,delta,oldValues,sigma2,phi);

  pot = - pot + pot2 + pot3 + pot4 - pot5 - pot6;

  cout << "UpdateDDeltaStar::Potential: " << pot << endl;
  */

  int k;
  for (k = 0; k < Q * G; k++)
    oldValues[k] = newValues[k];
  free(newValues);


  (*nAccept)++;

  *seed = ran.ChangeSeed(*seed);

  return;
}




void updateTau2RDDeltaStar_HyperInverseWishart(unsigned int *seed,
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
					      const double *r,
					      const double *sigma2,
					      const double *phi,
					      const double *b,
					      const vector<vector<vector<double> > > &Omega,
					      const vector<int> &oldClique,
					      const vector<vector<int> > &oldComponents) {
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
      // propose new values for DeltaStar from full conditionals
      //

      double *newDDelta = (double *) calloc(Q * G,sizeof(double));

      pot -= DeltaStarGibbs(oldClique,oldComponents,Q,G,S,newDDelta,r,
			    sigma2,phi,newValues,b,nu,delta,psi,x,
			    Omega,ran,1);

      pot += DeltaStarGibbs(oldClique,oldComponents,Q,G,S,Delta,r,
			    sigma2,phi,oldValues,b,nu,delta,psi,x,
			    Omega,ran,1);

      pot -= potentialTau2R();
      pot -= potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,oldValues,
						     r,Q,G,Omega,oldClique,
						     oldComponents);
      pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

      pot += potentialTau2R();
      pot += potentialDDeltaStar_HyperInverseWishart(newDDelta,b,sigma2,newValues,
						     r,Q,G,Omega,oldClique,
						     oldComponents);
      pot += potentialX(Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);


      if (ran.Unif01() <= exp(- pot)) {
	int q;
	for (q = 0; q < Q; q++)
	  tau2R[q] = newValues[q];

	int k;
	for (k = 0; k < Q * G; k++)
	  Delta[k] = newDDelta[k];

	(*nAccept)++;
      }

      free(newDDelta);
      free(oldValues);
      free(newValues);
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}






void updateBDDeltaStar_HyperInverseWishart(unsigned int *seed,
					   int nTry,
					   int *nAccept,
					   double epsilon,
					   double *b,
					   double *Delta,
					   int Q,
					   int G,
					   const int *S,
					   const double *x,
					   const int *psi,
					   const double *nu,
					   const int *delta,
					   const double *r,
					   const double *sigma2,
					   const double *phi,
					   const double *tau2R,
					   double pB0,
					   double pB1,
					   double alphaB,
					   double betaB,
					   const vector<vector<vector<double> > > &Omega,
					   const vector<int> &oldClique,
					   const vector<vector<int> > &oldComponents) {
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

    //
    // propose new values for DeltaStar from full conditionals
    //

    double *newDDelta = (double *) calloc(Q * G,sizeof(double));

    b[q] = newValue;
    pot -= DeltaStarGibbs(oldClique,oldComponents,Q,G,S,newDDelta,r,
			  sigma2,phi,tau2R,b,nu,delta,psi,x,
			  Omega,ran,1);

    b[q] = oldValue;
    pot += DeltaStarGibbs(oldClique,oldComponents,Q,G,S,Delta,r,
			  sigma2,phi,tau2R,b,nu,delta,psi,x,
			  Omega,ran,1);

    pot -= potentialB(Q,b,pB0,pB1,alphaB,betaB);
    pot -= potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,
						   r,Q,G,Omega,oldClique,
						   oldComponents);
    pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    b[q] = newValue;
    pot += potentialB(Q,b,pB0,pB1,alphaB,betaB);
    pot += potentialDDeltaStar_HyperInverseWishart(newDDelta,b,sigma2,tau2R,
						   r,Q,G,Omega,oldClique,
						   oldComponents);
    pot += potentialX(Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);
    b[q] = oldValue;

    if (ran.Unif01() <= exp(- pot)) {
      b[q] = newValue;

      int k;
      for (k = 0; k < Q * G; k++)
	Delta[k] = newDDelta[k];

      (*nAccept)++;
    }

    free(newDDelta);
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}





void updateRDDeltaStar_HyperInverseWishart(unsigned int *seed,
					   int nTry,
					   int *nAccept,
					   double epsilon,
					   double *r,
					   double *Delta,
					   int Q,
					   int G,
					   const int *S,
					   const double *x,
					   const int *psi,
					   const double *nu,
					   const int *delta,
					   const double *sigma2,
					   const double *phi,
					   const double *tau2R,
					   const double *b,
					   double nuR,
					   const vector<vector<vector<double> > > &Omega,
					   const vector<int> &oldClique,
					   const vector<vector<int> > &oldComponents) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    vector<vector<double> > oldR;
    vector<vector<double> > newR;
    oldR.resize(Q);
    newR.resize(Q);
    int p,q;
    for (p = 0; p < Q; p++) {
      oldR[p].resize(Q);
      newR[p].resize(Q);
    }
    for (p = 0; p < Q; p++) {
      oldR[p][p] = 1.0;
      newR[p][p] = 1.0;
      for (q = p + 1; q < Q; q++) {
	int kqq = qq2index(p,q,Q);
	oldR[p][q] = r[kqq];
	newR[p][q] = r[kqq];
	oldR[q][p] = r[kqq];
	newR[q][p] = r[kqq];
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

    newR[pp][qq] = oldR[pp][qq] * exp(u) /
      (1 - oldR[pp][qq] + oldR[pp][qq] * exp(u));
    newR[qq][pp] = newR[pp][qq];
    double *newRVector = (double *) calloc(Q * (Q - 1) / 2,sizeof(double));
    for (p = 0; p < Q; p++)
      for (q = p + 1; q < Q; q++) {
	int kqq = qq2index(p,q,Q);
	newRVector[kqq] = newR[p][q];
      }

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
	// propose new values for DeltaStar from full conditionals
	//

	double *newDDelta = (double *) calloc(Q * G,sizeof(double));

	pot -= DeltaStarGibbs(oldClique,oldComponents,Q,G,S,newDDelta,
			      newRVector,sigma2,phi,tau2R,b,nu,
			      delta,psi,x,Omega,ran,1);

	pot += DeltaStarGibbs(oldClique,oldComponents,Q,G,S,Delta,
			      r,sigma2,phi,tau2R,b,nu,
			      delta,psi,x,Omega,ran,1);

	pot -= potentialR(Q,r,nuR);
	pot -= potentialDDeltaStar_HyperInverseWishart(Delta,b,sigma2,tau2R,
						       r,Q,G,Omega,oldClique,
						       oldComponents);
	pot -= potentialX(Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

	pot += potentialR(Q,newRVector,nuR);
	pot += potentialDDeltaStar_HyperInverseWishart(newDDelta,b,sigma2,tau2R,
						       newRVector,Q,G,Omega,oldClique,
						       oldComponents);
	pot += potentialX(Q,G,S,x,psi,nu,delta,newDDelta,sigma2,phi);


	if (ran.Unif01() <= exp(- pot)) {
	  int kk;
	  for (kk = 0; kk < Q * (Q - 1) / 2; kk++)
	    r[kk] = newRVector[kk];

	  for (kk = 0; kk < Q * G; kk++)
	    Delta[kk] = newDDelta[kk];

	  (*nAccept)++;
	}

      free(newDDelta);
      }
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}




void updateDelta_HyperInverseWishart(unsigned int *seed,
				     int nTry,
				     int *nAccept,
				     int *delta,
				     int Q,
				     int G,
				     const int *S,
				     const double *x,
				     const int *psi,
				     const double *nu,
				     const double *Delta,
				     const double *r,
				     const double *sigma2,
				     const double *phi,
				     const double *xi,
				     const double *b) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    int q = (int) (ran.Unif01() * Q);
    int g = (int) (ran.Unif01() * G);

    int kqg = qg2index(q,g,Q,G);
    int oldDelta = delta[kqg];
    int newDelta = 1 - oldDelta;

    pot -= potentialDeltag(g,Q,G,delta,xi);
    pot -= potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    delta[kqg] = newDelta;
    pot += potentialDeltag(g,Q,G,delta,xi);
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    delta[kqg] = oldDelta;

    if (ran.Unif01() <= exp(- pot)) {
      delta[kqg] = newDelta;

      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}




void updateDelta_HyperInverseWishart_onedelta(unsigned int *seed,
					      int nTry,
					      int *nAccept,
					      int *delta,
					      int Q,
					      int G,
					      const int *S,
					      const double *x,
					      const int *psi,
					      const double *nu,
					      const double *Delta,
					      const double *r,
					      const double *sigma2,
					      const double *phi,
					      const double *xi,
					      const double *b) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    int g = (int) (ran.Unif01() * G);

    int nOn = 0;
    int q;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      nOn += delta[kqg];
    }
    if (nOn != 0 && nOn != Q) {
      cout << "Error found in function \"updateDeltaDDelta_onedelta\":" << endl;
      cout << "All delta's for any gene must be equal." << endl;
      cout << "For gene \"" << g << "\" this requirement is not fulfilled." <<
	endl;
      cout << "Aborting." << endl;
      exit(-1);
    }

    int kqg = qg2index(0,g,Q,G);
    int oldDelta = delta[kqg];
    int newDelta = 1 - oldDelta;

    pot -= potentialDeltag_onedelta(g,Q,G,delta,xi);
    pot -= potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = newDelta;
    }
    pot += potentialDeltag_onedelta(g,Q,G,delta,xi);
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }

    if (ran.Unif01() <= exp(- pot)) {
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	delta[kqg] = newDelta;
      }

      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}




void updateDelta_HyperInverseWishart_MRF2(unsigned int *seed,
					  int nTry,
					  int *nAccept,
					  int *delta,
					  int Q,
					  int G,
					  const int *S,
					  const double *x,
					  const int *psi,
					  const double *nu,
					  const double *Delta,
					  const double *r,
					  const double *sigma2,
					  const double *phi,
					  const double *b,
					  const vector<vector<int> > &neighbour,
					  double alpha,
					  double beta,
					  double betag) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    int q = (int) (ran.Unif01() * Q);
    int g = (int) (ran.Unif01() * G);

    int kqg = qg2index(q,g,Q,G);
    int oldDelta = delta[kqg];
    int newDelta = 1 - oldDelta;

    pot -= potentialDelta_MRF2(Q,G,delta,neighbour,alpha,beta,betag);
    pot -= potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    delta[kqg] = newDelta;
    pot += potentialDelta_MRF2(Q,G,delta,neighbour,alpha,beta,betag);
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    delta[kqg] = oldDelta;

    if (ran.Unif01() <= exp(- pot)) {
      delta[kqg] = newDelta;

      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}



void updateDelta_HyperInverseWishart_MRF2_onedelta(unsigned int *seed,
						   int nTry,
						   int *nAccept,
						   int *delta,
						   int Q,
						   int G,
						   const int *S,
						   const double *x,
						   const int *psi,
						   const double *nu,
						   const double *Delta,
						   const double *r,
						   const double *sigma2,
						   const double *phi,
						   const double *b,
						   const vector<vector<int> > &neighbour,
						   double alpha,
						   double beta) {
  Random ran(*seed);

  int k;
  for (k = 0; k < nTry; k++) {
    double pot = 0.0;

    int g = (int) (ran.Unif01() * G);

    int nOn = 0;
    int q;
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      nOn += delta[kqg];
    }
    if (nOn != 0 && nOn != Q) {
      cout << "Error found in function \"updateDelta_HyperInverseWishart_MRF2_onedelta\":" << endl;
      cout << "All delta's for any gene must be equal." << endl;
      cout << "For gene \"" << g << "\" this requirement is not fulfilled." <<
	endl;
      cout << "Aborting." << endl;
      exit(-1);
    }

    int kqg = qg2index(0,g,Q,G);
    int oldDelta = delta[kqg];
    int newDelta = 1 - oldDelta;


    pot -= potentialDelta_MRF2_onedelta(Q,G,delta,neighbour,alpha,beta);
    pot -= potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);

    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = newDelta;
    }
    pot += potentialDelta_MRF2_onedelta(Q,G,delta,neighbour,alpha,beta);
    pot += potentialXg(g,Q,G,S,x,psi,nu,delta,Delta,sigma2,phi);
    for (q = 0; q < Q; q++) {
      int kqg = qg2index(q,g,Q,G);
      delta[kqg] = oldDelta;
    }

    if (ran.Unif01() <= exp(- pot)) {
      for (q = 0; q < Q; q++) {
	int kqg = qg2index(q,g,Q,G);
	delta[kqg] = newDelta;
      }

      (*nAccept)++;
    }
  }

  *seed = ran.ChangeSeed(*seed);

  return;
}




