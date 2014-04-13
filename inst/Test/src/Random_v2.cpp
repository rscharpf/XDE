#include <iostream>
#include <limits.h>
#include <stdio.h>

#include "Matrix_v2.h"
#include "Random_v2.h"
#include "Cholesky.h"

Random::Random(unsigned int seed)
{
  unsigned int temp = UINT_MAX - 1;
  temp /= 2;
  temp += 1;
  temp /= FACTOR4;

  if (temp < FACTOR1*FACTOR2*FACTOR3)
    {
      cout << "Error: Implemented random number generator requires UINT_MAX to be at least " << 
	FACTOR1 << "*" << FACTOR2 << "*" << FACTOR3 << "*" << FACTOR4 << "*2-1. Your machine has UNIT_MAX=" << UINT_MAX << ". Aborting!\n";
      exit(-1);
    }

  modulus = FACTOR1;
  modulus *= FACTOR2;
  modulus *= FACTOR3;
  modulus *= FACTOR4;

  seedValue = seed;
  
  haveNorm01 = 0;

  return;
}



Random::~Random(void)
{
  return;
}




unsigned int Random::ChangeSeed(unsigned int seed)
{
  unsigned int old = seedValue;
  seedValue = seed;

  return old;
}



double Random::Unif01(void)
{
  seedValue = MULTIPLIER * seedValue + SHIFT;
  if (seedValue == 0) seedValue = MULTIPLIER * seedValue + SHIFT;
      
  if (seedValue > modulus*2-1)
    {
      double x = ((double) (seedValue - 1)) / ((double) modulus);
      x /= 2.0;
      int nn = (int) x;
      seedValue = seedValue - nn * modulus*2;
    }
  double r = ((double) seedValue) / ((double) modulus);
  r /= 2.0;
  
  return r;
}



double Random::Norm01(void)
{
  if (haveNorm01 == 1)
    {
      haveNorm01 = 0;
      return norm;
    }
  else
    {
      double u1,u2;
      double x,y;
      
      u1 = Unif01();
      u2 = Unif01();
      
      x = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
      y = sqrt(-2.0 * log(u1)) * sin(2.0 * PI * u2);
      
      haveNorm01 = 1;
      norm = x;
      
      return y;
    }
}




double Random::Exponential(double lambda)
{
  double u = Unif01();

  double v = - log(u) / lambda;

  return v;
}




int Random::Poisson(double lambda)
{
  //
  // Count the number of events before 1 in a
  // Poisson process with intensity lambda.
  // Note: this is very inefficient if lambda is large!!!
  //

  int v = 0;
  double sum = Exponential(lambda);
  while (sum <= 1.0)
    {
      v++;
      sum += Exponential(lambda);
    }

  return v;
}




int Random::Binomial(int n,double p)
{
  //
  // Count the number of successes.
  // Note: this is very inefficient if n is large
  //
  
  int v = 0;
  int k;
  for (k = 0; k < n; k++)
    v += (Unif01() <= p);

  return v;
}



int Random::Discrete(const vector<double> &prob)
{
  double sum = 0.0;
  int i;
  for (i = 0; i < prob.size(); i++)
    sum += prob[i];

  double u = sum * Unif01();

  int v;
  sum = prob[0];
  for (v = 0; u > sum; v++) sum += prob[v+1];

  return v;
}



//
// FUNCTION: Random::Gamma
//
// PURPOSE: Draw a sample from a gamma distribution (by rejection sampling)
// Algorithm is taken from Ripley (1987). Stochastic simulation.
// [Parameters: p = alpha, lambda = 1/beta]
//
double Random::Gamma(double p,double lambda)
{
  double x = 0.0;

  double alpha = p;
  if (alpha == 1.0) // this is an exponential distribution
    x = Exponential(1.0);
  else if (alpha < 1.0)
    {
      int accepted = 0;
      while (accepted == 0)
	{
	  double u0 = Unif01();
	  double u1 = Unif01();
	  if (u0 > exp(1.0)/(alpha + exp(1.0))) 
	    {
	      x = - log((alpha + exp(1.0))*(1-u0)/(alpha*exp(1.0)));
	      if (log(u1) > (alpha - 1.0)*log(x))
		accepted = 0;
	      else
		accepted = 1;
	    }
	  else
	    {
	      x = (1.0 / alpha) * log((alpha + exp(1.0))*u0/exp(1.0));
	      x = exp(x);
	      if (u1 > exp(-x))
		accepted = 0;
	      else
		accepted = 1;
	    }
	}
    }
  else // alpha > 1.0
    {
      double c = 2.0 / (alpha - 1);
      double c1 = alpha - 1;
      double c2 = (alpha - 1.0/(6.0*alpha)) / c1;
      double c3 = 2.0 / c1;
      double c4 = c + 2.0;
      double c5 = 1.0 / sqrt(alpha);

      int accepted1 = 0;
      while (accepted1 == 0)
	{
	  double u1 = 0.0;
	  double u2 = 0.0;

	  int accepted2 = 0;
	  while (accepted2 == 0)
	    {
	      u1 = Unif01();
	      u2 = Unif01();
	      if (alpha > 2.5) u1 = u2 + c5 * (1.0 - 1.86 * u1);
	      if (u1 > 0.0 && u1 < 1.0)
		accepted2 = 1;
	      else
		accepted2 = 0;
	    }
	  double w = c2 * u2 / u1;
	  if (c3 * u1 + w + 1.0 / w <= c4)
	    accepted1 = 1;
	  else if (c3 * log(u1) - log(w) + w >= 1.0)
	    accepted1 = 0;
	  else
	    accepted1 = 1;
	  x = c1 * w;
	}
    }

  x /= lambda;

  return x;

    /*
  double prob;
  double u,am,e,s,v1,v2,xx,yy;
  double fp,temp;
  int ip,j;
  double x = 0.0;

  fp = modf(p,&temp);
  ip = (int) (temp + 0.1);
  if (fp > 0.0)
    {
      do
	{
	  u = Unif01();
	  if (u <= exp(-1.0) / (exp(-1.0) + 1.0 / fp))
	    {
	      x = 1.0 - log(Unif01());
	      prob = exp((fp - 1.0) * log(x));
	    }
	  else
	    {
	      x = exp( log(Unif01()) / fp);
	      prob = exp( - x);
	    }
	}
      while (Unif01() > prob || x == 0.0);
    }
  
  if (ip > 0 && ip < 6)
    {
      xx = 1.0;
      for (j = 1; j <= ip; j++) xx *= Unif01();
      xx = - log(xx);
      
      x += xx;
    }
  else if (ip >= 6)
    {
      do
	{
	  do
	    {
	      do
		{
		  v1 = 2.0 * Unif01() - 1.0;
		  v2 = 2.0 * Unif01() - 1.0;
		}
	      while (v1 * v1 + v2 * v2 > 1.0);
	      yy = v2 / v1;
	      am = ip - 1.0;
	      s = sqrt(2.0 * am + 1.0);
	      xx = s * yy + am;
	    }
	  while (xx <= 0.0);
	  e = (1.0 + yy * yy) * exp( am * log(xx / am) - s * yy);
	}
      while (Unif01() > e);
      
      x += xx;
    }
  
  x /= lambda;

  return x;
  */
}





double Random::ChiSquared(double nu)
{
  double p = nu / 2.0;
  double lambda = 0.5;

  return Gamma(p,lambda);
}




double Random::InverseGamma(double p,double lambda)
{
  double x = Gamma(p,lambda);

  x = 1.0 / x;

  return x;
}



//
// FUNCTION: Random::Beta
//
// PURPOSE: Draw a sample from a beta distribution as X1/(X1+X2) where
//          X1 and X2 are independent gamma distributed variates
//
double Random::Beta(double alpha,double beta)
{
  double lambda = 1.0;
  double x1 = Gamma(alpha,lambda);
  double x2 = Gamma(beta,lambda);
  
  double b = x1 / (x1 + x2);
  
  return b;
}




vector<double> Random::MultiGaussian(const vector<vector<double> > &Sigma,
				     const vector<double> &mean)
{
  int n = mean.size();

  int err = 0;
  Cholesky chol(Sigma,err);
  if (err != 0) {
    cout << "Error in Cholesky!!\n";
    exit(-1);
  }

  int i;
  vector<double> vec(n);
  for (i = 0; i < n; i++)
    vec[i] = Norm01();

  vector<double> z;
  matrixMult(chol.q_L(),vec,z);
  
  vector<double> x(n);
  for (i = 0; i < n; i++)
    x[i] = z[i] + mean[i];

  return x;  
}




vector<vector<double> > Random::WishartAlternativeParam(double nu,const vector<vector<double> > &V)
{
  int err = 0;
  Cholesky chol(V,err);
  if (err != 0) {
    cout << "Error in Cholesky!!\n";
    exit(-1);
  }


  int p = V.size();
  vector<vector<double> > T;
  T.resize(p);
  int i,j;
  for (i = 0; i < p; i++) T[i].resize(p);
  for (i = 0; i < p; i++)
    for (j = 0; j < p; j++)
      T[i][j] = 0.0;
  for (i = 0; i < p; i++)
    T[i][i] = sqrt(ChiSquared(nu-i+1.0));
      
  for (i = 0; i < p; i++)
    for (j = 0; j < i; j++)
      T[i][j] = Norm01();
  
  vector<vector<double> > LT;
  matrixMult(chol.q_L(),T,LT);
  
  vector<vector<double> > SS;
  outerProduct(LT,SS);

  return SS;
}
    




vector<vector<double> > Random::InverseWishartAlternativeParam(double nu,const vector<vector<double> > &V)
{
  vector<vector<double> > Vinv;
  inverse(V,Vinv);

  vector<vector<double> > w(WishartAlternativeParam(nu,Vinv));

  vector<vector<double> > ww;
  inverse(w,ww);

  return ww;
}






vector<vector<double> >Random::StandardWishartAlternativeParam(int dim,double nu)
{
  int p = dim;
  vector<vector<double> > T;
  T.resize(p);
  int i;
  for (i = 0; i < p; i++)
    T[i].resize(p);
  int j;
  for (i = 0; i < p; i++)
    for (j = 0; j < p; j++)
      T[i][j] = 0.0;
  for (i = 0; i < p; i++)
    T[i][i] = sqrt(ChiSquared(nu-i+1.0));
      
  for (i = 0; i < p; i++)
    for (j = 0; j < i; j++)
      T[i][j] = Norm01();
  
  vector<vector<double> > AA;
  outerProduct(T,AA);

  return AA;
}
    






vector<vector<double> > Random::StandardInverseWishartAlternativeParam(int dim,double nu)
{
  vector<vector<double> > w(dim);
  int i;
  for (i = 0; i < dim; i++)
    w[i].resize(dim);
  w = StandardWishart(dim,nu);

  vector<vector<double> > ww;
  inverse(w,ww);

  return w;
}





vector<vector<double> > Random::CorrelationStandardInverseWishartAlternativeParam(int dim,double nu)
{
  vector<vector<double> > w(dim);
  int i;
  for (i = 0; i < dim; i++)
    w[i].resize(dim);
  w = StandardInverseWishart(dim,nu);

  vector<vector<double> > ww(w.size());
  int j;
  for (i = 0; i < w.size(); i++)
    {
      ww[i].resize(w[i].size());
      for (j = 0; j < w[i].size(); j++)
	ww[i][j] = w[i][j] / sqrt(w[i][i] * w[j][j]);
    }

  return ww;
}



vector<vector<vector<double> > > Random::HyperInverseWishart(double df,const vector<vector<vector<double> > > &D,
							     const vector<int> &oldClique,
							     const vector<vector<int> > &oldComponents) {
  /*
  int kk;
  for (kk = 0; kk < oldClique.size(); kk++) 
    cout << "D[" << kk << "].size() = " << D[kk].size() << endl;
  cout << endl;

  for (kk = 0; kk < oldClique.size(); kk++) 
    cout << "oldClique[" << kk << "]=" << oldClique[kk] << endl;
  cout << endl;

  int ll;
  for (kk = 0; kk < oldComponents.size(); kk++)
    for (ll = 0; ll < oldComponents[kk].size(); ll++)
      cout << "oldComponents[" << kk << "][" << ll << "]=" << oldComponents[kk][ll] << endl;
  */



  vector<vector<vector<double> > > Sigma;
  Sigma.resize(D.size());
  int k;
  for (k = 0; k < D.size(); k++) {
    Sigma[k].resize(D[k].size());
    int i;
    for (i = 0; i < Sigma[k].size(); i++)
      Sigma[k][i].resize(D[k].size());
  }

  int i,j;
  vector<vector<double> > temp(InverseWishart(df,D[0]));
  for (i = 0; i < temp.size(); i++)
    for (j = 0; j < temp[i].size(); j++)
      Sigma[0][i][j] = temp[i][j];
  
  /*
  cout << "Sigma[0]: " << endl;
  for (i = 0; i < Sigma[0].size(); i++) {
    for (j = 0; j < Sigma[0][i].size(); j++)
      cout << Sigma[0][i][j] << " ";
    cout << endl;
  }
  cout << endl;
  */

  for (k = 1; k < D.size(); k++) {

    if (oldComponents[k].size() > 0) {
      // generate DRR, DSS, DRS and DSR from D[k]
      
      vector<vector<double> > DRR;
      vector<vector<double> > DSS;
      vector<vector<double> > DRS;
      vector<vector<double> > DSR;
      DRR.resize(D[k].size() - oldComponents[k].size());
      DSS.resize(oldComponents[k].size());
      DRS.resize(DRR.size());
      DSR.resize(DSS.size());
      int i,j;
      for (i = 0; i < DRR.size(); i++) DRR[i].resize(DRR.size());
      for (i = 0; i < DSS.size(); i++) DSS[i].resize(DSS.size());
      for (i = 0; i < DRS.size(); i++) DRS[i].resize(DSS.size());
      for (i = 0; i < DSR.size(); i++) DSR[i].resize(DRR.size());
      
      for (i = 0; i < DSS.size(); i++)
	for (j = 0; j < DSS[i].size(); j++)
	  DSS[i][j] = D[k][i][j];
      
      for (i = 0; i < DRR.size(); i++)
	for (j = 0; j < DRR[i].size(); j++)
	  DRR[i][j] = D[k][i + DSS.size()][j + DSS.size()];
      
      for (i = 0; i < DSR.size(); i++)
	for (j = 0; j < DSR[i].size(); j++) {
	  DSR[i][j] = D[k][i][j + DSS.size()];
	  DRS[j][i] = D[k][j + DSS.size()][i];
	}
      
      // pick out SigmaSS, which is already simulated
      
      vector<vector<double> > SigmaSS;
      SigmaSS.resize(oldComponents[k].size());
      for (i = 0; i < SigmaSS.size(); i++) {
	SigmaSS[i].resize(oldComponents[k].size());
	for (j = 0; j < SigmaSS[i].size(); j++)
	  SigmaSS[i][j] = Sigma[oldClique[k]][oldComponents[k][i]][oldComponents[k][j]];
      }
      
      // simulate SigmaRGivenS
      
      vector<vector<double> > DSSInverse;
      inverse(DSS,DSSInverse);
      vector<vector<double> > temp1;
      matrixMult(DRS,DSSInverse,temp1);
      
      vector<vector<double> > temp2;
      matrixMult(temp1,DSR,temp2);
      vector<vector<double> > DRGivenS;
      DRGivenS.resize(DRR.size());
      for (i = 0; i < DRR.size(); i++) {
	DRGivenS[i].resize(DRR.size());
	for (j = 0; j < DRR[i].size(); j++)
	  DRGivenS[i][j] = DRR[i][j] - temp2[i][j];
      }
      
      vector<vector<double> > SigmaRGivenS(InverseWishart(df + DRR.size(),DRGivenS));
      
      // simulate U matrix
      
      vector<vector<double> > mean;
      matrixMult(DRS,DSSInverse,mean);
      /*
      cout << "mean: " << endl;
      for (i = 0; i < mean.size(); i++) {
	for (j = 0; j < mean[i].size(); j++)
	  cout << mean[i][j] << " ";
	cout << endl;
      }
      cout << endl;
      */

      int err = 0;
      Cholesky cholRGivenS(SigmaRGivenS,err);
      if (err != 0) {
	cout << "Error in Cholesky!!\n";
	exit(-1);
      }
      
      vector<vector<double> > LRGivenS(cholRGivenS.q_L());
      
      /*
      cout << "LRGivenS: " << endl;
      for (i = 0; i < LRGivenS.size(); i++) {
	for (j = 0; j < LRGivenS[i].size(); j++)
	  cout << LRGivenS[i][j] << " ";
	cout << endl;
      }
      cout << endl;
      */
      
      err = 0;
      Cholesky cholDSSInverse(DSSInverse,err);
      if (err != 0) {
	cout << "Error in Cholesky!!\n";
	exit(-1);
      }
      vector<vector<double> > LDSSInverse(cholDSSInverse.q_L());
      vector<vector<double> > LDSSInverseT;
      LDSSInverseT.resize(LDSSInverse[0].size());
      for (i = 0; i < LDSSInverseT.size(); i++) {
	LDSSInverseT[i].resize(LDSSInverse.size());
	for (j = 0; j < LDSSInverseT[i].size(); j++)
	  LDSSInverseT[i][j] = LDSSInverse[j][i];
      }
      
      /*
      cout << "LDSSInverseT: " << endl;
      for (i = 0; i < LDSSInverseT.size(); i++) {
	for (j = 0; j < LDSSInverseT[i].size(); j++)
	  cout << LDSSInverseT[i][j] << " ";
	cout << endl;
      }
      cout << endl;
      */

      vector<vector<double> > U;
      U.resize(DRR.size());
      for (i = 0; i < DRR.size(); i++) {
	U[i].resize(DSS.size());
	for (j = 0; j < DSS.size(); j++)
	  U[i][j] = Norm01();
      }
      
      vector<vector<double> > temp3(U);
      matrixMult(LRGivenS,temp3,U);
      
      vector<vector<double> > temp4(U);
      matrixMult(temp4,LDSSInverseT,U);

      for (i = 0; i < mean.size(); i++)
	for (j = 0; j < mean[i].size(); j++)
	  U[i][j] += mean[i][j];
      
      /*
      cout << "U: " << endl;
      for (i = 0; i < U.size(); i++) {
	for (j = 0; j < U[i].size(); j++)
	  cout << U[i][j] << " ";
	cout << endl;
      }
      cout << endl;
      */

      // compute the two remaining parts of Sigma
      
      vector<vector<double> > SigmaRS;
      matrixMult(U,SigmaSS,SigmaRS);
      
      vector<vector<double> > SigmaSR;
      SigmaSR.resize(SigmaRS[0].size());
      for (i = 0; i < SigmaSR.size(); i++) {
	SigmaSR[i].resize(SigmaRS.size());
	for (j = 0; j < SigmaSR[i].size(); j++)
	  SigmaSR[i][j] = SigmaRS[j][i];
      }
      
      vector<vector<double> > SigmaRR(SigmaRGivenS);
      vector<vector<double> > SigmaSSInverse;
      inverse(SigmaSS,SigmaSSInverse);
      vector<vector<double> > temp5;
      matrixMult(SigmaRS,SigmaSSInverse,temp5);
      vector<vector<double> > temp6;
      matrixMult(temp5,SigmaSR,temp6);
      
      for (i = 0; i < SigmaRGivenS.size(); i++)
	for (j = 0; j < SigmaRGivenS[i].size(); j++)
	  SigmaRR[i][j] += temp6[i][j];
      
      // put the four parts into one matrix
      
      for (i = 0; i < SigmaSS.size(); i++)
	for (j = 0; j < SigmaSS[i].size(); j++)
	  Sigma[k][i][j] = SigmaSS[i][j];
      for (i = 0; i < SigmaRR.size(); i++)
	for (j = 0; j < SigmaRR[i].size(); j++)
	  Sigma[k][i + SigmaSS.size()][j + SigmaSS.size()] = SigmaRR[i][j];
      for (i = 0; i < SigmaSS.size(); i++)
	for (j = 0; j < SigmaRR.size(); j++) {
	  Sigma[k][i][j + SigmaSS.size()] = SigmaSR[i][j];
	  Sigma[k][j + SigmaSS.size()][i] = SigmaRS[j][i];
	}    
    }
    else {
      vector<vector<double> > temp7(InverseWishart(df,D[k]));
      for (i = 0; i < temp7.size(); i++)
	for (j = 0; j < temp7[i].size(); j++)
	  Sigma[k][i][j] = temp7[i][j];
    }

    /*    
    cout << "Sigma[" << k << "]: " << endl;
    for (i = 0; i < Sigma[k].size(); i++) {
      for (j = 0; j < Sigma[k][i].size(); j++)
	cout << Sigma[k][i][j] << " ";
      cout << endl;
    }
    cout << endl;
    */
  }

  return Sigma;
}
    


vector<double> Random::GaussianGraphicalModel(const vector<double> &mean,
					      const vector<vector<vector<double> > > &Cov,
					      const vector<int> &oldClique,
					      const vector<vector<int> > &oldComponents) {
  // allocate space for result

  vector<vector<double> > UBlocks;
  UBlocks.resize(Cov.size());
  int k;
  for (k = 0; k < UBlocks[k].size(); k++)
    UBlocks[k].resize(Cov[k].size());

  vector<double> U(mean);
  
  // generate realisation from the first component

  int err = 0;
  Cholesky chol(Cov[0],err);
  if (err != 0) {
    cout << "Error in Cholesky!!\n";
    exit(-1);
  }

  int first = 0;
  vector<double> vec(Cov[0].size(),0.0);
  int i;
  for (i = 0; i < vec.size(); i++)
    vec[i] = Norm01();

  vector<double> z;
  matrixMult(chol.q_L(),vec,z);
  for (i = 0; i < vec.size(); i++) {
    U[i] = vec[i];
    UBlocks[0][i] = vec[i];
  }

  first += z.size();

  // generate realisations from the remaining components

  for (k = 1; k < Cov.size(); k++) {
    
    if (oldComponents[k].size() > 0) {
      // establish covariance submatrices

      vector<vector<double> > CovSS;
      vector<vector<double> > CovRR;
      vector<vector<double> > CovSR;
      vector<vector<double> > CovRS;
      CovSS.resize(oldComponents[k].size());
      CovRR.resize(Cov[k].size() - oldComponents[k].size());
      CovSR.resize(oldComponents[k].size());
      CovRS.resize(Cov[k].size() - oldComponents[k].size());
      for (i = 0; i < CovSS.size(); i++) CovSS[i].resize(CovSS.size());
      for (i = 0; i < CovRR.size(); i++) CovRR[i].resize(CovRR.size());
      for (i = 0; i < CovSR.size(); i++) CovSR[i].resize(CovRS.size());
      for (i = 0; i < CovRS.size(); i++) CovRS[i].resize(CovSR.size());

      int j;
      for (i = 0; i < CovSS.size(); i++)
	for (j = 0; j < CovSS[i].size(); j++)
	  CovSS[i][j] = Cov[k][i][j];

      for (i = 0; i < CovRR.size(); i++)
	for (j = 0; j < CovRR[i].size(); j++)
	  CovRR[i][j] = Cov[k][i + CovSS.size()][j + CovSS.size()];

      for (i = 0; i < CovSR.size(); i++)
	for (j = 0; j < CovSR[i].size(); j++) {
	  CovSR[i][j] = Cov[k][i][j + CovSS.size()];
	  CovRS[j][i] = Cov[k][j + CovSS.size()][i];
	}

      // establish conditional covariance matrix, and corresponding cholesky decomposition

      vector<vector<double> > CovSSInverse;
      inverse(CovSS,CovSSInverse);
      vector<vector<double> > temp1;
      matrixMult(CovRS,CovSSInverse,temp1);
      vector<vector<double> > temp2;
      matrixMult(temp1,CovSR,temp2);
      vector<vector<double> > CovRGivenS(CovRR);
      for (i = 0; i < CovRGivenS.size(); i++)
	for (j = 0; j < CovRGivenS[i].size(); j++)
	  CovRGivenS[i][j] -= temp2[i][j];
      
      int err = 0;
      Cholesky chol(CovRGivenS,err);
      if (err != 0) {
	cout << "Error in Cholesky!!\n";
	exit(-1);
      }

      // generate realisation

      vector<double> vec(CovRGivenS.size(),0.0);
      for (i = 0; i < vec.size(); i++)
	vec[i] = Norm01();
      
      vector<double> z;
      matrixMult(chol.q_L(),vec,z);

      vector<double> obs(CovSS.size(),0.0);
      for (i = 0; i < obs.size(); i++)
	obs[i] = UBlocks[oldClique[k]][oldComponents[k][i]];

      vector<double> mean;
      matrixMult(temp1,obs,mean);
      for (i = 0; i < mean.size(); i++)
	z[i] += mean[i];

      // insert generated values in data structure

      for (i = 0; i < z.size(); i++) 
	U[i + first] = z[i];
      
      first += z.size();

      for (i = 0; i < obs.size(); i++)
	UBlocks[k][i] = obs[i];
      for (i = 0; i < z.size(); i++)
	UBlocks[k][i + obs.size()] = z[i];
    }
    else {
      int err = 0;
      Cholesky chol(Cov[k],err);
      if (err != 0) {
	cout << "Error in Cholesky!!\n";
	exit(-1);
      }
      
      vector<double> vec(Cov[k].size(),0.0);
      for (i = 0; i < vec.size(); i++)
	vec[i] = Norm01();

      vector<double> z;
      matrixMult(chol.q_L(),vec,z);
      for (i = 0; i < z.size(); i++)
	U[i + first] = z[i];

      first += z.size();

      for (i = 0; i < z.size(); i++)
	UBlocks[k][i] = z[i];
    }
  }

  for (i = 0; i < mean.size(); i++)
    U[i] += mean[i];

  return U;
}





vector<vector<double> > Random::MatrixVariateNormal(const vector<vector<double> > &mean,
						    const vector<vector<double> > &R,
						    const vector<vector<vector<double> > > &Omega,
						    const vector<int> &oldClique,
						    const vector<vector<int> > &oldComponents) {
  // allocate space for result

  vector<vector<vector<double> > > UBlocks;
  UBlocks.resize(Omega.size());
  int i,j,k;
  for (k = 0; k < UBlocks.size(); k++) {
    UBlocks[k].resize(Omega[k].size());
    for (i = 0; i < UBlocks[k].size(); i++)
      UBlocks[k][i].resize(R.size());
  }

  vector<vector<double> > U(mean);

  // generate realisations from first component

  int err = 0;
  Cholesky chol(Omega[0],err);
  if (err != 0) {
    cout << "Error in Cholesky!!\n";
    exit(-1);
  }

  int first = 0;
  for (j = 0; j < U[0].size(); j++) {
    vector<double> vec(Omega[0].size(),0.0);
    for (i = 0; i < vec.size(); i++)
      vec[i] = Norm01();
    
    vector<double> z;
    matrixMult(chol.q_L(),vec,z);
    for (i = 0; i < vec.size(); i++)
      U[i][j] = z[i];

    if (j == U[0].size() - 1) first += z.size();

    for (i = 0; i < z.size(); i++)
      UBlocks[0][i][j] = z[i];
  }

  /*
  cout << "U:" << endl;
  for (i = 0; i < U.size(); i++) {
    for (j = 0; j < U[i].size(); j++)
      cout << U[i][j] << " ";
    cout << endl;
  }
  cout << endl;
  */

  // generate realisations from the remaining components
  
  for (k = 1; k < Omega.size(); k++) {
    
    if (oldComponents[k].size() > 0) {
      // establish covariance submatrices
      
      vector<vector<double> > OmegaSS;
      vector<vector<double> > OmegaRR;
      vector<vector<double> > OmegaSR;
      vector<vector<double> > OmegaRS;
      OmegaSS.resize(oldComponents[k].size());
      OmegaRR.resize(Omega[k].size() - oldComponents[k].size());
      OmegaSR.resize(oldComponents[k].size());
      OmegaRS.resize(Omega[k].size() - oldComponents[k].size());
      for (i = 0; i < OmegaSS.size(); i++) OmegaSS[i].resize(OmegaSS.size());
      for (i = 0; i < OmegaRR.size(); i++) OmegaRR[i].resize(OmegaRR.size());
      for (i = 0; i < OmegaSR.size(); i++) OmegaSR[i].resize(OmegaRS.size());
      for (i = 0; i < OmegaRS.size(); i++) OmegaRS[i].resize(OmegaSR.size());
      
      for (i = 0; i < OmegaSS.size(); i++)
	for (j = 0; j < OmegaSS[i].size(); j++)
	  OmegaSS[i][j] = Omega[k][i][j];
      
      for (i = 0; i < OmegaRR.size(); i++)
	for (j = 0; j < OmegaRR[i].size(); j++)
	  OmegaRR[i][j] = Omega[k][i + OmegaSS.size()][j + OmegaSS.size()];
      
      for (i = 0; i < OmegaSR.size(); i++)
	for (j = 0; j < OmegaSR[i].size(); j++) {
	  OmegaSR[i][j] = Omega[k][i][j + OmegaSS.size()];
	  OmegaRS[j][i] = Omega[k][j + OmegaSS.size()][i];
	}
      
      // establish conditional covariance matrix, and corresponding cholesky decomposition
      
      vector<vector<double> > OmegaSSInverse;
      inverse(OmegaSS,OmegaSSInverse);
      vector<vector<double> > temp1;
      matrixMult(OmegaRS,OmegaSSInverse,temp1);
      vector<vector<double> > temp2;
      matrixMult(temp1,OmegaSR,temp2);
      vector<vector<double> > OmegaRGivenS(OmegaRR);
      for (i = 0; i < OmegaRGivenS.size(); i++)
	for (j = 0; j < OmegaRGivenS[i].size(); j++)
	  OmegaRGivenS[i][j] -= temp2[i][j];
      
      int err = 0;
      Cholesky chol(OmegaRGivenS,err);
      if (err != 0) {
	cout << "Error in Cholesky!!\n";
	exit(-1);
      }
      
      // generate realisations
      
      for (j = 0; j < U[0].size(); j++) {
	vector<double> vec(OmegaRR.size(),0.0);
	for (i = 0; i < vec.size(); i++)
	  vec[i] = Norm01();
	
	vector<double> z;
	matrixMult(chol.q_L(),vec,z);
	
	vector<double> obs(OmegaSS.size(),0.0);
	for (i = 0; i < obs.size(); i++)
	  obs[i] = UBlocks[oldClique[k]][oldComponents[k][i]][j];
	
	vector<double> mean;
	matrixMult(temp1,obs,mean);
	for (i = 0; i < mean.size(); i++)
	  z[i] += mean[i];
	
	// insert generated values in data structure
	
	for (i = 0; i < z.size(); i++)
	  U[i + first][j] = z[i];
	
	if (j == U[0].size() - 1) first += z.size();
	
	for (i = 0; i < obs.size(); i++)
	  UBlocks[k][i][j] = obs[i];
	for (i = 0; i < z.size(); i++)
	  UBlocks[k][i + obs.size()][j] = z[i];
      }
    }
    else {
      int err = 0;
      Cholesky chol(Omega[k],err);
      if (err != 0) {
	cout << "Error in Cholesky!!\n";
	exit(-1);
      }

      for (j = 0; j < U[k].size(); j++) {
	vector<double> vec(Omega[k].size(),0.0);
	for (i = 0; i < vec.size(); i++)
	  vec[i] = Norm01();

	vector<double> z;
	matrixMult(chol.q_L(),vec,z);
	for (i = 0; i < vec.size(); i++)
	  U[i + first][j] = z[i];

	if (j == U[k].size() - 1) first += z.size();

	for (i = 0; i < z.size(); i++)
	  UBlocks[k][i][j] = z[i];
      }
    }

    /*
    cout << "U:" << endl;
    for (i = 0; i < U.size(); i++) {
      for (j = 0; j < U[i].size(); j++)
	cout << U[i][j] << " ";
      cout << endl;
    }
    cout << endl;
    */
  }

  /*
  cout << "U:" << endl;
  for (i = 0; i < U.size(); i++) {
    for (j = 0; j < U[i].size(); j++)
      cout << U[i][j] << " ";
    cout << endl;
  }
  cout << endl;
  */

  err = 0;
  Cholesky cholR(R,err);
  if (err != 0) {
    cout << "Error in Cholesky!!\n";
    exit(-1);
  }
  vector<vector<double> > L(cholR.q_L());
  vector<vector<double> > LT;
  LT.resize(L[0].size());
  for (i = 0; i < LT.size(); i++) {
    LT[i].resize(L.size());
    for (j = 0; j < LT[i].size(); j++)
      LT[i][j] = L[j][i];
  }
  
  vector<vector<double> > temp2(U);
  matrixMult(temp2,LT,U);
  
  
  for (i = 0; i < U.size(); i++)
    for (j = 0; j < U[i].size(); j++)
      U[i][j] += mean[i][j];
  
  /*
  cout << "U:" << endl;
  for (i = 0; i < U.size(); i++) {
    for (j = 0; j < U[i].size(); j++)
      cout << U[i][j] << " ";
    cout << endl;
  }
  cout << endl;
  */

  /*
  for (k = 0; k < UBlocks.size(); k++) {
    cout << "UBlocks[" << k << "]:" << endl;
    for (i = 0; i < UBlocks[k].size(); i++) {
      for (j = 0; j < UBlocks[k][i].size(); j++)
	cout << UBlocks[k][i][j] << " ";
      cout << endl;
    }
    cout << endl;
  }
  */
  
  return U;
}






double Random::PotentialHyperInverseWishart(double df,const vector<vector<vector<double> > > &D,
					    const vector<int> &oldClique,const vector<vector<int> > &oldComponents,
					    const vector<vector<vector<double> > > &Sigma) {
  double pot = 0.0;

  int k;
  for (k = 0; k < D.size(); k++) 
    pot += PotentialInverseWishart(df,D[k],Sigma[k]);
  
  // subtract potential for separators

  for (k = 1; k < D.size(); k++) 
    if (oldComponents[k].size() > 0) {
      vector<vector<double> > DSS;
      vector<vector<double> > SigmaSS;
      
      DSS.resize(oldComponents[k].size());
      SigmaSS.resize(oldComponents[k].size());
      int i,j;
      for (i = 0; i < DSS.size(); i++) DSS[i].resize(DSS.size());
      for (i = 0; i < SigmaSS.size(); i++) SigmaSS[i].resize(SigmaSS.size());
      
      for (i = 0; i < DSS.size(); i++)
	for (j = 0; j < DSS[i].size(); j++) {
	  DSS[i][j] = D[oldClique[k]][oldComponents[k][i]][oldComponents[k][j]];
	  SigmaSS[i][j] = Sigma[oldClique[k]][oldComponents[k][i]][oldComponents[k][j]];
	}
      
      pot -= PotentialInverseWishart(df,DSS,SigmaSS);
    }
  
  

  return pot;
}


double Random::PotentialGaussianGraphicalModel(const vector<double> &mean,
					       const vector<vector<vector<double> > > &Cov,
					       const vector<int> &oldClique,
					       const vector<vector<int> > &oldComponents,
					       const vector<double> &U) {
  double pot = 0.0;

  // subtract mean values

  vector<double> UU(U);
  int i,k;
  for (i = 0; i < UU.size(); i++)
    UU[i] -= mean[i];

  // allocate space and initialise temporal strage of U in blocks

  vector<vector<double> > UBlocks;
  UBlocks.resize(Cov.size());
  for (k = 0; k < UBlocks.size(); k++) 
    UBlocks[k].resize(Cov[k].size());

  int first = 0;
  for (i = 0; i < Cov[0].size(); i++)
    UBlocks[0][i] = UU[first + i];
  first += Cov[0].size();

  for (k = 1; k < Cov.size(); k++) {
    for (i = 0; i < oldComponents[k].size(); i++)
      UBlocks[k][i] = UBlocks[oldClique[k]][oldComponents[k][i]];
    
    for (i = 0; i < Cov[k].size() - oldComponents[k].size(); i++)
      UBlocks[k][i + oldComponents[k].size()] = UU[first + i];
    first += Cov[k].size() - oldComponents[k].size();
  }

  // add potential for each clique
  
  for (k = 0; k < Cov.size(); k++) {
    vector<double> zero(UBlocks[k].size(),0.0);
    pot += PotentialMultiGaussian(Cov[k],zero,UBlocks[k]);
  }

  // subtract potential for each seperator

  for (k = 1; k < Cov.size(); k++) {
    if (oldComponents[k].size() > 0) {
      vector<vector<double> > CovSub;
      vector<double> USub;
      CovSub.resize(oldComponents[k].size());
      USub.resize(oldComponents[k].size());
      int j;
      for (i = 0; i < CovSub.size(); i++) {
	CovSub[i].resize(oldComponents[k].size());
	for (j = 0; j < CovSub[i].size(); j++)
	  CovSub[i][j] = Cov[k][i][j];
      }

      for (i = 0; i < USub.size(); i++)
	USub[i] = UBlocks[k][i];
      
      vector<double> zero(USub.size(),0.0);
      pot -= PotentialMultiGaussian(CovSub,zero,USub);
    }
  }

  return pot;
}


double Random::PotentialMatrixVariateNormal(const vector<vector<double> > &Omega,
					    const vector<vector<double> > &R,
					    const vector<vector<double> > &X) {  
  double pot = 0.0;

  vector<vector<double> > OmegaInverse;
  double lndetOmega = inverseLnDeterminant(Omega,OmegaInverse);
  
  vector<vector<double> > RInverse;
  double lndetR = inverseLnDeterminant(R,RInverse);

  vector<vector<double> > XT;
  XT.resize(X[0].size());
  int i,j;
  for (i = 0; i < XT.size(); i++) {
    XT[i].resize(X.size());
    for (j = 0; j < XT[i].size(); j++)
      XT[i][j] = X[j][i];
  }

  vector<vector<double> > temp1;
  matrixMult(XT,OmegaInverse,temp1);
  vector<vector<double> > temp2;
  matrixMult(temp1,X,temp2);
  vector<vector<double> > temp3;
  matrixMult(temp2,RInverse,temp3);

  for (i = 0; i < temp3.size(); i++)
    pot += 0.5 * temp3[i][i];

  pot += 0.5 * ((double) (Omega.size() * R.size())) * log(2.0 * 3.14159265);
  pot += 0.5 * ((double) R.size()) * lndetOmega;
  pot += 0.5 * ((double) Omega.size()) * lndetR;
  
  return pot;
}




double Random::PotentialMatrixVariateNormal(const vector<vector<double> > &mean,
					    const vector<vector<double> > &R,
					    const vector<vector<vector<double> > > &Omega,
					    const vector<int> &oldClique,
					    const vector<vector<int> > &oldComponents,
					    const vector<vector<double> > &U) {
  double pot = 0.0;

  // subtract mean values

  vector<vector<double> > UU(U);
  int i,j,k;
  for (i = 0; i < UU.size(); i++)
    for (j = 0; j < UU[i].size(); j++)
      UU[i][j] -= mean[i][j];

  // allocate space and initialise temporal storage of U in blocks

  vector<vector<vector<double> > > UBlocks;
  UBlocks.resize(Omega.size());
  for (k = 0; k < UBlocks.size(); k++) {
    UBlocks[k].resize(Omega[k].size());
    for (i = 0; i < UBlocks[k].size(); i++)
      UBlocks[k][i].resize(R.size());
  }

  int first = 0;
  for (i = 0; i < Omega[0].size(); i++)
    for (j = 0; j < UU[first].size(); j++)
      UBlocks[0][i][j] = UU[first + i][j];
  first += Omega[0].size();

  for (k = 1; k < Omega.size(); k++) {
    for (i = 0; i < oldComponents[k].size(); i++)
      for (j = 0; j < UU[first].size(); j++)
	UBlocks[k][i][j] = UBlocks[oldClique[k]][oldComponents[k][i]][j];
    
    for (i = 0; i < Omega[k].size() - oldComponents[k].size(); i++)
      for (j = 0; j < UU[first].size(); j++)
	UBlocks[k][i + oldComponents[k].size()][j] = UU[first + i][j];
    first += Omega[k].size() - oldComponents[k].size();
  }
  
  /*
  for (k = 0; k < UBlocks.size(); k++) {
    cout << "UBlocks[" << k << "]:" << endl;
    for (i = 0; i < UBlocks[k].size(); i++) {
      for (j = 0; j < UBlocks[k][i].size(); j++)
	cout << UBlocks[k][i][j] << " ";
      cout << endl;
    }
    cout << endl;
  }
  */
  
  // add potential for each clique
  
  for (k = 0; k < Omega.size(); k++)
    pot += PotentialMatrixVariateNormal(Omega[k],R,UBlocks[k]);
  
  // subtract potential for each separator

  for (k = 1; k < Omega.size(); k++) 
    if (oldComponents[k].size() > 0) {
      vector<vector<double> > OmegaSub;
      vector<vector<double> > USub;
      OmegaSub.resize(oldComponents[k].size());
      USub.resize(oldComponents[k].size());
      for (i = 0; i < OmegaSub.size(); i++) {
	OmegaSub[i].resize(oldComponents[k].size());
	for (j = 0; j < OmegaSub[i].size(); j++)
	  OmegaSub[i][j] = Omega[k][i][j];
      }
      for (i = 0; i < USub.size(); i++) {
	USub[i].resize(UBlocks[k][i].size());
	for (j = 0; j < USub[i].size(); j++)
	  USub[i][j] = UBlocks[k][i][j];
      }
      
      pot -= PotentialMatrixVariateNormal(OmegaSub,R,USub);
    }
  

  return pot;
}
			 





double Random::PotentialGaussian(double variance,double mean,double x)
{
  double diff = x - mean;
  double pot = diff * diff / variance;
  pot += log(2.0 * PI);
  pot += log(variance);
  pot *= 0.5;

  return pot;
}





double Random::PotentialPoisson(double lambda,int x)
{
  double pot;
  
  pot = - x * log(lambda);
  pot += lambda;
  
  int k;
  for (k = 2; k <= x; k++)
    pot += log((double) k);

  return pot;
}





double Random::PotentialBinomial(int n,double p,int x)
{
  double pot = - x * log(p) - (n - x) * log(1.0 - p);

  if (x > 0)
    {
      int k;
      for (k = 1; k <= x; k++)
	pot += - log((double) (n - k + 1)) + log((double) k);
    }

  return pot;
}




double Random::PotentialGamma(double p,double lambda,double x)
{
  double pot = - p * log(lambda) + lnGamma(p);
  pot += - (p - 1.0) * log(x);
  pot += lambda * x;
  
  return pot;
}




double Random::PotentialInverseGamma(double p,double lambda,double x)
{
  double pot = - p * log(lambda) + lnGamma(p);
  pot += (p + 1.0) * log(x);
  pot += lambda / x;

  return pot;
}



double Random::PotentialBeta(double alpha,double beta,double x)
{
  double pot = 0.0;
  
  pot -= lnGamma(alpha + beta);
  pot += lnGamma(alpha);
  pot += lnGamma(beta);

  pot -= (alpha - 1) * log(x);
  pot -= (beta - 1) * log(1.0 - x);

  return pot;
}



double Random::PotentialMultiGaussian(const vector<vector<double> > &Sigma,
				      const vector<double> &mean,
				      const vector<double> &x)
{
  int n = x.size();

  vector<double> diff(n);
  vector<vector<double> > SigmaInverse;
  double determinant = inverse(Sigma,SigmaInverse);

  int i;
  for (i = 0; i < n; i++)
    diff[i] = x[i] - mean[i];

  double pot = 0.5 * quadratic(SigmaInverse,diff);
  pot += 0.5 * log(determinant);

  pot += ((double) n) * log(2.0 * PI) / 2.0;


  return pot;
}




double Random::PotentialMultiGaussian(const vector<vector<double> > &SigmaInv,
				      double determinant,const vector<double> &mean,
				      const vector<double> &x)
{
  int n = x.size();

  vector<double> diff(mean.size());
  int i;
  for (i = 0; i < n; i++)
    diff[i] = x[i] - mean[i];

  double sum = 0.0;
  for (i = 0; i < n; i++)
    sum += SigmaInv[i][i] * diff[i] * diff[i];

  int j;
  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++)
      sum += 2.0 * SigmaInv[i][j] * diff[i] * diff[j];

  double pot = 0.5 * sum;
  pot += 0.5 * log(determinant);

  pot += ((double) n) * log(2.0 * PI) / 2.0;

  return pot;
}




double Random::PotentialMultiGaussian(const vector<vector<double> > &Sigma,
				      const vector<double> &x)
{
  int n = x.size();

  vector<double> diff(n);
  vector<vector<double> > SigmaInverse;
  double determinant = inverse(Sigma,SigmaInverse);

  int i;
  for (i = 0; i < n; i++)
    diff[i] = x[i];

  double pot = 0.5 * quadratic(SigmaInverse,diff);
  pot += 0.5 * log(determinant);

  pot += ((double) n) * log(2.0 * PI) / 2.0;


  return pot;
}




double Random::PotentialMultiGaussian(const vector<vector<double> > &SigmaInv,
				      double determinant,const vector<double> &x)
{
  int n = x.size();

  double sum = 0.0;
  int i;
  for (i = 0; i < n; i++)
    sum += SigmaInv[i][i] * x[i] * x[i];

  int j;
  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++)
      sum += 2.0 * SigmaInv[i][j] * x[i] * x[j];

  double pot = 0.5 * sum;
  pot += 0.5 * log(determinant);

  pot += ((double) n) * log(2.0 * PI) / 2.0;


  return pot;
}




double Random::PotentialCauchy(double var,double mean,double x)
{
  double sigma = sqrt(var);

  x = (x - mean) / sigma;

  double pot = log(PI) + log(sigma) + log(1.0 + x * x);

  return pot;
}





double Random::PotentialTScaled(double var,double mean,double df,double x)
{
  double sigma = sqrt(var);
  x = (x - mean) / sigma;

  double pot = PotentialT(df,x);
  pot += log(sigma);

  return pot;
}




double Random::PotentialT(double df,double x)
{
  double pot = 0.0;

  pot -= lnGamma((df + 1.0) / 2.0);
  pot += lnGamma(df / 2.0);
  pot += 0.5 * log(PI * df);
  pot += 0.5 * (df + 1.0) * log(1.0 + x * x / df);

  return pot;
}











double Random::PotentialWishartAlternativeParam(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x)
{
  int m = x.size();

  vector<vector<double> > temp1;
  vector<vector<double> > temp2;
  double detV = inverse(V,temp1);
  double detX = inverse(x,temp2);

  vector<vector<double> > prod;
  matrixMult(temp1,x,prod);

  double trace = 0.0;
  int i;
  for (i = 0; i < m; i++)
    trace += prod[i][i];

  double pot = - (nu - ((double) (m + 1))) * log(detX) / 2.0;
  pot += trace / 2.0;
  pot += nu * log(detV) / 2.0;
  pot +=  nu  * ((double) m) * log(2.0) / 2.0;
  pot += m * (m - 1) * log(PI) / 4.0;
  for (i = 1; i <= m; i++)
    pot += lnGamma((nu - ((double) (i - 1))) / 2.0);

  return pot;
}






double Random::PotentialStandardWishartAlternativeParam(double nu,const vector<vector<double> > &x)
{
  int m = x.size();
  vector<vector<double> > prod;
  double detX = inverse(x,prod);

  double trace = 0.0;
  int i;
  for (i = 0; i < m; i++)
    trace += prod[i][i];

  double pot = - (nu - ((double) (m + 1))) * log(detX) / 2.0;
  pot += trace / 2.0;
  pot +=  nu  * ((double) m) * log(2.0) / 2.0;
  pot += m * (m - 1) * log(PI) / 4.0;
  for (i = 1; i <= m; i++)
    pot += lnGamma((nu - ((double) (i - 1))) / 2.0);

  return pot;
}






double Random::PotentialInverseWishartAlternativeParam(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x)
{
  int m = x.size();
  
  vector<vector<double> > temp1;
  vector<vector<double> > xInverse;
  double detV = inverse(V,temp1);
  double lndetX = inverseLnDeterminant(x,xInverse);

  int pp = 0;
  if (pp == 1) {
    FILE *out = fopen("x.txt","w");
    int i,j;
    for (i = 0; i < m; i++) {
      for (j = 0; j < m; j++)
	fprintf(out,"%20.18e ",x[i][j]);
      fprintf(out,"\n");
    }
    fclose(out);
  }

  vector<vector<double> > prod;
  matrixMult(V,xInverse,prod);

  int i;
  double trace = 0.0;
  for (i = 0; i < m; i++)
    trace += prod[i][i];

  double pot = (nu + ((double) (m + 1))) * lndetX / 2.0;
  pot += trace / 2.0;
  pot += - nu * log(detV) / 2.0;
  pot +=  nu  * ((double) m) * log(2.0) / 2.0;
  pot += m * (m - 1) * log(PI) / 4.0;
  for (i = 1; i <= m; i++)
    pot += lnGamma(nu - ((double) (i - 1)) / 2.0);

  return pot;
}






double Random::PotentialStandardInverseWishartAlternativeParam(double nu,const vector<vector<double> > &x)
{
  int m = x.size();

  vector<vector<double> > xInverse;
  double detX = inverse(x,xInverse);

  int i;
  double trace = 0.0;
  for (i = 0; i < m; i++)
    trace += xInverse[i][i];

  double pot = (nu + ((double) (m + 1))) * log(detX) / 2.0;
  pot += trace / 2.0;
  pot +=  nu  * ((double) m) * log(2.0) / 2.0;
  pot += m * (m - 1) * log(PI) / 4.0;
  for (i = 1; i <= m; i++)
    pot += lnGamma((nu - ((double) (i - 1))) / 2.0);

  return pot;
}




double Random::PotentialCorrelationStandardInverseWishartAlternativeParam(double nu,const vector<vector<double> > &x)
{
  int m = x.size();
  
  vector<vector<double> > xInverse;
  double detX = inverse(x,xInverse);

  double pot = (nu + ((double) (m + 1))) * log(detX) / 2.0;
  int i;
  for (i = 0; i < m; i++)
    pot += (nu / 2.0) * log(xInverse[i][i]);

  pot -=  m  * log(2.0);
  pot -= m * lnGamma(nu / 2.0);
  pot += m * (m - 1) * log(PI) / 4.0;
  for (i = 1; i <= m; i++)
    pot += lnGamma((nu - ((double) (i - 1))) / 2.0);

  return pot;
}






vector<int> Random::Permutation(int n)
{
  vector<int> perm(n,0);
  int k;
  for (k = 0; k < perm.size(); k++)
    perm[k] = k;

  for (k = perm.size() - 1; k >= 0; k--)
    {
      int kk = (int) (Unif01() * (k + 1));
      int temp = perm[kk];
      perm[kk] = perm[k];
      perm[k] = temp;
    }

  return perm;
}






double Random::lnGamma(double xx)
{
  //
  // return ln(Gamma(xx)). 
  //

  double x,y,tt,sum;
  static double coeff[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
			  -0.5395239384953e-5};
  int j;
  
  y = xx;
  tt = xx + 5.5 - (xx + 0.5) * log(xx + 5.5);
  sum = 1.000000000190015;
  for (j = 0; j <= 5; j++) sum += coeff[j]/++y;
  
  double answer = - tt + log(2.5066282746310005 * sum / xx);

  return answer;
}



vector<vector<double> > Random::Wishart(double nu,const vector<vector<double> > &V) {
  double nuAlt = nu  + ((double) V.size()) - 1.0;

  return WishartAlternativeParam(nuAlt,V);
}

vector<vector<double> > Random::StandardWishart(int dim,double nu) {
  double nuAlt = nu  + ((double) dim) - 1.0;
  
  return StandardWishartAlternativeParam(dim,nuAlt);
}

vector<vector<double> > Random::InverseWishart(double nu,const vector<vector<double> > &V) {
  double nuAlt = nu  + ((double) V.size()) - 1.0;

  return InverseWishartAlternativeParam(nuAlt,V);
}

vector<vector<double> > Random::StandardInverseWishart(int dim,double nu) {
  double nuAlt = nu  + ((double) dim) - 1.0;
  
  return StandardInverseWishartAlternativeParam(dim,nuAlt);
}

vector<vector<double> > Random::CorrelationStandardInverseWishart(int dim,double nu) {
  double nuAlt = nu  + ((double) dim) - 1.0;
  
  return CorrelationStandardInverseWishartAlternativeParam(dim,nuAlt);
}


double Random::PotentialWishart(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x) {
  double nuAlt = nu  + ((double) V.size()) - 1.0;
  
  return PotentialWishartAlternativeParam(nuAlt,V,x);
}

double Random::PotentialStandardWishart(double nu,const vector<vector<double> > &x) {
  double nuAlt = nu  + ((double) x.size()) - 1.0;
  
  return PotentialStandardWishartAlternativeParam(nuAlt,x);
}

double Random::PotentialInverseWishart(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x) {
  double nuAlt = nu  + ((double) V.size()) - 1.0;
  
  return PotentialInverseWishartAlternativeParam(nuAlt,V,x);
}

double Random::PotentialStandardInverseWishart(double nu,const vector<vector<double> > &x) {
  double nuAlt = nu  + ((double) x.size()) - 1.0;
  
  return PotentialStandardInverseWishartAlternativeParam(nuAlt,x);
}

double Random::PotentialCorrelationStandardInverseWishart(double nu,const vector<vector<double> > &x) {
  double nuAlt = nu  + ((double) x.size()) - 1.0;
  
  return PotentialCorrelationStandardInverseWishartAlternativeParam(nuAlt,x);
}

 

