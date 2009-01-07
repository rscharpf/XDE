#include <iostream>
#include <limits.h>

#include "Random.h"
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
  if (err != 0)
    {
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




vector<vector<double> > Random::Wishart(double nu,const vector<vector<double> > &V)
{
  int err = 0;
  Cholesky chol(V,err);
  //  assert(err == 0);

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
    




vector<vector<double> > Random::InverseWishart(double nu,const vector<vector<double> > &V)
{
  vector<vector<double> > w(V.size());
  int i;
  for (i = 0; i < V.size(); i++)
    w[i].resize(V.size());
  w = Wishart(nu,V);

  vector<vector<double> > ww;
  inverse(w,ww);

  return ww;
}






vector<vector<double> >Random::StandardWishart(int dim,double nu)
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
    






vector<vector<double> > Random::StandardInverseWishart(int dim,double nu)
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





vector<vector<double> > Random::CorrelationStandardInverseWishart(int dim,double nu)
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











double Random::PotentialWishart(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x)
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






double Random::PotentialStandardWishart(double nu,const vector<vector<double> > &x)
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






double Random::PotentialInverseWishart(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x)
{
  int m = x.size();
  
  vector<vector<double> > temp1;
  vector<vector<double> > xInverse;
  double detV = inverse(V,temp1);
  double detX = inverse(x,xInverse);

  vector<vector<double> > prod;
  matrixMult(temp1,xInverse,prod);

  int i;
  double trace = 0.0;
  for (i = 0; i < m; i++)
    trace += prod[i][i];

  double pot = (nu + ((double) (m + 1))) * log(detX) / 2.0;
  pot += trace / 2.0;
  pot += nu * log(detV) / 2.0;
  pot +=  nu  * ((double) m) * log(2.0) / 2.0;
  pot += m * (m - 1) * log(PI) / 4.0;
  for (i = 1; i <= m; i++)
    pot += lnGamma((nu - ((double) (i - 1))) / 2.0);

  return pot;
}






double Random::PotentialStandardInverseWishart(double nu,const vector<vector<double> > &x)
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




double Random::PotentialCorrelationStandardInverseWishart(double nu,const vector<vector<double> > &x)
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



