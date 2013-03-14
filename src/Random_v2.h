#ifndef RANDOM_H
#define RANDOM_H

#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

using namespace std;



#define MULTIPLIER 69069
#define SHIFT          1
#define FACTOR1    256
#define FACTOR2    256
#define FACTOR3    256
#define FACTOR4    128
#define PI         3.14159265


class Random
{
 public:
  Random(unsigned int seed);
  ~Random(void);

  unsigned int ChangeSeed(unsigned int seed);
  double Unif01(void);
  double Norm01(void);
  double Exponential(double lambda);
  int Poisson(double lambda);
  int Binomial(int n,double p);
  int Discrete(const vector<double> &prob);
  double Gamma(double p,double lambda);
  double ChiSquared(double nu);
  double InverseGamma(double p,double lambda);
  double Beta(double alpha,double beta);
  vector<double> MultiGaussian(const vector<vector<double> > &Sigma,
			       const vector<double> &mean);

  vector<vector<double> > Wishart(double nu,const vector<vector<double> > &V);
  vector<vector<double> > StandardWishart(int dim,double nu);
  vector<vector<double> > InverseWishart(double nu,const vector<vector<double> > &V);
  vector<vector<double> > StandardInverseWishart(int dim,double nu);
  vector<vector<double> > CorrelationStandardInverseWishart(int dim,double nu);

  vector<vector<vector<double> > > HyperInverseWishart(double df,const vector<vector<vector<double> > > &D,
						       const vector<int> &oldClique,
						       const vector<vector<int> > &oldComponents);
  vector<double> GaussianGraphicalModel(const vector<double> &mean,
					const vector<vector<vector<double> > > &Cov,
					const vector<int> &oldClique,
					const vector<vector<int> > &oldComponents);
  vector<vector<double> > MatrixVariateNormal(const vector<vector<double> > &mean,
					      const vector<vector<double> > &R,
					      const vector<vector<vector<double> > > &Omega,
					      const vector<int> &oldClique,
					      const vector<vector<int> > &oldComponents);



  double PotentialGaussian(double variance,double mean,double x);
  double PotentialPoisson(double lambda,int x);
  double PotentialBinomial(int n,double p,int x);
  double PotentialGamma(double p,double lambda,double x);
  double PotentialInverseGamma(double p,double lambda,double x);
  double PotentialBeta(double alpha,double beta,double x);
  double PotentialMultiGaussian(const vector<vector<double> > &Sigma,
				const vector<double> &mean,
				const vector<double> &x);
  double PotentialMultiGaussian(const vector<vector<double> > &SigmaInv,
				double determinant,const vector<double> &mean,
				const vector<double> &x);
  double PotentialMultiGaussian(const vector<vector<double> > &Sigma,
				const vector<double> &x);
  double PotentialMultiGaussian(const vector<vector<double> > &SigmaInv,
				double determinant,const vector<double> &x);
  double PotentialTScaled(double var,double mean,double df,double x);
  double PotentialT(double df,double x);
  double PotentialCauchy(double var,double mean,double x);

  double PotentialWishart(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x);
  double PotentialStandardWishart(double nu,const vector<vector<double> > &x);
  double PotentialInverseWishart(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x);
  double PotentialStandardInverseWishart(double nu,const vector<vector<double> > &x);
  double PotentialCorrelationStandardInverseWishart(double nu,const vector<vector<double> > &x);

  double PotentialHyperInverseWishart(double df,const vector<vector<vector<double> > > &D,
				      const vector<int> &oldClique,const vector<vector<int> > &oldComponents,
				      const vector<vector<vector<double> > > &Sigma);
  double PotentialGaussianGraphicalModel(const vector<double> &mean,
					 const vector<vector<vector<double> > > &Cov,
					 const vector<int> &oldClique,
					 const vector<vector<int> > &oldComponents,
					 const vector<double> &U);
  double PotentialMatrixVariateNormal(const vector<vector<double> > &Omega,
				      const vector<vector<double> > &R,
				      const vector<vector<double> > &X);
  double PotentialMatrixVariateNormal(const vector<vector<double> > &mean,
				      const vector<vector<double> > &R,
				      const vector<vector<vector<double> > > &Omega,
				      const vector<int> &oldClique,
				      const vector<vector<int> > &oldComponents,
				      const vector<vector<double> > &U);

  vector<int> Permutation(int n);

  double lnGamma(double x);

 private:
  vector<vector<double> > WishartAlternativeParam(double nu,const vector<vector<double> > &V);
  vector<vector<double> > StandardWishartAlternativeParam(int dim,double nu);
  vector<vector<double> > InverseWishartAlternativeParam(double nu,const vector<vector<double> > &V);
  vector<vector<double> > StandardInverseWishartAlternativeParam(int dim,double nu);
  vector<vector<double> > CorrelationStandardInverseWishartAlternativeParam(int dim,double nu);

  double PotentialWishartAlternativeParam(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x);
  double PotentialStandardWishartAlternativeParam(double nu,const vector<vector<double> > &x);
  double PotentialInverseWishartAlternativeParam(double nu,const vector<vector<double> > &V,const vector<vector<double> > &x);
  double PotentialStandardInverseWishartAlternativeParam(double nu,const vector<vector<double> > &x);
  double PotentialCorrelationStandardInverseWishartAlternativeParam(double nu,const vector<vector<double> > &x);

  unsigned int modulus;

  unsigned int seedValue;
  void setseed();
  int haveNorm01;
  double norm;
};

#endif
