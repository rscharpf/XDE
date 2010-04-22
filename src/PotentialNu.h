#ifndef POTENTIALNU_H
#define POTENTIALNU_H

#include "Potential.h"
#include "PotentialNug.h"
#include "Structure.h"
#include "Random.h"
#include "Matrix.h"


class PotentialNu : public Potential
{
 public:

  PotentialNu(const Structure *str);
  virtual ~PotentialNu(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
};



inline PotentialNu::PotentialNu(const Structure *str) : Potential()
{
  this->str = str;

  return;
}



inline PotentialNu::~PotentialNu(void)
{
  return;
}


inline Potential *PotentialNu::copy(void) const
{
  Potential *pp = new PotentialNu(str);

  return pp;
}


inline double PotentialNu::potential(Random &ran) const
{
  double pot = 0.0;
  int Q = str->Q;
  int G = str->G;

  int p,q;
  vector<vector<double> > Rho(str->Q);
  vector<vector<double> > RhoInv;
  for (p = 0; p < Q; p++)
    Rho[p].resize(Q);
  
  for (p = 0; p < Q; p++)
    for (q = 0; q < Q; q++)
      Rho[p][q] = str->rho[p][q];
  
  double determinant = inverse(Rho,RhoInv);

  //  Matrix Cov(Q,Q);
  //  int p,q;
  //  for (p = 0; p < Q; p++)
  //    for (q = 0; q < Q; q++)
  //      Cov(p+1,q+1) = str->rho[p][q];
  //
  //  double determinant = 0.0;
  //  Cov.invert(&determinant);
  //
  //  vector<vector<double> > RhoInv(str->Q);
  //  for (p = 0; p < Q; p++)
  //    RhoInv[p].resize(Q);
  //
  //  for (p = 0; p < Q; p++)
  //    for (q = 0; q < Q; q++)
  //      RhoInv[p][q] = Cov(p+1,q+1);

  vector<double> vv(str->Q);
  for (p = 0; p < Q; p++)
    vv[p] = str->gamma2 * str->tau2Rho[p];

  vector<double> value(Q);

  int g;
  for (g = 0; g < G; g++)
    {
      double det = determinant;

      for (p = 0; p < Q; p++)
	{
	  double variance = vv[p] * exp(str->a[p] * log(str->sigma2[p][g]));
	  det *= variance;
	  value[p] = str->nu[p][g] / sqrt(variance);
	}
      
      pot += ran.PotentialMultiGaussian(RhoInv,det,value);
    }

  return pot;
}



#endif
