#ifndef POTENTIALNUG_H
#define POTENTIALNUG_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"
#include "Matrix.h"

class PotentialNug : public Potential
{
 public:

  PotentialNug(int g,const Structure *str);
  virtual ~PotentialNug(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  int g;
  const Structure *str;
};



inline PotentialNug::PotentialNug(int g,const Structure *str) : Potential()
{
  this->g = g;
  this->str = str;

  return;
}



inline PotentialNug::~PotentialNug(void)
{
  return;
}


inline Potential *PotentialNug::copy(void) const
{
  Potential *pp = new PotentialNug(g,str);

  return pp;
}


inline double PotentialNug::potential(Random &ran) const
{
  double pot = 0.0;

  vector<vector<double> > Sigma;
  Sigma.resize(str->Q);
  int q;
  for (q = 0; q < str->Q; q++)
    Sigma[q].resize(str->Q);
  
  int p;
  for (p = 0; p < str->Q; p++)
    {
      Sigma[p][p] = str->gamma2 * str->tau2Rho[p];
      Sigma[p][p] *= exp(str->a[p] * log(str->sigma2[p][g]));
    }

  for (p = 0; p < str->Q; p++)
    for (q = p + 1; q < str->Q; q++)
      {
	Sigma[p][q] = str->gamma2;
	Sigma[p][q] *= str->rho[p][q];
	Sigma[p][q] *= sqrt(str->tau2Rho[p] * str->tau2Rho[q]);
	Sigma[p][q] *= exp(0.5 * (str->a[q] * log(str->sigma2[q][g]) + str->a[p] * log(str->sigma2[p][g])));

	Sigma[q][p] = Sigma[p][q];
      }
  vector<double> value(str->Q,0);
  for (q = 0; q < str->Q; q++)
    value[q] = str->nu[q][g];
  
  vector<vector<double> > SSigma;
  double determinant = inverse(Sigma,SSigma);

//  Matrix Cov(Sigma.size(),Sigma.size());
//  for (p = 0; p < str->Q; p++)
//    for (q = 0; q < str->Q; q++)
//      Cov(p+1,q+1) = Sigma[p][q];
//
//  double determinant = 0.0;
//  Cov.invert(&determinant);
//  for (p = 0; p < str->Q; p++)
//    for (q = 0; q < str->Q; q++)
//      Sigma[p][q] = Cov(p+1,q+1);
  
  pot += ran.PotentialMultiGaussian(SSigma,determinant,value);
  
  return pot;
}



#endif
