#ifndef POTENTIALDDELTA_H
#define POTENTIALDDELTA_H

#include "Potential.h"
#include "PotentialDDeltag.h"
#include "Structure.h"
#include "Random.h"
#include "Matrix.h"


class PotentialDDelta : public Potential
{
 public:

  PotentialDDelta(const Structure *str);
  virtual ~PotentialDDelta(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
};



inline PotentialDDelta::PotentialDDelta(const Structure *str) : Potential()
{
  this->str = str;

  return;
}



inline PotentialDDelta::~PotentialDDelta(void)
{
  return;
}


inline Potential *PotentialDDelta::copy(void) const
{
  Potential *pp = new PotentialDDelta(str);

  return pp;
}


inline double PotentialDDelta::potential(Random &ran) const
{
  double pot = 0.0;
  int Q = str->Q;
  int G = str->G;

  int p,q;
  vector<vector<double> > Cov;
  vector<vector<double> > RInv;
  Cov.resize(Q);
  for (p = 0; p < Q; p++)
    Cov[p].resize(Q);

  for (p = 0; p < Q; p++)
    for (q = 0; q < Q; q++)
      Cov[p][q] = str->r[p][q];

  double determinant = inverse(Cov,RInv);

  //  Matrix Cov(Q,Q);
  //  int p,q;
  //  for (p = 0; p < Q; p++)
  //    for (q = 0; q < Q; q++)
  //      Cov(p+1,q+1) = str->r[p][q];
  //
  //  double determinant = 0.0;
  //  Cov.invert(&determinant);
  //
  //  vector<vector<double> > RInv(Q);
  //  for (p = 0; p < Q; p++)
  //    RInv[p].resize(Q);
  //
  //  for (p = 0; p < Q; p++)
  //    for (q = 0; q < Q; q++)
  //      RInv[p][q] = Cov(p+1,q+1);


  vector<double> vv(Q);
  for (p = 0; p < Q; p++)
    vv[p] = str->c2 * str->tau2R[p];

  vector<double> value(Q);

  int g;
  for (g = 0; g < G; g++)
    {
      double det = determinant;

      for (p = 0; p < Q; p++)
	{
	  double variance = vv[p] * exp(str->b[p] * log(str->sigma2[p][g]));
	  det *= variance;
	  value[p] = str->Delta[p][g] / sqrt(variance);
	}
      
      pot += ran.PotentialMultiGaussian(RInv,det,value);
    }

  return pot;
} 



#endif
