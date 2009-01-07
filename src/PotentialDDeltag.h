#ifndef POTENTIALDDELTAG_H
#define POTENTIALDDELTAG_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"
#include "Matrix.h"


class PotentialDDeltag : public Potential
{
 public:

  PotentialDDeltag(int g,const Structure *str);
  virtual ~PotentialDDeltag(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  int g;
  const Structure *str;
};



inline PotentialDDeltag::PotentialDDeltag(int g,const Structure *str) : Potential()
{
  this->g = g;
  this->str = str;

  return;
}



inline PotentialDDeltag::~PotentialDDeltag(void)
{
  return;
}


inline Potential *PotentialDDeltag::copy(void) const
{
  Potential *pp = new PotentialDDeltag(g,str);

  return pp;
}


inline double PotentialDDeltag::potential(Random &ran) const
{
  double pot = 0.0;

  vector<vector<double> > R;
  R.resize(str->Q);
  int q;
  for (q = 0; q < str->Q; q++)
    R[q].resize(str->Q);
  
  int p;
  for (p = 0; p < str->Q; p++)
    {
      R[p][p] = str->c2 * str->tau2R[p];
      R[p][p] *= exp(str->b[p] * log(str->sigma2[p][g]));
    }

  for (p = 0; p < str->Q; p++)
    for (q = p + 1; q < str->Q; q++)
      {
	R[p][q] = str->c2;
	R[p][q] *= str->r[p][q];
	R[p][q] *= sqrt(str->tau2R[p] * str->tau2R[q]);
	R[p][q] *= exp(0.5 * (str->b[q] * log(str->sigma2[q][g]) + str->b[p] * log(str->sigma2[p][g])));

	R[q][p] = R[p][q];
      }
  vector<double> value(str->Q,0);
  for (q = 0; q < str->Q; q++)
    value[q] = str->Delta[q][g];
  
  vector<vector<double> > RInv;
  double determinant = inverse(R,RInv);

  //  Matrix Cov(R.size(),R.size());
  //  for (p = 0; p < str->Q; p++)
  //    for (q = 0; q < str->Q; q++)
  //      Cov(p+1,q+1) = R[p][q];
  //
  //  double determinant = 0.0;
  //  Cov.invert(&determinant);
  //
  //  for (p = 0; p < str->Q; p++)
  //    for (q = 0; q < str->Q; q++)
  //      R[p][q] = Cov(p+1,q+1);

  pot += ran.PotentialMultiGaussian(RInv,determinant,value);
  
  return pot;
}



#endif
