#ifndef POTENTIALPHIQG_H
#define POTENTIALPHIQG_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialPhiqg : public Potential
{
 public:

  PotentialPhiqg(int q,int g,const Structure *str);
  virtual ~PotentialPhiqg(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  int q,g;
  const Structure *str;
};



inline PotentialPhiqg::PotentialPhiqg(int q,int g,const Structure *str) : Potential()
{
  this->q = q;
  this->g = g;
  this->str = str;

  return;
}



inline PotentialPhiqg::~PotentialPhiqg(void)
{
  return;
}


inline Potential *PotentialPhiqg::copy(void) const
{
  Potential *pp = new PotentialPhiqg(q,g,str);

  return pp;
}


inline double PotentialPhiqg::potential(Random &ran) const
{
  double pot = 0.0;

  double param2 = str->lambda[q] / str->theta[q];
  double param1 = str->lambda[q] * param2;

  pot += ran.PotentialGamma(param1,param2,str->phi[q][g]);
  
  return pot;
}



#endif
