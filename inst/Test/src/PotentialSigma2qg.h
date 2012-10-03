#ifndef POTENTIALSIGMA2QG_H
#define POTENTIALSIGMA2QG_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialSigma2qg : public Potential
{
 public:

  PotentialSigma2qg(int q,int g,const Structure *str);
  virtual ~PotentialSigma2qg(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  int q,g;
  const Structure *str;
};



inline PotentialSigma2qg::PotentialSigma2qg(int q,int g,const Structure *str) : Potential()
{
  this->q = q;
  this->g = g;
  this->str = str;

  return;
}



inline PotentialSigma2qg::~PotentialSigma2qg(void)
{
  return;
}


inline Potential *PotentialSigma2qg::copy(void) const
{
  Potential *pp = new PotentialSigma2qg(q,g,str);

  return pp;
}


inline double PotentialSigma2qg::potential(Random &ran) const
{
  double param2 = str->l[q] / str->t[q];
  double param1 = str->l[q] * param2;
  
  double pot = ran.PotentialGamma(param1,param2,str->sigma2[q][g]);
  
  return pot;
}



#endif
