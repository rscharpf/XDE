#ifndef POTENTIALTAU2_H
#define POTENTIALTAU2_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialTau2 : public Potential
{
 public:

  PotentialTau2(const Structure *str);
  virtual ~PotentialTau2(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
};



inline PotentialTau2::PotentialTau2(const Structure *str) : Potential()
{
  this->str = str;

  return;
};



inline PotentialTau2::~PotentialTau2(void)
{
  return;
};


inline Potential *PotentialTau2::copy(void) const
{
  Potential *pp = new PotentialTau2(str);

  return pp;
};


inline double PotentialTau2::potential(Random &ran) const
{
  double pot = 0.0;

  int q;
  for (q = 0; q < str->Q; q++)
    pot += ran.PotentialGamma(str->sT,str->sT,str->tau2[q]);

  return pot;
};



#endif
