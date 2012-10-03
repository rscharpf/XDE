#ifndef POTENTIALR_H
#define POTENTIALR_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialR : public Potential
{
 public:

  PotentialR(const Structure *str);
  virtual ~PotentialR(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
};



inline PotentialR::PotentialR(const Structure *str) : Potential()
{
  this->str = str;

  return;
}



inline PotentialR::~PotentialR(void)
{
  return;
}


inline Potential *PotentialR::copy(void) const
{
  Potential *pp = new PotentialR(str);

  return pp;
}


inline double PotentialR::potential(Random &ran) const
{
  double pot = ran.PotentialCorrelationStandardInverseWishart(str->nuR,str->r);

  return pot;
} 



#endif
