#ifndef POTENTIALRHO_H
#define POTENTIALRHO_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialRho : public Potential
{
 public:

  PotentialRho(const Structure *str);
  virtual ~PotentialRho(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
};



inline PotentialRho::PotentialRho(const Structure *str) : Potential()
{
  this->str = str;

  return;
}



inline PotentialRho::~PotentialRho(void)
{
  return;
}


inline Potential *PotentialRho::copy(void) const
{
  Potential *pp = new PotentialRho(str);

  return pp;
}


inline double PotentialRho::potential(Random &ran) const
{
  double pot = ran.PotentialCorrelationStandardInverseWishart(str->nuRho,str->rho);

  return pot;
}



#endif
