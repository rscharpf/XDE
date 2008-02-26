#ifndef POTENTIALXI_H
#define POTENTIALXI_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialXi : public Potential
{
 public:

  PotentialXi(const Structure *str);
  virtual ~PotentialXi(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
};



inline PotentialXi::PotentialXi(const Structure *str) : Potential()
{
  this->str = str;

  return;
}



inline PotentialXi::~PotentialXi(void)
{
  return;
}


inline Potential *PotentialXi::copy(void) const
{
  Potential *pp = new PotentialXi(str);

  return pp;
}


inline double PotentialXi::potential(Random &ran) const
{
  double pot = ran.PotentialBeta(str->alphaXi,str->betaXi,str->xi);

  return pot;
}



#endif
