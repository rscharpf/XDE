#ifndef POTENTIALZERO_H
#define POTENTIALZERO_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialZero : public Potential
{
 public:

  PotentialZero(void);
  virtual ~PotentialZero(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
};



inline PotentialZero::PotentialZero(void) : Potential()
{
  return;
}



inline PotentialZero::~PotentialZero(void)
{
  return;
}


inline Potential *PotentialZero::copy(void) const
{
  Potential *pp = new PotentialZero;

  return pp;
}


inline double PotentialZero::potential(Random &ran) const
{
  return 0.0;
}



#endif
