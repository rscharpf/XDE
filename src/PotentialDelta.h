#ifndef POTENTIALDELTA_H
#define POTENTIALDELTA_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialDelta : public Potential
{
 public:

  PotentialDelta(const Structure *str);
  virtual ~PotentialDelta(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
};



inline PotentialDelta::PotentialDelta(const Structure *str) : Potential()
{
  this->str = str;

  return;
}



inline PotentialDelta::~PotentialDelta(void)
{
  return;
}


inline Potential *PotentialDelta::copy(void) const
{
  Potential *pp = new PotentialDelta(str);

  return pp;
}


inline double PotentialDelta::potential(Random &ran) const
{
  double pot = 0.0;

  int g;
  for (g = 0; g < str->G; g++)
    {
      if (str->delta[g] == 1)
	pot += - log(str->xi);
      else
	pot += - log(1.0 - str->xi);
    }
  
  return pot;
}



#endif
