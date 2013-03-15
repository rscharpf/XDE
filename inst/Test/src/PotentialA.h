#ifndef POTENTIALA_H
#define POTENTIALA_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialA : public Potential
{
 public:

  PotentialA(const Structure *str);
  virtual ~PotentialA(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
};



inline PotentialA::PotentialA(const Structure *str) : Potential()
{
  this->str = str;

  return;
}



inline PotentialA::~PotentialA(void)
{
  return;
}


inline Potential *PotentialA::copy(void) const
{
  Potential *pp = new PotentialA(str);

  return pp;
}


inline double PotentialA::potential(Random &ran) const
{
  double pot = 0.0;

  int q;
  for (q = 0; q < str->Q; q++)
    {
      if (str->a[q] == 0.0)
	pot += - log(str->pA0);
      else if (str->a[q] == 1.0)
	pot += - log(str->pA1);
      else
	{
	  pot += - log(1.0 - str->pA0 - str->pA1);
	  pot += ran.PotentialBeta(str->alphaA,str->betaA,str->a[q]);
	}
    }

  return pot;
}



#endif
