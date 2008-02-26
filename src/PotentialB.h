#ifndef POTENTIALB_H
#define POTENTIALB_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialB : public Potential
{
 public:

  PotentialB(const Structure *str);
  virtual ~PotentialB(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
};



inline PotentialB::PotentialB(const Structure *str) : Potential()
{
  this->str = str;

  return;
}



inline PotentialB::~PotentialB(void)
{
  return;
}


inline Potential *PotentialB::copy(void) const
{
  Potential *pp = new PotentialB(str);

  return pp;
}


inline double PotentialB::potential(Random &ran) const
{
  double pot = 0.0;

  int q;
  for (q = 0; q < str->Q; q++)
    {
      if (str->b[q] == 0.0)
	pot += - log(str->pB0);
      else if (str->b[q] == 1.0)
	pot += - log(str->pB1);
      else
	{
	  pot += - log(1.0 - str->pB0 - str->pB1);
	  pot += ran.PotentialBeta(str->alphaB,str->betaB,str->b[q]);
	}
    }

  return pot;
}



#endif
