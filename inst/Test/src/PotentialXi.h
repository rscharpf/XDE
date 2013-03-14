#ifndef POTENTIALXI_H
#define POTENTIALXI_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialXi : public Potential
{
 public:

  PotentialXi(const Structure *str,int oneDelta);
  virtual ~PotentialXi(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
  int oneDelta;
};



inline PotentialXi::PotentialXi(const Structure *str,int oneDelta) : Potential()
{
  this->str = str;
  this->oneDelta = oneDelta;

  return;
}



inline PotentialXi::~PotentialXi(void)
{
  return;
}


inline Potential *PotentialXi::copy(void) const
{
  Potential *pp = new PotentialXi(str,oneDelta);

  return pp;
}


inline double PotentialXi::potential(Random &ran) const
{
  double pot = 0.0;

  if (oneDelta == 0)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	pot += ran.PotentialBeta(str->alphaXi,str->betaXi,str->xi[q]);
    }
  else
    {
      pot += ran.PotentialBeta(str->alphaXi,str->betaXi,str->xi[0]);
    }
  
  return pot;
}



#endif
