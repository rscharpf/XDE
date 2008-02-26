#ifndef POTENTIALPHI_H
#define POTENTIALPHI_H

#include "Potential.h"
#include "PotentialPhiqg.h"
#include "Structure.h"
#include "Random.h"


class PotentialPhi : public Potential
{
 public:

  PotentialPhi(const Structure *str);
  virtual ~PotentialPhi(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  const Structure *str;
  vector<Potential *> model;
};



inline PotentialPhi::PotentialPhi(const Structure *str) : Potential()
{
  this->str = str;

  int q,g;
  for (q = 0; q < str->Q; q++)
    for (g = 0; g < str->G; g++)
      model.push_back(new PotentialPhiqg(q,g,str));

  return;
}



inline PotentialPhi::~PotentialPhi(void)
{
  int i;
  for (i = 0; i < model.size(); i++)
    delete model[i];

  return;
}


inline Potential *PotentialPhi::copy(void) const
{
  Potential *pp = new PotentialPhi(str);

  return pp;
}


inline double PotentialPhi::potential(Random &ran) const
{
  double pot = 0.0;

  int i;
  for (i = 0; i < model.size(); i++)
    pot += model[i]->potential(ran);

  return pot;
}



#endif
