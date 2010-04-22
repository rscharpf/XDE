#ifndef POTENTIALSUM_H
#define POTENTIALSUM_H

#include "Potential.h"
#include "Structure.h"
#include "Random.h"


class PotentialSum : public Potential
{
 public:

  PotentialSum(const vector<Potential *> &term);
  virtual ~PotentialSum(void);

  double potential(Random &ran) const;
  Potential *copy(void) const;

 private:
  vector<Potential *> term;
};



inline PotentialSum::PotentialSum(const vector<Potential *> &term) : Potential()
{
  this->term.resize(term.size());

  int i;
  for (i = 0; i < term.size(); i++)
    this->term[i] = term[i]->copy();

  return;
}



inline PotentialSum::~PotentialSum(void)
{
  int i;
  for (i = 0; i < term.size(); i++)
    delete term[i];
  
  return;
}


inline Potential *PotentialSum::copy(void) const
{
  Potential *pp = new PotentialSum(term);

  return pp;
}


inline double PotentialSum::potential(Random &ran) const
{
  double pot = 0.0;

  int i;
  for (i = 0; i < term.size(); i++)
    pot += term[i]->potential(ran);

  return pot;
}



#endif
