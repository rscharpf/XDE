#ifndef UPDATELMH_H
#define UPDATELMH_H

#include "UpdateMultiplicativePositive.h"

#include "PotentialSigma2qg.h"


class UpdateLMH : public Update
{
 public:

  UpdateLMH(Structure *str,const Potential *modelL,double epsilon);
  ~UpdateLMH(void);

  int update(Random &ran);
  Update *copy(void) const;
  void setEpsilon(double epsilon);

 private:
  Potential *modelL;

  Structure *str;
  
  vector<Update *> up;
};



inline UpdateLMH::UpdateLMH(Structure *str,const Potential *modelL,double epsilon) : Update(epsilon)
{
  this->modelL = modelL->copy();

  this->str = str;
  
  int q;
  for (q = 0; q < str->Q; q++)
    {
      vector<Potential *> term;
      term.push_back(modelL->copy());
      int g;
      for (g = 0; g < str->G; g++)
	term.push_back(new PotentialSigma2qg(q,g,str));
      PotentialSum model(term);
	
      up.push_back(new UpdateMultiplicativePositive(&model,&(str->l[q]),epsilon));

      int i;
      for (i = 0; i < term.size(); i++)
	delete term[i];
    }
  
  return;
}



inline UpdateLMH::~UpdateLMH(void)
{
  int i;
  for (i = 0; i < up.size(); i++)
    delete up[i];

  delete modelL;

  return;
}


inline void UpdateLMH::setEpsilon(double epsilon)
{
  int i;
  for (i = 0; i < up.size(); i++)
    up[i]->setEpsilon(epsilon);
  
  this->epsilon = epsilon;

  return;
}




inline Update *UpdateLMH::copy(void) const
{
  Update *u = new UpdateLMH(str,modelL,epsilon);

  return u;
}




inline int UpdateLMH::update(Random &ran)
{
  int nAccept = 0;

  int i;
  for (i = 0; i < up.size(); i++)
    {
      addTry();
      int acc = up[i]->update(ran);
      if (acc) addAccept();
      nAccept += acc;
    }

  return nAccept;
}


#endif
