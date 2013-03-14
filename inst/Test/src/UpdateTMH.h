#ifndef UPDATETMH_H
#define UPDATETMH_H

#include "UpdateMultiplicativePositive.h"

#include "PotentialSigma2qg.h"


class UpdateTMH : public Update
{
 public:

  UpdateTMH(Structure *str,const Potential *modelT,double epsilon);
  ~UpdateTMH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Potential *modelT;

  Structure *str;
  
  vector<Update *> up;
};



inline UpdateTMH::UpdateTMH(Structure *str,const Potential *modelT,double epsilon) : Update(epsilon)
{
  this->modelT = modelT->copy();

  this->str = str;
  
  int q;
  for (q = 0; q < str->Q; q++)
    {
      vector<Potential *> term;
      term.push_back(modelT->copy());
      int g;
      for (g = 0; g < str->G; g++)
	term.push_back(new PotentialSigma2qg(q,g,str));
      PotentialSum model(term);
	
      up.push_back(new UpdateMultiplicativePositive(&model,&(str->t[q]),epsilon));

      int i;
      for (i = 0; i < term.size(); i++)
	delete term[i];
    }
  
  return;
}



inline UpdateTMH::~UpdateTMH(void)
{
  int i;
  for (i = 0; i < up.size(); i++)
    delete up[i];

  delete modelT;

  return;
}


inline Update *UpdateTMH::copy(void) const
{
  Update *u = new UpdateTMH(str,modelT,epsilon);

  return u;
}




inline int UpdateTMH::update(Random &ran)
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
