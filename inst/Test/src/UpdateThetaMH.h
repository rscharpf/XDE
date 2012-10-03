#ifndef UPDATETHETAMH_H
#define UPDATETHETAMH_H

#include "UpdateMultiplicativePositive.h"

#include "PotentialPhiqg.h"


class UpdateThetaMH : public Update
{
 public:

  UpdateThetaMH(Structure *str,const Potential *modelTheta,double epsilon);
  ~UpdateThetaMH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Potential *modelTheta;

  Structure *str;
  
  vector<Update *> up;
};



inline UpdateThetaMH::UpdateThetaMH(Structure *str,const Potential *modelTheta,
				    double epsilon) : Update(epsilon)
{
  this->modelTheta = modelTheta->copy();

  this->str = str;
  
  int q;
  for (q = 0; q < str->Q; q++)
    {
      vector<Potential *> term;
      term.push_back(modelTheta->copy());
      int g;
      for (g = 0; g < str->G; g++)
	term.push_back(new PotentialPhiqg(q,g,str));
      PotentialSum model(term);
	
      up.push_back(new UpdateMultiplicativePositive(&model,&(str->theta[q]),epsilon));

      int i;
      for (i = 0; i < term.size(); i++)
	delete term[i];
    }
  
  return;
}



inline UpdateThetaMH::~UpdateThetaMH(void)
{
  int i;
  for (i = 0; i < up.size(); i++)
    delete up[i];

  delete modelTheta;

  return;
}


inline Update *UpdateThetaMH::copy(void) const
{
  Update *u = new UpdateThetaMH(str,modelTheta,epsilon);

  return u;
}




inline int UpdateThetaMH::update(Random &ran)
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
