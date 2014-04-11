#ifndef UPDATELAMBDAMH_H
#define UPDATELAMBDAMH_H

#include "UpdateMultiplicativePositive.h"

#include "PotentialPhiqg.h"

class UpdateLambdaMH : public Update
{
 public:

  UpdateLambdaMH(Structure *str,const Potential *modelLambda,double epsilon);
  ~UpdateLambdaMH(void);

  int update(Random &ran);
  Update *copy(void) const;

 private:
  Potential *modelLambda;

  Structure *str;
  
  vector<Update *> up;
};



inline UpdateLambdaMH::UpdateLambdaMH(Structure *str,const Potential *modelLambda,double epsilon) : Update(epsilon)
{
  this->modelLambda = modelLambda->copy();

  this->str = str;
  
  int q;
  for (q = 0; q < str->Q; q++)
    {
      vector<Potential *> term;
      term.push_back(modelLambda->copy());
      int g;
      for (g = 0; g < str->G; g++)
	term.push_back(new PotentialPhiqg(q,g,str));
      PotentialSum model(term);
	
      up.push_back(new UpdateMultiplicativePositive(&model,&(str->lambda[q]),epsilon));

      int i;
      for (i = 0; i < term.size(); i++)
	delete term[i];
    }
  
  return;
}



inline UpdateLambdaMH::~UpdateLambdaMH(void)
{
  int i;
  for (i = 0; i < up.size(); i++)
    delete up[i];

  delete modelLambda;

  return;
}


inline Update *UpdateLambdaMH::copy(void) const
{
  Update *u = new UpdateLambdaMH(str,modelLambda,epsilon);

  return u;
}




inline int UpdateLambdaMH::update(Random &ran)
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
