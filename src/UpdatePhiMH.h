#ifndef UPDATEPHIMH_H
#define UPDATEPHIMH_H

#include "UpdateMultiplicativePositive.h"

#include "PotentialPhiqg.h"
#include "PotentialXqg.h"


class UpdatePhiMH : public Update
{
 public:

  UpdatePhiMH(Structure *str,double epsilon);
  ~UpdatePhiMH(void);

  int update(Random &ran);
  Update *copy(void) const;
  void setEpsilon(double epsilon);

 private:
  Structure *str;
  
  vector<Update *> up;
};



inline UpdatePhiMH::UpdatePhiMH(Structure *str,double epsilon) : Update(epsilon)
{
  this->str = str;
  
  int q,g;
  for (q = 0; q < str->Q; q++)
    for (g = 0; g < str->G; g++)
      {
  	vector<Potential *> term;
	term.push_back(new PotentialPhiqg(q,g,str));
	term.push_back(new PotentialXqg(q,g,str));
	PotentialSum model(term);
	
	up.push_back(new UpdateMultiplicativePositive(&model,&(str->phi[q][g]),epsilon));

	int i;
	for (i = 0; i < term.size(); i++)
	  delete term[i];
      }
  
  return;
}



inline UpdatePhiMH::~UpdatePhiMH(void)
{
  int i;
  for (i = 0; i < up.size(); i++)
    delete up[i];

  return;
}




inline void UpdatePhiMH::setEpsilon(double epsilon)
{
  int i;
  for (i = 0; i < up.size(); i++)
    up[i]->setEpsilon(epsilon);

  this->epsilon = epsilon;

  return;
}





inline Update *UpdatePhiMH::copy(void) const
{
  Update *u = new UpdatePhiMH(str,epsilon);

  return u;
}




inline int UpdatePhiMH::update(Random &ran)
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
