#ifndef UPDATESIGMA2MH_H
#define UPDATESIGMA2MH_H

#include "UpdateMultiplicativePositive.h"

#include "PotentialSigma2qg.h"
#include "PotentialXqg.h"
#include "PotentialNug.h"
#include "PotentialDDeltag.h"


class UpdateSigma2MH : public Update
{
 public:

  UpdateSigma2MH(Structure *str,double epsilon);
  ~UpdateSigma2MH(void);

  int update(Random &ran);
  Update *copy(void) const;
  void setEpsilon(double epsilon);

 private:
  Structure *str;
  
  vector<Update *> up;
};



inline UpdateSigma2MH::UpdateSigma2MH(Structure *str,double epsilon) : Update(epsilon)
{
  this->str = str;
  
  int q,g;
  for (q = 0; q < str->Q; q++)
    for (g = 0; g < str->G; g++)
      {
  	vector<Potential *> term;
	term.push_back(new PotentialSigma2qg(q,g,str));
	term.push_back(new PotentialXqg(q,g,str));
	term.push_back(new PotentialNug(g,str));
	term.push_back(new PotentialDDeltag(g,str));
	PotentialSum model(term);
	
	up.push_back(new UpdateMultiplicativePositive(&model,&(str->sigma2[q][g]),epsilon));

	int i;
	for (i = 0; i < term.size(); i++)
	  delete term[i];
      }
  
  return;
}



inline UpdateSigma2MH::~UpdateSigma2MH(void)
{
  int i;
  for (i = 0; i < up.size(); i++)
    delete up[i];

  return;
}


inline void UpdateSigma2MH::setEpsilon(double epsilon)
{
  int i;
  for (i = 0; i < up.size(); i++)
    up[i]->setEpsilon(epsilon);

  this->epsilon = epsilon;

  return;
}



inline Update *UpdateSigma2MH::copy(void) const
{
  Update *u = new UpdateSigma2MH(str,epsilon);

  return u;
}




inline int UpdateSigma2MH::update(Random &ran)
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
