#ifndef REPORTPROBDELTA_H
#define REPORTPROBDELTA_H

#include "Structure.h"
#include "Report.h"
#include "PotentialDelta.h"
#include "PotentialDDeltag.h"
#include "PotentialXqg.h"
#include "PotentialSum.h"

class ReportProbDelta : public Report
{
 public:

  ReportProbDelta(const string &filename,Structure *str,int oneDelta);
  ReportProbDelta(double *value,Structure *str,int oneDelta);
  ~ReportProbDelta(void);

  void report(const Structure *str);

 private:
  int writeToFile;
  
  Structure *str;
  vector<Potential *> model;

  double *value;
  int nr;
};



inline ReportProbDelta::ReportProbDelta(const string &filename,Structure *str,int oneDelta) : Report(filename)
{
  writeToFile = 1;

  this->str = str;

  model.resize(0);
  int g;
  for (g = 0; g < str->G; g++)
    {
      vector <Potential *> term;
      term.push_back(new PotentialDelta(str,oneDelta));
      term.push_back(new PotentialDDeltag(g,str));
      int q;
      for (q = 0; q < str->Q; q++)
	term.push_back(new PotentialXqg(q,g,str));

      model.push_back(new PotentialSum(term));

      int i;
      for (i = 0; i < term.size(); i++)
	delete term[i];
    }
  
  return;
}



inline ReportProbDelta::ReportProbDelta(double *value,Structure *str,int oneDelta) : Report()
{
  writeToFile = 0;

  this->str = str;

  this->value = value;
  nr = 0;

  model.resize(0);
  int g;
  for (g = 0; g < str->G; g++)
    {
      vector <Potential *> term;
      term.push_back(new PotentialDelta(str,oneDelta));
      term.push_back(new PotentialDDeltag(g,str));
      int q;
      for (q = 0; q < str->Q; q++)
	term.push_back(new PotentialXqg(q,g,str));

      model.push_back(new PotentialSum(term));

      int i;
      for (i = 0; i < term.size(); i++)
	delete term[i];
    }

  return;
}


inline ReportProbDelta::~ReportProbDelta(void)
{
  int i;
  for (i = 0; i < model.size(); i++)
    delete model[i];

  return;
}


inline void ReportProbDelta::report(const Structure *str)
{
  double prob = 0.0;
  double pot0,pot1;

  Random ran(1);

  int q,g;
  for (g = 0; g < str->G; g++)
    for (q = 0; q < str->Q; q++)
      {
	int oldValue = str->delta[q][g];
	
	this->str->delta[q][g] = 0;
	pot0 = model[g]->potential(ran);
	this->str->delta[q][g] = 1;
	pot1 = model[g]->potential(ran);
	
	double minPot = pot0 < pot1 ? pot0 : pot1;
	pot0 -= minPot;
	pot1 -= minPot;
	
	prob = exp(- pot1) / (exp(- pot0) + exp(- pot1));
	
	Structure *ss = (Structure *) str;
	ss->delta[q][g] = oldValue;
	
	if (writeToFile)
	  out << prob << " ";
	else
	  {
	    value[nr] = prob;
	    nr++;
	  }
      }
  
  if (writeToFile)
    {
      out << "\n";
      out.flush();
    }
  
  
  return;
}


#endif
