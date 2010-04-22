#ifndef REPORTPOTENTIAL_H
#define REPORTPOTENTIAL_H

#include "Report.h"
#include "Potential.h"

class ReportPotential : public Report
{
 public:

  ReportPotential(const string &filename,const vector<Potential *> &model);
  ReportPotential(double *value,const vector<Potential *> &model);
  ~ReportPotential(void);

  void report(const Structure *str);

 private:
  int writeToFile;
  double *value;
  int nr;
  vector<Potential *> model;
};



inline ReportPotential::ReportPotential(const string &filename,const vector<Potential *> &model) : Report(filename)
{
  writeToFile = 1;

  this->model.resize(model.size());
  int i;
  for (i = 0; i < model.size(); i++)
    this->model[i] = model[i];

  return;
}



inline ReportPotential::ReportPotential(double *value,const vector<Potential *> &model) : Report()
{
  writeToFile = 0;
  this->value = value;
  nr = 0;

  this->model.resize(model.size());
  int i;
  for (i = 0; i < model.size(); i++)
    this->model[i] = model[i];

  return;
}



inline ReportPotential::~ReportPotential(void)
{
  return;
}


inline void ReportPotential::report(const Structure *str)
{
  Random ran(1);

  if (writeToFile)
    {
      int i;
      for (i = 0; i < model.size(); i++)
	out << model[i]->potential(ran) << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int i;
      for (i = 0; i < model.size(); i++)
	{
	  value[nr] = model[i]->potential(ran);
	  nr++;
	}
    }

  return;
}


#endif
