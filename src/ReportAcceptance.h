#ifndef REPORTACCEPTANCE_H
#define REPORTACCEPTANCE_H

#include "Report.h"
#include "Update.h"

class ReportAcceptance : public Report
{
 public:

  ReportAcceptance(const string &filename,const vector<Update *> &update);
  ReportAcceptance(double *value,const vector<Update *> &update);
  ~ReportAcceptance(void);

  void report(const Structure *str);

 private:
  int writeToFile;
  double *value;
  int nr;
  vector<Update *> update;
};



inline ReportAcceptance::ReportAcceptance(const string &filename,const vector<Update *> &update) : Report(filename)
{
  writeToFile = 1;
  
  this->update.resize(update.size());
  int i;
  for (i = 0; i < update.size(); i++)
    this->update[i] = update[i];

  return;
}



inline ReportAcceptance::ReportAcceptance(double *value,const vector<Update *> &update) : Report()
{
  writeToFile = 0;
  this->value = value;
  nr = 0;

  this->update.resize(update.size());
  int i;
  for (i = 0; i < update.size(); i++)
    this->update[i] = update[i];

  return;
}



inline ReportAcceptance::~ReportAcceptance(void)
{
  return;
}


inline void ReportAcceptance::report(const Structure *str)
{
  if (writeToFile)
    {
      int i;
      for (i = 0; i < update.size(); i++)
	out << update[i]->acceptRate() << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int i;
      for (i = 0; i < update.size(); i++)
	{
	  value[nr] = update[i]->acceptRate();
	  nr++;
	}
    }



  return;
}


#endif
