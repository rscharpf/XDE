#ifndef REPORTPHI_H
#define REPORTPHI_H

#include "Report.h"


class ReportPhi : public Report
{
 public:

  ReportPhi(const string &filename);
  ReportPhi(double *value);
  ~ReportPhi(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportPhi::ReportPhi(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportPhi::ReportPhi(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportPhi::~ReportPhi(void)
{
  return;
}


inline void ReportPhi::report(const Structure *str)
{
  if (writeToFile)
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  out << str->phi[q][g] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  {
	    value[nr] = str->phi[q][g];
	    nr++;
	  }
    }

  return;
}


#endif
