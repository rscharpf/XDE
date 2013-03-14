#ifndef REPORTNU_H
#define REPORTNU_H

#include "Report.h"


class ReportNu : public Report
{
 public:

  ReportNu(const string &filename);
  ReportNu(double *value);
  ~ReportNu(void);

  void report(const Structure *str);

 private:
  int writeToFile;
  double *value;
  int nr;
};



inline ReportNu::ReportNu(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportNu::ReportNu(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportNu::~ReportNu(void)
{
  return;
}


inline void ReportNu::report(const Structure *str)
{
  if (writeToFile)
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  out << str->nu[q][g] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  {
	    value[nr] = str->nu[q][g];
	    nr++;
	  }
    }

      
  return;
}


#endif
