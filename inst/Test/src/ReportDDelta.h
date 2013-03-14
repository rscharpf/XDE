#ifndef REPORTDDELTA_H
#define REPORTDDELTA_H

#include "Report.h"


class ReportDDelta : public Report
{
 public:

  ReportDDelta(const string &filename);
  ReportDDelta(double *value);
  ~ReportDDelta(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportDDelta::ReportDDelta(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportDDelta::ReportDDelta(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportDDelta::~ReportDDelta(void)
{
  return;
}


inline void ReportDDelta::report(const Structure *str)
{
  if (writeToFile)
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  out << str->Delta[q][g] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  {
	    value[nr] = str->Delta[q][g];
	    nr++;
	  }
    }

  return;
}


#endif
