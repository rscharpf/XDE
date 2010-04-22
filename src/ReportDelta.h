#ifndef REPORTDELTA_H
#define REPORTDELTA_H

#include "Report.h"


class ReportDelta : public Report
{
 public:

  ReportDelta(const string &filename);
  ReportDelta(int *value);
  ~ReportDelta(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  int *value;
  int nr;
};



inline ReportDelta::ReportDelta(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportDelta::ReportDelta(int *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportDelta::~ReportDelta(void)
{
  return;
}


inline void ReportDelta::report(const Structure *str)
{
  if (writeToFile)
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  out << str->delta[q][g] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q,g;
      for (g = 0; g < str->G; g++)
	for (q = 0; q < str->Q; q++)
	  {
	    value[nr] = str->delta[q][g];
	    nr++;
	  }      
    }

  return;
}


#endif
