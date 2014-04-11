#ifndef REPORTT_H
#define REPORTT_H

#include "Report.h"


class ReportT : public Report
{
 public:

  ReportT(const string &filename);
  ReportT(double *value);
  ~ReportT(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportT::ReportT(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportT::ReportT(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportT::~ReportT(void)
{
  return;
}


inline void ReportT::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->t[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->t[q];
	  nr++;
	}
    }

  return;
}


#endif
