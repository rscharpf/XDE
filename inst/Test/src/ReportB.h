#ifndef REPORTB_H
#define REPORTB_H

#include "Report.h"


class ReportB : public Report
{
 public:

  ReportB(const string &filename);
  ReportB(double *value);
  ~ReportB(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportB::ReportB(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportB::ReportB(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportB::~ReportB(void)
{
  return;
}


inline void ReportB::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->b[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->b[q];
	  nr++;
	}
    }

  return;
}


#endif
