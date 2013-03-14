#ifndef REPORTXI_H
#define REPORTXI_H

#include "Report.h"


class ReportXi : public Report
{
 public:

  ReportXi(const string &filename);
  ReportXi(double *value);
  ~ReportXi(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportXi::ReportXi(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}



inline ReportXi::ReportXi(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}


inline ReportXi::~ReportXi(void)
{
  return;
}


inline void ReportXi::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->xi[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->xi[q];
	  nr++;
	}
    }

  return;
}


#endif
