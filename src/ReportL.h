#ifndef REPORTL_H
#define REPORTL_H

#include "Report.h"


class ReportL : public Report
{
 public:

  ReportL(const string &filename);
  ReportL(double *value);
  ~ReportL(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportL::ReportL(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportL::ReportL(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportL::~ReportL(void)
{
  return;
}


inline void ReportL::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->l[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->l[q];
	  nr++;
	}
    }

  return;
}


#endif
