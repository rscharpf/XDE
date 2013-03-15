#ifndef REPORTLAMBDA_H
#define REPORTLAMBDA_H

#include "Report.h"


class ReportLambda : public Report
{
 public:

  ReportLambda(const string &filename);
  ReportLambda(double *value);
  ~ReportLambda(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportLambda::ReportLambda(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportLambda::ReportLambda(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportLambda::~ReportLambda(void)
{
  return;
}


inline void ReportLambda::report(const Structure *str)
{
  if (writeToFile)
    {
      int q;
      for (q = 0; q < str->Q; q++)
	out << str->lambda[q] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int q;
      for (q = 0; q < str->Q; q++)
	{
	  value[nr] = str->lambda[q];
	  nr++;
	}
    }

  return;
}


#endif
