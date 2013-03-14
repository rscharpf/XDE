#ifndef REPORTR_H
#define REPORTR_H

#include "Report.h"


class ReportR : public Report
{
 public:

  ReportR(const string &filename);
  ReportR(double *value);
  ~ReportR(void);

  void report(const Structure *str);

 private:
  int writeToFile;

  double *value;
  int nr;
};



inline ReportR::ReportR(const string &filename) : Report(filename)
{
  writeToFile = 1;

  return;
}


inline ReportR::ReportR(double *value) : Report()
{
  writeToFile = 0;

  this->value = value;
  nr = 0;

  return;
}



inline ReportR::~ReportR(void)
{
  return;
}


inline void ReportR::report(const Structure *str)
{
  if (writeToFile)
    {
      int i,j;
      for (i = 0; i < str->Q; i++)
	for (j = i+1; j < str->Q; j++)
	  out << str->r[i][j] << " ";
      
      out << "\n";
      out.flush();
    }
  else
    {
      int i,j;
      for (i = 0; i < str->Q; i++)
	for (j = i+1; j < str->Q; j++)
	  {
	    value[nr] = str->r[i][j];
	    nr++;
	  }
    }

  return;
}


#endif
