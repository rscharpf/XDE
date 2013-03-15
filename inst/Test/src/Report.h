#ifndef REPORT_H
#define REPORT_H

#include "Structure.h"


class Report
{
 public:
  
  Report(const string &filename);
  Report(void) {file = 0; return;};
  virtual ~Report(void);

  virtual void report(const Structure *str) = 0;

 protected:
  int file;
  ofstream out;
};



inline Report::Report(const string &filename)
{
  file = 1;
  
  out.open(filename.c_str());
  if (out.fail())
    {
      cout << "ERROR: Unable to open file " << filename << ". Aborting.\n\n";
      exit(-1);
    }

  return;
}



inline Report::~Report(void)
{
  if (file) out.close();

  return;
}



#endif
