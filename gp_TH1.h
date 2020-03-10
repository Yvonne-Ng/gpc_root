#ifndef GP_TH1_H
#define GP_TH1_H

//including all dependence oof gp.h
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "ndlexceptions.h"
#include "ndlstrutil.h"
#include "CMatrix.h"
#include "CKern.h"
#include "CGp.h"
//including also gp_config
#include "gp_config.h"
#include "TH1.h"

//Do not need to be an executable. just need to be a libray with functions
//int main(int argc, char* argv[]);

class gp_TH1 {
 public:
  gp_TH1(gp_config* config);
  void readFromTH1(TH1* hist);
  void learn();
  void relearn();

  void draw_root();
  //void display();
  //void gnuplot();

  //void helpInfo();
  //void helpHeader();
};

#else /* GP_H */
#endif /* GP_H */

