#ifndef GeomAcc_h
#define GeomAcc_h

#include <sstream>
#include <iostream>
#include "TString.h"
#include "TH3.h"
#include "TH2.h"

class GeomAcc{


 public :
 
  void Calculate(TString, TString);
  bool Cut(double, double);
  TH2D Eff(TH2D* , TH2D*, int, int, TString );
  void Format(TString, TString);
 private:

};


#endif
