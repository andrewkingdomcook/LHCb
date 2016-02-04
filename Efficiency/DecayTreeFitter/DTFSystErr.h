#ifndef DTFSystErr_h
#define DTFSystErr_h

#include <sstream>
#include <iostream>
#include "TString.h"
#include "TH3.h"
#include "TTree.h"

class DTFSystErr{
 public :

  enum Quarks{up = 1,
	      down,
	      strange,
	      charm,
	      bottom,
	      top};


  double CalcErrMC(double, double);
  double CalcErrData(double, double);
  void Fit(TString,TString,TString,TString);
  void Optimise();
  void CalculateEfficiency();
  double FitDCB(TTree *, TString, double);
  bool MCTruthTopAncestor(Int_t );
  void CalculateEfficiencyMC(TString);
  void CompareDoubleDTF(TString, TString);
  void CalculateEfficiencyDoubleJpsiMC(TString);
 private:

};


#endif
