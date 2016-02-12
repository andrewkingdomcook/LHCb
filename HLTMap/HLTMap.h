#ifndef HLTMap_h
#define HLTMap_h

#include <sstream>
#include <iostream>
#include "TString.h"
#include "TH3.h"
#include "TH2.h"



class HLTMap{



 public :
  void MakePhotonMap();
  void TrackPhotonDistance();
  void MakeTrackMap();
  void AreaMap();
  void ScanMap();
 
  

 private:
};


#endif
