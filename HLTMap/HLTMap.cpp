
#define HLTMap_cpp
#include "HLTMap.h"
#include "LHCbStyle.h"

#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TKey.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TObject.h"
#include <TMath.h>
#include <cmath>
#include <sstream>
#include "RooFit.h"
#include "TString.h"
#include "TObjArray.h"
#include "RooRealVar.h"

#include "TROOT.h"
#include "TEnv.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TBrowser.h"
#include "THStack.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TKey.h"
#include "TRandom.h"
#include "TLegend.h"
#include "TMath.h"
#include "TCanvas.h"
#include <boost/lexical_cast.hpp>
#include "limits.h"
#include "math.h"

using             boost::lexical_cast;

int main(){
  HLTMap HLTClass;

  HLTClass.ScanMap();
  // HLTClass.AreaMap();
  // HLTClass.MakeTrackMap();
  // HLTClass.TrackPhotonDistance(); 
 
  return 0;
}



//*************************************************************************************************************
//
//Scans over Mirror area with different binning, creates efficiency plot for photon hits in a specific mirror pair using
//Track position in X, Y and Phi distribution efficiency.
//
//*************************************************************************************************************
void HLTMap::ScanMap()
{
  gStyle->SetOptStat(0);
  TFile *f =
    TFile::Open("/afs/cern.ch/work/p/phadc/private/gangadir/workspace/phadc/LocalXML/123/output/rich.tuples.root");
  TTree * treePhoton = (TTree*)f->Get("RICH/RichAlignMoniR2Gas/RichAlign");

  UInt_t primeMirrNum ;
  UInt_t secMirrNum ;
  double primeMirrX, primeMirrY;
  double trackX, trackY;
  double deltaTheta, phi , momentum;
  Bool_t unAmbigious;
  int eventCount, segmentCount;

  TObjArray* histList = new TObjArray(0);

  //photon variables
  treePhoton->SetBranchAddress("sphMirror",&primeMirrNum); //Prime mirror photon hits
  treePhoton->SetBranchAddress("secMirror", &secMirrNum); // Secondary mirror photon hits
  treePhoton->SetBranchAddress("SphericalMirrX", &primeMirrX); //photon X co-ordinate in prime mirror plane
  treePhoton->SetBranchAddress("SphericalMirrY", &primeMirrY); //photon Y co-ordinate in prime mirror plane
  treePhoton->SetBranchAddress("unAmbiguousPhoton", &unAmbigious); 
  treePhoton->SetBranchAddress("deltaTheta",&deltaTheta); 
  treePhoton->SetBranchAddress("phiRec",&phi);
  //track variables
  treePhoton->SetBranchAddress("xExitpoint",&trackX); //track X co-ordinate in prime mirror plane
  treePhoton->SetBranchAddress("yExitpoint",&trackY); //track Y co-ordinate in prime mirror plane
  treePhoton->SetBranchAddress("eventCount",&eventCount); //event number 
  treePhoton->SetBranchAddress("segmentKey",&segmentCount);
  treePhoton->SetBranchAddress("momentum",&momentum); // track momentum

  bool fromTuple = false; //If get data from the tuple or from an already run .root file
  int numbMirrPair = 47; //number of mirror pairs to iterate over;default should be 47
  //photons in large tuple = 398027636 (aprox 4x10^8) (small tuple 1.2x10^7))
  //new tuple (cut at 30GeV):4.3 x10^8  Number over 40GeV: 2.36x10^8
  int tuplePhotons = treePhoton->GetEntries(); 
  int rebinXMin = 5; //set the minimum rebinning value that is *saved* - it still runs over the lowest
  int rebinYMin = 5; 
  double maxAllowedEvents = 90000000;//300000000; //max events to filter through that will be considered.
  int maxBinMerge = 50; // highest number of bin mergings 


  //
  //Setting number of bins in TH3 plots - this sets the lower limit for the bin merging later. 
  //
  int xBin = 120;
  int yBin = 120;

  //
  //Array to hold final results
  //
  double results[9][numbMirrPair];
  double highestEventsToFilter = 0; //highest number of events to filter to normalise rest.
 
  //
  //Left hand side Mirror pair variables (should have *either* left or right mirror pairs switched on)
  //
  
  //dimensions of mirror plots (xstart,xend,ystart,yend)
  int p00[4] ={-50,550,-2150,-850};//
  int p01[4] ={150,750,-2000,-920};
  int p02[4] ={580,1180,-2000,-920};//
  int p03[4] ={1000,2300,-1880,-720};//

  int p04[4] ={-50,550,-1080,480};
  int p05[4] ={400,1000,-1080,-480};
  int p06[4] ={800,1400,-1050,-450};
  int p07[4] ={1200,2650,-1600,-400};//

  int p08[4] ={-50,550,-680,-120};
  int p09[4] ={100,700,-680,-120};
  int p10[4] ={550,1150,-680,-120};
  int p11[4] ={1040,1680,-680,-40};//

  int p12[4] ={0,600,-300,300};
  int p13[4] ={400,1000,-350,250};
  int p14[4] ={800,1400,-320,280};
  int p15[4] ={1100,2700,-450,450};//

  int p16[4] ={-50,550,0,600};
  int p17[4] ={100,700,50,650};
  int p18[4] ={550,1150,100,700};
  int p19[4] ={1000,1700,0,700};//

  int p20[4] ={0,600,450,1050};
  int p21[4] ={400,1000,450,1050};
  int p22[4] ={800,1400,450,1050};
  int p23[4] ={1200,2500,270,1300};//

  int p24[4] ={-100,500,810,1770};//
  int p25[4] ={150,750,880,1800}; //
  int p26[4] ={600,1200,880,1740};//
  int p27[4] ={1000,2300,880,1740};//

  //primary mirror numbers
  unsigned int primeMirrArray[47]=
    {12,16,16,20,20,24, 8, 8, 4, 4, 0,12,17,17,21,21,25, 9, 9, 5, 5, 1,13,13,18,18,22,22,26,22,27,10,10, 6, 6, 2, 6, 3,14,14,19,19
     ,23,11,11, 7,15};
  //secondary mirror numbers  
  unsigned int secMirrArray[47]=
    { 8, 8,12,12,16,16, 8, 4, 4, 0, 0, 9, 9,13,13,17,17, 9, 5, 5, 1, 1, 9,10,10,14,14,18,18,19,19,10, 6, 6, 2, 2, 3, 3,10,11,11,15
     ,15,11, 7, 7,11};
 
  //plot dimensions for each mirror pair
  int xScanStart[47]={p12[0],p16[0],p16[0],p20[0],p20[0],p24[0],p08[0],p08[0],p04[0],p04[0],p00[0],p12[0],p17[0],p17[0],p21[0],
                       p21[0],p25[0],p09[0],p09[0],p05[0],p05[0],p01[0],p13[0],p13[0],p18[0],p18[0],p22[0],p22[0],p26[0],p22[0],
                       p27[0],p10[0],p10[0],p06[0],p06[0],p02[0],p06[0],p03[0],p14[0],p14[0],p19[0],p19[0],p23[0],p11[0],p11[0],
                       p07[0],p15[0]};

  int xScanEnd[47]={p12[1],p16[1],p16[1],p20[1],p20[1],p24[1],p08[1],p08[1],p04[1],p04[1],p00[1],p12[1],p17[1],p17[1],p21[1],
                     p21[1],p25[1],p09[1],p09[1],p05[1],p05[1],p01[1],p13[1],p13[1],p18[1],p18[1],p22[1],p22[1],p26[1],p22[1],
                     p27[1],p10[1],p10[1],p06[1],p06[1],p02[1],p06[1],p03[1],p14[1],p14[1],p19[1],p19[1],p23[1],p11[1],p11[1],
                     p07[1],p15[1]};

  int yScanStart[47]={p12[2],p16[2],p16[2],p20[2],p20[2],p24[2],p08[2],p08[2],p04[2],p04[2],p00[2],p12[2],p17[2],p17[2],p21[2],
                       p21[2],p25[2],p09[2],p09[2],p05[2],p05[2],p01[2],p13[2],p13[2],p18[2],p18[2],p22[2],p22[2],p26[2],p22[2],
                       p27[2],p10[2],p10[2],p06[2],p06[2],p02[2],p06[2],p03[2],p14[2],p14[2],p19[2],p19[2],p23[2],p11[2],p11[2],
                       p07[2],p15[2]};

  int yScanEnd[47]={p12[3],p16[3],p16[3],p20[3],p20[3],p24[3],p08[3],p08[3],p04[3],p04[3],p00[3],p12[3],p17[3],p17[3],p21[3],
                     p21[3],p25[3],p09[3],p09[3],p05[3],p05[3],p01[3],p13[3],p13[3],p18[3],p18[3],p22[3],p22[3],p26[3],p22[3],
                     p27[3],p10[3],p10[3],p06[3],p06[3],p02[3],p06[3],p03[3],p14[3],p14[3],p19[3],p19[3],p23[3],p11[3],p11[3],
                     p07[3],p15[3]};
                    
  //
  //Right hand side Mirror pair variables (should have *either* left or right mirror pairs switched on)
  //
  /*
  int p31[4] ={-550,50,-2150,-850};//
  int p30[4] ={-750,-150,-2000,-920};
  int p29[4] ={-1180,-580,-2000,-920};//
  int p28[4] ={-2300,-1000,-1880,-720};//

  int p35[4] ={-550,50,-1080,480};
  int p34[4] ={-1000,-400,-1080,-480};
  int p33[4] ={-1400,-800,-1050,-450};
  int p32[4] ={-2650,-1200,-1600,-400};//

  int p39[4] ={-550,50,-680,-120};
  int p38[4] ={-700,-100,-680,-120};
  int p37[4] ={-1150,-550,-680,-120};
  int p36[4] ={-1680,-1040,-680,-40};//

  int p43[4] ={-600,0,-300,300};
  int p42[4] ={-1000,-400,-350,250};
  int p41[4] ={-1400,-800,-320,280};
  int p40[4] ={-2700,-1100,-450,450};//

  int p47[4] ={-550,50,0,600};
  int p46[4] ={-700,-100,50,650};
  int p45[4] ={-1150,-550,100,700};
  int p44[4] ={-1700,-1000,0,700};//

  int p51[4] ={-600,0,450,1050};
  int p50[4] ={-1000,-400,450,1050};
  int p49[4] ={-1400,-800,450,1050};
  int p48[4] ={-2500,-1200,270,1300};//

  int p55[4] ={-500,100,810,1770};//
  int p54[4] ={-750,-150,880,1800}; //
  int p53[4] ={-1200,-600,880,1740};//
  int p52[4] ={-2300,-1000,880,1740};//

 //primary mirror numbers
  unsigned int primeMirrArray[47]=
    {43,47,47,51,51,55,39,39,35,35,31,43,46,46,50,50,54,38,38,34,34,30,42,42,45,45,49,49,53,49,52,37,37,33,33,29,33,28,41,41,44,44
     ,48,36,36,32,40};
  //secondary mirror numbers  
  unsigned int secMirrArray[47]=
    {31,31,35,35,39,39,31,27,27,23,23,30,30,34,34,38,38,30,26,26,22,22,30,29,29,33,33,37,37,36,36,29,25,25,21,21,20,20,29,28,28,32
     ,32,28,24,24,28};
 
  //plot dimensions for each mirror pair

  int xScanStart[47]={p43[0],p47[0],p47[0],p51[0],p51[0],p55[0],p39[0],p39[0],p35[0],p35[0],p31[0],p43[0],p46[0],p46[0],p50[0],
                      p50[0],p54[0],p38[0],p38[0],p34[0],p34[0],p30[0],p42[0],p42[0],p45[0],p45[0],p49[0],p49[0],p53[0],p49[0],
                      p52[0],p37[0],p37[0],p33[0],p33[0],p29[0],p33[0],p28[0],p41[0],p41[0],p44[0],p44[0],p48[0],p36[0],p36[0],
                      p32[0],p40[0]};
  
  int xScanEnd[47]={p43[1],p47[1],p47[1],p51[1],p51[1],p55[1],p39[1],p39[1],p35[1],p35[1],p31[1],p43[1],p46[1],p46[1],p50[1],
                    p50[1],p54[1],p38[1],p38[1],p34[1],p34[1],p30[1],p42[1],p42[1],p45[1],p45[1],p49[1],p49[1],p53[1],p49[1],
                    p52[1],p37[1],p37[1],p33[1],p33[1],p29[1],p33[1],p28[1],p41[1],p41[1],p44[1],p44[1],p48[1],p36[1],p36[1],
                    p32[1],p40[1]};

  int yScanStart[47]={p43[2],p47[2],p47[2],p51[2],p51[2],p55[2],p39[2],p39[2],p35[2],p35[2],p31[2],p43[2],p46[2],p46[2],p50[2],
                      p50[2],p54[2],p38[2],p38[2],p34[2],p34[2],p30[2],p42[2],p42[2],p45[2],p45[2],p49[2],p49[2],p53[2],p49[2],
                      p52[2],p37[2],p37[2],p33[2],p33[2],p29[2],p33[2],p28[2],p41[2],p41[2],p44[2],p44[2],p48[2],p36[2],p36[2],
                      p32[2],p40[2]};

  int yScanEnd[47]={p43[3],p47[3],p47[3],p51[3],p51[3],p55[3],p39[3],p39[3],p35[3],p35[3],p31[3],p43[3],p46[3],p46[3],p50[3],
                    p50[3],p54[3],p38[3],p38[3],p34[3],p34[3],p30[3],p42[3],p42[3],p45[3],p45[3],p49[3],p49[3],p53[3],p49[3],
                    p52[3],p37[3],p37[3],p33[3],p33[3],p29[3],p33[3],p28[3],p41[3],p41[3],p44[3],p44[3],p48[3],p36[3],p36[3],
                    p32[3],p40[3]};
    
  */ 

 
  //
  //Loop over mirror pairs
  //
  for(int j=0; j<numbMirrPair; j++){
  //for(int j=36; j<38; j++){
    std::cout  << "// PrimeMirror " << primeMirrArray[j]  << " secondaryMirror " << secMirrArray[j] <<  std::endl;
    double photonCount = 0;

    //
    //Variables to store highest efficiency bin paramenters when loop over different merged bins
    //
    double maxEff = 0;
    double maxEvents = 0;
    double maxXBinLow = 0;
    double maxXBinHigh = 0;
    double maxYBinLow = 0;
    double maxYBinHigh = 0;
    TH2D * maxTotalEfficiencyPlot;
    TH2D * maxTotalEventsToFilterPlot;

    double maxEffBin1 = 0;
    double maxEventsBin1 = 0; 
    double maxXBinLowBin1 = 0;
    double maxXBinHighBin1 = 0;
    double maxYBinLowBin1 = 0;
    double maxYBinHighBin1 = 0;
    TH2D * maxTotalEfficiencyPlotBin1;
    TH2D * maxTotalEventsToFilterPlotBin1;

    //
    //Set names and initialize TH3 plots which will hold track x,y and phi data
    //
    std::ostringstream scanAreaName, scanAreaMirrPairName;  
    scanAreaName << "scanArea" <<  primeMirrArray[j] << "Sec" <<secMirrArray[j];
    scanAreaMirrPairName << "scanAreaMirrPair" <<  primeMirrArray[j] << "Sec" <<secMirrArray[j];
   
    //should probably change these to TH3D if run entire tuple again!!
    TH3D * scanAreaPlotUnMerged  = new TH3D(scanAreaName.str().c_str(),scanAreaName.str().c_str(),
                       xBin,xScanStart[j],xScanEnd[j],yBin,yScanStart[j],yScanEnd[j],20, 0.0, 6.2827226);

    TH3D * scanAreaMirrPairPlotUnMerged  = new TH3D(scanAreaMirrPairName.str().c_str(),scanAreaMirrPairName.str().c_str(),
                       xBin,xScanStart[j],xScanEnd[j],yBin,yScanStart[j],yScanEnd[j],20, 0.0, 6.2827226);
  
    //
    //get data from a tuple or if ran this program once already can get TH3 files from an already processed root file.
    //
    if(fromTuple)
    {
      //
      //loop over photons in tuple to load up TH3 plots
      //
      for(int i=0; i<tuplePhotons; i++)
      {
        treePhoton->GetEntry(i);
        if(momentum <= 40.0)continue; //minumum momentum cut  (tuple is 30 GeV)     
        
        photonCount++;
       
        if((trackX < xScanStart[j]) ||  (trackX > xScanEnd[j])  )continue;
        if((trackY < yScanStart[j]) ||  (trackY > yScanEnd[j])  )continue;
      
        scanAreaPlotUnMerged->Fill(trackX,trackY,phi); //fill plot for denominator of efficiency

        if(primeMirrNum != primeMirrArray[j]) continue; //check photon in correct primary mirror
        if(secMirrNum != secMirrArray[j]) continue;  // check photon in correct secondary mirror combination
        if(unAmbigious == 0) continue; //remove ambigious photons 
        
        scanAreaMirrPairPlotUnMerged->Fill(trackX,trackY,phi); // fill plot for numerator of efficiency
      }
      
      histList->Add(scanAreaMirrPairPlotUnMerged);
      histList->Add(scanAreaPlotUnMerged);
    }
    else
    {
      photonCount = 236010000; //bit of a hack due to p cut being at 30 GeV, otherwise could just do tuple->numberofentries
      TFile *runFile = TFile::Open("FINAL/ScanMap2julyLeft.root"); 
      // TFile *runFile = TFile::Open("testFloatProblemLeft.root");
      scanAreaPlotUnMerged = (TH3D*) runFile->Get(scanAreaName.str().c_str()); 
      scanAreaMirrPairPlotUnMerged = (TH3D*) runFile->Get(scanAreaMirrPairName.str().c_str()); 
      // histList->Add(scanAreaMirrPairPlotUnMerged); //don't really need this if run from pre-made TH3D
      //histList->Add(scanAreaPlotUnMerged);
    }
    
    //
    //loop to vary size and shape of bins with variable bin merging
    //
    for(int rebinx=1; rebinx<=maxBinMerge; rebinx++)
    {
      for(int rebiny=1; rebiny<=maxBinMerge; rebiny++)
      {
        int xScanBin = xBin/rebinx;
        int yScanBin = yBin/rebiny;

        std::ostringstream effBinScanName, eventsBinScanName;  
        effBinScanName << "TotalEfficiency" <<primeMirrArray[j]<< "Sec" <<secMirrArray[j] << "("<< rebinx<<","<<rebiny<<")";
        eventsBinScanName << "EventstoFilter" <<primeMirrArray[j] << "Sec" <<secMirrArray[j] << "("<<rebinx<<","<<rebiny<<")";
        
        TH2D * recoEfficiencyPlot  = new TH2D(effBinScanName.str().c_str() , effBinScanName.str().c_str(),
                                               xScanBin,xScanStart[j],xScanEnd[j],yScanBin,yScanStart[j],yScanEnd[j]);
        
        TH2D * totalEventsToFilterPlot  = new TH2D(eventsBinScanName.str().c_str(), eventsBinScanName.str().c_str(),
                                                   xScanBin,xScanStart[j],xScanEnd[j],yScanBin,yScanStart[j],yScanEnd[j]);

        //
        //Merge bins of TH3 plots holding x,y and phi info to create efficiency maps
        //
        TH3D * scanAreaPlotMergedX = (TH3D*) scanAreaPlotUnMerged->RebinX(rebinx, "scanAreaPlotMergedX");
        TH3D * scanAreaPlot =  (TH3D*)scanAreaPlotMergedX->RebinY(rebiny, "scanAreaPlot");
        // scanAreaPlotMergedX->Delete();

        TH3D * scanAreaMirrPairPlotMergedX=(TH3D*) scanAreaMirrPairPlotUnMerged->RebinX(rebinx,"scanAreaMirrPairPlotMergedX");
        TH3D * scanAreaMirrPairPlot =  (TH3D*)scanAreaMirrPairPlotMergedX->RebinY(rebiny, "scanAreaMirrPairPlot");
        //scanAreaMirrPairPlotMergedX->Delete();
        
        //
        //loop over merged bins to find highest efficiency bin
        //
        for(int x=1; x<=xScanBin;x++)
        {
          for(int y=1; y<=yScanBin;y++)
          {

            //constants
            double photonsPerTrack = 60.0; //average number of photons per track
            double tracksPerEvent = 21.0; //average number of tracks per event -  NEED TO CHECK THIS IN NEW SAMPLE!
            double photonsPhiBin = 300.0; // number of photons needed in each phi bin
            double mometumCutFactor = 0.1; //factor to account for loss of events due to min momtum cut - at 10% > 40 Gev 

            double areaCount = 0;// photons hitting Primary mirror in specific bin
            double finalCount = 0; //Photons hitting mirror pair in specific bin
            int binSize[20];//Phi distribution in specific bin

            //
            //For each bin in X, Y: load phi distribution(Z) into 'binSize[]' and get number of photons hitting primary mirror bin
            //and photons hitting primary mirror pair bin.
            //
            for(int i=1; i<21;i++)
            {
              binSize[i-1] = scanAreaMirrPairPlot->GetBinContent(x,y,i); 
              areaCount = areaCount + scanAreaPlot->GetBinContent(x,y,i); 
              finalCount = finalCount + scanAreaMirrPairPlot->GetBinContent(x,y,i);
            }
            //
            //Extract 16th lowest phi bin, and total entries in phi dist, in order to work out how skewed it is. 
            //
            int elements = sizeof(binSize) / sizeof(binSize[0]); 
            std::sort(binSize, binSize + elements);

            double photonHits;
            double phiEfficiency;

            if(binSize[4] != 0)
            {
            photonHits =(photonsPhiBin/binSize[4])*finalCount;//photons required so 16 phi bins have at least 300 entries
            phiEfficiency = (photonsPhiBin*16)/photonHits; //efficiency caused by uneven phi distribution
            }
            else
            {
              photonHits = 0;
              phiEfficiency = 0;
            }
            
            //
            //Reconstruction efficiency
            //
            double areaEfficiency;
            if(areaCount != 0) areaEfficiency = finalCount/areaCount; //efficiency of mirror pair events hitting primary mirror
            else areaEfficiency = 0;            
            double recoEfficiency = areaEfficiency*phiEfficiency; //Reconstruction efficiency

            double proportion = finalCount/photonCount; //efficiency to filter correct mirror pair in selected area from total

            double sortingTracks; 
            double sortingEvents;            
            if((recoEfficiency)!= 0 && !(std::isnan(recoEfficiency))&& !(isinf(float (recoEfficiency))) && proportion != 0) 
            {
              sortingTracks = (photonsPhiBin*16)/(proportion*recoEfficiency*photonsPerTrack*mometumCutFactor);
              sortingEvents = sortingTracks/tracksPerEvent; //events to sort through
            }
            else
            {
              sortingTracks = -1; 
              sortingEvents = -1;
            }
            
            
            //
            //load up efficiency values and events to filter
            //
            recoEfficiencyPlot->SetBinContent(x,y,recoEfficiency);
            totalEventsToFilterPlot->SetBinContent(x,y,sortingEvents);
          } 
        }
        //clean up memory
        scanAreaPlot->Delete();
        scanAreaMirrPairPlot->Delete();
        scanAreaPlotMergedX->Delete();
        scanAreaMirrPairPlotMergedX->Delete();
        
        recoEfficiencyPlot->GetXaxis()->SetTitle("X coordinate(mm)");
        recoEfficiencyPlot->GetYaxis()->SetTitle("Y coordinate(mm)");
        totalEventsToFilterPlot->GetXaxis()->SetTitle("X coordinate(mm)");
        totalEventsToFilterPlot->GetYaxis()->SetTitle("Y coordinate(mm)");
     
        Int_t binGlobal = recoEfficiencyPlot->GetMaximumBin() ;
        Int_t nx = recoEfficiencyPlot->GetNbinsX() + 2;
        Int_t ny = recoEfficiencyPlot->GetNbinsY() + 2;
        Int_t binCoordinateX = binGlobal%nx;
        Int_t binCoordinateY = ((binGlobal-binCoordinateX)/nx)%ny;
        
        //get number of events
        double numberofEvents = totalEventsToFilterPlot->GetBinContent(binCoordinateX,binCoordinateY); 
        //get x and y measurement rather than just bins
        double xBinLow =recoEfficiencyPlot->GetXaxis()->GetBinLowEdge(binCoordinateX);
        double xBinHigh=xBinLow+ recoEfficiencyPlot->GetXaxis()->GetBinWidth(binCoordinateX);
        double yBinLow = recoEfficiencyPlot->GetYaxis()->GetBinLowEdge(binCoordinateY);
        double yBinHigh= yBinLow+ recoEfficiencyPlot->GetYaxis()->GetBinWidth(binCoordinateY); ;

        //load up variables again if have a higher value
        if(maxEff <= recoEfficiencyPlot->GetMaximum()){
          //don't load up max value for smallest bins but print values and load up plots. High overbinning entries could stop 
          //'legitimate' entries being loaded - although, they may be fine, depending on number of entries
            if((rebinx >= rebinXMin) && (rebiny >= rebinYMin) && (numberofEvents < maxAllowedEvents) && (numberofEvents > 0)) 
          { 
            maxEff = recoEfficiencyPlot->GetMaximum();
            maxEvents = numberofEvents ;
            maxXBinLow = xBinLow;
            maxXBinHigh = xBinHigh;
            maxYBinLow = yBinLow;
            maxYBinHigh = yBinHigh;
            maxTotalEfficiencyPlot = recoEfficiencyPlot;
            maxTotalEventsToFilterPlot = totalEventsToFilterPlot;
          }
          else if (maxEffBin1 < recoEfficiencyPlot->GetMaximum())
          { maxEffBin1 =  recoEfficiencyPlot->GetMaximum();
            maxEventsBin1 = numberofEvents;
            maxXBinLowBin1 = xBinLow;
            maxXBinHighBin1 = xBinHigh;
            maxYBinLowBin1 = yBinLow;
            maxYBinHighBin1 = yBinHigh;
            maxTotalEfficiencyPlotBin1 = recoEfficiencyPlot;
            maxTotalEventsToFilterPlotBin1 = totalEventsToFilterPlot;
          }
        }

        //load up compariable example for each one (but don't store as highest variables)
        if((rebinx == 6 && rebiny == 6 && fromTuple) || (rebinx == 4 && rebiny == 4 && fromTuple)){
          histList->Add(recoEfficiencyPlot);
          histList->Add(totalEventsToFilterPlot);
          std::cout << "- Reference efficiency value "   <<  recoEfficiencyPlot->GetMaximum()  << std::endl;
          std::cout << "- Events to filter "  << numberofEvents << std::endl;
          std::cout << "- Area co-ordinates ("<<xBinLow<<","<<xBinHigh<<","<<yBinLow<<","<<yBinHigh<<")";
          std::cout << " Rebin " << rebinx  << " , "  << rebiny  << std::endl;
          std::cout << " " << std::endl;
        }

      }
    }
    //
    //add highest efficiency plot and corresponding number of events to filter plot
    //
    std::cout << "***************************************************"  << std::endl;
    if(maxEff> maxEffBin1)std::cout << "* Best Value (rebin>1)" 
                                       << " P" << primeMirrArray[j]<< " S" <<secMirrArray[j]<<std::endl;
    else std::cout << "* Best Value with rebinx>" <<rebinXMin <<" rinbiny>" <<rebinYMin << 
           " P" << primeMirrArray[j]<< " S" <<secMirrArray[j]<<std::endl;;
    if(maxEvents > 5000000) std::cout << "******* Over 5 million events required *******"  << std::endl;
    if(maxEff<0.1) std::cout<<"************** Efficiency LESS THAN 0.1 ************"<<std::endl;
    else if(maxEff<0.3) std::cout << " **** less than 0.3 ****  "<< std::endl;
    std::cout << "* Max efficiency value "   << maxEff  << std::endl;
    std::cout << "* Events to filter "  << maxEvents << std::endl;
    std::cout << "* Area co-ordinates ("<<maxXBinLow<<","<<maxXBinHigh<<","<<maxYBinLow<<","<<maxYBinHigh<<")" << std::endl;
   
    histList->Add(maxTotalEfficiencyPlot);
    histList->Add(maxTotalEventsToFilterPlot);
    if(maxEff < maxEffBin1)
    { std::cout << " " << std::endl;
      std::cout << "* Best Value (rebinx <)" <<rebinXMin <<" rinbiny<" <<rebinYMin << std::endl;
      std::cout << "* Max efficiency value "   << maxEffBin1  << std::endl;
      std::cout << "* Events to filter "  << maxEventsBin1 << std::endl;
      std::cout << "* Area co-ordinates ("<<maxXBinLowBin1<<","<<maxXBinHighBin1
                <<","<<maxYBinLowBin1<<","<<maxYBinHighBin1<<")" << std::endl;
      histList->Add(maxTotalEfficiencyPlotBin1);
      histList->Add(maxTotalEventsToFilterPlotBin1);
    }
    std::cout << "***************************************************"  << std::endl;
    std::cout << "              " << std::endl;

    //load up highest number of events to normalise with
    if( highestEventsToFilter < maxEvents) highestEventsToFilter = maxEvents; 
    
    //load up results
    results[0][j] = (double) primeMirrArray[j];
    results[1][j] = (double) secMirrArray[j];
    results[2][j] = maxXBinLow;
    results[3][j] = maxXBinHigh;
    results[4][j] = maxYBinLow;
    results[5][j] = maxYBinHigh;
    results[6][j] = maxEff;
    results[7][j] = maxEvents;

    //cleanup memory
    
    // maxTotalEfficiencyPlot->Delete();
    //maxTotalEventsToFilterPlot->Delete();
    //maxTotalEfficiencyPlotBin1->Delete();
    //maxTotalEventsToFilterPlotBin1->Delete();
    scanAreaPlotUnMerged->Delete();
    scanAreaMirrPairPlotUnMerged->Delete();
  }
  //
  //Print out to screen final results
  //
  std::cout << "*Highest events to filter " << highestEventsToFilter  << std::endl;
  std::cout << "*Final results:   "  << std::endl;
  for(int j=0; j<numbMirrPair; j++){
    results[8][j] = (results[7][j]/highestEventsToFilter); // normalise to highest amount of events
    std::cout << "Pri"<< results[0][j] << " Sec" << results[1][j];
    std::cout << " (" << results[2][j] << "," << results[3][j] << "|";
    std::cout <<  results[4][j] << "," << results[5][j] << ") ";
    std::cout <<"Normalised Weight: "  << results[8][j] << " Efficiency: " << results[6][j] << std::endl;
  }

  std::cout << std::endl;  
  std::cout << "*Mirror pairs with efficiency lower than 0.3 greater than 0.1:" << std::endl;
  for(int j=0; j<numbMirrPair; j++){
    if(results[6][j] < 0.3 && results[6][j] > 0.1 )
    {
      std::cout << "Pri"<< results[0][j] << " Sec" << results[1][j];
      std::cout << " (" << results[2][j] << "," << results[3][j] << "|";
      std::cout << results[4][j] << "," << results[5][j] << ") ";
      std::cout <<"Normalised Weight: " << results[8][j] << " Efficiency: " << results[6][j] << std::endl;
    }
  }

  std::cout << std::endl;  
  std::cout << "*Mirror pairs with efficiency lower than 0.1:" << std::endl;
  for(int j=0; j<numbMirrPair; j++){
    if(results[6][j] < 0.1)
    {
      std::cout << "Pri"<< results[0][j] << " Sec" << results[1][j];
      std::cout << " (" << results[2][j] << "," << results[3][j] << "|";
      std::cout << results[4][j] << "," << results[5][j] << ") ";
      std::cout <<"Normalised Weight: " << results[8][j] << " Efficiency: " << results[6][j] << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "*Mirror pairs requiring greater than 5 million events to filter through:" << std::endl;
  for(int j=0; j<numbMirrPair; j++){
    if(results[7][j] > 5000000)
    {
      std::cout << "Pri"<< results[0][j] << " Sec" << results[1][j];
      std::cout << " (" << results[2][j] << "," << results[3][j] << "|";
      std::cout <<  results[4][j] << "," << results[5][j] << ") ";
      std::cout << "Normalised Weight: " << results[8][j] << " Efficiency: " << results[6][j];
      std::cout << "Events to filter: "<< results[7][j] << std::endl;
    }
  }
  
  
  
  TFile* andrewsFile =
new 
TFile("~/cmtuser/Brunel_v45r0/Rich/RichAlignment/options/HLTMap/MapMaker/ScanMap2julyLeft.root","recreate");
  andrewsFile->cd();
  histList->Write();
  andrewsFile->Close();
  histList->Delete();
  f->Close();
  
  
  //
  //Create a root file to hold final variables
  //
  TFile fRoot("finalVariables2julyBatchLeft.root","recreate");
  TTree t1("t1","HLT Map Variables");

  //variables for loading up root file
  Int_t primeMirrRoot, secMirrRoot;
  Double_t xBinLowRoot, xBinHighRoot, yBinLowRoot, yBinHighRoot;
  Double_t maxEffRoot, maxEventsRoot, normalisedWeightingRoot;
  t1.Branch("primeMirr",&primeMirrRoot,"primeMirrRoot/I");
  t1.Branch("secMirr",&secMirrRoot,"secMirrRoot/I");
  t1.Branch("xBinLow",&xBinLowRoot,"xBinLowRoot/D");
  t1.Branch("xBinHigh",&xBinHighRoot,"xBinHighRoot/D");
  t1.Branch("yBinLow",&yBinLowRoot,"yBinLowRoot/D");
  t1.Branch("yBinHigh",&yBinHighRoot,"yBinHighRoot/D");
  t1.Branch("maxEff",&maxEffRoot,"maxEffRoot/D");
  t1.Branch("maxEvents",&maxEventsRoot,"maxEventsRoot/D");
  t1.Branch("normalisedWeighting",&normalisedWeightingRoot,"normalisedWeightingRoot/D");
  
  std::cout <<  " point 1 "  << std::endl;

  //fill the tree
  for (Int_t i=0; i<numbMirrPair; i++) {
    primeMirrRoot =  primeMirrArray[i]; ;//(Int_t) results[0][i];
    secMirrRoot =  secMirrArray[i];//(Int_t) results[1][i];
    xBinLowRoot = (Double_t) results[2][i] ;
    xBinHighRoot =(Double_t) results[3][i];
    yBinLowRoot =(Double_t) results[4][i];
    yBinHighRoot =(Double_t) results[5][i];
    maxEffRoot =(Double_t) results[6][i];
    maxEventsRoot =(Double_t) results[7][i];
    normalisedWeightingRoot =(Double_t) results[8][i];
    t1.Fill();
  }
  t1.Write();
  fRoot.Close();
 
   


  return;
}



//*************************************************************************************************************
//
//For a potential HLT-selection area, creates efficiency plot for photon hits in a specific mirror pair using
//Track position in X, Y and Phi distribution efficiency.
//
//*************************************************************************************************************
void HLTMap::AreaMap()
{
  gStyle->SetOptStat(0);

  TFile *f = TFile::Open("/afs/cern.ch/work/p/phadc/private/gangadir/workspace/phadc/LocalXML/120/output/rich.tuples.root");
  TTree * treePhoton = (TTree*)f->Get("RICH/RichAlignMoniR2Gas/RichAlign");
	
  UInt_t primeMirrNum ;
  UInt_t secMirrNum ;
  double primeMirrX, primeMirrY;
  double trackX, trackY;
  double deltaTheta, phi;
  int eventCount, segmentCount;

  //photon
  treePhoton->SetBranchAddress("sphMirror",&primeMirrNum);
  treePhoton->SetBranchAddress("secMirror", &secMirrNum);
  treePhoton->SetBranchAddress("SphericalMirrX", &primeMirrX);
  treePhoton->SetBranchAddress("SphericalMirrY", &primeMirrY);
  //track
  treePhoton->SetBranchAddress("xExitpoint",&trackX);
  treePhoton->SetBranchAddress("yExitpoint",&trackY);
  treePhoton->SetBranchAddress("eventCount",&eventCount);
  treePhoton->SetBranchAddress("segmentKey",&segmentCount);
  treePhoton->SetBranchAddress("deltaTheta",&deltaTheta);
  treePhoton->SetBranchAddress("phiRec",&phi);
  //
  //mirror pair arrays
  //
  unsigned int primeMirrArray[47]=
  {12,12,16,13,17,13,14,18,14,15,19,16,20,17,21,18,22,19,23,20,24,21,25,22,26,22,27,11,10,9,8,11,7,10,6,9,5,8,4,6,3,6,2,5,1,4,0};
  
  unsigned int secMirrArray[47]=
  { 9, 8, 8, 9, 9,10,10,10,11,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,11,10,9,8, 7,7, 6,6,5,5,4,4,3,3,2,2,1,1,0,0};
  //
  //control boundaries of primary mirror plots (not needed at moment in this program)
  //
  /*
  int xStart[47]={0,0,-750,0,0,0,0,0,0,1200,300,-750,-750,0,0,0,0,300,1000,-750,
                 -500,0,0,500,0,500,1000,1000,0,0,-750,1000,1200,0,500,0,0,-750,-750,0,1000,0,0,0,0,-750,-750};

  int xEnd[47]={1500,1500,750,1500,1500,1500,1500,1500,1500,2700,1800,750,750,1500,1500,1500,1500,1800,2500,7500,1000,1500,
     1500,2000,1500,2000,2500,2500,1500,1500,750,2500,2700,1500,2000,1500,1500,750,750,1500,2500,1500,1500,1500,1500,750,750};

  int yStart[47]={-750,-750,-750,-750,-750,-750,-750,-750,-750,-750,-750,-750,0,-750,0,0,0,0,0,0,500,0,500,0,750,0,500,-750,-750,
                 -750,-750,-750,-1700,-750,-1500,-750,-1500,-750,-1500,-1500,-2000,-1500,-2200,-1500,-2100,-1500,-2200};

  int yEnd[47]={750,750,750,750,750,750,750,750,750,750,750,750,1500,750,1500,1500,1500,1500,1500,1500,2000,1500,2000,1500,2250,
               1500,2000,750,750,750,750,750,-200,750,0,750,0,750,0,0,-500,0,-700,0,-600,0,-700};
  */
  //
  //Holds details for test areas to select tracks by in HLT
  //
  double xAreaStart[5] = {380.0,120.0,30.0,420.0,390.0}; 
  double xAreaEnd[5] = {410.0, 150.0,60,450.0, 420.0};
  double yAreaStart[5] = {-140.0, -120.0,360.0,-100.0, 122.0};
  double yAreaEnd[5] = {100.0, -60,390.0,120.0, 148.0};

  TObjArray* histList = new TObjArray(0);  

  //
  //loop over mirror pairs
  //
  for(int j=0; j<47; j++){
    std::cout << "j: "<<j << " matrix " << primeMirrArray[j]  << "  " << secMirrArray[j] <<  std::endl;
    
    std::ostringstream thetaPhiName,phiName,thetaName;
    std::ostringstream thetaPhiAreaName,phiAreaName,thetaAreaName;
    std::ostringstream thetaPhiAreaName6000,phiAreaName6000,thetaAreaName6000;
    
    thetaPhiName << "thetaVsPhi" <<  primeMirrArray[j] << "Sec" <<secMirrArray[j];
    phiName << "phi" << primeMirrArray[j] << "Sec" << secMirrArray[j];
    thetaName << "theta" << primeMirrArray[j] << "Sec" << secMirrArray[j];
    
    thetaPhiAreaName << "AreathetaVsPhi" <<  primeMirrArray[j] << "Sec" <<secMirrArray[j] << "("<< xAreaStart[j]<<","
                     << xAreaEnd[j]<<")("<<yAreaStart[j]<<","<<yAreaEnd[j]<<")" ;
    phiAreaName << "Areaphi" << primeMirrArray[j] << "Sec" << secMirrArray[j] << "("<< xAreaStart[j]<<","
                     << xAreaEnd[j]<<")("<<yAreaStart[j]<<","<<yAreaEnd[j]<<")";
    thetaAreaName << "Areatheta" << primeMirrArray[j] << "Sec" << secMirrArray[j] << "("<< xAreaStart[j]<<","
                     << xAreaEnd[j]<<")("<<yAreaStart[j]<<","<<yAreaEnd[j]<<")";

    thetaPhiAreaName6000 << "6000AreathetaVsPhi" <<  primeMirrArray[j] << "Sec" <<secMirrArray[j] << "("<< xAreaStart[j]<<","
                     << xAreaEnd[j]<<")("<<yAreaStart[j]<<","<<yAreaEnd[j]<<")" ;
    phiAreaName6000 << "6000Areaphi" << primeMirrArray[j] << "Sec" << secMirrArray[j] << "("<< xAreaStart[j]<<","
                     << xAreaEnd[j]<<")("<<yAreaStart[j]<<","<<yAreaEnd[j]<<")";
    thetaAreaName6000 << "6000Areatheta" << primeMirrArray[j] << "Sec" << secMirrArray[j] << "("<< xAreaStart[j]<<","
                     << xAreaEnd[j]<<")("<<yAreaStart[j]<<","<<yAreaEnd[j]<<")";

  
    TH2F * thetaPhiPlot = 
      new TH2F(thetaPhiName.str().c_str(),thetaPhiName.str().c_str(),20.0,0.0,6.2827226,20.0,-0.006,0.006);
    TH1F * phiPlot = new TH1F(phiName.str().c_str(),phiName.str().c_str(),20.0, 0.0, 6.2827226 );
    TH1F * thetaPlot =new TH1F(thetaName.str().c_str(),thetaName.str().c_str(),20.0, -0.006,0.006);

    TH2F * thetaPhiAreaPlot = 
      new TH2F(thetaPhiAreaName.str().c_str(),thetaPhiAreaName.str().c_str(),20.0,0.0,6.2827226,20.0,-0.006,0.006);
    TH1F * phiAreaPlot = new TH1F(phiAreaName.str().c_str(),phiAreaName.str().c_str(),20.0, 0.0, 6.2827226 );
    TH1F * thetaAreaPlot =new TH1F(thetaAreaName.str().c_str(),thetaAreaName.str().c_str(),20.0, -0.006,0.006);
  
    double photonCount = 0;
    double areaCount = 0;
    double finalCount = 0;
   
    
   
    //total number of photons in large tuple = 398027636 (aprox 4x10^8) (small tuple 1.2x10^7)

    //  for(int i=0; i<treePhoton->GetEntries(); i++)
    for(int i=0; i<1000000; i++)
    {
      treePhoton->GetEntry(i);
      photonCount++;
      //select area
      if(trackX < xAreaStart[j] || trackX > xAreaEnd[j]) continue;
      if(trackY < yAreaStart[j] || trackY > yAreaEnd[j]) continue;
      areaCount++;
      // check photons in correct mirror combination
      if(primeMirrNum != primeMirrArray[j]) continue;
      if(secMirrNum != secMirrArray[j]) continue;
      finalCount++;
    }

    //for(int i=0; i<treePhoton->GetEntries(); i++)
    for(int i=0; i<1000000; i++)
    {
      treePhoton->GetEntry(i);
      
      // check photons in correct mirror combination
      if(primeMirrNum != primeMirrArray[j]) continue;
      if(secMirrNum != secMirrArray[j]) continue;
      
      phiPlot->Fill(phi);
      thetaPlot->Fill(deltaTheta);
      thetaPhiPlot->Fill(phi,deltaTheta);
      
      //select area
      if(trackX < xAreaStart[j] || trackX > xAreaEnd[j]) continue;
      if(trackY < yAreaStart[j] || trackY > yAreaEnd[j]) continue;
      phiAreaPlot->Fill(phi);
      thetaAreaPlot->Fill(deltaTheta);
      thetaPhiAreaPlot->Fill(phi,deltaTheta);
    }

    phiPlot->GetXaxis()->SetTitle("#phi [rad]");
    phiAreaPlot->GetXaxis()->SetTitle("#phi [rad]");
    thetaPlot->GetXaxis()->SetTitle("#Delta#theta [rad]");
    thetaAreaPlot->GetXaxis()->SetTitle("#Delta#theta [rad]");
    thetaPhiPlot->GetXaxis()->SetTitle("#phi [rad]");
    thetaPhiPlot->GetYaxis()->SetTitle("#Delta#theta [rad]");
    thetaPhiAreaPlot->GetXaxis()->SetTitle("#phi [rad]");
    thetaPhiAreaPlot->GetYaxis()->SetTitle("#Delta#theta [rad]");

    histList->Add(phiPlot);
    histList->Add(thetaPlot);
    histList->Add(thetaPhiPlot);
    histList->Add(phiAreaPlot);
    histList->Add(thetaAreaPlot);
    histList->Add(thetaPhiAreaPlot);
   
    //
    //calculate effect of phi on efficiency by finding size of 16th bin and comparing to total distribution
    //
    int binSize[20];
    //get phi bin entries and order by size
    for(int j=1; j<21; j++)
    {
      binSize[j-1] =  phiAreaPlot->GetBinContent(j);  
    }
    int elements = sizeof(binSize) / sizeof(binSize[0]); 
    std::sort(binSize, binSize + elements);

    double photonsPerTrack = 60.0; //average number of photons per track
    double tracksPerEvent = 21.0; //average number of tracks per event 
    double photonsPhiBin = 300; // number of photons needed in each phi bin
    double photonHits =(photonsPhiBin/binSize[4])*finalCount;//photons required so 16 phi bins have at least 300 entries
    double phiEfficiency = (photonsPhiBin*16)/photonHits; //efficiency caused by uneven phi distribution
    double areaEfficiency = finalCount/areaCount; //efficiency of correct mirror-pair-selection within selected-area
    double proportion = finalCount/photonCount;
    double sortingTracks=((photonsPhiBin*16/proportion)/photonsPerTrack)/(phiEfficiency*areaEfficiency);//tracks to filter 
    double sortingEvents = sortingTracks/tracksPerEvent; //events to sort through
    double tracksReconstruct = ((photonsPhiBin*16.0)/photonsPerTrack)/(areaEfficiency*phiEfficiency);

    std::cout <<"* Number of entries in 16th least populated theta bin: "<< binSize[4] << std::endl;
    std::cout << "* Total number of photons filtered: " << photonCount << std::endl;
    std::cout << "* Number of filtered photons with tracks in selected area: " << finalCount << std::endl;
    std::cout << "* Proportion of total photons in selected area: " << proportion << std::endl;
    std::cout << "* Correct mirror-pair photon efficiency in selected-area: " << areaEfficiency << std::endl;
    std::cout << "* Phi distribution efficiency: " << phiEfficiency << std::endl;
    std::cout << "*** Total Efficiency: " << areaEfficiency*phiEfficiency  << std::endl;
    std::cout << "*** Tracks required to sort through: "<<sortingTracks<<" Events: "<< sortingEvents
              <<" photons "<< sortingTracks*60 <<std::endl;
    std::cout << "*** Tracks to reconstruct " << tracksReconstruct << "(cf 100% efficient =  80 tracks)" << std::endl;


  }
  
  TFile* andrewsFile = new TFile("AreaMap.root", "recreate");
  histList->Write();
  andrewsFile->Close();

  return;
}


//*************************************************************************************************************
//
// Creates efficiency plot for photon hits in a specific mirror pair using Track position in RICH2 primary mirror plane
//
//*************************************************************************************************************

void HLTMap::MakeTrackMap()
{
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open("/afs/cern.ch/work/p/phadc/private/gangadir/workspace/phadc/LocalXML/120/output/rich.tuples.root");
  TTree * treePhoton = (TTree*)f->Get("RICH/RichAlignMoniR2Gas/RichAlign");
	
  UInt_t primeMirrNum ;
  UInt_t secMirrNum ;
  double primeMirrX, primeMirrY;
  double trackX, trackY;
  int eventCount, segmentCount;

  //photon variables
  treePhoton->SetBranchAddress("sphMirror",&primeMirrNum);
  treePhoton->SetBranchAddress("secMirror", &secMirrNum);
  treePhoton->SetBranchAddress("SphericalMirrX", &primeMirrX);
  treePhoton->SetBranchAddress("SphericalMirrY", &primeMirrY);
  //track variables
  treePhoton->SetBranchAddress("xExitpoint",&trackX);
  treePhoton->SetBranchAddress("yExitpoint",&trackY);
  treePhoton->SetBranchAddress("eventCount",&eventCount);
  treePhoton->SetBranchAddress("segmentKey",&segmentCount);  



  unsigned int primeMirrArray[47]=
  {12,12,16,13,17,13,14,18,14,15,19,16,20,17,21,18,22,19,23,20,24,21,25,22,26,22,27,11,10,9,8,11,7,10,6,9,5,8,4,6,3,6,2,5,1,4,0};
  
  unsigned int secMirrArray[47]=
  { 9, 8, 8, 9, 9,10,10,10,11,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,11,10,9,8, 7,7, 6,6,5,5,4,4,3,3,2,2,1,1,0,0};

  int xStart[47]={0,0,-750,0,0,0,0,0,0,1200,300,-750,-750,0,0,0,0,300,1000,-750,
                 -500,0,0,500,0,500,1000,1000,0,0,-750,1000,1200,0,500,0,0,-750,-750,0,1000,0,0,0,0,-750,-750};

  int xEnd[47]={1500,1500,750,1500,1500,1500,1500,1500,1500,2700,1800,750,750,1500,1500,1500,1500,1800,2500,7500,1000,1500,
     1500,2000,1500,2000,2500,2500,1500,1500,750,2500,2700,1500,2000,1500,1500,750,750,1500,2500,1500,1500,1500,1500,750,750};

  int yStart[47]={-750,-750,-750,-750,-750,-750,-750,-750,-750,-750,-750,-750,0,-750,0,0,0,0,0,0,500,0,500,0,750,0,500,-750,-750,
                 -750,-750,-750,-1700,-750,-1500,-750,-1500,-750,-1500,-1500,-2000,-1500,-2200,-1500,-2100,-1500,-2200};

  int yEnd[47]={750,750,750,750,750,750,750,750,750,750,750,750,1500,750,1500,1500,1500,1500,1500,1500,2000,1500,2000,1500,2250,
               1500,2000,750,750,750,750,750,-200,750,0,750,0,750,0,0,-500,0,-700,0,-600,0,-700};
  int nbins = 50.0;

  TObjArray* histList = new TObjArray(0);  
  
  for(int j=0; j<47; j++){
    
    std::cout << "j: "<<j << " matrix " << primeMirrArray[j]  << "  " << secMirrArray[j] <<  std::endl;
    
    std::ostringstream specificName, primaryName, purityName;
    specificName <<"specificPri" << primeMirrArray[j] << "Sec" <<secMirrArray[j];
    primaryName <<"primary" << primeMirrArray[j];
    purityName <<"efficiencyPri" << primeMirrArray[j] << "Sec" <<secMirrArray[j];

    TH2F * primaryXYPlotSpecific = 
      new TH2F(specificName.str().c_str(),specificName.str().c_str(),nbins,xStart[j],xEnd[j],nbins,yStart[j],yEnd[j]);
    TH2F * primaryXYPlot = 
      new TH2F(primaryName.str().c_str(),primaryName.str().c_str(),nbins,xStart[j],xEnd[j],nbins,yStart[j],yEnd[j]);
    TH2F * purityXYPlot = 
      new TH2F(purityName.str().c_str(),purityName.str().c_str(),nbins,xStart[j],xEnd[j],nbins,yStart[j],yEnd[j]);

    int correctPriMirr = 0;
    int wrongPriMirr = 0;
    int correctSecMirr = 0;
    int wrongSecMirr = 0;
    int trackNumber = 0;
    double prevTrackX = 0;
    double prevTrackY = 0;
    int prevEventNumber = 0;
 
    
    //total number of photons in large tuple = 398027636 (aprox 4x10^8) (small tuple 1.2x10^7)
    
     for(int i=0; i<treePhoton->GetEntries(); i++)
    // for(int i=0; i<100000000; i++)
    {
      treePhoton->GetEntry(i);

      //initialise track number for first loop
      if(i == 0){
        trackNumber = segmentCount;
        prevTrackX = trackX;
        prevTrackY = trackY;
      }

      //potentially write plot when get new track / on last entry/ when change event in case new track  has same track number
      if(segmentCount != trackNumber || prevEventNumber != eventCount || i == (treePhoton->GetEntries() - 1))
      { 
        //tracks with any hits in primary mirror
        if (correctPriMirr >= 1) primaryXYPlot->Fill(prevTrackX, prevTrackY);
        
        // hits in specific pair
        if (wrongPriMirr == 0 && wrongSecMirr == 0) primaryXYPlotSpecific->Fill(prevTrackX, prevTrackY);
        
        //reset numbers
        trackNumber = segmentCount; 
        correctPriMirr = 0;
        wrongPriMirr = 0;
        correctSecMirr = 0;
        wrongSecMirr = 0;
      }
      
      //load up X and Y for next loop
      prevTrackX = trackX;
      prevTrackY = trackY;
      prevEventNumber = eventCount;
      if (primeMirrNum == primeMirrArray[j]) correctPriMirr++;
      else wrongPriMirr++;
      if(secMirrNum == secMirrArray[j]) correctSecMirr++;
      else wrongSecMirr++;
    }

    //
    //Calculate efficiency for specific mirror pair. 
    //
    for (int xbin =1; xbin <=100; xbin++){ 
      for (int ybin=1;ybin<=100; ybin++){
        if((primaryXYPlot->GetBinContent(xbin,ybin))==0)continue;
        if((primaryXYPlotSpecific->GetBinContent(xbin,ybin))==0)
        {
          purityXYPlot->SetBinContent(xbin,ybin,0.00001); //bit of a fudge
        }
        else
        {
          purityXYPlot->SetBinContent(
          xbin,ybin,((primaryXYPlotSpecific->GetBinContent(xbin,ybin))/(primaryXYPlot->GetBinContent(xbin,ybin))));
        }
        
        //
        //error
        //
        double kyield = (primaryXYPlotSpecific->GetBinContent(xbin,ybin));
        double nyield = (primaryXYPlot->GetBinContent(xbin,ybin));
        double var =fabs((((kyield+1)*(kyield+2))/((nyield+2)*(nyield+3)))-(((kyield+1)*(kyield+1))/((nyield+2)*(nyield+2))));
        double error = TMath::Sqrt(var);
        if(nyield !=0) purityXYPlot->SetBinError(xbin,ybin,error);    
        
      }
    } 
    primaryXYPlotSpecific->GetXaxis()->SetTitle("X coordinate(mm)");
    primaryXYPlotSpecific->GetYaxis()->SetTitle("Y coordinate(mm)");
    primaryXYPlot->GetXaxis()->SetTitle("X coordinate(mm)");
    primaryXYPlot->GetYaxis()->SetTitle("Y coordinate(mm)");
    purityXYPlot->GetXaxis()->SetTitle("X coordinate(mm)");
    purityXYPlot ->GetYaxis()->SetTitle("Y coordinate(mm)");
    purityXYPlot->SetMaximum(1.0);//setting max for colz scale
    histList->Add(primaryXYPlotSpecific);
    histList->Add(primaryXYPlot);
    histList->Add(purityXYPlot);

  }

  TFile* andrewsFile = new TFile("TrackMapforHLTBATCH.root", "recreate");
  histList->Write();
  andrewsFile->Close();
  

  return;
}




//*************************************************************************************************************
//
// Creates efficiency map for photon hits in a specific mirror pair in the RICH2 mirror plane
//
//*************************************************************************************************************

void HLTMap::MakePhotonMap()
{
  gStyle->SetOptStat(0);//
	TFile *f = TFile::Open("rich.tuples.root");
  TTree * treePhoton = (TTree*)f->Get("RICH/RichAlignMoniR2Gas/RichAlign");
	
  UInt_t primeMirrNum ;
  UInt_t secMirrNum ;
  double primeMirrX, primeMirrY;

  treePhoton->SetBranchAddress("sphMirror",&primeMirrNum);
  treePhoton->SetBranchAddress("secMirror", &secMirrNum);
  treePhoton->SetBranchAddress("SphericalMirrX", &primeMirrX);
  treePhoton->SetBranchAddress("SphericalMirrY", &primeMirrY);
  //
  //primary and secondary mirror numbers
  //
  unsigned int primeMirrArray[47]=
  {12,12,16,13,17,13,14,18,14,15,19,16,20,17,21,18,22,19,23,20,24,21,25,22,26,22,27,11,10,9,8,11,7,10,6,9,5,8,4,6,3,6,2,5,1,4,0};
  
  unsigned int secMirrArray[47]=
  { 9, 8, 8, 9, 9,10,10,10,11,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,19,19,11,10,9,8, 7,7, 6,6,5,5,4,4,3,3,2,2,1,1,0,0};
  //
  //plot co-ordinates for relevant mirror pair
  //
  int xStart[47]={0,0,-750,0,0,0,0,0,0,1200,300,-750,-750,0,0,0,0,300,1000,-750,
                 -500,0,0,500,0,500,1000,1000,0,0,-750,1000,1200,0,500,0,0,-750,-750,0,1000,0,0,0,0,-750,-750};

  int xEnd[47]={1500,1500,750,1500,1500,1500,1500,1500,1500,2700,1800,750,750,1500,1500,1500,1500,1800,2500,7500,1000,1500,
     1500,2000,1500,2000,2500,2500,1500,1500,750,2500,2700,1500,2000,1500,1500,750,750,1500,2500,1500,1500,1500,1500,750,750};

  int yStart[47]={-750,-750,-750,-750,-750,-750,-750,-750,-750,-750,-750,-750,0,-750,0,0,0,0,0,0,500,0,500,0,750,0,500,-750,-750,
                 -750,-750,-750,-1700,-750,-1500,-750,-1500,-750,-1500,-1500,-2000,-1500,-2200,-1500,-2100,-1500,-2200};

  int yEnd[47]={750,750,750,750,750,750,750,750,750,750,750,750,1500,750,1500,1500,1500,1500,1500,1500,2000,1500,2000,1500,2250,
               1500,2000,750,750,750,750,750,-200,750,0,750,0,750,0,0,-500,0,-700,0,-600,0,-700};
  int nbins = 50.0;

  TObjArray* histList = new TObjArray(0);  
  //
  //loop over mirror pairs
  //
  for(int j=0; j<47; j++){

    std::cout << "j: "<<j << " matrix " << primeMirrArray[j]  << "  " << secMirrArray[j] <<  std::endl;
    
    std::ostringstream specificName, primaryName, purityName;
    specificName <<"specificPri" << primeMirrArray[j] << "Sec" <<secMirrArray[j];
    primaryName <<"primary" << primeMirrArray[j];
    purityName <<"efficiencyPri" << primeMirrArray[j] << "Sec" <<secMirrArray[j];

    TH2F * primaryXYPlotSpecific = 
      new TH2F(specificName.str().c_str(),specificName.str().c_str(),nbins,xStart[j],xEnd[j],nbins,yStart[j],yEnd[j]);
    TH2F * primaryXYPlot = 
      new TH2F(primaryName.str().c_str(),primaryName.str().c_str(),nbins,xStart[j],xEnd[j],nbins,yStart[j],yEnd[j]);
    TH2F * purityXYPlot = 
      new TH2F(purityName.str().c_str(),purityName.str().c_str(),nbins,xStart[j],xEnd[j],nbins,yStart[j],yEnd[j]);

    //
    // Fill x,y componets of all primary mirror photon hits of specific mirror pair and just specific primary mirror
    //
    for(int i=0; i<treePhoton->GetEntries(); i++)
    {
      treePhoton->GetEntry(i);
      if (primeMirrNum != primeMirrArray[j]) continue;
      primaryXYPlot->Fill(primeMirrX, primeMirrY); 
      if (secMirrNum != secMirrArray[j]) continue;
      primaryXYPlotSpecific->Fill(primeMirrX, primeMirrY);
    }
    //
    //Calculate efficiency for specific mirror pair. 
    //
    for (int xbin =1; xbin <=100; xbin++){ 
      for (int ybin=1;ybin<=100; ybin++){
        if((primaryXYPlot->GetBinContent(xbin,ybin))==0)continue;
        if((primaryXYPlotSpecific->GetBinContent(xbin,ybin))==0)
        {
          purityXYPlot->SetBinContent(xbin,ybin,0.00001); //bit of a hack
        }
        else
        {
          purityXYPlot->SetBinContent(
          xbin,ybin,((primaryXYPlotSpecific->GetBinContent(xbin,ybin))/(primaryXYPlot->GetBinContent(xbin,ybin))));
        }
        
        //
        //error
        //
        double kyield = (primaryXYPlotSpecific->GetBinContent(xbin,ybin));
        double nyield = (primaryXYPlot->GetBinContent(xbin,ybin));
        double var =fabs((((kyield+1)*(kyield+2))/((nyield+2)*(nyield+3)))-(((kyield+1)*(kyield+1))/((nyield+2)*(nyield+2))));
        double error = TMath::Sqrt(var);
        if(nyield !=0) purityXYPlot->SetBinError(xbin,ybin,error);    
        
      }
    } 
    primaryXYPlotSpecific->GetXaxis()->SetTitle("X coordinate(mm)");
    primaryXYPlotSpecific->GetYaxis()->SetTitle("Y coordinate(mm)");
    primaryXYPlot->GetXaxis()->SetTitle("X coordinate(mm)");
    primaryXYPlot->GetYaxis()->SetTitle("Y coordinate(mm)");
    purityXYPlot->GetXaxis()->SetTitle("X coordinate(mm)");
    purityXYPlot ->GetYaxis()->SetTitle("Y coordinate(mm)");
    purityXYPlot->SetMaximum(1.0);//setting max for colz scale
    histList->Add(primaryXYPlotSpecific);
    histList->Add(primaryXYPlot);
    histList->Add(purityXYPlot);

  }

  TFile* andrewsFile = new TFile("PhotonMapforHLT.root", "recreate");
  histList->Write();
  andrewsFile->Close();
  

  return;
}


//*************************************************************************************************************
//
// Creates plot of average track to photon distance, plots all photons X,Y posn on Rich2 Mirror plane.
// Also plots photon hits from sample events.
//
//*************************************************************************************************************
void HLTMap::TrackPhotonDistance()
{

 TFile *f =
    TFile::Open("/afs/cern.ch/work/p/phadc/private/gangadir/workspace/phadc/LocalXML/123/output/rich.tuples.root");
  //TFile *f = 
 //TFile::Open("/afs/cern.ch/work/p/phadc/private/gangadir/workspace/phadc/LocalXML/120/output/rich.tuples.root");
  TTree * treePhoton = (TTree*)f->Get("RICH/RichAlignMoniR2Gas/RichAlign");
	
  UInt_t primeMirrNum ;
  UInt_t secMirrNum ;
  double primeMirrX, primeMirrY;
  double trackX, trackY;
  int eventCount, segmentCount;
  bool unAmbigious;
  
  //photon
  treePhoton->SetBranchAddress("sphMirror",&primeMirrNum);
  treePhoton->SetBranchAddress("secMirror", &secMirrNum);
  treePhoton->SetBranchAddress("SphericalMirrX", &primeMirrX);
  treePhoton->SetBranchAddress("SphericalMirrY", &primeMirrY);
  treePhoton->SetBranchAddress("unAmbiguousPhoton", &unAmbigious);
  //track
  treePhoton->SetBranchAddress("xExitpoint",&trackX);
  treePhoton->SetBranchAddress("yExitpoint",&trackY);
  treePhoton->SetBranchAddress("eventCount",&eventCount);
  treePhoton->SetBranchAddress("segmentKey",&segmentCount);  

  TH1F * distancePlot = new TH1F("photontrackdistance","photontrackdistance" ,400,0,400);
  TH2F * photonsTotal = new TH2F("photonsTotalXY","photonsTotalXY",100,-3000,3000.,100,-3000.0,3000.0);
  TH2F * photonsEvent1 = new TH2F("photonsEvent1","photonsEvent1",100,-3000,3000.,100,-3000.0,3000.0);
  TH2F * photonsEvent2 = new TH2F("photonsEvent2","photonsEvent2",100,-3000,3000.,100,-3000.0,3000.0);
  

  for(int i=0; i<treePhoton->GetEntries(); i++)
  // for(int i=0; i<12000000; i++)
    {
      treePhoton->GetEntry(i);
      double distance  = TMath::Sqrt((trackX - primeMirrX)*(trackX - primeMirrX)+ (trackY - primeMirrY)*(trackY - primeMirrY)) ;

      if(!unAmbigious) continue;
      distancePlot->Fill(distance);
      photonsTotal->Fill(primeMirrX, primeMirrY);  
      if(eventCount == 1)
      {
        photonsEvent1->Fill(primeMirrX,primeMirrY);
      }
      if(eventCount == 2)
      {
        photonsEvent2->Fill(primeMirrX, primeMirrY);
      }
    }

  TObjArray* histList = new TObjArray(0);

  distancePlot->GetXaxis()->SetTitle("distance (mm)");
  photonsEvent1->GetXaxis()->SetTitle("X coordinate(mm)");
  photonsEvent1->GetYaxis()->SetTitle("Y coordinate(mm)");
  photonsEvent2->GetXaxis()->SetTitle("X coordinate(mm)");
  photonsEvent2->GetYaxis()->SetTitle("Y coordinate(mm)");
 
  histList->Add(distancePlot);
  histList->Add(photonsEvent1);
  histList->Add(photonsEvent2);
  photonsTotal->GetXaxis()->SetTitle("X coordinate(mm)");
  photonsTotal->GetYaxis()->SetTitle("Y coordinate(mm)");
  histList->Add(photonsTotal);

  TFile* distanceFile = new TFile("DistanceforHLTFULLNEWTUPLEunambig.root", "recreate");

  histList->Write();
  distanceFile->Close();

  return;
}
