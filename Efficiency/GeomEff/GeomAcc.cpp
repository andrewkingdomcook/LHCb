#define GeomAcc_cpp
#include "GeomAcc.h"
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
#include "TMath.h"


#define MAX_SIN_THETA_DAU_IN_LHCB   0.389418342308651
#define MIN_SIN_THETA_DAU_IN_LHCB   0.009999833334167


int main(){
  GeomAcc geomClass;

  TString saveName = "Up.root";
  TString geomTuple = "/Volumes/TimeMachineBackups/LHCb/Analysis/DoubleJpsi/ntuples/2011/Analysis/Stripping20r1/MicroDSTDiMuonDiMuonIncLine/Geometric/Up/TupleUp.root";
  geomClass.Calculate(geomTuple,saveName);
 
  saveName = "Down.root";
  geomTuple = "/Volumes/TimeMachineBackups/LHCb/Analysis/DoubleJpsi/ntuples/2011/Analysis/Stripping20r1/MicroDSTDiMuonDiMuonIncLine/Geometric/Down/TupleDown.root";
  geomClass.Calculate(geomTuple,saveName);


  return 0;
}



/* *************************************************************************** *
 *                                                                             *
 * Calculate:                                                                  *
 *                                                                             *
 * Main program to calculate and save the Geometric Acceptance Efficiency      *
 *                                                                             *
 * *************************************************************************** */

void GeomAcc::Calculate(TString geomTuple, TString saveName){

  TFile* file = new TFile(geomTuple);
  TTree* tree=(TTree*)file->Get("SelectJpsi/JpsiGen"); 

  double JpsiPT, JpsiRap, CosTheta;
  double MuNegDau, MuPosDau;

  tree->SetBranchAddress("JPsi_PT",&JpsiPT);
  tree->SetBranchAddress("JPsi_RAPIDITY",&JpsiRap);
  tree->SetBranchAddress("JPsi_COSTHETA",&CosTheta);
  tree->SetBranchAddress("MuMinus_THETA",&MuNegDau);
  tree->SetBranchAddress("MuPlus_THETA",&MuPosDau);


  int xbin = 20;
  int ybin = 5;

  TH2D *Num = new TH2D("Num","Num",xbin,0.0,10000.0, ybin,2.0,4.2);
  TH2D *Den = new TH2D("Den","Den",xbin,0.0,10000.0, ybin,2.0,4.2);
  
  TH2D *NumLong = new TH2D("NumLong","NumLong",xbin,0.0,10000.0, ybin,2.0,4.2);
  TH2D *DenLong = new TH2D("DenLong","DenLong",xbin,0.0,10000.0, ybin,2.0,4.2);
 
  TH2D *NumTrans = new TH2D("NumTrans","NumTrans",xbin,0.0,10000.0, ybin,2.0,4.2);
  TH2D *DenTrans = new TH2D("DenTrans","DenTrans",xbin,0.0,10000.0, ybin,2.0,4.2);

  for(int i=0; i<tree->GetEntries(); i++)
  {
    tree->GetEntry(i);

    Den->Fill(JpsiPT,JpsiRap);
    DenLong->Fill(JpsiPT,JpsiRap,(1-CosTheta*CosTheta));
    DenTrans->Fill(JpsiPT,JpsiRap,(1+CosTheta*CosTheta));

    //
    //Calculate Daughters in acceptance criteria
    //
    bool inAcc = Cut(MuPosDau, MuNegDau);

    if(inAcc) {
      Num->Fill(JpsiPT,JpsiRap);
      NumLong->Fill(JpsiPT,JpsiRap,(1-CosTheta*CosTheta));
      NumTrans->Fill(JpsiPT,JpsiRap,(1+CosTheta*CosTheta));
    }
  }

  //
  //Calculate efficiency
  //

  TH2D GeoEff =  (Eff(Num, Den, xbin, ybin, "GeoEff"));
  TH2D GeoLongEff =  (Eff(NumLong, DenLong, xbin, ybin,"GeoLongEff"));
  TH2D GeoTransEff =  (Eff(NumTrans, DenTrans, xbin, ybin,"GeoTransEff"));


  TObjArray* histList = new TObjArray(0);
  histList->Add(Num);
  histList->Add(Den);
  histList->Add(&GeoEff);
  histList->Add(NumLong);
  histList->Add(DenLong);
  histList->Add(&GeoLongEff);
  histList->Add(NumTrans);
  histList->Add(DenTrans);
  histList->Add(&GeoTransEff);
  TFile* saveFile = new TFile("Output/"+saveName, "recreate");
  histList->Write();
  saveFile->Close();

  //
  // Save into final usable format in Efficiency folder
  //

  Format(saveName, "NoPol/");
  Format(saveName, "Long/");
  Format(saveName, "Trans/");

  return;
}


/* *************************************************************************** *
 *                                                                             *
 * Cut:                                                                        *
 *                                                                             *
 * Calculate daughters in acceptance cut on Muons                              *
 *                                                                             *
 * *************************************************************************** */
bool GeomAcc::Cut(double MuPosDau, double MuNegDau){

  double absSinPos = TMath::Abs(TMath::Sin(MuPosDau));
  double absSinNeg = TMath::Abs(TMath::Sin(MuNegDau));
  if(absSinPos <  MIN_SIN_THETA_DAU_IN_LHCB || absSinPos >  MAX_SIN_THETA_DAU_IN_LHCB ) return false;
  if(absSinNeg <  MIN_SIN_THETA_DAU_IN_LHCB || absSinNeg >  MAX_SIN_THETA_DAU_IN_LHCB) return false;
   

  return true;
}



/* *************************************************************************** *
 *                                                                             *
 * Eff:                                                                        *
 *                                                                             *
 * Return Geometric Efficency TH2D                                             *
 *                                                                             *
 * *************************************************************************** */

TH2D GeomAcc::Eff(TH2D* Num, TH2D* Den, int xbin, int ybin, TString plotName){

 TH2D GeoEff(plotName,plotName,xbin,0.0,10000.0, ybin,2.0,4.2);
 GeoEff.GetXaxis()->SetTitle("P_{T}");
 GeoEff.GetYaxis()->SetTitle("Rapidity");
 


 for (int y =1; y<=ybin; y++){
    for (int x=1; x<= xbin; x++){
      
      if((Den->GetBinContent(x,y))==0)continue;
      if((Num->GetBinContent(x,y))==0)continue;
      GeoEff.SetBinContent(x,y,((Num->GetBinContent(x,y))/(Den->GetBinContent(x,y))));
      //
      //error
      //
      double kyield = (Num->GetBinContent(xbin,ybin));
      double nyield = (Den->GetBinContent(xbin,ybin));
      double var =fabs((((kyield+1)*(kyield+2))/((nyield+2)*(nyield+3)))-(((kyield+1)*(kyield+1))/((nyield+2)*(nyield+2))));
      double error = TMath::Sqrt(var);
      if(nyield !=0) GeoEff.SetBinError(xbin,ybin,error);
    }
  }


 return GeoEff;

}



/* *************************************************************************** *
 *                                                                             *
 * Format:                                                                     *
 *                                                                             *
 * Save into a format to be used in main code                                  *
 *                                                                             *
 * *************************************************************************** */
void GeomAcc::Format(TString saveName, TString pol){

  TString plotName;
  if(pol == "NoPol/") plotName = "GeoEff";
  else if (pol == "Long/") plotName = "GeoLongEff";
  else if (pol == "Trans/") plotName = "GeoTransEff"; 

  TFile *f1 = TFile::Open("Output/" + saveName);
  TH2D *GeoEff2 = (TH2D*) f1->Get(plotName);
  GeoEff2->SetName("Geometric Acceptance");

  //
  //hack to use existing description plot
  //

  TFile *f2 = TFile::Open("Description/" + saveName);
  TH1 * Descrip = (TH1*) f2->Get("Description");

  TObjArray * histList = new TObjArray(0);
  histList->Add(GeoEff2);
  histList->Add(Descrip);

  TFile* saveFile = new TFile("Efficiency/"+ pol+ saveName, "recreate");
  histList->Write();
  saveFile->Close();
 

  return;
}





