
#define DTFSystErr_cpp
#include "DTFSystErr.h"
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
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooPlot.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "TStyle.h"
#include "TPaveLabel.h"
#include <TAttFill.h>
#include <TColor.h>
#include <TAttLine.h>
#include <Rtypes.h>
#include <TString.h>
#include <THistPainter.h>
#include <TStyle.h>


using namespace RooFit;




int main(){
  DTFSystErr DTFSystClass;
  //DTFSystClass.CalculateEfficiency();
  TString MCTree = "/Volumes/TimeMachineBackups/LHCb/Analysis/DoubleJpsi/ntuples/2011/Analysis/Stripping20r1/MicroDSTDiMuonDiMuonIncLine/MC/MySingleMC/Up/Tuple.root";
  DTFSystClass.CalculateEfficiencyMC(MCTree);
  
  TString  filename1 =  "Output/DoubleDTF";
  TString  savename =  "Output/Comparision_DoubleDTF_TopAncestor_withNonPrompt.eps";
  //DTFSystClass.CalculateEfficiencyDoubleJpsiMC(filename1);
  //DTFSystClass.CompareDoubleDTF(filename1, savename);

  return 0;
  
}



//=========================================================================================================
//
// Plot two existing Double DTF efficency plots on same plot
//
//=========================================================================================================
void DTFSystErr::CompareDoubleDTF(TString filename1, TString savename){

  TFile *f1 = TFile::Open("/Users/phadc/LHCb/Analysis/DoubleJpsi/Output/2011/Stripping20r1/DTFEfficiency/DoubleDTFdfg.root") ;
  TH1D *efficiencyData = (TH1D*) f1->Get("Cut efficiency fitter/Cut Efficiency/Efficiency DTF #chi^{2}");
  TFile *f2 = TFile::Open(filename1+".root") ;
  TH1D *efficiencyMC = (TH1D*) f2->Get("DTFefficiency");
  TH1D *efficiencyMCNonPrompt = (TH1D*) f2->Get("backgroundEfficiency");

  // TFile *f4 = TFile::Open("DoubleDTF_1PV.root") ;
  //TH1D *efficiencyNPV1 = (TH1D*) f4->Get("DTFefficiency");

  TCanvas *effCan = new TCanvas("effCan");
  
  efficiencyData->SetMinimum(0.0);
  efficiencyData->SetStats(0);
  efficiencyData->GetXaxis()->SetTitle("DTF(Double)");
  efficiencyData->GetYaxis()->SetTitle("Efficiency");
  efficiencyData->SetLineColor(kRed);
  efficiencyData->SetMarkerSize(0.4);
  efficiencyData->Draw("E0");

  efficiencyMC->SetLineColor(kBlue);
  efficiencyMC->SetMarkerSize(0.4);
  efficiencyMC->Draw("same");
  
  efficiencyMCNonPrompt->SetLineColor(kBlack);
  efficiencyMCNonPrompt->SetMarkerSize(0.4);
  efficiencyMCNonPrompt->Draw("same");
   
  //efficiencyNPV1->SetLineColor(kBlack);
  //efficiencyNPV1->SetMarkerSize(0.4);
  //efficiencyNPV1->Draw("same");

  effCan->SaveAs(savename); 

  return;

}


//=========================================================================================================
//
// Calculate DTF efficiency from data.
//
//=========================================================================================================
void DTFSystErr::CalculateEfficiency()
{
 gStyle->SetOptStat(0);

 TString treeName;
 treeName = "/Volumes/TimeMachineBackups/LHCb/Analysis/DoubleJpsi/ntuples/2011/Analysis/Stripping20r1/MicroDSTDiMuonDiMuonIncLine/Down/Tuple_0x5A0032_20r1_Down.root";
 
 LHCbStyle();
 gStyle->SetPaintTextFormat("0.2f");//limits decimal places COLZText
 gStyle->SetOptStat(0);//
 
 TFile* file = new TFile(treeName);//MC
 TTree* tree=(TTree*)file->Get("Jpsi2MuMu/AllJpsi"); //MC filter tree
 
 double MCeff, MCeff2,MCeff3,MCeff4,MCeff5,MCeff6,MCeff7,MCeff8,MCeff9,MCeff10,MCeff11,MCeff12;
 double MCeff13,MCeff14,MCeff15,MCeff16,MCeff17,MCeff18,MCeff19,MCeff20,MCeff21,MCeff22,MCeff23;
 double MCeff24,MCeff25;
 Int_t MCeff26;
 double Theta, deltaPhi;


 tree->SetBranchAddress("Jpsi_PT",&MCeff);
 tree->SetBranchAddress("Jpsi_RAPIDITY",&MCeff2);
 tree->SetBranchAddress("Jpsi_COSTHETA",&MCeff3);
 tree->SetBranchAddress("MuPlus_TRACK_GHOSTPROB",&MCeff4);
 tree->SetBranchAddress("MuMinus_TRACK_GHOSTPROB",&MCeff5);
 tree->SetBranchAddress("MuPlus_PSEUDORAPIDITY",&MCeff6);
 tree->SetBranchAddress("MuMinus_PSEUDORAPIDITY",&MCeff7);
 tree->SetBranchAddress("MuPlus_ProbNNMu",&MCeff8);
 tree->SetBranchAddress("MuMinus_ProbNNMu",&MCeff9);
 tree->SetBranchAddress("Jpsi_DTF_CHI2PERNDOF",&MCeff13);
 tree->SetBranchAddress("Jpsi_VCHI2",&MCeff14);
 tree->SetBranchAddress("MuPlus_PT",&MCeff15);
 tree->SetBranchAddress("MuMinus_PT",&MCeff16);
 tree->SetBranchAddress("MuPlus_KL",&MCeff17);
 tree->SetBranchAddress("MuMinus_KL",&MCeff18);
 tree->SetBranchAddress("MuPlus_CHI2PERNDOF",&MCeff19);
 tree->SetBranchAddress("MuMinus_CHI2PERNDOF",&MCeff20);
 tree->SetBranchAddress("Jpsi_MASS",&MCeff22);
 tree->SetBranchAddress("Jpsi_P",&MCeff23);
 tree->SetBranchAddress("MuPlus_P",&MCeff24);
 tree->SetBranchAddress("MuMinus_P",&MCeff25);

 
  std::ostringstream name, filenameroot, filenameFit, plotAxis;
  name << "DTF";
  plotAxis << "DTF";
 
  filenameroot  <<  name.str().c_str() << "Significance.root";
  filenameFit << name.str().c_str() << "MassFit.root"  ;

  double binStart;
  double binEnd ;
  int binTotal ;
  double itStart, itEnd, itStep;
    
  binStart = 0.0;
  binEnd = 20.0;
  binTotal = 40;
  itStart = 0.5;
  itEnd = 20.05;
  itStep = 0.5;
  
  TH1D *efficiency = new TH1D("efficiency","efficiency",binTotal,binStart,binEnd);
 
  //
  //loop over cut of interest
  //
  TObjArray * histListFit = new TObjArray(0);

  for(double j=itStart; j<itEnd; j = j+itStep){

    std::ostringstream fitPlotName ;
    fitPlotName << "massPlot" << j ;  
    
    TFile *fileHold = new TFile("fileHold.root","recreate");
    TTree *treeHold = new TTree("treeHold","treeHold");
    treeHold->Branch("holdMass",&MCeff22);

    TFile *fileHold2 = new TFile("fileHold2.root","recreate");
    TTree *treeHold2 = new TTree("treeHold2","treeHold2");
    treeHold2->Branch("holdMass",&MCeff22);

    double countBefore = 0;
    double countAfter = 0;
 
    for(int i=0; i<1000000; i++)
      {
	tree->GetEntry(i);

	if(MCeff8 < 0.5) continue; //ProbNNMu positive 
	if(MCeff9 < 0.5) continue; //ProbNNMu negative 
	if(MCeff14 >20) continue; // Jpsi Vertex  	
	if(MCeff22 < 3000.0) continue; // Jpsi Mass
	if(MCeff22 > 3200.0) continue; // Jpsi Mass
	if(MCeff6 < 2.0) continue; //pseudorapMuPlus 
	if(MCeff6 > 5.0) continue; //pseudorapMuPlus
	if(MCeff7 < 2.0) continue; //pseudorapMuMinus  
	if(MCeff7 > 5.0) continue; //pseudoapMuMinus
	if(MCeff2 < 2.0)continue; //Jpsi rapidity  // 
	if(MCeff2 > 4.2)continue; //Jpsi rapdidity // NOTE!!! 
	if(MCeff >10000.0)continue; //Jpsi Pt
	if(MCeff15 <650) continue; // MuPlus Pt (default 650)
	if(MCeff16 <650) continue; // MuMinus Pt (default 650)
	if(MCeff19 >1.7) continue; // MuPlus Track Chi2 
	if(MCeff20 >1.7) continue; // MuMinus Track Chi2 
	if(MCeff17 != -1) continue; // MuPlus KL
	if(MCeff18 != -1) continue; // MuMinus KL
	if(MCeff4 > 0.05) continue; //ghostMuPlus 
	if(MCeff5 > 0.05) continue; //ghostMuMinus
	if(MCeff24 < 6000.0)continue; //MuPlus momentum cut 
	if(MCeff25 < 6000.0)continue; //MuMinus momentum cut
	if(MCeff24 > 200000.0)continue; //MuPlus momentum cut 
	if(MCeff25 > 200000.0)continue; //MuMinus momentum cut
	
	countBefore++;
	
	treeHold->Fill();
	if(MCeff13 > j ) continue; // Jpsi DTF  - default 11.0
	treeHold2->Fill();
	countAfter++;	
      }
          
    //
    // Fit 'before' and 'after'
    //
    TString fitPlotName2 = fitPlotName.str().c_str(); 
    double before = FitDCB(treeHold,fitPlotName2, j);
    double after = FitDCB(treeHold2,fitPlotName2, j);

    double DTFEffVal = after/before;  
    std::cout << "count before " << countBefore  << "signal before** " << before  << std::endl;
    std::cout << "count after " << countAfter << "signal after** " << after  << std::endl;
    std::cout << " efficiency is " << DTFEffVal << std::endl;
    double jplot =1;
    jplot =  j*2.0;   
    int jplotInt = (int) round( jplot); 
    efficiency->SetBinContent(jplotInt, DTFEffVal);
    efficiency->SetBinError(jplotInt, CalcErrData(countBefore, countAfter)); //wrong method - needs changing
  }

  TObjArray * histList2 = new TObjArray(0);
  TFile* fitFile = new TFile(filenameroot.str().c_str(), "recreate");
  histList2->Add(efficiency);
  histList2->Write();
  fitFile->Close();  
  file->Close();

  return;
}




//=========================================================================================================
//
// Double Crystal Ball + exp fitting function.
//
//=========================================================================================================


double DTFSystErr::FitDCB(TTree *treeHold,  TString fitPlotName, double j )
{
    RooRealVar holdMass("holdMass", "holdMass", 3000.0,3200.0); 
    RooDataSet * data = new RooDataSet("data","data",treeHold, holdMass);

    RooRealVar cbmean("cbmean", "cbmean" ,3094.45, 3090.0 , 3099.45 ) ; 
    RooRealVar cbsigma("cbsigma", "cbsigma",13.48 , 00.48, 16.48) ;
    RooRealVar alpha("alpha","alpha", 1.7,1.7, 1.7); //Default
    RooRealVar n("n","n", 2.7, 2.7,2.7); //default
    RooRealVar alpha2("alpha2","alpha2", -1.85,-1.85, -1.85); // default
    RooRealVar n2("n2","n2", 3.0, 3.0,3.0); //default
      
    RooCBShape cball("cball", "crystal ball", holdMass, cbmean, cbsigma, alpha, n);
    RooCBShape cball2("cball2", "crystal ball2", holdMass, cbmean, cbsigma, alpha2, n2);  
    RooRealVar c1("c1","backgroundslope",0.0,-100.0,100.0); //default
    RooExponential p0("p0","exp function for background",holdMass,c1);
 
    RooRealVar sig_yield("N_{Signal}", "yield signal peak", 200, 0, 1000000);
    RooRealVar bkgd_yield("N_{BackGnd}", "yield of background", 10, 0, 1000000);
 

    //
    //add together the two crystal balls
    //
    RooArgList crystalshape;
    crystalshape.add(cball);
    crystalshape.add(cball2);
    RooRealVar Ratio("ratio","ratio",0.5,0.0,1.0);
    RooAddPdf cballsum("cballsum","cballsum",cball,cball2,Ratio);
     
    gStyle->SetOptStat(0);//
    
    RooArgList pdfshapes;
    RooArgList yields; 
    pdfshapes.add(cballsum); 
    pdfshapes.add(p0);
    yields.add(sig_yield);
    yields.add(bkgd_yield);
     std::cout << "point 5" << std::endl;
    RooAddPdf sum("sum","sum of signal and background",pdfshapes, yields);

    TCanvas* Optimiseplot = new TCanvas(fitPlotName);
    std::cout << "point 4" << std::endl;

    RooPlot * xframe;
    xframe = holdMass.frame();   
    sum.fitTo(*data);
    sum.plotOn(xframe); 
      
    RooArgSet plotStuff("plotStuff");
    plotStuff.add(cbmean);
    plotStuff.add(cbsigma);
    plotStuff.add(sig_yield);
    plotStuff.add(bkgd_yield);
    plotStuff.add(alpha);
    plotStuff.add(alpha2);
    plotStuff.add(n);
    plotStuff.add(n2);
    plotStuff.add(c1);
 
    gStyle->SetPaintTextFormat(".3f");
    
    sum.plotOn(xframe, RooFit::Components(cballsum), RooFit::LineColor(kBlue));
    sum.plotOn(xframe, RooFit::Components(p0), RooFit::LineColor(kCyan));
    sum.plotOn(xframe, RooFit::LineColor(kRed));
 
    xframe->Draw();   
    
    std::cout << "Fit DCB signal yield " << sig_yield.getVal() << "Fit DCB background " << bkgd_yield.getVal() 
	      << " total " <<  bkgd_yield.getVal()+sig_yield.getVal() << std::endl;
    std::cout << "data entries " << data->numEntries() << std::endl;
  
    double significance = sig_yield.getVal()/std::sqrt(sig_yield.getVal()+bkgd_yield.getVal());
    double purity = sig_yield.getVal()/(sig_yield.getVal()+bkgd_yield.getVal());
    std::cout << "Cut value is " <<j << std::endl;
 
    Optimiseplot->SaveAs("lastmassplot.root");
    double signal = sig_yield.getVal() ;
    return signal;
}



//=========================================================================================================
//
// Takes single jpsi MC and calculates a DTF efficiency (both weighted and unweighted).
//
//=========================================================================================================

void DTFSystErr::CalculateEfficiencyMC(TString MCTree )
{

  LHCbStyle();
  gStyle->SetPaintTextFormat("0.2f");//limits decimal places COLZText
  gStyle->SetOptStat(0);//
  double colzScaleMax = 1.0;
 
  TFile* file = new TFile(MCTree);//MC
  TTree* tree=(TTree*)file->Get("Jpsi2MuMu/AllJpsi"); //MC filter tree
 

  double MCeff, MCeff2,MCeff3,MCeff4,MCeff5,MCeff6,MCeff7,MCeff8,MCeff9,MCeff10,MCeff11,MCeff12;
  double MCeff13,MCeff14,MCeff15,MCeff16,MCeff17,MCeff18,MCeff19,MCeff20,MCeff21,MCeff22,MCeff23;
  double MCeff24,MCeff25;
  Int_t MCeff26, MCeff27;
  double Theta, deltaPhi;
  tree->SetBranchAddress("Jpsi_Theta",&Theta);
  tree->SetBranchAddress("Jpsi_Phi",&deltaPhi);
  tree->SetBranchAddress("Jpsi_PT",&MCeff);
  tree->SetBranchAddress("Jpsi_RAPIDITY",&MCeff2);
  tree->SetBranchAddress("Jpsi_COSTHETA",&MCeff3);
  tree->SetBranchAddress("MuPlus_TRACK_GHOSTPROB",&MCeff4);
  tree->SetBranchAddress("MuMinus_TRACK_GHOSTPROB",&MCeff5);
  tree->SetBranchAddress("MuPlus_Pseudorapidity",&MCeff6);
  tree->SetBranchAddress("MuMinus_Pseudorapidity",&MCeff7);
  tree->SetBranchAddress("MuPlus_ProbNNMu",&MCeff8);
  tree->SetBranchAddress("MuMinus_ProbNNMu",&MCeff9);
  tree->SetBranchAddress("Jpsi_DTF_CHI2PERNDOF",&MCeff13);
  tree->SetBranchAddress("Jpsi_VCHI2",&MCeff14);
  tree->SetBranchAddress("MuPlus_PT",&MCeff15);
  tree->SetBranchAddress("MuMinus_PT",&MCeff16);
  tree->SetBranchAddress("MuPlus_KL",&MCeff17);
  tree->SetBranchAddress("MuMinus_KL",&MCeff18);
  tree->SetBranchAddress("MuPlus_CHI2PERNDOF",&MCeff19);
  tree->SetBranchAddress("MuMinus_CHI2PERNDOF",&MCeff20);
  tree->SetBranchAddress("Jpsi_MASS",&MCeff22);
  tree->SetBranchAddress("Jpsi_P",&MCeff23);
  tree->SetBranchAddress("MuPlus_P",&MCeff24);
  tree->SetBranchAddress("MuMinus_P",&MCeff25);
  tree->SetBranchAddress("Event_BestTracksMult",&MCeff26);
  tree->SetBranchAddress("Jpsi_TOP_ANCESTOR_PID",&MCeff27);

  //
  //Get Weighting plots 
  //
  TFile *f5 = TFile::Open("/Users/phadc/LHCb/Analysis/DoubleJpsi/Input/2011/Stripping20r1/MicroDSTDiMuonDiMuonIncLine/0x5A0032-PT-Weights.root");
  TFile *f6 = TFile::Open("/Users/phadc/LHCb/Analysis/DoubleJpsi/Input/2011/Stripping20r1/MicroDSTDiMuonDiMuonIncLine/0x5A0032-Rapidity-Weights.root");
  TFile *f7 = TFile::Open("/Users/phadc/LHCb/Analysis/DoubleJpsi/Input/2011/Stripping20r1/MicroDSTDiMuonDiMuonIncLine/0x5A0032-BestTracksMultiplicity-Weights.root");
  TH1D *ptWeightPlot = (TH1D*)f5->Get("P_{T} Weight");
  TH1D *rapWeightPlot =(TH1D*)f6->Get("Rapidity Weight");
  TH1D *multWeightPlot =(TH1D*)f7->Get("Multiplicity Weight");
  
  std::cout <<"number of entries " << tree->GetEntries()   << std::endl;
  
  double binStart;
  double binEnd ;
  int binTotal ;
  double itStart, itEnd, itStep;
 
  binStart = 0.0;
  binEnd = 20.0;
  binTotal = 40;
  itStart = 0.5;
  itEnd = 20.05;
  itStep = 0.5;    

  TH1D *efficiencyPlotMC = new TH1D("efficiencyPlotMC","efficiencyPlotMC",binTotal,binStart,binEnd);
  TH1D *efficiencyPlotMCWeighted = new TH1D("efficiencyPlotMCWeighted","efficiencyPlotMCWeighted",binTotal,binStart,binEnd);

  for(double j=itStart; j<itEnd; j = j+itStep){

    double countBefore = 0;
    double countAfter = 0;
    double count_minus1 = 0;
    //
    //hack to get around problem with bin filling
    //
    double jplot =  j*2.0;   
    int jplotInt = (int) round( jplot); 
    // std::cout << "jplotInt "<< jplotInt  << " j "<< j  << std::endl;

    double totalWeightBefore = 0;
    double totalWeightAfter = 0;
    double countNonPrompt = 0;
    double countPrompt = 0;

      for(int i=0; i<tree->GetEntries(); i++)
      {
	
	if(MCeff27 == -1) count_minus1++;
	tree->GetEntry(i);

	if(MCeff8 < 0.5) continue; //ProbNNMu positive 
	if(MCeff9 < 0.5) continue; //ProbNNMu negative 
	if(MCeff14 >20) continue; // Jpsi Vertex  	
	if(MCeff22 < 3000.0) continue; // Jpsi Mass
	if(MCeff22 > 3200.0) continue; // Jpsi Mass
	if(MCeff6 < 2.0) continue; //pseudorapMuPlus 
	if(MCeff6 > 5.0) continue; //pseudorapMuPlus
	if(MCeff7 < 2.0) continue; //pseudorapMuMinus  
	if(MCeff7 > 5.0) continue; //pseudoapMuMinus
	if(MCeff2 < 2.0)continue; //Jpsi rapidity  // 
	if(MCeff2 > 4.2)continue; //Jpsi rapdidity // NOTE!!! 
	if(MCeff >10000.0)continue; //Jpsi Pt
	if(MCeff15 <650) continue; // MuPlus Pt (default 650)
	if(MCeff16 <650) continue; // MuMinus Pt (default 650)
	if(MCeff19 >1.7) continue; // MuPlus Track Chi2 
	if(MCeff20 >1.7) continue; // MuMinus Track Chi2 
	if(MCeff17 != -1) continue; // MuPlus KL
	if(MCeff18 != -1) continue; // MuMinus KL
	if(MCeff4 > 0.05) continue; //ghostMuPlus 
	if(MCeff5 > 0.05) continue; //ghostMuMinus
	if(MCeff24 < 6000.0)continue; //MuPlus momentum cut 
	if(MCeff25 < 6000.0)continue; //MuMinus momentum cut
	if(MCeff24 > 200000.0)continue; //MuPlus momentum cut 
	if(MCeff25 > 200000.0)continue; //MuMinus momentum cut

	bool isSignalMC = MCTruthTopAncestor(MCeff27); // truth info
	if(!isSignalMC)continue;
	if(!isSignalMC) countNonPrompt++;
	else countPrompt++;

	//
	//Info to weight by pt, rapdity, best track multiplicity
	//
	int binPt = (ptWeightPlot->GetXaxis())->FindBin(MCeff);
	int binRap = (rapWeightPlot->GetXaxis())->FindBin(MCeff2);
	int binMult = (multWeightPlot->GetXaxis())->FindBin(MCeff26); 
	double reWeighting =
        (ptWeightPlot->GetBinContent(binPt))*(rapWeightPlot->GetBinContent(binRap))*(multWeightPlot->GetBinContent(binMult));
	totalWeightBefore = totalWeightBefore + reWeighting;
	countBefore++;
	if(MCeff13 > j ) continue; // Jpsi DTF, default 11

	totalWeightAfter = totalWeightAfter + reWeighting;
       	countAfter++;

      }
  
    double efficiency  =  countAfter/countBefore ;  
    std::cout << "DTF cut val " << j << " unweighted efficiency "   << efficiency << " (raw numbers " <<countBefore << " "<< countAfter << " ) " << std::endl;
    double weightedEfficiency  = totalWeightAfter / totalWeightBefore;
    std::cout << "DTF cut val " << j << " weighted efficiency "   << weightedEfficiency << std::endl;
    std::cout << "number nonprompt " << countNonPrompt  << " prompt " << countPrompt << std::endl;
  
    efficiencyPlotMC->SetBinContent(jplotInt, efficiency);
    efficiencyPlotMC->SetBinError(jplotInt, CalcErrMC(countBefore, countAfter));
    efficiencyPlotMCWeighted->SetBinContent(jplotInt, weightedEfficiency);
    efficiencyPlotMCWeighted->SetBinError(jplotInt, CalcErrMC(totalWeightBefore, totalWeightAfter));
  }

  TObjArray * histList2 = new TObjArray(0);
  TFile* fitFile = new TFile("DTFMC.root", "recreate"); //change
  histList2->Add(efficiencyPlotMC);
  histList2->Add(efficiencyPlotMCWeighted);
  histList2->Write();
  fitFile->Close();

  TFile *f4 = TFile::Open("DTFSignificance.root") ;
  TH1D *efficiencyData = (TH1D*) f4->Get("efficiency");
  TCanvas *effCan = new TCanvas("effCan");
  efficiencyData->SetMinimum(0.0);
  efficiencyData->SetStats(0);
  efficiencyData->GetXaxis()->SetTitle("DTF(single)");
  efficiencyData->GetYaxis()->SetTitle("Efficiency");
  efficiencyData->SetLineColor(kRed);
  efficiencyData->SetMarkerSize(0.4);
  efficiencyData->Draw("E0");
  //efficiencyPlotMC->SetLineColor(kGreen);
  //efficiencyPlotMC->SetMarkerSize(0.4);
  //efficiencyPlotMC->Draw("same");
  efficiencyPlotMCWeighted->SetLineColor(kBlue);
  efficiencyPlotMCWeighted->SetMarkerSize(0.4);
  efficiencyPlotMCWeighted->Draw("same");

  effCan->SaveAs("comparison.eps"); //change
  return;
  
  
}



//=========================================================================================================
//
// Takes Double jpsi MC and calculates a DTF efficiency.
//
//=========================================================================================================
void DTFSystErr::CalculateEfficiencyDoubleJpsiMC(TString filename1)
{
 TString filename, twoDplot, MCTree, MCTruth, MCTreeDown, MCTruthDown, folder;
 twoDplot = "MagUp2D.eps";
  filename = "20r1RecSel2011MagUpNoPol";
  MCTree = "/Volumes/TimeMachineBackups/LHCb/Analysis/DoubleJpsi/ntuples/2011/Analysis/Stripping20r1/MicroDSTDiMuonDiMuonIncLine/MC/MyDoubleJpsi/New/Up/Tuple.root";
  MCTruth = "/Volumes/TimeMachineBackups/LHCb/Analysis/DoubleJpsi/ntuples/2011/Analysis/Stripping20r1/MicroDSTDiMuonDiMuonIncLine/MC/MyDoubleJpsi/Old/Up/TruthTuple.root";

  int  threeDybin = 15; 

  LHCbStyle();
  gStyle->SetPaintTextFormat("0.2f");//limits decimal places COLZText
  gStyle->SetOptStat(0);//

  
  TFile* file = new TFile(MCTree);//MC
  TFile* file2 = new TFile(MCTruth);//True
  TTree* tree=(TTree*)file->Get("Jpsi2MuMuTruth/AllDoubleJpsi"); //MC filter tree
  TTree* tree2 = (TTree*)file2->Get("MyJpsi2MuMu/AllDoubleMCTrueJpsi"); //MC Truth filter tree
  std::cout << "point1" << std::endl;
  
  double MCeff, MCeff2,MCeff3,MCeff4,MCeff5,MCeff6,MCeff7,MCeff8,MCeff9,MCeff10,MCeff11,MCeff12;
  double MCeff13,MCeff14,MCeff15,MCeff16,MCeff17,MCeff18,MCeff19,MCeff20,MCeff21,MCeff22,MCeff23;
  double MCeff24,MCeff25, MCeff27;
  Int_t jpsi1_topAnc;
  Int_t MCeff26;
  double Theta, deltaPhi;
  tree->SetBranchAddress("Jpsi_1_Theta",&Theta); //----
  tree->SetBranchAddress("Jpsi_1_Phi",&deltaPhi); //-----
  tree->SetBranchAddress("JPsi_1_PT",&MCeff);
  tree->SetBranchAddress("JPsi_1_RAPIDITY",&MCeff2);
  tree->SetBranchAddress("JPsi_1_COSTHETA",&MCeff3);
  tree->SetBranchAddress("MuPlus1_TRACK_GHOSTPROB",&MCeff4);
  tree->SetBranchAddress("MuMinus1_TRACK_GHOSTPROB",&MCeff5);
  tree->SetBranchAddress("MuPlus1_Pseudorapidity",&MCeff6);
  tree->SetBranchAddress("MuMinus1_Pseudorapidity",&MCeff7);
  tree->SetBranchAddress("MuPlus1_ProbNNMu",&MCeff8);
  tree->SetBranchAddress("MuMinus1_ProbNNMu",&MCeff9);
  tree->SetBranchAddress("JPsi_1_DTF_CHI2PERNDOF",&MCeff13); //----
  tree->SetBranchAddress("JPsi_1_VCHI2",&MCeff14);
  tree->SetBranchAddress("MuPlus1_PT",&MCeff15);
  tree->SetBranchAddress("MuMinus1_PT",&MCeff16);
  tree->SetBranchAddress("MuPlus1_KL",&MCeff17);
  tree->SetBranchAddress("MuMinus1_KL",&MCeff18);
  tree->SetBranchAddress("MuPlus1_CHI2PERNDOF",&MCeff19);
  tree->SetBranchAddress("MuMinus1_CHI2PERNDOF",&MCeff20);
  tree->SetBranchAddress("JPsi_1_MASS",&MCeff22);
  tree->SetBranchAddress("JPsi_1_P",&MCeff23);
  tree->SetBranchAddress("MuPlus1_P",&MCeff24);
  tree->SetBranchAddress("MuMinus1_P",&MCeff25);
  tree->SetBranchAddress("JPsi_1_TOP_ANCESTOR_PID",&jpsi1_topAnc);
  tree->SetBranchAddress("Event_BestTracksMult",&MCeff26); //----
  tree->SetBranchAddress("DoubleJpsi_DTF_CHI2PERNDOF",&MCeff27); //Double Jpsi DTF


  double MC2eff, MC2eff2,MC2eff3,MC2eff4,MC2eff5,MC2eff6,MC2eff7,MC2eff8,MC2eff9,MC2eff10,MC2eff11,MC2eff12;
  double MC2eff13,MC2eff14,MC2eff15,MC2eff16,MC2eff17,MC2eff18,MC2eff19,MC2eff20,MC2eff21,MC2eff22,MC2eff23;
  double MC2eff24,MC2eff25;
  Int_t jpsi2_topAnc;
  Int_t MC2eff26;
  double Theta2, deltaPhi2;
  tree->SetBranchAddress("Jpsi_2_Theta",&Theta2); //-----
  tree->SetBranchAddress("Jpsi_2_Phi",&deltaPhi2); //-----
  tree->SetBranchAddress("JPsi_2_PT",&MC2eff);
  tree->SetBranchAddress("JPsi_2_RAPIDITY",&MC2eff2);
  tree->SetBranchAddress("JPsi_2_COSTHETA",&MC2eff3);
  tree->SetBranchAddress("MuPlus2_TRACK_GHOSTPROB",&MC2eff4);
  tree->SetBranchAddress("MuMinus2_TRACK_GHOSTPROB",&MC2eff5);
  tree->SetBranchAddress("MuPlus2_Pseudorapidity",&MC2eff6);
  tree->SetBranchAddress("MuMinus2_Pseudorapidity",&MC2eff7);
  tree->SetBranchAddress("MuPlus2_ProbNNMu",&MC2eff8);
  tree->SetBranchAddress("MuMinus2_ProbNNMu",&MC2eff9);
  tree->SetBranchAddress("JPsi_2_DTF_CHI2PERNDOF",&MC2eff13); //------
  tree->SetBranchAddress("JPsi_2_VCHI2",&MC2eff14);
  tree->SetBranchAddress("MuPlus2_PT",&MC2eff15);
  tree->SetBranchAddress("MuMinus2_PT",&MC2eff16);
  tree->SetBranchAddress("MuPlus2_KL",&MC2eff17);
  tree->SetBranchAddress("MuMinus2_KL",&MC2eff18);
  tree->SetBranchAddress("MuPlus2_CHI2PERNDOF",&MC2eff19);
  tree->SetBranchAddress("MuMinus2_CHI2PERNDOF",&MC2eff20);
  tree->SetBranchAddress("JPsi_2_MASS",&MC2eff22);
  tree->SetBranchAddress("JPsi_2_P",&MC2eff23);
  tree->SetBranchAddress("MuPlus2_P",&MC2eff24);
  tree->SetBranchAddress("MuMinus2_P",&MC2eff25);
  tree->SetBranchAddress("JPsi_2_TOP_ANCESTOR_PID",&jpsi2_topAnc);

  Int_t MC_npv;

  tree->SetBranchAddress("Event_nPV",&MC_npv);

  double effTruth, eff2Truth,eff3Truth, eff4Truth, eff5Truth, eff6Truth, eff8Truth,eff9Truth,eff10Truth ;
  double TruthTheta, TruthPhi;
  Int_t eff7Truth;
  //  tree2->SetBranchAddress("Jpsi1_Theta",&TruthTheta); // NEEDS ADDING
  // tree2->SetBranchAddress("Jpsi_1_Phi",&TruthPhi); //NEEDS Adding to tuple
  tree2->SetBranchAddress("JPsi_1_PT",&effTruth); 
  tree2->SetBranchAddress("JPsi_1_RAPIDITY",&eff2Truth);
  tree2->SetBranchAddress("JPsi_1_COSTHETA",&eff3Truth);
  tree2->SetBranchAddress("JPsi_1_P",&eff4Truth);
  tree2->SetBranchAddress("MuPlus1_P",&eff5Truth);
  tree2->SetBranchAddress("MuMinus1_P",&eff6Truth);
  tree2->SetBranchAddress("Event_BestTracksMult",&eff7Truth);
  tree2->SetBranchAddress("MuMinus1_PSEUDORAPIDITY",&eff8Truth);
  tree2->SetBranchAddress("MuPlus1_PSEUDORAPIDITY",&eff9Truth);
  tree2->SetBranchAddress("JPsi_1_MASS",&eff10Truth);

 
  double effTruth2, eff2Truth2,eff3Truth2, eff4Truth2, eff5Truth2, eff6Truth2, eff8Truth2,eff9Truth2,eff10Truth2 ;
  double TruthTheta2, TruthPhi2;
  Int_t eff7Truth2;
  
  tree2->SetBranchAddress("JPsi_2_PT",&effTruth2); 
  tree2->SetBranchAddress("JPsi_2_RAPIDITY",&eff2Truth2);
  tree2->SetBranchAddress("JPsi_2_COSTHETA",&eff3Truth2);
  tree2->SetBranchAddress("JPsi_2_P",&eff4Truth2);
  tree2->SetBranchAddress("MuPlus2_P",&eff5Truth2);
  tree2->SetBranchAddress("MuMinus2_P",&eff6Truth2);
  tree2->SetBranchAddress("MuMinus2_PSEUDORAPIDITY",&eff8Truth2);
  tree2->SetBranchAddress("MuPlus2_PSEUDORAPIDITY",&eff9Truth2);
  tree2->SetBranchAddress("JPsi_2_MASS",&eff10Truth2);


  double binStart = 0.0; //DTF = 0.0 // Mu-track = 0 // ghostProb = 0.05 
  double binEnd = 20.0; //DTF = 30.0 //Mu-track = 5.5// ghostProb = 1.0
  //int binTotal = 160; //DTF = 30 // Mu-track = 18// ghostProb = 20
  int binTotal = 100; //DTF = 30 // Mu-track = 18// ghostProb = 20


  TH1D *sig = new TH1D("significance","significance",binTotal,binStart,binEnd);
  TH1D *signal = new TH1D("signal","signal",binTotal,binStart,binEnd);
  TH1D *back = new TH1D("background","background",binTotal,binStart,binEnd);
  TH1D *purity = new TH1D("purity","purity",binTotal,binStart,binEnd);
  TH1D *DTFefficiency = new TH1D("DTFefficiency","DTFefficiency",binTotal,binStart,binEnd);
  TH1D *totalEfficiency = new TH1D("totalEfficiency","totalEfficiency",binTotal,binStart,binEnd);
  TH1D *backgroundEfficiency = new TH1D("backgroundEfficiency","backgroundEfficiency",binTotal,binStart,binEnd);
  double countTruth = 0.0;
   
 for(int j=1; j<tree2->GetEntries(); j++)
  {
    tree2->GetEntry(j);
    //
    //Jpsi1
    //
    if(eff2Truth < 2.0)continue; //Jpsi rapidity  
    if(eff2Truth > 4.2)continue; //Jpsi rapdidity
    if(effTruth > 10000.0)continue; //jpsi pt
    if(eff10Truth < 3000.0) continue; // Jpsi Mass
    if(eff10Truth > 3200.0) continue; // Jpsi Mass
  
    //
    //Jpsi2
    //
    if(eff2Truth2 < 2.0)continue; //Jpsi rapidity  
    if(eff2Truth2 > 4.2)continue; //Jpsi rapdidity
    if(effTruth2 > 10000.0)continue; //jpsi pt
    if(eff10Truth2 < 3000.0) continue; // Jpsi Mass
    if(eff10Truth2 > 3200.0) continue; // Jpsi Mass
    countTruth++;
    
  }
 
 int count_npv1 = 0;

  //
  //loop over cut of interest
  //
 double jplot = 0;

  //for(double j=0.5; j<5.25; j = j+0.25){ //Mu track
  // for(double j=0.05; j<1.05; j = j+0.05){//ghostProb
 //  for(double j=1.0; j<31; j++){   //DTF
 for(double j=0.0; j<20; j= j+ 0.2){   //DTF
    double countSignalPreCut = 0;
    double countBackgroundPreCut = 0;

    double countSignal = 0;
    double countBackground = 0;
    double countTest1 = 0;
    double countTest2 = 0;
    double countTest3 = 0;

    double count_minus1 = 0;
    double count_4 = 0;
    double count_5 = 0;

    jplot++;
    //for(int i=0; i<tree->GetEntries(); i++)
    for(int i=0; i<100000; i++)
    {
      tree->GetEntry(i);
      // if(jpsi1_topAnc == 5 || jpsi1_topAnc == 4 ){
      //	MCTruthTopAncestor(jpsi1_topAnc);
      // }
      //
      //Cut on Number of pv for test - TAKE OUT WHEN NOT TESTING!
      //
      // if(MC_npv != 1) { 
      //	count_npv1++;
      //	continue;
      //}

      //
      //Jpsi1
      //
      if(jpsi1_topAnc == -1) count_minus1++;
      if(abs(jpsi1_topAnc) == 4) count_4++;
      if(abs(jpsi1_topAnc) == 5) count_5++;

      if(MCeff21 == 0.0) countTest1++;
      if(MC2eff21 == 0.0) countTest2++;
      if(MCeff21 == 0.0 || MC2eff21 == 0.0)  countTest3++;

      if(MCeff14 >20) continue; // Jpsi Vertex
      if(MCeff13 > 11.0 ) continue; // Jpsi DTF: CHANGE BACK TO 11!
      if(MCeff22 < 3000.0) continue; // Jpsi Mass
      if(MCeff22 > 3200.0) continue; // Jpsi Mass
      if(MCeff6 < 2.0) continue; //pseudorapMuPlus 
      if(MCeff6 > 5.0) continue; //pseudorapMuPlus
      if(MCeff7 < 2.0) continue; //pseudorapMuMinus  
      if(MCeff7 > 5.0) continue; //pseudoapMuMinus
      if(MCeff2 < 2.0)continue; //Jpsi rapidity 
      if(MCeff2 > 4.2)continue; //Jpsi rapdidity 
      if(MCeff >10000.0)continue; //Jpsi Pt
      if(MCeff15 <650) continue; // MuPlus Pt
      if(MCeff16 <650) continue; // MuMinus Pt
      if(MCeff19 >1.7) continue; // MuPlus Track Chi2
      if(MCeff20 >1.7) continue; // MuMinus Track Chi2
      if(MCeff17 != -1) continue; // MuPlus KL
      if(MCeff18 != -1) continue; // MuMinus KL
      if(MCeff4 > 0.05) continue; //ghostMuPlus
      if(MCeff5 > 0.05) continue; //ghostMuMinus
      // if(MCeff21 == 0.0 ) continue; // Truth Info
      if(MCeff24 < 6000.0)continue; //MuPlus momentum cut
      if(MCeff25 < 6000.0)continue; //MuMinus momentum cut
      if(MCeff24 > 200000.0)continue; //MuPlus momentum cut
      if(MCeff25 > 200000.0)continue; //MuMinus momentum cut

      bool isSignalMC_1 = MCTruthTopAncestor(jpsi1_topAnc); // truth info
     

      //
      //Jpsi2
      //
      
      if(MC2eff14 >20) continue; // Jpsi Vertex
      if(MC2eff13 > 11.0 ) continue; // Jpsi DTF: CHANGE BACK TO 11!
      if(MC2eff22 < 3000.0) continue; // Jpsi Mass
      if(MC2eff22 > 3200.0) continue; // Jpsi Mass
      if(MC2eff6 < 2.0) continue; //pseudorapMuPlus 
      if(MC2eff6 > 5.0) continue; //pseudorapMuPlus
      if(MC2eff7 < 2.0) continue; //pseudorapMuMinus  
      if(MC2eff7 > 5.0) continue; //pseudoapMuMinus
      if(MC2eff2 < 2.0)continue; //Jpsi rapidity  
      if(MC2eff2 > 4.2)continue; //Jpsi rapdidity
      if(MC2eff >10000.0)continue; //Jpsi Pt
      if(MC2eff15 <650) continue; // MuPlus Pt
      if(MC2eff16 <650) continue; // MuMinus Pt
      if(MC2eff19 >1.7) continue; // MuPlus Track Chi2
      if(MC2eff20 >1.7) continue; // MuMinus Track Chi2
      if(MC2eff17 != -1) continue; // MuPlus KL
      if(MC2eff18 != -1) continue; // MuMinus KL
      if(MC2eff4 > 0.05) continue; //ghostMuPlus
      if(MC2eff5 > 0.05) continue; //ghostMuMinus
      // if(MC2eff21 == 0.0 ) continue; // Truth Info
      if(MC2eff24 < 6000.0)continue; //MuPlus momentum cut
      if(MC2eff25 < 6000.0)continue; //MuMinus momentum cut
      if(MC2eff24 > 200000.0)continue; //MuPlus momentum cut
      if(MC2eff25 > 200000.0)continue; //MuMinus momentum cut
      bool isSignalMC_2 = MCTruthTopAncestor(jpsi2_topAnc); // truth info
      // if(MCeff21 == 0.0 || MC2eff21 == 0.0 )  countBackgroundPreCut++;
      if(!isSignalMC_1 || !isSignalMC_2 )  countBackgroundPreCut++;
      else countSignalPreCut++;
      if(MCeff27 >j) continue; // Double Jpsi Vertex 

      //
      //get background and signal count from truth info
      //
      // if(MCeff21 == 0.0 || MC2eff21 == 0.0)  countBackground++;
      if(!isSignalMC_1 || !isSignalMC_2 ) countBackground++;
      else countSignal++;
    }

    std::cout <<"*********** count -1 " << count_minus1 << " frac  " << count_minus1/100000 << " count 4 " << count_4 << " count 5 " << count_5 << std::endl;
    std::cout << "number of n_pv = 1 events " << count_npv1 << std::endl;
    std::cout << " either jpsi1 or 2 background  " << countTest3 << std::endl;
    std::cout << "jpsi1 background " << countTest1 <<" jpsi2 background  "<<countTest2<<" total "<< tree->GetEntries()<<std::endl;
    std::cout <<"background1 " << countBackgroundPreCut <<  " back2 " << countBackground << " endSignal "<<countSignal<<std::endl;    
    std::cout << "DTF cut value is " <<j << std::endl;
    
    //  double jplot = (j+0.2)*5.0 ; //DTF = j, mutrack = (j*4)-2 //ghost = j*20 
    sig->SetBinContent(jplot,(countSignal/std::sqrt(countSignal+(countBackground)))); 
    signal->SetBinContent(jplot,countSignal); 
    back->SetBinContent(jplot,countBackground); 
    purity->SetBinContent(jplot,countSignal/(countSignal+countBackground)); 

    DTFefficiency->SetBinContent(jplot,countSignal/countSignalPreCut); 
    DTFefficiency->SetBinError(jplot,CalcErrMC(countSignalPreCut, countSignal)); 
    DTFefficiency->SetMarkerSize(0.5);

    totalEfficiency->SetBinContent(jplot,countSignal/countTruth);
    if(countBackgroundPreCut =! 0){
    backgroundEfficiency->SetBinContent(jplot, countBackground/countBackgroundPreCut);
    backgroundEfficiency->SetBinError(jplot, CalcErrMC(countBackgroundPreCut, countBackground));
    std::cout << "total efficiency " << countSignal/countTruth << " DTF Eff " << countSignal/countSignalPreCut << std::endl;
    }
  }
  
 //TString  filename1 =  "Output/DoubleDTF";

  TObjArray* histList = new TObjArray(0);
  LHCbStyle();
  gStyle->SetOptStat(0);//
  DTFefficiency->GetXaxis()->SetTitle("#Chi^{2}_{DTF}/nDoF");
  DTFefficiency->GetYaxis()->SetTitle("#epsilon^{#Chi^{2}_{DTF}}");
  DTFefficiency->SetTitle("");
  TCanvas *Canvas1 = new TCanvas("DTFEff");
  DTFefficiency->Draw();
  Canvas1->SaveAs(filename1 + ".eps");
  histList->Add(sig);
  histList->Add(signal);
  histList->Add(back);
  histList->Add(purity);
  histList->Add(totalEfficiency);
  histList->Add(DTFefficiency);
  histList->Add(backgroundEfficiency);
  
  TFile* andrewsFile = new TFile(filename1 + ".root", "recreate");
  histList->Write();
  andrewsFile->Close();

  return;

}





//=========================================================================================================
//
// Calculate error in MC - Note: wrong method for now 
//
//=========================================================================================================
double DTFSystErr::CalcErrMC(double countBefore, double countAfter){

  double kyield = (countAfter);
  double nyield = (countBefore);
  double var =fabs((((kyield+1)*(kyield+2))/((nyield+2)*(nyield+3)))-(((kyield+1)*(kyield+1))/((nyield+2)*(nyield+2))));
  double error = TMath::Sqrt(var);

  return error;
}


//=========================================================================================================
//
// Calculate error in Data - Note: Need to add in proper method 
//
//=========================================================================================================

double DTFSystErr::CalcErrData(double countBefore, double countAfter){

  double error  = CalcErrMC(countBefore, countAfter);

  return error;
}


//=========================================================================================================
//
// Provides MC-Truth info for MC program.
//
//=========================================================================================================

bool DTFSystErr::MCTruthTopAncestor(Int_t topAncestor)
{

  //
  //Check if prompt
  //
  bool isPrompt  = (topAncestor == 0) ? true : false;

  //
  //has bottom
  //
  
  int q = bottom;
  
  bool retVal = false;
  
  int count = 0;
  int pid   = abs(topAncestor) / 10 ;
  
  while (count < 4){
    if ((pid % 10) == q){
      retVal = true;
      break;
    }
    pid = pid / 10;
    ++count;
  }
  //  if (retVal) continue;
  
  //std::cout << " has bottom  " << retVal ;


  //
  //has Charm
  //
  
  int q2 = charm;
  
  bool retVal2 = false;
  
  int count2 = 0;
  int pid2   = abs(topAncestor) / 10 ;
  
  while (count2 < 4){
    if ((pid2 % 10) == q2){
      retVal2 = true;
      break;
    }
    pid2 = pid2 / 10;
    ++count2;
  }
  
  bool isFeedDown = false;
  if(retVal) isFeedDown = false;
  else if(retVal2) isFeedDown = true;
 
  //  if(!(isPrompt || isFeedDown)) continue;
  bool returnDecision = true;
  if(!(isPrompt || isFeedDown)) returnDecision = false;
  

  return returnDecision;
}




