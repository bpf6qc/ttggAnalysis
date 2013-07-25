#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRandom3.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <algorithm>
#include <utility>
#include <vector>
#include <stdio.h>
#include <stdarg.h>
#include <exception>

#include "analyze_mc.h"

using namespace std;

const TString ffColor = "kOrange+10";
const TString eeColor = "kBlue";
const TString egColor = "kGreen";

void analyze(TString input, bool addMC, int channel, int intLumi_int) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);

  const int nChannels = 7;
  TString channels[nChannels] = {
    "nojet",
    "j", "b",
    "bj",
    "muJets", // gg+mu+bj + X (dilep veto)
    "eleJets",
    "hadronic" // gg+5j1b + X (lep veto)
  };

  prep_signal(channels[channel]);

  TFile * in = new TFile(input, "READ");

  // pixel veto
  Float_t fakeRate = 0.0199;
  Float_t fakeRate_err = 0.0001;

  // electron conversion veto
  //Float_t fakeRate = 0.08;
  //Float_t fakeRate_err = 0.02;

  Float_t fakeRate_sys = 0.0006;

  // Scale the backgrounds
  Float_t egScale = fakeRate/(1. - fakeRate);
  Float_t egScaleErr = fakeRate_err/(1. - fakeRate)/(1. - fakeRate);

  TTree * ggTree = (TTree*)in->Get("gg_"+channels[channel]+"_EvtTree");
  TTree * egTree = (TTree*)in->Get("eg_"+channels[channel]+"_EvtTree");

  TFile * fTTGJets = new TFile("inputs/signal_contamination_ttgjets.root", "READ");
  TTree * ttgjetsTree = (TTree*)fTTGJets->Get("gg_"+channels[channel]+"_EvtTree_ttgjets");

  TFile * fQCD30to40 = new TFile("inputs/signal_contamination_qcd30to40.root", "READ");
  TTree * qcd30to40Tree = (TTree*)fQCD30to40->Get("gg_"+channels[channel]+"_EvtTree_qcd30to40");

  TFile * fQCD40 = new TFile("inputs/signal_contamination_qcd40.root", "READ");
  TTree * qcd40Tree = (TTree*)fQCD40->Get("gg_"+channels[channel]+"_EvtTree_qcd40");

  TFile * fGJet20to40 = new TFile("inputs/signal_contamination_GJet20to40.root", "READ");
  TTree * gjet20to40Tree = (TTree*)fGJet20to40->Get("gg_"+channels[channel]+"_EvtTree_GJet20to40");

  TFile * fGJet40 = new TFile("inputs/signal_contamination_GJet40.root", "READ");
  TTree * gjet40Tree = (TTree*)fGJet40->Get("gg_"+channels[channel]+"_EvtTree_GJet40");

  TFile * fTTHadronic = new TFile("inputs/signal_contamination_ttJetsHadronic.root", "READ");
  TTree * ttHadronicTree = (TTree*)fTTHadronic->Get("gg_"+channels[channel]+"_EvtTree_ttJetsHadronic");
  
  TFile * fTTSemiLep = new TFile("inputs/signal_contamination_ttJetsSemiLep.root", "READ");
  TTree * ttSemiLepTree = (TTree*)fTTSemiLep->Get("gg_"+channels[channel]+"_EvtTree_ttJetsSemiLep");

  TFile * fSigA = new TFile("../acceptance/signal_contamination_mst_460_m1_175.root", "READ");
  TTree * sigaTree = (TTree*)fSigA->Get("gg_"+channels[channel]+"_EvtTree_mst_460_m1_175");

  TFile * fSigB = new TFile("../acceptance/signal_contamination_mst_560_m1_325.root", "READ");
  //TTree * sigbTree = (TTree*)fSigB->Get("gg_"+channels[channel]+"_EvtTree_mst_560_m1_325");
  TTree * sigbTree = (TTree*)fSigB->Get("gg_"+channels[channel]+"_EvtTree_mst_460_m1_175"); // DURP

  // Now save the met plots out to file -- use these later for the limit-setting
  TFile * out = new TFile("met_reweighted_"+channels[channel]+".root", "RECREATE");

  // durp

  out->Write();
  out->Close();

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);

  // Make the correlation plot for MET filters
  TH2D * metFilter = (TH2D*)in->Get("metFilter");
  if(channel == 0) {
    metFilter->GetXaxis()->SetBinLabel(1, "CSCBeamHalo");
    metFilter->GetYaxis()->SetBinLabel(1, "CSCBeamHalo");
    metFilter->GetXaxis()->SetBinLabel(2, "HcalNoise");
    metFilter->GetYaxis()->SetBinLabel(2, "HcalNoise");
    metFilter->GetXaxis()->SetBinLabel(3, "EcalDeadCellTP");
    metFilter->GetYaxis()->SetBinLabel(3, "EcalDeadCellTP");
    metFilter->GetXaxis()->SetBinLabel(4, "EcalDeadCellBE");
    metFilter->GetYaxis()->SetBinLabel(4, "EcalDeadCellBE");
    metFilter->GetXaxis()->SetBinLabel(5, "HcalLaser");
    metFilter->GetYaxis()->SetBinLabel(5, "HcalLaser");
    metFilter->GetXaxis()->SetBinLabel(6, "TrackingFailure");
    metFilter->GetYaxis()->SetBinLabel(6, "TrackingFailure");
    metFilter->GetXaxis()->SetBinLabel(7, "EEBadSC");
    metFilter->GetYaxis()->SetBinLabel(7, "EEBadSC");
    metFilter->GetXaxis()->SetBinLabel(8, "EERingOfFire");
    metFilter->GetYaxis()->SetBinLabel(8, "EERingOfFire");
    metFilter->GetXaxis()->SetBinLabel(9, "InconsistentMuon");
    metFilter->GetYaxis()->SetBinLabel(9, "InconsistentMuon");
    metFilter->GetXaxis()->SetBinLabel(10, "GreedyMuon");
    metFilter->GetYaxis()->SetBinLabel(10, "GreedyMuon");

    metFilter->GetXaxis()->SetLabelSize(0.035);
    metFilter->GetYaxis()->SetLabelSize(0.015);
    metFilter->GetZaxis()->SetLabelSize(0.02);

    metFilter->Draw("colz");
    metFilter->SetMarkerColor(kWhite);
    metFilter->Draw("text same");
    can->SetLogz(true);
    can->SaveAs("metFilter"+gifOrPdf);

    can->SetLogz(false);
  }

  TH2D * DR_jet_gg = (TH2D*)in->Get("DR_jet_gg");
  if(channel == 0) {
    DR_jet_gg->Draw("colz");
    TLine * vertline = new TLine(0.5, 0, 0.5, 5);
    vertline->SetLineColor(kRed);
    vertline->SetLineWidth(3);
    vertline->Draw("same");
    TLine * horiline = new TLine(0, 0.5, 5, 0.5);
    horiline->SetLineColor(kRed);
    horiline->SetLineWidth(3);
    horiline->Draw("same");
    can->SetLogz(true);
    can->SaveAs("DR_jet_gg_"+channels[channel]+gifOrPdf);

    can->SetLogz(false);
  }

  PlotMaker * pMaker = new PlotMaker(intLumi_int, egScale, egScaleErr, channels[channel]);
  pMaker->SetTrees(ggTree, egTree,
		   qcd30to40Tree, qcd40Tree,
		   gjet20to40Tree, gjet40Tree,
		   ttHadronicTree, ttSemiLepTree,
		   ttgjetsTree,
		   sigaTree, sigbTree);

  const int nMetBins = 16;
  Double_t xbins_met[nMetBins+1] = {
    0,
    5,
    10,
    15,
    20,
    25,
    30,
    35,
    40,
    45,
    50,
    60,
    70,
    80,
    100,
    150,
    300};
  //650};

  pMaker->CreateMETPlot(nMetBins, xbins_met,
			xbins_met[0], xbins_met[nMetBins],
			7.e-4, 25000.,
			0., 9.1,
			true, true, true);
    
  in->Close();
  fTTGJets->Close();
  fQCD30to40->Close();
  fQCD40->Close();
  fGJet20to40->Close();
  fGJet40->Close();
  fTTHadronic->Close();
  fTTSemiLep->Close();
  fSigA->Close();
  fSigB->Close();

}
