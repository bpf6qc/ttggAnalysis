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

#include "analyze.h"
#include "rootRoutines.h"

using namespace std;

const TString ffColor = "kOrange+10";
const TString eeColor = "kBlue";
const TString egColor = "kGreen";

void analyze(TString input, bool addMC, int channel, TString intLumi, int intLumi_int, bool useFF, bool useDifferenceSystematic, bool useTTGJets, bool useMCforQCD) {

  gROOT->Reset();
  gROOT->SetBatch(true);
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0000);
  gStyle->SetOptTitle(0);
  //gSystem->Load("libRooFit");

  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

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

  // Start grabbing objects from input file
  TTree * gfTree = (TTree*)in->Get("gf_"+channels[channel]+"_EvtTree");
  TTree * ffTree = (TTree*)in->Get("ff_"+channels[channel]+"_EvtTree");
  TTree * eeTree = (TTree*)in->Get("ee_"+channels[channel]+"_EvtTree");

  TString emOrJet = "jet";

  TH2D * diempt_gg =        (TH2D*)in->Get("di"+emOrJet+"pt_gg_"+channels[channel]); diempt_gg->Sumw2();
  TH2D * diempt_gf =        (TH2D*)in->Get("di"+emOrJet+"pt_gf_"+channels[channel]); diempt_gf->Sumw2();
  TH2D * diempt_ff =        (TH2D*)in->Get("di"+emOrJet+"pt_ff_"+channels[channel]); diempt_ff->Sumw2();
  TH2D * diempt_ee =        (TH2D*)in->Get("di"+emOrJet+"pt_ee_"+channels[channel]); diempt_ee->Sumw2();
  TH2D * diempt_ee_lo =     (TH2D*)in->Get("di"+emOrJet+"pt_ee_loMass_"+channels[channel]); diempt_ee_lo->Sumw2();
  TH2D * diempt_ee_hi =     (TH2D*)in->Get("di"+emOrJet+"pt_ee_hiMass_"+channels[channel]); diempt_ee_hi->Sumw2();
  TH2D * diempt_ee_signal = (TH2D*)in->Get("di"+emOrJet+"pt_ee_onMass_"+channels[channel]); diempt_ee_signal->Sumw2();
  
  // Make ratios
  TH1D * diempt_gg_0 = (TH1D*)diempt_gg->ProjectionX("gg_px0", 1, 1, "e");
  TH1D * diempt_gg_1 = (TH1D*)diempt_gg->ProjectionX("gg_px1", 2, 2, "e");
  TH1D * diempt_gg_2 = (TH1D*)diempt_gg->ProjectionX("gg_px2", 3, -1, "e");

  Float_t gg_test = diempt_gg->Integral();
  Float_t gf_test = diempt_gf->Integral();
  Float_t ff_test = diempt_ff->Integral();
  Float_t ee_test = diempt_ee->Integral();
  Float_t ee_lo_test = diempt_ee_lo->Integral();
  Float_t ee_hi_test = diempt_ee_hi->Integral();
  Float_t ee_signal_test = diempt_ee_signal->Integral();

  TH1D * ratio_gf_0 = GetWeights(diempt_gg_0, diempt_gf->ProjectionX("gf_px0", 1, 1, "e"), gg_test, gf_test);
  TH1D * ratio_gf_1 = GetWeights(diempt_gg_1, diempt_gf->ProjectionX("gf_px1", 2, 2, "e"), gg_test, gf_test);
  TH1D * ratio_gf_2 = GetWeights(diempt_gg_2, diempt_gf->ProjectionX("gf_px2", 3, -1, "e"), gg_test, gf_test);

  TH1D * ratio_ff_0 = GetWeights(diempt_gg_0, diempt_ff->ProjectionX("ff_px0", 1, 1, "e"), gg_test, ff_test);
  TH1D * ratio_ff_1 = GetWeights(diempt_gg_1, diempt_ff->ProjectionX("ff_px1", 2, 2, "e"), gg_test, ff_test);
  TH1D * ratio_ff_2 = GetWeights(diempt_gg_2, diempt_ff->ProjectionX("ff_px2", 3, -1, "e"), gg_test, ff_test);

  TH1D * ratio_ee_0 = GetWeights(diempt_gg_0, diempt_ee->ProjectionX("ee_px0", 1, 1, "e"), gg_test, ee_test);
  TH1D * ratio_ee_1 = GetWeights(diempt_gg_1, diempt_ee->ProjectionX("ee_px1", 2, 2, "e"), gg_test, ee_test);
  TH1D * ratio_ee_2 = GetWeights(diempt_gg_2, diempt_ee->ProjectionX("ee_px2", 3, -1, "e"), gg_test, ee_test);

  TH1D * ratio_ee_lo_0 = GetWeights(diempt_gg_0, diempt_ee_lo->ProjectionX("ee_lo_px0", 1, 1, "e"), gg_test, ee_lo_test);
  TH1D * ratio_ee_lo_1 = GetWeights(diempt_gg_1, diempt_ee_lo->ProjectionX("ee_lo_px1", 2, 2, "e"), gg_test, ee_lo_test);
  TH1D * ratio_ee_lo_2 = GetWeights(diempt_gg_2, diempt_ee_lo->ProjectionX("ee_lo_px2", 3, -1, "e"), gg_test, ee_lo_test);

  TH1D * ratio_ee_hi_0 = GetWeights(diempt_gg_0, diempt_ee_hi->ProjectionX("ee_hi_px0", 1, 1, "e"), gg_test, ee_hi_test);
  TH1D * ratio_ee_hi_1 = GetWeights(diempt_gg_1, diempt_ee_hi->ProjectionX("ee_hi_px1", 2, 2, "e"), gg_test, ee_hi_test);
  TH1D * ratio_ee_hi_2 = GetWeights(diempt_gg_2, diempt_ee_hi->ProjectionX("ee_hi_px2", 3, -1, "e"), gg_test, ee_hi_test);

  TH1D * ratio_ee_signal_0 = GetWeights(diempt_gg_0, diempt_ee_signal->ProjectionX("ee_signal_px0", 1, 1, "e"), gg_test, ee_signal_test);
  TH1D * ratio_ee_signal_1 = GetWeights(diempt_gg_1, diempt_ee_signal->ProjectionX("ee_signal_px1", 2, 2, "e"), gg_test, ee_signal_test);
  TH1D * ratio_ee_signal_2 = GetWeights(diempt_gg_2, diempt_ee_signal->ProjectionX("ee_signal_px2", 3, -1, "e"), gg_test, ee_signal_test);

  // Reweight MET spectra
  TH1D * met_default = (TH1D*)in->Get("met_ff_"+channels[channel]); met_default->Sumw2();

  TH1D * met_ff_Rew = Reweight(true, met_default, "pfMET", ratio_ff_0, ratio_ff_1, ratio_ff_2, ffTree, 0, true, "ff_"+channels[channel], false);
  TH1D * met_gf_Rew = Reweight(true, met_default, "pfMET", ratio_gf_0, ratio_gf_1, ratio_gf_2, gfTree, 0, true, "gf_"+channels[channel], false);
  TH1D * met_ee_Rew = Reweight(true, met_default, "pfMET", ratio_ee_signal_0, ratio_ee_signal_1, ratio_ee_signal_2, eeTree, 1, true, "ee_"+channels[channel], false);
  TH1D * met_ee_lo_Rew = Reweight(true, met_default, "pfMET", ratio_ee_lo_0, ratio_ee_lo_1, ratio_ee_lo_2, eeTree, 2, true, "ee_"+channels[channel]+"_lo", false);
  TH1D * met_ee_hi_Rew = Reweight(true, met_default, "pfMET", ratio_ee_hi_0, ratio_ee_hi_1, ratio_ee_hi_2, eeTree, 3, true, "ee_"+channels[channel]+"_hi", false);

  // pixel veto
  Float_t fakeRate = 0.0199;
  Float_t fakeRate_err = 0.0001;

  // electron conversion veto
  //Float_t fakeRate = 0.08;
  //Float_t fakeRate_err = 0.02;

  Float_t fakeRate_sys = 0.0006;

  // Do the sideband subtraction on ee
  for(int i = 0; i < met_ee_Rew->GetNbinsX(); i++) {
    Double_t e1 = met_ee_Rew->GetBinError(i+1);
    Double_t e2 = met_ee_lo_Rew->GetBinError(i+1);
    Double_t e3 = met_ee_hi_Rew->GetBinError(i+1);

    met_ee_Rew->SetBinContent(i+1, met_ee_Rew->GetBinContent(i+1) - met_ee_lo_Rew->GetBinContent(i+1) - met_ee_hi_Rew->GetBinContent(i+1));
    met_ee_Rew->SetBinError(i+1, sqrt(e1*e1 + e2*e2 + e3*e3));
  }

  // Get eg and gg MET
  TH1D * met_gg = (TH1D*)in->Get("met_gg_"+channels[channel]); met_gg->Sumw2();
  TH1D * met_eg = (TH1D*)in->Get("met_eg_"+channels[channel]); met_eg->Sumw2();

  met_gg->SetName("met_gg_"+channels[channel]);
  met_eg->SetName("met_eg_"+channels[channel]);

  // Scale the backgrounds
  Float_t egScale = fakeRate/(1. - fakeRate);
  Float_t egScaleErr = fakeRate_err/(1. - fakeRate)/(1. - fakeRate);
  TH1D * met_eg_noNorm = (TH1D*)met_eg->Clone();
  met_eg->Scale(egScale);

  TH1D * ewk_normErr2 = (TH1D*)met_eg->Clone();

  for(int i = 0; i < met_eg->GetNbinsX(); i++) {
    Float_t normerr = egScaleErr*(met_eg_noNorm->GetBinContent(i+1));
    ewk_normErr2->SetBinContent(i+1, normerr*normerr);

    Float_t staterr = met_eg->GetBinError(i+1);

    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);

    met_eg->SetBinError(i+1, new_err);
  }

  TTree * ggTree = (TTree*)in->Get("gg_"+channels[channel]+"_EvtTree");
  TTree * egTree = (TTree*)in->Get("eg_"+channels[channel]+"_EvtTree");

  TFile * fTTGJets = new TFile("/uscms_data/d2/bfrancis/btagRA3/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/signal_contamination_ttgjets.root", "READ");
  TTree * ttgjetsTree = (TTree*)fTTGJets->Get("gg_"+channels[channel]+"_EvtTree_ttgjets");
  Float_t ttgjetsScale = intLumi_int * 1.444 / 71598.;

  TFile * fQCD30to40 = new TFile("/uscms_data/d2/bfrancis/btagRA3/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/signal_contamination_qcd30to40.root", "READ");
  TTree * qcd30to40Tree = (TTree*)fQCD30to40->Get("gg_"+channels[channel]+"_EvtTree_qcd30to40");
  Float_t qcd30to40Scale = intLumi_int * 5.195E7 * 2.35E-4 * 1.019 * 1.019 / 6061407.;

  TFile * fQCD40 = new TFile("/uscms_data/d2/bfrancis/btagRA3/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/signal_contamination_qcd40.root", "READ");
  TTree * qcd40Tree = (TTree*)fQCD40->Get("gg_"+channels[channel]+"_EvtTree_qcd40");
  Float_t qcd40Scale = intLumi_int * 5.195E7 * 0.002175 * 1.019 * 1.019 / 9782735.;

  TH1D * met_ttgjets = SignalHistoFromTree(ttgjetsScale, true, "pfMET", ttgjetsTree, "met_ttgjets", "met_ttgjets", 400, 0., 2000.);
  met_ttgjets->Scale(ttgjetsScale);

  int maxBinToNorm = met_ee_Rew->FindBin(20.0);
  if(met_ee_Rew->FindBin(20.001) == maxBinToNorm) maxBinToNorm--;

  Float_t eeScale, eeScaleErr_num, eeScaleErr_den, eeScaleErr;

  if(useTTGJets) eeScale = (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm) - met_ttgjets->Integral(1, maxBinToNorm)) / met_ee_Rew->Integral(1, maxBinToNorm);
  else eeScale = (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm)) / met_ee_Rew->Integral(1, maxBinToNorm);
  eeScaleErr_num = 0.;
  for(int i = 0; i < maxBinToNorm; i++) {
    eeScaleErr_num += met_eg->GetBinError(i+1) * met_eg->GetBinError(i+1);
  }
  if(useTTGJets) {
    Float_t eeScaleErr_num2 = 0.;
    for(int i = 0; i < maxBinToNorm; i++) {
      eeScaleErr_num2 += met_ttgjets->GetBinError(i+1) * met_ttgjets->GetBinError(i+1);
    }
    eeScaleErr_num = (sqrt(met_gg->Integral(1, maxBinToNorm)) - sqrt(eeScaleErr_num) - sqrt(eeScaleErr_num2)) / (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm) - met_ttgjets->Integral(1, maxBinToNorm));
  }
  else eeScaleErr_num = (sqrt(met_gg->Integral(1, maxBinToNorm)) - sqrt(eeScaleErr_num)) / (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm));

  eeScaleErr_den = 0.;
  for(int i = 0; i < maxBinToNorm; i++) {
    eeScaleErr_den += met_ee_Rew->GetBinError(i+1) * met_ee_Rew->GetBinError(i+1);
  }
  eeScaleErr_den = sqrt(eeScaleErr_den) / (met_ee_Rew->Integral(1, maxBinToNorm));
  
  eeScaleErr = sqrt(eeScaleErr_num*eeScaleErr_num + eeScaleErr_den*eeScaleErr_den) * eeScale;
  
  bool ee_failure = false;
  if(met_ee_Rew->Integral(1, maxBinToNorm) == 0.) {
    ee_failure = true;
    eeScale = 1.;
    eeScaleErr = 0.;
  }
  
  TH1D * met_ee_Rew_noNorm = (TH1D*)met_ee_Rew->Clone();
  met_ee_Rew->Scale(eeScale);

  TH1D * ee_norm2 = (TH1D*)met_ee_Rew->Clone();

  for(int i = 0; i < met_ee_Rew->GetNbinsX(); i++) {
    Float_t normerr = eeScaleErr*(met_ee_Rew_noNorm->GetBinContent(i+1));
    ee_norm2->SetBinContent(i+1, normerr*normerr);

    Float_t staterr = met_ee_Rew->GetBinError(i+1);

    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);
    
    met_ee_Rew->SetBinError(i+1, new_err);
  }

  maxBinToNorm = met_ff_Rew->FindBin(20.0);
  if(met_ff_Rew->FindBin(20.001) == maxBinToNorm) maxBinToNorm--;

  Float_t gfScale, gfScaleErr_num, gfScaleErr_den, gfScaleErr;

  if(useTTGJets) gfScale = (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm) - met_ttgjets->Integral(1, maxBinToNorm))/met_gf_Rew->Integral(1, maxBinToNorm);
  else gfScale = (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm))/met_gf_Rew->Integral(1, maxBinToNorm);
  gfScaleErr_num = 0.;
  for(int i = 0; i < maxBinToNorm; i++)	{
    gfScaleErr_num += met_eg->GetBinError(i+1) * met_eg->GetBinError(i+1);
  }
  if(useTTGJets) {
    Float_t gfScaleErr_num2 = 0.;
    for(int i = 0; i < maxBinToNorm; i++) {
      gfScaleErr_num2 += met_ttgjets->GetBinError(i+1) * met_ttgjets->GetBinError(i+1);
    }
    gfScaleErr_num = (sqrt(met_gg->Integral(1, maxBinToNorm)) - sqrt(gfScaleErr_num) - sqrt(gfScaleErr_num2)) / (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm) - met_ttgjets->Integral(1, maxBinToNorm));
  }
  else gfScaleErr_num = (sqrt(met_gg->Integral(1, maxBinToNorm)) - sqrt(gfScaleErr_num)) / (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm));

  gfScaleErr_den = 0.;
  for(int i = 0; i < maxBinToNorm; i++) {
    gfScaleErr_den += met_gf_Rew->GetBinError(i+1) * met_gf_Rew->GetBinError(i+1);
  }
  gfScaleErr_den = sqrt(gfScaleErr_den) / met_gf_Rew->Integral(1, maxBinToNorm);
  
  gfScaleErr = sqrt(gfScaleErr_num*gfScaleErr_num + gfScaleErr_den*gfScaleErr_den) * gfScale;

  bool gf_failure = false;
  if(met_gf_Rew->Integral(1, maxBinToNorm) == 0.) {
    gf_failure = true;
    gfScale = 1.;
    gfScaleErr = 0.;
  }

  TH1D * met_gf_Rew_noNorm = (TH1D*)met_gf_Rew->Clone();
  met_gf_Rew->Scale(gfScale);

  TH1D * gf_norm2 = (TH1D*)met_gf_Rew->Clone();

  for(int i = 0; i < met_gf_Rew->GetNbinsX(); i++) {
    Float_t normerr = gfScaleErr*(met_gf_Rew_noNorm->GetBinContent(i+1));
    gf_norm2->SetBinContent(i+1, normerr*normerr);

    Float_t staterr = met_gf_Rew->GetBinError(i+1);
    
    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);

    met_gf_Rew->SetBinError(i+1, new_err);
  }

  Float_t ffScale, ffScaleErr_num, ffScaleErr_den, ffScaleErr;

  if(useTTGJets) ffScale = (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm) - met_ttgjets->Integral(1, maxBinToNorm))/met_ff_Rew->Integral(1, maxBinToNorm);
  else ffScale = (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm))/met_ff_Rew->Integral(1, maxBinToNorm);
  ffScaleErr_num = 0.;
  for(int i = 0; i < maxBinToNorm; i++)	{
    ffScaleErr_num += met_eg->GetBinError(i+1) * met_eg->GetBinError(i+1);
  }
  if(useTTGJets) {
    Float_t ffScaleErr_num2 = 0.;
    for(int i = 0; i < maxBinToNorm; i++) {
      ffScaleErr_num2 += met_ttgjets->GetBinError(i+1) * met_ttgjets->GetBinError(i+1);
    }
    ffScaleErr_num = (sqrt(met_gg->Integral(1, maxBinToNorm)) - sqrt(ffScaleErr_num) - sqrt(ffScaleErr_num2)) / (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm) - met_ttgjets->Integral(1, maxBinToNorm));
  }
  else ffScaleErr_num = (sqrt(met_gg->Integral(1, maxBinToNorm)) - sqrt(ffScaleErr_num)) / (met_gg->Integral(1, maxBinToNorm) - met_eg->Integral(1, maxBinToNorm));

  ffScaleErr_den = 0.;
  for(int i = 0; i < maxBinToNorm; i++) {
    ffScaleErr_den += met_ff_Rew->GetBinError(i+1) * met_ff_Rew->GetBinError(i+1);
  }
  ffScaleErr_den = sqrt(ffScaleErr_den) / met_ff_Rew->Integral(1, maxBinToNorm);
  
  ffScaleErr = sqrt(ffScaleErr_num*ffScaleErr_num + ffScaleErr_den*ffScaleErr_den) * ffScale;

  bool ff_failure = false;
  if(met_ff_Rew->Integral(1, maxBinToNorm) == 0.) {
    ff_failure = true;
    ffScale = 1.;
    ffScaleErr = 0.;
  }

  TH1D * met_ff_Rew_noNorm = (TH1D*)met_ff_Rew->Clone();
  met_ff_Rew->Scale(ffScale);

  TH1D * ff_norm2 = (TH1D*)met_ff_Rew->Clone();

  for(int i = 0; i < met_ff_Rew->GetNbinsX(); i++) {
    Float_t normerr = ffScaleErr*(met_ff_Rew_noNorm->GetBinContent(i+1));
    ff_norm2->SetBinContent(i+1, normerr*normerr);

    Float_t staterr = met_ff_Rew->GetBinError(i+1);
    
    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);

    met_ff_Rew->SetBinError(i+1, new_err);
  }

  formatTable(met_gg,
	      met_eg_noNorm, ewk_normErr2, fakeRate, fakeRate_sys, egScale,
	      met_ff_Rew_noNorm, ff_norm2, ffScale,
	      met_gf_Rew_noNorm, gf_norm2, gfScale,
	      met_ee_Rew_noNorm, ee_norm2, eeScale,
	      met_ttgjets, useTTGJets,
	      channels[channel]);

  // Now save the met plots out to file -- use these later for the limit-setting
  TFile * out = new TFile("met_reweighted_"+channels[channel]+".root", "RECREATE");

  met_gg->Write();
  met_eg->Write();
  met_gf_Rew->Write();
  met_ff_Rew->Write();
  met_ee_Rew->Write();
  met_ttgjets->Write();

ratio_gf_0->Write();
ratio_gf_1->Write();
ratio_gf_2->Write();

ratio_ff_0->Write();
ratio_ff_1->Write();
ratio_ff_2->Write();

  out->Write();
  out->Close();

  // Make plots
  const int nBins = 16;
  Double_t xbins[nBins+1] = {
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

  TH1D * h_gg = (TH1D*)met_gg->Rebin(nBins, "h_gg", xbins);
  TH1D * h_ewk = (TH1D*)met_eg->Rebin(nBins, "h_ewk", xbins);
  TH1D * h_qcd_ee = (TH1D*)met_ee_Rew->Rebin(nBins, "h_qcd_ee", xbins);
  TH1D * h_qcd_ff = (TH1D*)met_ff_Rew->Rebin(nBins, "h_qcd_ff", xbins);
  TH1D * h_qcd_gf = (TH1D*)met_gf_Rew->Rebin(nBins, "h_qcd_gf", xbins);

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);

  // Make the correlation plot for MET filters
  TH2D * metFilter = (TH2D*)in->Get("metFilter");
  if(channel == 0) {
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

  //if(channel == 0) analyzeReweighting(ggTree, ffTree, gfTree);

  TFile * fSig460 = new TFile("/uscms_data/d2/bfrancis/btagRA3/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/acceptance/signal_contamination_mst_460_m1_175.root", "READ");
  TTree * ggTree_460 = (TTree*)fSig460->Get("gg_"+channels[channel]+"_EvtTree_mst_460_m1_175");

  TFile * fSig800 = new TFile("/uscms_data/d2/bfrancis/btagRA3/CMSSW_5_3_8_patch3/src/SusyAnalysis/SusyNtuplizer/macro/acceptance/signal_contamination_mst_560_m1_325.root", "READ");
  //TTree * ggTree_800 = (TTree*)fSig800->Get("gg_"+channels[channel]+"_EvtTree_mst_560_m1_325");
  TTree * ggTree_800 = (TTree*)fSig800->Get("gg_"+channels[channel]+"_EvtTree_mst_460_m1_175"); // DURP

  Float_t signalScale_460 = 0.147492 * intLumi_int * 1.019 * 1.019 / 15000.;
  Float_t signalScale_800 = 0.0399591 * intLumi_int * 1.019 * 1.019 / 15000.;

  vector<TH1D*> leadMatchedJetPt = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
						ratio_ff_0, ratio_ff_1, ratio_ff_2,
						ratio_gf_0, ratio_gf_1, ratio_gf_2,
						egScale, egScaleErr,
						ffScale, ffScaleErr,
						gfScale, gfScaleErr, ttgjetsScale,
						"leadMatchedJetPt", true,
						200, 0., 2000.,
						channels[channel],
						1,
						true);

  TH1D * leadMatchedJetPt_460 = SignalHistoFromTree(signalScale_460, true, "leadMatchedJetPt", ggTree_460, "leadMatchedJetPt_460", "leadMatchedJetPt_460", 200, 0., 2000.);
  TH1D * leadMatchedJetPt_800 = SignalHistoFromTree(signalScale_800, true, "leadMatchedJetPt", ggTree_800, "leadMatchedJetPt_800", "leadMatchedJetPt_800", 200, 0., 2000.);
  TH1D * leadMatchedJetPt_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "leadMatchedJetPt", qcd30to40Tree, "leadMatchedJetPt_qcd30to40", "leadMatchedJetPt_qcd30to40", 200, 0., 2000.);
  TH1D * leadMatchedJetPt_qcd40 = SignalHistoFromTree(qcd40Scale, true, "leadMatchedJetPt", qcd40Tree, "leadMatchedJetPt_qcd40", "leadMatchedJetPt_qcd40", 200, 0., 2000.);

  prettyPlot(leadMatchedJetPt[0], leadMatchedJetPt[1], leadMatchedJetPt[2], leadMatchedJetPt[3], leadMatchedJetPt[4],
	     useDifferenceSystematic, useTTGJets, useFF, false,
	     leadMatchedJetPt_460, leadMatchedJetPt_800,
	     leadMatchedJetPt_qcd30to40, leadMatchedJetPt_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "leadMatchedJetPt", "lead #gamma's matched jet Pt ", "Number of Events",
	     0., 1400., 1.e-2, 3.e6,
	     0., 1.9,
	     false, true, true);

  vector<TH1D*> trailMatchedJetPt = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
						 ratio_ff_0, ratio_ff_1, ratio_ff_2,
						 ratio_gf_0, ratio_gf_1, ratio_gf_2,
						 egScale, egScaleErr,
						 ffScale, ffScaleErr,
						 gfScale, gfScaleErr, ttgjetsScale,
						 "trailMatchedJetPt", true,
						 200, 0., 2000.,
						 channels[channel],
						 1,
						 true);
  TH1D * trailMatchedJetPt_460 = SignalHistoFromTree(signalScale_460, true, "trailMatchedJetpt", ggTree_460, "trailMatchedJetPt_460", "trailMatchedJetPt_460", 200, 0., 2000.);
  TH1D * trailMatchedJetPt_800 = SignalHistoFromTree(signalScale_800, true, "trailMatchedJetPt", ggTree_800, "trailMatchedJetPt_800", "trailMatchedJetPt_800", 200, 0., 2000.);
  TH1D * trailMatchedJetPt_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "trailMatchedJetPt", qcd30to40Tree, "trailMatchedJetPt_qcd30to40", "trailMatchedJetPt_qcd30to40", 200, 0., 2000.);
  TH1D * trailMatchedJetPt_qcd40 = SignalHistoFromTree(qcd40Scale, true, "trailMatchedJetPt", qcd40Tree, "trailMatchedJetPt_qcd40", "trailMatchedJetPt_qcd40", 200, 0., 2000.);

  prettyPlot(trailMatchedJetPt[0], trailMatchedJetPt[1], trailMatchedJetPt[2], trailMatchedJetPt[3], trailMatchedJetPt[4],
	     useDifferenceSystematic, useTTGJets, useFF, false,
	     trailMatchedJetPt_460, trailMatchedJetPt_800,
	     trailMatchedJetPt_qcd30to40, trailMatchedJetPt_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "trailMatchedJetPt", "trail #gamma's matched jet Pt ", "Number of Events",
	     0., 1400., 1.e-2, 3.e6,
	     0., 1.9,
	     false, true, true);

  vector<TH1D*> leadPtOverInvmass = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
						 ratio_ff_0, ratio_ff_1, ratio_ff_2,
						 ratio_gf_0, ratio_gf_1, ratio_gf_2,
						 egScale, egScaleErr,
						 ffScale, ffScaleErr,
						 gfScale, gfScaleErr, ttgjetsScale,
						 "leadptOverInvmass", true,
						 300, 0., 6.,
						 channels[channel],
						 3,
						 true);

  TH1D * leadPtOverInvmass_460 = SignalHistoFromTree(signalScale_460, true, "leadptOverInvmass", ggTree_460, "leadPtOverInvmass_460", "leadPtOverInvmass_460", 100, 0., 6.);
  TH1D * leadPtOverInvmass_800 = SignalHistoFromTree(signalScale_800, true, "leadptOverInvmass", ggTree_800, "leadPtOverInvmass_800", "leadPtOverInvmass_800", 100, 0., 6.);
  TH1D * leadPtOverInvmass_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "leadPtOverInvmass", qcd30to40Tree, "leadPtOverInvmass_qcd30to40", "leadPtOverInvmass_qcd30to40", 100, 0., 6.);
  TH1D * leadPtOverInvmass_qcd40 = SignalHistoFromTree(qcd40Scale, true, "leadPtOverInvmass", qcd40Tree, "leadPtOverInvmass_qcd40", "leadPtOverInvmass_qcd40", 100, 0., 6.);

  prettyPlot(leadPtOverInvmass[0], leadPtOverInvmass[1], leadPtOverInvmass[2], leadPtOverInvmass[3], leadPtOverInvmass[4],
	     useDifferenceSystematic, useTTGJets, useFF, false,
	     leadPtOverInvmass_460, leadPtOverInvmass_800,
	     leadPtOverInvmass_qcd30to40, leadPtOverInvmass_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "leadptOverInvmass", "lead Photon Pt / m(gg)", "Number of Events",
	     0., 6., 1.e-2, 3.e6,
	     0., 1.9,
	     true, true, true);

  vector<TH1D*> trailPtOverInvmass = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
						  ratio_ff_0, ratio_ff_1, ratio_ff_2,
						  ratio_gf_0, ratio_gf_1, ratio_gf_2,
						  egScale, egScaleErr,
						  ffScale, ffScaleErr,
						  gfScale, gfScaleErr, ttgjetsScale,
						  "trailptOverInvmass", true,
						  300, 0., 2.,
						  channels[channel],
						  3,
						  true);

  TH1D * trailPtOverInvmass_460 = SignalHistoFromTree(signalScale_460, true, "trailptOverInvmass", ggTree_460, "trailPtOverInvmass_460", "trailPtOverInvmass_460", 100, 0., 2.);
  TH1D * trailPtOverInvmass_800 = SignalHistoFromTree(signalScale_800, true, "trailptOverInvmass", ggTree_800, "trailPtOverInvmass_800", "trailPtOverInvmass_800", 100, 0., 2.);
  TH1D * trailPtOverInvmass_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "trailPtOverInvmass", qcd30to40Tree, "trailPtOverInvmass_qcd30to40", "trailPtOverInvmass_qcd30to40", 100, 0., 2.);
  TH1D * trailPtOverInvmass_qcd40 = SignalHistoFromTree(qcd40Scale, true, "trailPtOverInvmass", qcd40Tree, "trailPtOverInvmass_qcd40", "trailPtOverInvmass_qcd40", 100, 0., 2.);

  prettyPlot(trailPtOverInvmass[0], trailPtOverInvmass[1], trailPtOverInvmass[2], trailPtOverInvmass[3], trailPtOverInvmass[4],
	     useDifferenceSystematic, useTTGJets, useFF, false,
	     trailPtOverInvmass_460, trailPtOverInvmass_800,
	     trailPtOverInvmass_qcd30to40, trailPtOverInvmass_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "trailptOverInvmass", "trail Photon Pt / m(gg)", "Number of Events",
	     0., 2., 1.e-2, 3.e6,
	     0., 1.9,
	     true, true, true);

  vector<TH1D*> minDPhi_gMET = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
					    ratio_ff_0, ratio_ff_1, ratio_ff_2,
					    ratio_gf_0, ratio_gf_1, ratio_gf_2,
					    egScale, egScaleErr,
					    ffScale, ffScaleErr,
					    gfScale, gfScaleErr, ttgjetsScale,
					    "minDPhi_gMET", true,
					    170, 0., 3.2,
					    channels[channel],
					    2,
					    false);

  TH1D * minDPhi_gMET_460 = SignalHistoFromTree(signalScale_460, true, "minDPhi_gMET", ggTree_460, "minDPhi_gMET_460", "minDPhi_gMET_460", 85, 0., 3.2);
  TH1D * minDPhi_gMET_800 = SignalHistoFromTree(signalScale_800, true, "minDPhi_gMET", ggTree_800, "minDPhi_gMET_800", "minDPhi_gMET_800", 85, 0., 3.2);
  TH1D * minDPhi_gMET_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "minDPhi_gMET", qcd30to40Tree, "minDPhi_gMET_qcd30to40", "minDPhi_gMET_qcd30to40", 85, 0., 3.2);
  TH1D * minDPhi_gMET_qcd40 = SignalHistoFromTree(qcd40Scale, true, "minDPhi_gMET", qcd40Tree, "minDPhi_gMET_qcd40", "minDPhi_gMET_qcd40", 85, 0., 3.2);

  prettyPlot(minDPhi_gMET[0], minDPhi_gMET[1], minDPhi_gMET[2], minDPhi_gMET[3], minDPhi_gMET[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     minDPhi_gMET_460, minDPhi_gMET_800,
	     minDPhi_gMET_qcd30to40, minDPhi_gMET_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "minDPhi_gMET", "min dPhi between photons and MET", "Number of Events",
	     0., 3.2, 1.e-2, 3.e6,
	     0., 4.5,
	     true, false, false);

  vector<TH1D*> minDPhi_jMET = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
					    ratio_ff_0, ratio_ff_1, ratio_ff_2,
					    ratio_gf_0, ratio_gf_1, ratio_gf_2,
					    egScale, egScaleErr,
					    ffScale, ffScaleErr,
					    gfScale, gfScaleErr, ttgjetsScale,
					    "minDPhi_jMET", true,
					    170, 0., 3.2,
					    channels[channel],
					    2,
					    false);

  TH1D * minDPhi_jMET_460 = SignalHistoFromTree(signalScale_460, true, "minDPhi_jMET", ggTree_460, "minDPhi_jMET_460", "minDPhi_jMET_460", 85, 0., 3.2);
  TH1D * minDPhi_jMET_800 = SignalHistoFromTree(signalScale_800, true, "minDPhi_jMET", ggTree_800, "minDPhi_jMET_800", "minDPhi_jMET_800", 85, 0., 3.2);
  TH1D * minDPhi_jMET_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "minDPhi_jMET", qcd30to40Tree, "minDPhi_jMET_qcd30to40", "minDPhi_jMET_qcd30to40", 85, 0., 3.2);
  TH1D * minDPhi_jMET_qcd40 = SignalHistoFromTree(qcd40Scale, true, "minDPhi_jMET", qcd40Tree, "minDPhi_jMET_qcd40", "minDPhi_jMET_qcd40", 85, 0., 3.2);

  prettyPlot(minDPhi_jMET[0], minDPhi_jMET[1], minDPhi_jMET[2], minDPhi_jMET[3], minDPhi_jMET[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     minDPhi_jMET_460, minDPhi_jMET_800,
	     minDPhi_jMET_qcd30to40, minDPhi_jMET_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "minDPhi_jMET", "min dPhi between jets and MET", "Number of Events",
	     0., 3.2, 1.e-2, 3.e6,
	     0., 4.5,
	     true, false, false);

  vector<TH1D*> minDR_leadPhoton_jets = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
						     ratio_ff_0, ratio_ff_1, ratio_ff_2,
						     ratio_gf_0, ratio_gf_1, ratio_gf_2,
						     egScale, egScaleErr,
						     ffScale, ffScaleErr,
						     gfScale, gfScaleErr, ttgjetsScale,
						     "minDR_leadPhoton_jets", true,
						     250, 0., 5.,
						     channels[channel],
						     2,
						     false);

  TH1D * minDR_leadPhoton_jets_460 = SignalHistoFromTree(signalScale_460, true, "minDR_leadPhoton_jets", ggTree_460, "minDR_leadPhoton_jets_460", "minDR_leadPhoton_jets_460", 125, 0., 5.);
  TH1D * minDR_leadPhoton_jets_800 = SignalHistoFromTree(signalScale_800, true, "minDR_leadPhoton_jets", ggTree_800, "minDR_leadPhoton_jets_800", "minDR_leadPhoton_jets_800", 125, 0., 5.);
  TH1D * minDR_leadPhoton_jets_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "minDR_leadPhoton_jets", qcd30to40Tree, "minDR_leadPhoton_jets_qcd30to40", "minDR_leadPhoton_jets_qcd30to40", 125, 0., 5.);
  TH1D * minDR_leadPhoton_jets_qcd40 = SignalHistoFromTree(qcd40Scale, true, "minDR_leadPhoton_jets", qcd40Tree, "minDR_leadPhoton_jets_qcd40", "minDR_leadPhoton_jets_qcd40", 125, 0., 5.);

  prettyPlot(minDR_leadPhoton_jets[0], minDR_leadPhoton_jets[1], minDR_leadPhoton_jets[2], minDR_leadPhoton_jets[3], minDR_leadPhoton_jets[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     minDR_leadPhoton_jets_460, minDR_leadPhoton_jets_800,
	     minDR_leadPhoton_jets_qcd30to40, minDR_leadPhoton_jets_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "minDR_leadPhoton_jets", "min dR between lead and jets", "Number of Events",
	     0, 5, 1.e-2, 3.e6,
	     0., 4.5,
	     true, false, false);

  vector<TH1D*> minDR_trailPhoton_jets = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
						      ratio_ff_0, ratio_ff_1, ratio_ff_2,
						      ratio_gf_0, ratio_gf_1, ratio_gf_2,
						      egScale, egScaleErr,
						      ffScale, ffScaleErr,
						      gfScale, gfScaleErr, ttgjetsScale,
						      "minDR_trailPhoton_jets", true,
						      250, 0., 5.,
						      channels[channel],
						      2,
						      false);

  TH1D * minDR_trailPhoton_jets_460 = SignalHistoFromTree(signalScale_460, true, "minDR_trailPhoton_jets", ggTree_460, "minDR_trailPhoton_jets_460", "minDR_trailPhoton_jets_460", 125, 0., 5.);
  TH1D * minDR_trailPhoton_jets_800 = SignalHistoFromTree(signalScale_800, true, "minDR_trailPhoton_jets", ggTree_800, "minDR_trailPhoton_jets_800", "minDR_trailPhoton_jets_800", 125, 0., 5.);
  TH1D * minDR_trailPhoton_jets_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "minDR_trailPhoton_jets", qcd30to40Tree, "minDR_trailPhoton_jets_qcd30to40", "minDR_trailPhoton_jets_qcd30to40", 125, 0., 5.);
  TH1D * minDR_trailPhoton_jets_qcd40 = SignalHistoFromTree(qcd40Scale, true, "minDR_trailPhoton_jets", qcd40Tree, "minDR_trailPhoton_jets_qcd40", "minDR_trailPhoton_jets_qcd40", 125, 0., 5.);

  prettyPlot(minDR_trailPhoton_jets[0], minDR_trailPhoton_jets[1], minDR_trailPhoton_jets[2], minDR_trailPhoton_jets[3], minDR_trailPhoton_jets[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     minDR_trailPhoton_jets_460, minDR_trailPhoton_jets_800,
	     minDR_trailPhoton_jets_qcd30to40, minDR_trailPhoton_jets_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "minDR_trailPhoton_jets", "min dR between trail and jets", "Number of Events",
	     0, 5, 1.e-2, 3.e6,
	     0., 4.5,
	     true, false, false);

  vector<TH1D*> jet1pt = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				      ratio_ff_0, ratio_ff_1, ratio_ff_2,
				      ratio_gf_0, ratio_gf_1, ratio_gf_2,
				      egScale, egScaleErr,
				      ffScale, ffScaleErr,
				      gfScale, gfScaleErr, ttgjetsScale,
				      "jet1_pt", true,
				      400, 0., 2000.,
				      channels[channel],
				      2,
				      false);

  TH1D * jet1pt_460 = SignalHistoFromTree(signalScale_460, true, "jet1_pt", ggTree_460, "jet1pt_460", "jet1pt_460", 200, 0., 2000.);
  TH1D * jet1pt_800 = SignalHistoFromTree(signalScale_800, true, "jet1_pt", ggTree_800, "jet1pt_800", "jet1pt_800", 200, 0., 2000.);
  TH1D * jet1pt_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "jet1_pt", qcd30to40Tree, "jet1pt_qcd30to40", "jet1pt_qcd30to40", 200, 0., 2000.);
  TH1D * jet1pt_qcd40 = SignalHistoFromTree(qcd40Scale, true, "jet1_pt", qcd40Tree, "jet1pt_qcd40", "jet1pt_qcd40", 200, 0., 2000.);

  prettyPlot(jet1pt[0], jet1pt[1], jet1pt[2], jet1pt[3], jet1pt[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     jet1pt_460, jet1pt_800,
	     jet1pt_qcd30to40, jet1pt_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "jet1pt", "Pt of leading jet", "Number of Events",
	     0, 1400, 1.e-2, 3.e6,
	     0., 4.5,
	     true, true, true);

  vector<TH1D*> jet2pt = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				      ratio_ff_0, ratio_ff_1, ratio_ff_2,
				      ratio_gf_0, ratio_gf_1, ratio_gf_2,
				      egScale, egScaleErr,
				      ffScale, ffScaleErr,
				      gfScale, gfScaleErr, ttgjetsScale,
				      "jet2_pt", true,
				      400, 0., 2000.,
				      channels[channel],
				      2,
				      false);

  TH1D * jet2pt_460 = SignalHistoFromTree(signalScale_460, true, "jet2_pt", ggTree_460, "jet2pt_460", "jet2pt_460", 200, 0., 2000.);
  TH1D * jet2pt_800 = SignalHistoFromTree(signalScale_800, true, "jet2_pt", ggTree_800, "jet2pt_800", "jet2pt_800", 200, 0., 2000.);
  TH1D * jet2pt_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "jet2_pt", qcd30to40Tree, "jet2pt_qcd30to40", "jet2pt_qcd30to40", 200, 0., 2000.);
  TH1D * jet2pt_qcd40 = SignalHistoFromTree(qcd40Scale, true, "jet2_pt", qcd40Tree, "jet2pt_qcd40", "jet2pt_qcd40", 200, 0., 2000.);

  prettyPlot(jet2pt[0], jet2pt[1], jet2pt[2], jet2pt[3], jet2pt[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     jet2pt_460, jet2pt_800,
	     jet2pt_qcd30to40, jet2pt_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "jet2pt", "Pt of sub-leading jet", "Number of Events",
	     0, 1400, 1.e-2, 3.e6,
	     0., 4.5,
	     true, true, true);

  vector<TH1D*> leadEta = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				       ratio_ff_0, ratio_ff_1, ratio_ff_2,
				       ratio_gf_0, ratio_gf_1, ratio_gf_2,
				       egScale, egScaleErr,
				       ffScale, ffScaleErr,
				       gfScale, gfScaleErr, ttgjetsScale,
				       "leadPhotonEta", true,
				       800, -1.5, 1.5,
				       channels[channel],
				       5,
				       false);

  TH1D * leadEta_460 = SignalHistoFromTree(signalScale_460, true, "leadPhotonEta", ggTree_460, "leadEta_460", "leadEta_460", 160, -1.5, 1.5);
  TH1D * leadEta_800 = SignalHistoFromTree(signalScale_800, true, "leadPhotonEta", ggTree_800, "leadEta_800", "leadEta_800", 160, -1.5, 1.5);
  TH1D * leadEta_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "leadPhotonEta", qcd30to40Tree, "leadEta_qcd30to40", "leadEta_qcd30to40", 160, -1.5, 1.5);
  TH1D * leadEta_qcd40 = SignalHistoFromTree(qcd40Scale, true, "leadPhotonEta", qcd40Tree, "leadEta_qcd40", "leadEta_qcd40", 160, -1.5, 1.5);

  prettyPlot(leadEta[0], leadEta[1], leadEta[2], leadEta[3], leadEta[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     leadEta_460, leadEta_800,
	     leadEta_qcd30to40, leadEta_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "leadEta", "#eta of leading #gamma", "Number of Events",
	     -1.5, 1.5, 1.e-2, 3.e6,
	     0., 2.1,
	     true, false, false);

  vector<TH1D*> trailEta = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
					ratio_ff_0, ratio_ff_1, ratio_ff_2,
					ratio_gf_0, ratio_gf_1, ratio_gf_2,
					egScale, egScaleErr,
					ffScale, ffScaleErr,
					gfScale, gfScaleErr, ttgjetsScale,
					"trailPhotonEta", true,
					800, -1.5, 1.5,
					channels[channel],
					5,
					false);

  TH1D * trailEta_460 = SignalHistoFromTree(signalScale_460, true, "trailPhotonEta", ggTree_460, "trailEta_460", "trailEta_460", 160, -1.5, 1.5);
  TH1D * trailEta_800 = SignalHistoFromTree(signalScale_800, true, "trailPhotonEta", ggTree_800, "trailEta_800", "trailEta_800", 160, -1.5, 1.5);
  TH1D * trailEta_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "trailPhotonEta", qcd30to40Tree, "trailEta_qcd30to40", "trailEta_qcd30to40", 160, -1.5, 1.5);
  TH1D * trailEta_qcd40 = SignalHistoFromTree(qcd40Scale, true, "trailPhotonEta", qcd40Tree, "trailEta_qcd40", "trailEta_qcd40", 160, -1.5, 1.5);

  prettyPlot(trailEta[0], trailEta[1], trailEta[2], trailEta[3], trailEta[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     trailEta_460, trailEta_800,
	     trailEta_qcd30to40, trailEta_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "trailEta", "#eta of trailing #gamma", "Number of Events",
	     -1.5, 1.5, 1.e-2, 3.e6,
	     0., 2.1,
	     true, false, false);

  vector<TH1D*> leadPhi = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				       ratio_ff_0, ratio_ff_1, ratio_ff_2,
				       ratio_gf_0, ratio_gf_1, ratio_gf_2,
				       egScale, egScaleErr,
				       ffScale, ffScaleErr,
				       gfScale, gfScaleErr, ttgjetsScale,
				       "leadPhotonPhi", true,
				       630, -3.14159, 3.14159,
				       channels[channel],
				       5,
				       false);

  TH1D * leadPhi_460 = SignalHistoFromTree(signalScale_460, true, "leadPhotonPhi", ggTree_460, "leadPhi_460", "leadPhi_460", 126, -3.14159, 3.14159);
  TH1D * leadPhi_800 = SignalHistoFromTree(signalScale_800, true, "leadPhotonPhi", ggTree_800, "leadPhi_800", "leadPhi_800", 126, -3.14159, 3.14159);
  TH1D * leadPhi_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "leadPhotonPhi", qcd30to40Tree, "leadPhi_qcd30to40", "leadPhi_qcd30to40", 126, -3.14159, 3.14159);
  TH1D * leadPhi_qcd40 = SignalHistoFromTree(qcd40Scale, true, "leadPhotonPhi", qcd40Tree, "leadPhi_qcd40", "leadPhi_qcd40", 126, -3.14159, 3.14159);

  prettyPlot(leadPhi[0], leadPhi[1], leadPhi[2], leadPhi[3], leadPhi[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     leadPhi_460, leadPhi_800,
	     leadPhi_qcd30to40, leadPhi_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "leadPhi", "#phi of leading #gamma", "Number of Events",
	     -3.2, 3.2, 1.e-2, 3.e6,
	     0., 2.1,
	     false, false, false);

  vector<TH1D*> trailPhi = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
					ratio_ff_0, ratio_ff_1, ratio_ff_2,
					ratio_gf_0, ratio_gf_1, ratio_gf_2,
					egScale, egScaleErr,
					ffScale, ffScaleErr,
					gfScale, gfScaleErr, ttgjetsScale,
					"trailPhotonPhi", true,
					630, -3.14159, 3.14159,
					channels[channel],
					5,
					false);
  
  TH1D * trailPhi_460 = SignalHistoFromTree(signalScale_460, true, "trailPhotonPhi", ggTree_460, "trailPhi_460", "trailPhi_460", 126, -3.14159, 3.14159);
  TH1D * trailPhi_800 = SignalHistoFromTree(signalScale_800, true, "trailPhotonPhi", ggTree_800, "trailPhi_800", "trailPhi_800", 126, -3.14159, 3.14159);
  TH1D * trailPhi_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "trailPhotonPhi", qcd30to40Tree, "trailPhi_qcd30to40", "trailPhi_qcd30to40", 126, -3.14159, 3.14159);
  TH1D * trailPhi_qcd40 = SignalHistoFromTree(qcd40Scale, true, "trailPhotonPhi", qcd40Tree, "trailPhi_qcd40", "trailPhi_qcd40", 126, -3.14159, 3.14159);

  prettyPlot(trailPhi[0], trailPhi[1], trailPhi[2], trailPhi[3], trailPhi[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     trailPhi_460, trailPhi_800,
	     trailPhi_qcd30to40, trailPhi_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "trailPhi", "#phi of trailing #gamma", "Number of Events",
	     -3.2, 3.2, 1.e-2, 3.e6,
	     0., 2.1,
	     true, false, false);

  vector<TH1D*> leadEt = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				      ratio_ff_0, ratio_ff_1, ratio_ff_2,
				      ratio_gf_0, ratio_gf_1, ratio_gf_2,
				      egScale, egScaleErr,
				      ffScale, ffScaleErr,
				      gfScale, gfScaleErr, ttgjetsScale,
				      "leadPhotonEt", true,
				      400, 0., 2000.,
				      channels[channel],
				      2,
				      true);

  TH1D * leadEt_460 = SignalHistoFromTree(signalScale_460, true, "leadPhotonEt", ggTree_460, "leadEt_460", "leadEt_460", 200, 0., 2000.);
  TH1D * leadEt_800 = SignalHistoFromTree(signalScale_800, true, "leadPhotonEt", ggTree_800, "leadEt_800", "leadEt_800", 200, 0., 2000.);
  TH1D * leadEt_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "leadPhotonEt", qcd30to40Tree, "leadEt_qcd30to40", "leadEt_qcd30to40", 200, 0., 2000.);
  TH1D * leadEt_qcd40 = SignalHistoFromTree(qcd40Scale, true, "leadPhotonEt", qcd40Tree, "leadEt_qcd40", "leadEt_qcd40", 200, 0., 2000.);

  prettyPlot(leadEt[0], leadEt[1], leadEt[2], leadEt[3], leadEt[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     leadEt_460, leadEt_800,
	     leadEt_qcd30to40, leadEt_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "leadEt", "Et of leading #gamma", "Number of Events",
	     0, 1200, 1.e-2, 3.e6,
	     0., 5.1,
	     true, true, true);

  vector<TH1D*> trailEt = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				       ratio_ff_0, ratio_ff_1, ratio_ff_2,
				       ratio_gf_0, ratio_gf_1, ratio_gf_2,
				       egScale, egScaleErr,
				       ffScale, ffScaleErr,
				       gfScale, gfScaleErr, ttgjetsScale,
				       "trailPhotonEt", true,
				       400, 0., 2000.,
				       channels[channel],
				       2,
				       true);

  TH1D * trailEt_460 = SignalHistoFromTree(signalScale_460, true, "trailPhotonEt", ggTree_460, "trailEt_460", "trailEt_460", 200, 0., 2000.);
  TH1D * trailEt_800 = SignalHistoFromTree(signalScale_800, true, "trailPhotonEt", ggTree_800, "trailEt_800", "trailEt_800", 200, 0., 2000.);
  TH1D * trailEt_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "trailPhotonEt", qcd30to40Tree, "trailEt_qcd30to40", "trailEt_qcd30to40", 200, 0., 2000.);
  TH1D * trailEt_qcd40 = SignalHistoFromTree(qcd40Scale, true, "trailPhotonEt", qcd40Tree, "trailEt_qcd40", "trailEt_qcd40", 200, 0., 2000.);

  prettyPlot(trailEt[0], trailEt[1], trailEt[2], trailEt[3], trailEt[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     trailEt_460, trailEt_800,
	     trailEt_qcd30to40, trailEt_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "trailEt", "Et of trailing #gamma", "Number of Events",
	     0, 800, 1.e-2, 3.e6,
	     0., 9.9,
	     true, true, true);

  vector<TH1D*> nJets = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				     ratio_ff_0, ratio_ff_1, ratio_ff_2,
				     ratio_gf_0, ratio_gf_1, ratio_gf_2,
				     egScale, egScaleErr,
				     ffScale, ffScaleErr,
				     gfScale, gfScaleErr, ttgjetsScale,
				     "Njets", false,
				     20, 0., 20.,
				     channels[channel],
				     1,
				     true);

  TH1D * nJets_460 = SignalHistoFromTree(signalScale_460, false, "Njets", ggTree_460, "nJets_460", "nJets_460", 20, 0., 20.);
  TH1D * nJets_800 = SignalHistoFromTree(signalScale_800, false, "Njets", ggTree_800, "nJets_800", "nJets_800", 20, 0., 20.);
  TH1D * nJets_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, false, "Njets", qcd30to40Tree, "nJets_qcd30to40", "nJets_qcd30to40", 20, 0., 20.);
  TH1D * nJets_qcd40 = SignalHistoFromTree(qcd40Scale, false, "Njets", qcd40Tree, "nJets_qcd40", "nJets_qcd40", 20, 0., 20.);

  prettyPlot(nJets[0], nJets[1], nJets[2], nJets[3], nJets[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     nJets_460, nJets_800,
	     nJets_qcd30to40, nJets_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "nJets", "nJets", "Number of Events",
	     0, 9, 1.5e-2, 3.e6,
	     0., 2.1,
	     true, true, false);

  vector<TH1D*> nBtags = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				      ratio_ff_0, ratio_ff_1, ratio_ff_2,
				      ratio_gf_0, ratio_gf_1, ratio_gf_2,
				      egScale, egScaleErr,
				      ffScale, ffScaleErr,
				      gfScale, gfScaleErr, ttgjetsScale,
				      "Nbtags", false,
				      20, 0., 20.,
				      channels[channel],
				      1,
				      true);

  TH1D * nBtags_460 = SignalHistoFromTree(signalScale_460, false, "Nbtags", ggTree_460, "nBtags_460", "nBtags_460", 20, 0., 20.);
  TH1D * nBtags_800 = SignalHistoFromTree(signalScale_800, false, "Nbtags", ggTree_800, "nBtags_800", "nBtags_800", 20, 0., 20.);
  TH1D * nBtags_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, false, "Nbtags", qcd30to40Tree, "nBtags_qcd30to40", "nBtags_qcd30to40", 20, 0., 20.);
  TH1D * nBtags_qcd40 = SignalHistoFromTree(qcd40Scale, false, "Nbtags", qcd40Tree, "nBtags_qcd40", "nBtags_qcd40", 20, 0., 20.);

  prettyPlot(nBtags[0], nBtags[1], nBtags[2], nBtags[3], nBtags[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     nBtags_460, nBtags_800,
	     nBtags_qcd30to40, nBtags_qcd40,
	     intLumi, ff_failure, 
	     channels[channel], "nBtags", "nBtags", "Number of Events",
	     0, 4, 1.5e-2, 3.e6,
	     0., 2.1,
	     true, true, false);

  vector<TH1D*> ht = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				  ratio_ff_0, ratio_ff_1, ratio_ff_2,
				  ratio_gf_0, ratio_gf_1, ratio_gf_2,
				  egScale, egScaleErr,
				  ffScale, ffScaleErr,
				  gfScale, gfScaleErr, ttgjetsScale,
				  "HT", true,
				  400, 0., 2000.,
				  channels[channel],
				  4,
				  false);

  TH1D * ht_460 = SignalHistoFromTree(signalScale_460, true, "HT", ggTree_460, "ht_460", "ht_460", 100, 0., 2000.);
  TH1D * ht_800 = SignalHistoFromTree(signalScale_800, true, "HT", ggTree_800, "ht_800", "ht_800", 100, 0., 2000.);
  TH1D * ht_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "HT", qcd30to40Tree, "ht_qcd30to40", "ht_qcd30to40", 100, 0., 2000.);
  TH1D * ht_qcd40 = SignalHistoFromTree(qcd40Scale, true, "HT", qcd40Tree, "ht_qcd40", "ht_qcd40", 100, 0., 2000.);

  prettyPlot(ht[0], ht[1], ht[2], ht[3], ht[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     ht_460, ht_800,
	     ht_qcd30to40, ht_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "HT", "HT (GeV/c)", "Number of Events",
	     0, 2000, 1.e-2, 3.e6,
	     0., 5.1,
	     true, true, true);
  
  vector<TH1D*> ht_jets = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				       ratio_ff_0, ratio_ff_1, ratio_ff_2,
				       ratio_gf_0, ratio_gf_1, ratio_gf_2,
				       egScale, egScaleErr,
				       ffScale, ffScaleErr,
				       gfScale, gfScaleErr, ttgjetsScale,
				       "HT_jets", true,
				       400, 0., 2000.,
				       channels[channel],
				       4,
				       false);

  TH1D * ht_jets_460 = SignalHistoFromTree(signalScale_460, true, "HT_jets", ggTree_460, "ht_jets_460", "ht_jets_460", 100, 0., 2000.);
  TH1D * ht_jets_800 = SignalHistoFromTree(signalScale_800, true, "HT_jets", ggTree_800, "ht_jets_800", "ht_jets_800", 100, 0., 2000.);
  TH1D * ht_jets_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "HT_jets", qcd30to40Tree, "ht_jets_qcd30to40", "ht_jets_qcd30to40", 100, 0., 2000.);
  TH1D * ht_jets_qcd40 = SignalHistoFromTree(qcd40Scale, true, "HT_jets", qcd40Tree, "ht_jets_qcd40", "ht_jets_qcd40", 100, 0., 2000.);

  prettyPlot(ht_jets[0], ht_jets[1], ht_jets[2], ht_jets[3], ht_jets[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     ht_jets_460, ht_jets_800,
	     ht_jets_qcd30to40, ht_jets_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "HT_jets", "HT_jets (GeV/c)", "Number of Events",
	     0, 2000, 1.e-2, 3.e6,
	     0., 5.1,
	     true, true, true);

  vector<TH1D*> hadronicPt = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
					  ratio_ff_0, ratio_ff_1, ratio_ff_2,
					  ratio_gf_0, ratio_gf_1, ratio_gf_2,
					  egScale, egScaleErr,
					  ffScale, ffScaleErr,
					  gfScale, gfScaleErr, ttgjetsScale,
					  "hadronic_pt", true,
					  400, 0., 2000.,
					  channels[channel],
					  4,
					  true);

  TH1D * hadronicPt_460 = SignalHistoFromTree(signalScale_460, true, "hadronic_pt", ggTree_460, "hadronicPt_460", "hadronicPt_460", 100, 0., 2000.);
  TH1D * hadronicPt_800 = SignalHistoFromTree(signalScale_800, true, "hadronic_pt", ggTree_800, "hadronicPt_800", "hadronicPt_800", 100, 0., 2000.);
  TH1D * hadronicPt_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "hadronic_pt", qcd30to40Tree, "hadronicPt_qcd30to40", "hadronicPt_qcd30to40", 100, 0., 2000.);
  TH1D * hadronicPt_qcd40 = SignalHistoFromTree(qcd40Scale, true, "hadronic_pt", qcd40Tree, "hadronicPt_qcd40", "hadronicPt_qcd40", 100, 0., 2000.);

  prettyPlot(hadronicPt[0], hadronicPt[1], hadronicPt[2], hadronicPt[3], hadronicPt[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     hadronicPt_460, hadronicPt_800,
	     hadronicPt_qcd30to40, hadronicPt_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "hadronicPt", "hadronic Pt (GeV/c)", "Number of Events",
	     0, 1250, 1.e-2, 3.e6,
	     0., 5.1,
	     true, true, true);

  vector<TH1D*> dR = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				  ratio_ff_0, ratio_ff_1, ratio_ff_2,
				  ratio_gf_0, ratio_gf_1, ratio_gf_2,
				  egScale, egScaleErr,
				  ffScale, ffScaleErr,
				  gfScale, gfScaleErr, ttgjetsScale,
				  "photon_dR", true,
				  500, 0., 5.,
				  channels[channel],
				  10,
				  false);

  TH1D * dR_460 = SignalHistoFromTree(signalScale_460, true, "photon_dR", ggTree_460, "dR_460", "dR_460", 50, 0., 5.);
  TH1D * dR_800 = SignalHistoFromTree(signalScale_800, true, "photon_dR", ggTree_800, "dR_800", "dR_800", 50, 0., 5.);
  TH1D * dR_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "photon_dR", qcd30to40Tree, "dR_qcd30to40", "dR_qcd30to40", 50, 0., 5.);
  TH1D * dR_qcd40 = SignalHistoFromTree(qcd40Scale, true, "photon_dR", qcd40Tree, "dR_qcd40", "dR_qcd40", 50, 0., 5.);

  prettyPlot(dR[0], dR[1], dR[2], dR[3], dR[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     dR_460, dR_800,
	     dR_qcd30to40, dR_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "dR", "#DeltaR_{#gamma#gamma}", "Number of Events",
	     0.4, 5., 1.e-2, 3.e6,
	     0., 1.9,
	     true, false, false);

  vector<TH1D*> dPhi = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				    ratio_ff_0, ratio_ff_1, ratio_ff_2,
				    ratio_gf_0, ratio_gf_1, ratio_gf_2,
				    egScale, egScaleErr,
				    ffScale, ffScaleErr,
				    gfScale, gfScaleErr, ttgjetsScale,
				    "photon_dPhi", true,
				    315, 0., 3.14159,
				    channels[channel],
				    9,
				    false);
  
  TH1D * dPhi_460 = SignalHistoFromTree(signalScale_460, true, "photon_dPhi", ggTree_460, "dPhi_460", "dPhi_460", 35, 0., 3.14159);
  TH1D * dPhi_800 = SignalHistoFromTree(signalScale_800, true, "photon_dPhi", ggTree_800, "dPhi_800", "dPhi_800", 35, 0., 3.14159);
  TH1D * dPhi_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "photon_dPhi", qcd30to40Tree, "dPhi_qcd30to40", "dPhi_qcd30to40", 35, 0., 3.14159);
  TH1D * dPhi_qcd40 = SignalHistoFromTree(qcd40Scale, true, "photon_dPhi", qcd40Tree, "dPhi_qcd40", "dPhi_qcd40", 35, 0., 3.14159);

  prettyPlot(dPhi[0], dPhi[1], dPhi[2], dPhi[3], dPhi[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     dPhi_460, dPhi_800,
	     dPhi_qcd30to40, dPhi_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "dPhi", "#Delta#phi_{#gamma#gamma}", "Number of Events",
	     0, 3.14159, 1.e-2, 3.e6,
	     0., 1.9,
	     true, false, false);

  vector<TH1D*> invmass = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				       ratio_ff_0, ratio_ff_1, ratio_ff_2,
				       ratio_gf_0, ratio_gf_1, ratio_gf_2,
				       egScale, egScaleErr,
				       ffScale, ffScaleErr,
				       gfScale, gfScaleErr, ttgjetsScale,
				       "invmass", true,
				       400, 0., 2000.,
				       channels[channel],
				       4,
				       true);

  TH1D * invmass_460 = SignalHistoFromTree(signalScale_460, true, "invmass", ggTree_460, "invmass_460", "invmass_460", 100, 0., 2000.);
  TH1D * invmass_800 = SignalHistoFromTree(signalScale_800, true, "invmass", ggTree_800, "invmass_800", "invmass_800", 100, 0., 2000.);
  TH1D * invmass_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "invmass", qcd30to40Tree, "invmass_qcd30to40", "invmass_qcd30to40", 100, 0., 2000.);
  TH1D * invmass_qcd40 = SignalHistoFromTree(qcd40Scale, true, "invmass", qcd40Tree, "invmass_qcd40", "invmass_qcd40", 100, 0., 2000.);

  prettyPlot(invmass[0], invmass[1], invmass[2], invmass[3], invmass[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     invmass_460, invmass_800,
	     invmass_qcd30to40, invmass_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "invmass", "m_{#gamma#gamma} (GeV/c^{2})", "Number of Events",
	     0, 2000, 1.e-2, 3.e6,
	     0., 11.5,
	     true, true, true);

  vector<TH1D*> max_csv = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				       ratio_ff_0, ratio_ff_1, ratio_ff_2,
				       ratio_gf_0, ratio_gf_1, ratio_gf_2,
				       egScale, egScaleErr,
				       ffScale, ffScaleErr,
				       gfScale, gfScaleErr, ttgjetsScale,
				       "max_csv", true,
				       50, 0., 1.,
				       channels[channel],
				       1,
				       false);
  
  TH1D * max_csv_460 = SignalHistoFromTree(signalScale_460, true, "max_csv", ggTree_460, "max_csv_460", "max_csv_460", 50, 0., 1.);
  TH1D * max_csv_800 = SignalHistoFromTree(signalScale_800, true, "max_csv", ggTree_800, "max_csv_800", "max_csv_800", 50, 0., 1.);
  TH1D * max_csv_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "max_csv", qcd30to40Tree, "max_csv_qcd30to40", "max_csv_qcd30to40", 50, 0., 1.);
  TH1D * max_csv_qcd40 = SignalHistoFromTree(qcd40Scale, true, "max_csv", qcd40Tree, "max_csv_qcd40", "max_csv_qcd40", 50, 0., 1.);

  prettyPlot(max_csv[0], max_csv[1], max_csv[2], max_csv[3], max_csv[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     max_csv_460, max_csv_460,
	     max_csv_qcd30to40, max_csv_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "max_csv", "max csv", "Number of Events",
	     0, 1, 1.e-2, 3.e6,
	     0., 1.9,
	     true, true, false);

  vector<TH1D*> submax_csv = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
					  ratio_ff_0, ratio_ff_1, ratio_ff_2,
					  ratio_gf_0, ratio_gf_1, ratio_gf_2,
					  egScale, egScaleErr,
					  ffScale, ffScaleErr,
					  gfScale, gfScaleErr, ttgjetsScale,
					  "submax_csv", true,
					  50, 0., 1.,
					  channels[channel],
					  1,
					  false);
  
  TH1D * submax_csv_460 = SignalHistoFromTree(signalScale_460, true, "submax_csv", ggTree_460, "submax_csv_460", "submax_csv_460", 50, 0., 1.);
  TH1D * submax_csv_800 = SignalHistoFromTree(signalScale_800, true, "submax_csv", ggTree_800, "submax_csv_800", "submax_csv_800", 50, 0., 1.);
  TH1D * submax_csv_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "submax_csv", qcd30to40Tree, "submax_csv_qcd30to40", "submax_csv_qcd30to40", 50, 0., 1.);
  TH1D * submax_csv_qcd40 = SignalHistoFromTree(qcd40Scale, true, "submax_csv", qcd40Tree, "submax_csv_qcd40", "submax_csv_qcd40", 50, 0., 1.);

  prettyPlot(submax_csv[0], submax_csv[1], submax_csv[2], submax_csv[3], submax_csv[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     submax_csv_460, submax_csv_800,
	     submax_csv_qcd30to40, submax_csv_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "submax_csv", "second max csv", "Number of Events",
	     0, 1, 1.e-2, 3.e6,
	     0., 5.,
	     true, true, false);

  vector<TH1D*> nPV = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				   ratio_ff_0, ratio_ff_1, ratio_ff_2,
				   ratio_gf_0, ratio_gf_1, ratio_gf_2,
				   egScale, egScaleErr,
				   ffScale, ffScaleErr,
				   gfScale, gfScaleErr, ttgjetsScale,
				   "nPV", false,
				   50, 0., 50.,
				   channels[channel],
				   1,
				   false);

  TH1D * nPV_460 = SignalHistoFromTree(signalScale_460, false, "nPV", ggTree_460, "nPV_460", "nPV_460", 50, 0., 50.);
  TH1D * nPV_800 = SignalHistoFromTree(signalScale_800, false, "nPV", ggTree_800, "nPV_800", "nPV_800", 50, 0., 50.);
  TH1D * nPV_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, false, "nPV", qcd30to40Tree, "nPV_qcd30to40", "nPV_qcd30to40", 50, 0., 50.);
  TH1D * nPV_qcd40 = SignalHistoFromTree(qcd40Scale, false, "nPV", qcd40Tree, "nPV_qcd40", "nPV_qcd40", 50, 0., 50.);

  prettyPlot(nPV[0], nPV[1], nPV[2], nPV[3], nPV[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     nPV_460, nPV_800,
	     nPV_qcd30to40, nPV_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "nPV", "nPV", "Number of Events",
	     0, 50, 1.e-2, 3.e6,
	     0., 5.,
	     true, true, false);

  vector<TH1D*> nMuons = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
				      ratio_ff_0, ratio_ff_1, ratio_ff_2,
				      ratio_gf_0, ratio_gf_1, ratio_gf_2,
				      egScale, egScaleErr,
				      ffScale, ffScaleErr,
				      gfScale, gfScaleErr, ttgjetsScale,
				      "Nmuons", false,
				      4, 0., 4.,
				      channels[channel],
				      1,
				      false);
  
  TH1D * nMuons_460 = SignalHistoFromTree(signalScale_460, false, "Nmuons", ggTree_460, "nMuons_460", "nMuons_460", 4, 0., 4.);
  TH1D * nMuons_800 = SignalHistoFromTree(signalScale_800, false, "Nmuons", ggTree_800, "nMuons_800", "nMuons_800", 4, 0., 4.);
  TH1D * nMuons_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, false, "Nmuons", qcd30to40Tree, "nMuons_qcd30to40", "nMuons_qcd30to40", 4, 0., 4.);
  TH1D * nMuons_qcd40 = SignalHistoFromTree(qcd40Scale, false, "Nmuons", qcd40Tree, "nMuons_qcd40", "nMuons_qcd40", 4, 0., 4.);

  prettyPlot(nMuons[0], nMuons[1], nMuons[2], nMuons[3], nMuons[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     nMuons_460, nMuons_800,
	     nMuons_qcd30to40, nMuons_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "nMuons", "nMuons", "Number of Events",
	     0, 3, 1.5e-2, 3.e6,
	     0., 3.6,
	     true, true, false);

  vector<TH1D*> nElectrons = preparePlots(ggTree,  egTree, ffTree, gfTree, ttgjetsTree,
					  ratio_ff_0, ratio_ff_1, ratio_ff_2,
					  ratio_gf_0, ratio_gf_1, ratio_gf_2,
					  egScale, egScaleErr,
					  ffScale, ffScaleErr,
					  gfScale, gfScaleErr, ttgjetsScale,
					  "Nelectrons", false,
					  4, 0., 4.,
					  channels[channel],
					  1,
					  false);
  
  TH1D * nElectrons_460 = SignalHistoFromTree(signalScale_460, false, "Nelectrons", ggTree_460, "nElectrons_460", "nElectrons_460", 4, 0., 4.);
  TH1D * nElectrons_800 = SignalHistoFromTree(signalScale_800, false, "Nelectrons", ggTree_800, "nElectrons_800", "nElectrons_800", 4, 0., 4.);
  TH1D * nElectrons_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, false, "Nelectrons", qcd30to40Tree, "nElectrons_qcd30to40", "nElectrons_qcd30to40", 4, 0., 4.);
  TH1D * nElectrons_qcd40 = SignalHistoFromTree(qcd40Scale, false, "Nelectrons", qcd40Tree, "nElectrons_qcd40", "nElectrons_qcd40", 4, 0., 4.);

  prettyPlot(nElectrons[0], nElectrons[1], nElectrons[2], nElectrons[3], nElectrons[4],
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     nElectrons_460, nElectrons_800,
	     nElectrons_qcd30to40, nElectrons_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "nElectrons", "nElectrons", "Number of Events",
	     0, 3, 1.5e-2, 3.e6,
	     0., 1.9,
	     true, true, false);

  // Now create combined qcd+ewk plots
  // Divide by bin width

  h_gg = (TH1D*)DivideByBinWidth(h_gg);
  h_ewk = (TH1D*)DivideByBinWidth(h_ewk);
  h_qcd_gf = (TH1D*)DivideByBinWidth(h_qcd_gf);
  h_qcd_ff = (TH1D*)DivideByBinWidth(h_qcd_ff);
  h_qcd_ee = (TH1D*)DivideByBinWidth(h_qcd_ee);

  TH1D * h_sigMET_a = (TH1D*)fSig460->Get("met_mst_460_m1_175_gg_"+channels[channel]);
  Double_t xsec_a = 0.147492;
  Double_t acc_a = 1. / 15000.;
  h_sigMET_a->Scale(acc_a * xsec_a * intLumi_int * 1.019 * 1.019);
  TH1D * h_sigMET_aReb = (TH1D*)h_sigMET_a->Rebin(nBins, "h_sigMET_aReb", xbins);

  for(int i = 0; i < h_sigMET_aReb->GetNbinsX(); i++) h_sigMET_aReb->SetBinContent(i+1, h_sigMET_aReb->GetBinContent(i+1) / h_sigMET_aReb->GetBinWidth(i+1));

  //TH1D * h_sigMET_b = (TH1D*)fSig800->Get("met_mst_560_m1_325_gg_"+channels[channel]);
  TH1D * h_sigMET_b = (TH1D*)fSig800->Get("met_mst_460_m1_175_gg_"+channels[channel]); // DURP
  Double_t xsec_b = 0.0399591;
  Double_t acc_b = 1. / 15000.;
  h_sigMET_b->Scale(acc_b * xsec_b * intLumi_int * 1.019 * 1.019);
  TH1D * h_sigMET_bReb = (TH1D*)h_sigMET_b->Rebin(nBins, "h_sigMET_bReb", xbins);

  for(int i = 0; i < h_sigMET_bReb->GetNbinsX(); i++) h_sigMET_bReb->SetBinContent(i+1, h_sigMET_bReb->GetBinContent(i+1) / h_sigMET_bReb->GetBinWidth(i+1));

  met_ttgjets = (TH1D*)met_ttgjets->Rebin(nBins, "met_ttgjets_reb", xbins);
  for(int i = 0; i < met_ttgjets->GetNbinsX(); i++) {
    met_ttgjets->SetBinContent(i+1, met_ttgjets->GetBinContent(i+1) / met_ttgjets->GetBinWidth(i+1));
    met_ttgjets->SetBinError(i+1, met_ttgjets->GetBinError(i+1) / met_ttgjets->GetBinWidth(i+1));
  }

  TH1D * met_qcd30to40 = SignalHistoFromTree(qcd30to40Scale, true, "pfMET", qcd30to40Tree, "met_qcd30to40", "met_qcd30to40", 400, 0., 2000.);
  met_qcd30to40 = (TH1D*)met_qcd30to40->Rebin(nBins, "met_qcd30to40_reb", xbins);
  for(int i = 0; i < met_qcd30to40->GetNbinsX(); i++) {
    met_qcd30to40->SetBinContent(i+1, met_qcd30to40->GetBinContent(i+1) / met_qcd30to40->GetBinWidth(i+1));
    met_qcd30to40->SetBinError(i+1, met_qcd30to40->GetBinError(i+1) / met_qcd30to40->GetBinWidth(i+1));
  }

  TH1D * met_qcd40 = SignalHistoFromTree(qcd40Scale, true, "pfMET", qcd40Tree, "met_qcd40", "met_qcd40", 400, 0., 2000.);
  met_qcd40 = (TH1D*)met_qcd40->Rebin(nBins, "met_qcd40_reb", xbins);
  for(int i = 0; i < met_qcd40->GetNbinsX(); i++) {
    met_qcd40->SetBinContent(i+1, met_qcd40->GetBinContent(i+1) / met_qcd40->GetBinWidth(i+1));
    met_qcd40->SetBinError(i+1, met_qcd40->GetBinError(i+1) / met_qcd40->GetBinWidth(i+1));
  }

  prettyPlot(h_gg, h_ewk, h_qcd_ff, h_qcd_gf, met_ttgjets,
	     useDifferenceSystematic, useTTGJets, useFF, useMCforQCD,
	     h_sigMET_aReb, h_sigMET_bReb,
	     met_qcd30to40, met_qcd40,
	     intLumi, ff_failure,
	     channels[channel], "MET", "#slash{E}_{T} (GeV)", "Number of Events / GeV",
	     xbins[0], xbins[nBins], 7.e-4, 25000.,
	     0., 9.1,
	     true, true, true);

  // make a ff/ee ratio plot
  TCanvas * can2 = new TCanvas("can2", "Plot", 10, 10, 2000, 2000);

  TH1D * ffDIVee = (TH1D*)h_qcd_ff->Clone("ffDIVee");
  ffDIVee->Divide(h_qcd_ee);

  ffDIVee->Draw("e1");
  TLine * oneLine = new TLine(0, 1, ffDIVee->GetXaxis()->GetBinUpEdge(ffDIVee->GetNbinsX()), 1);
  oneLine->SetLineStyle(2);
  oneLine->Draw();

  can2->SaveAs("ffDIVee_"+channels[channel]+gifOrPdf);

  in->Close();
  fSig460->Close();
  fSig800->Close();
  fQCD30to40->Close();
  fQCD40->Close();

}
