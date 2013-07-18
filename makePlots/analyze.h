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

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooPolynomial.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "/uscms_data/d2/bfrancis/btagRA3/CMSSW_4_2_8/src/PhysicsTools/TagAndProbe/interface/RooCMSShape.h"

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

using namespace std;
using namespace RooFit;

const TString gifOrPdf = ".pdf";

void formatTable(TH1D * h_gg,
		 TH1D * h_ewk, TH1D * ewk_norm2, Float_t fakeRate, Float_t fakeRate_sys, Float_t egScale,
		 TH1D * qcd_ff, TH1D * ff_norm2, Float_t ffScale,
		 TH1D * qcd_gf, TH1D * gf_norm2, Float_t gfScale,
		 TH1D * qcd_ee, TH1D * ee_norm2, Float_t eeScale,
		 TH1D * ttgjets, bool useTTGJets,
		 TString fileType) {

  FILE* tableFile = fopen("errorTable_"+fileType+".temp", "w");

  Double_t rangeLow[5] = {0, 0, 50, 80, 100};
  Double_t rangeHigh[5] = {20, 50, -1, -1, -1};

  Double_t binLow[5], binHigh[5];
  for(int i = 0; i < 5; i++) {
    binLow[i] = h_gg->GetXaxis()->FindBin(rangeLow[i]);
    binHigh[i] = (rangeHigh[i] == -1) ? -1 : h_gg->GetXaxis()->FindBin(rangeHigh[i]) - 1;
  }

  for(int i = 0; i < 5; i++) {

    Double_t gg, ggerr;
    gg = h_gg->IntegralAndError(binLow[i], binHigh[i], ggerr);
    fprintf(tableFile, "ggval%dx:%.0f\n", i+1, gg);
    fprintf(tableFile, "ggstat%dx:%.1f\n", i+1, ggerr);

    Double_t ff, fferr;
    ff = qcd_ff->IntegralAndError(binLow[i], binHigh[i], fferr) * ffScale;
    fferr *= ffScale;
    fprintf(tableFile, "ffval%dx:%.1f\n", i+1, ff);
    fprintf(tableFile, "ffstat%dx:%.2f\n", i+1, fferr);
    Double_t ff_norm = sqrt(ff_norm2->Integral(binLow[i], binHigh[i])) * ffScale;
    fprintf(tableFile, "ffnorm%dx:%.2f\n", i+1, ff_norm);

    Double_t gf, gferr;
    gf = qcd_gf->IntegralAndError(binLow[i], binHigh[i], gferr) * gfScale;
    gferr *= gfScale;
    fprintf(tableFile, "gfval%dx:%.1f\n", i+1, gf);
    fprintf(tableFile, "gfstat%dx:%.2f\n", i+1, gferr);
    Double_t gf_norm = sqrt(gf_norm2->Integral(binLow[i], binHigh[i])) * gfScale;
    fprintf(tableFile, "gfnorm%dx:%.2f\n", i+1, gf_norm);

    Double_t ee, eeerr;
    ee = qcd_ee->IntegralAndError(binLow[i], binHigh[i], eeerr) * eeScale;
    eeerr *= eeScale;
    fprintf(tableFile, "eeval%dx:%.1f\n", i+1, ee);
    fprintf(tableFile, "eestat%dx:%.2f\n", i+1, eeerr);
    Double_t ee_norm = sqrt(ee_norm2->Integral(binLow[i], binHigh[i])) * eeScale;
    fprintf(tableFile, "eenorm%dx:%.2f\n", i+1, ee_norm);

    Double_t eg, egerr;
    eg = h_ewk->IntegralAndError(binLow[i], binHigh[i], egerr) * egScale;
    egerr *= egScale;
    fprintf(tableFile, "ewkval%dx:%.1f\n", i+1, eg);
    fprintf(tableFile, "ewkstat%dx:%.2f\n", i+1, egerr);
    Double_t ewk_norm = sqrt(ewk_norm2->Integral(binLow[i], binHigh[i])) * egScale;
    fprintf(tableFile, "ewknorm%dx:%.2f\n", i+1, ewk_norm);
    Double_t ewk_sys = eg * fakeRate_sys / fakeRate;
    fprintf(tableFile, "ewksyst%dx:%.2f\n", i+1, ewk_sys);
    
    Double_t ttgg, ttggerr;
    ttgg = ttgjets->IntegralAndError(binLow[i], binHigh[i], ttggerr);
    if(useTTGJets) fprintf(tableFile, "ttgval%dx:%.1f\nttgstat%dx:%.2f\n", i+1, ttgg, i+1, ttggerr);

    Double_t bkg = useTTGJets ? eg + ff + ttgg : eg + ff;
    Double_t bkgstat = useTTGJets ? sqrt(egerr*egerr + fferr*fferr + ttggerr*ttggerr) : sqrt(egerr*egerr + fferr*fferr);
    Double_t bkgnorm = sqrt(ewk_norm*ewk_norm + ff_norm2->Integral(binLow[i], binHigh[i]));
    fprintf(tableFile, "fftotalval%dx:%.1f\nfftotalstat%dx:%.2f\nfftotalnorm%dx:%.2f\nfftotalsyst%dx:%.2f\n", i+1, bkg, i+1, bkgstat, i+1, bkgnorm, i+1, ewk_sys);

    bkg = useTTGJets ? eg + gf + ttgg : eg + gf;
    bkgstat = useTTGJets ? sqrt(egerr*egerr + gferr*gferr + ttggerr*ttggerr) : sqrt(egerr*egerr + gferr*gferr);
    bkgnorm = sqrt(ewk_norm*ewk_norm + gf_norm2->Integral(binLow[i], binHigh[i]));
    fprintf(tableFile, "gftotalval%dx:%.1f\ngftotalstat%dx:%.2f\ngftotalnorm%dx:%.2f\ngftotalsyst%dx:%.2f\n", i+1, bkg, i+1, bkgstat, i+1, bkgnorm, i+1, ewk_sys);

    bkg = useTTGJets ? eg + ee + ttgg : eg + ee;
    bkgstat = useTTGJets ? sqrt(egerr*egerr + eeerr*eeerr + ttggerr*ttggerr) : sqrt(egerr*egerr + eeerr*eeerr);
    bkgnorm = sqrt(ewk_norm*ewk_norm + ee_norm2->Integral(binLow[i], binHigh[i]));
    fprintf(tableFile, "eetotalval%dx:%.1f\neetotalstat%dx:%.2f\neetotalnorm%dx:%.2f\neetotalsyst%dx:%.2f\n", i+1, bkg, i+1, bkgstat, i+1, bkgnorm, i+1, ewk_sys);
    
  }

  fclose(tableFile);
}

// Return the weights (ratio) of diempt
TH1D * GetWeights(TH1D* ggDiEMpT, TH1D* h_diempt, Float_t gg_test, Float_t other_test) {
  TH1D * ratio = (TH1D*)ggDiEMpT->Clone();
  TH1D * diempt = (TH1D*)h_diempt->Clone();

  const int ndiemptbins = 30;
  Double_t diemptbins[ndiemptbins+1] = {
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
    55,
    60,
    65,
    70,
    75,
    80,
    85,
    90,
    95,
    100,
    110,
    120,
    130,
    140,
    150,
    200,
    300,
    400,
    600,
    1000};

  TH1D * ratio2 = (TH1D*)ratio->Rebin(ndiemptbins, "ratio2", diemptbins);
  TH1D * diempt2 = (TH1D*)diempt->Rebin(ndiemptbins, "diempt2", diemptbins);

  ratio2->Scale(1./gg_test);
  diempt2->Scale(1./other_test);

  ratio2->Divide(diempt2);

  return ratio2;
}

TH1D * GetAlternativeWeights(TTree * ggtree, TTree * bkgtree, TString variable, bool isAFloat, Int_t nbins, Double_t* xbins, TString name) {

  TH1D * gg = new TH1D("gg_"+variable+"_"+name, "gg_"+variable+"_"+name, nbins, xbins); gg->Sumw2();
  TH1D * bkg = new TH1D("bkg_"+variable+"_"+name, "bkg_"+variable+"_"+name, nbins, xbins); bkg->Sumw2();

  float met;
  ggtree->SetBranchAddress("pfMET", &met);
  bkgtree->SetBranchAddress("pfMET", &met);

  float var;
  int var_int;

  if(isAFloat) {
    ggtree->SetBranchAddress(variable, &var);
    bkgtree->SetBranchAddress(variable, &var);
  }
  else {
    ggtree->SetBranchAddress(variable, &var_int);
    bkgtree->SetBranchAddress(variable, &var_int);
  }

  for(int i = 0; i < ggtree->GetEntries(); i++) {
    ggtree->GetEntry(i);
    if(met >= 50.) continue;
    if(isAFloat) gg->Fill(var);
    else gg->Fill(var_int);
  }

  for(int i = 0; i < bkgtree->GetEntries(); i++) {
    bkgtree->GetEntry(i);
    if(met >= 50.) continue;
    if(isAFloat) bkg->Fill(var);
    else bkg->Fill(var_int);
  }

  //gg->Scale(1./(float)ggtree->GetEntries());
  //bkg->Scale(1./(float)bkgtree->GetEntries());

  gg->Scale(1./(float)gg->Integral());
  bkg->Scale(1./(float)bkg->Integral());

  TH1D * weight = (TH1D*)gg->Clone("weights_"+variable+"_"+name);
  weight->Divide(bkg);

  ggtree->ResetBranchAddresses();
  bkgtree->ResetBranchAddresses();

  return weight;
}

// Return a di-em pt reweighted met histogram
TH1D * Reweight(bool isAFloat, TH1D * h, TString variable, TH1D* ratio_0, TH1D* ratio_1, TH1D* ratio_2, TTree* tree, int type, bool save, TString name, bool isMC) {
  
  // type 0 = ff or ee full spectrum
  //      1 = ee signal (81-101)
  //      2 = ee low    (71-81)
  //      3 = ee high   (101-111)

  TH1D * rew = (TH1D*)h->Clone(variable+"_"+name);
  rew->Reset();

  TH1D * unrew = (TH1D*)rew->Clone(variable+"_"+name+"_unrew");

  float var_, diempt_, invmass_, puWeight_, puWeightErr_;
  int nJets_, var_int;
  if(variable != "invmass" && 
     variable != "Njets" && 
     variable != "diJetPt") {
    if(isAFloat) tree->SetBranchAddress(variable, &var_);
    else tree->SetBranchAddress(variable, &var_int);
  }
  //tree->SetBranchAddress("diEMpT", &diempt_);
  tree->SetBranchAddress("diJetPt", &diempt_);
  tree->SetBranchAddress("invmass", &invmass_);
  tree->SetBranchAddress("Njets", &nJets_);
  if(isMC) {
    tree->SetBranchAddress("puWeight", &puWeight_);
    tree->SetBranchAddress("puWeightError", &puWeightErr_);
  }

  Double_t ratioWeight, ratioWeightErr;
  Double_t puWeight, puWeightErr;

  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    if(variable == "invmass") var_ = invmass_;
    if(variable == "diJetPt") var_ = diempt_;
    if(variable == "Njets") var_int = nJets_;

    Double_t eventWeight = 1.;
    Double_t eventWeightErr = 0.;
    if(isMC) {
      eventWeight = puWeight_;
      eventWeightErr = puWeightErr_;
    }

    if(eventWeight == 0.) continue;

    Double_t ratioWeight, ratioWeightErr;
    if(nJets_ == 0) {
      ratioWeight = ratio_0->GetBinContent(ratio_0->FindBin(diempt_));
      ratioWeightErr = ratio_0->GetBinError(ratio_0->FindBin(diempt_));
    }
    else if(nJets_ == 1) {
      ratioWeight = ratio_1->GetBinContent(ratio_1->FindBin(diempt_));
      ratioWeightErr = ratio_1->GetBinError(ratio_1->FindBin(diempt_));
    }
    else {
      ratioWeight = ratio_2->GetBinContent(ratio_2->FindBin(diempt_));
      ratioWeightErr = ratio_2->GetBinError(ratio_2->FindBin(diempt_));
    }

    if(ratioWeight == 0.) continue;

    Double_t oldError = (isAFloat) ? rew->GetBinError(rew->FindBin(var_)) : rew->GetBinError(rew->FindBin(var_int));

    if(
       (type == 1 && invmass_ > 81.0 && invmass_ < 101.0) ||
       (type == 2 && invmass_ > 71.0 && invmass_ < 81.0) ||
       (type == 3 && invmass_ > 101.0 && invmass_ < 111.0) ||
       (type == 0)
       ) {
      
      if(isAFloat) {
	rew->Fill(var_, 1.0*ratioWeight*eventWeight);
	unrew->Fill(var_);
      }
      else {
	rew->Fill(var_int, 1.0*ratioWeight*eventWeight);
	unrew->Fill(var_int);
      }

      Double_t newError = ratioWeightErr*ratioWeightErr / (ratioWeight*ratioWeight);
      newError += eventWeightErr*eventWeightErr / (eventWeight*eventWeight);
      newError = ratioWeight*eventWeight*sqrt(newError);

      if(isAFloat) rew->SetBinError(rew->FindBin(var_), sqrt(oldError*oldError + newError*newError));
      else rew->SetBinError(rew->FindBin(var_int), sqrt(oldError*oldError + newError*newError));
    }

  }  // loop over entries

  if(save) {
    TCanvas * canv = new TCanvas("canv", "Plot", 10, 10, 2000, 2000);

    if(ratio_0->Integral() > 0.) {
      ratio_0->Draw();
      canv->SaveAs("ratio_0_"+name+gifOrPdf);
    }

    if(ratio_1->Integral() > 0.) {
      ratio_1->Draw();
      canv->SaveAs("ratio_1_"+name+gifOrPdf);
    }

    if(ratio_2->Integral() > 0.) {
      ratio_2->Draw();
      canv->SaveAs("ratio_2_"+name+gifOrPdf);
    }

    //rew->Draw();
    //canv->SaveAs(name+"CenMETValRew"+gifOrPdf);

    delete canv;
  }

  tree->ResetBranchAddresses();

  return rew;
}

TH1D * AlternativeReweight(TString outVariable, bool outIsAFloat,
			   TString reweightVariable, bool rewByAFloat,
			   TH1D * h, // dummy to copy
			   TH1D * ratio, // actual weights
			   TTree * tree,
			   int type, bool save, TString name, bool isMC) {

  TH1D * rew = (TH1D*)h->Clone(name);
  rew->Reset();

  float outVar, rewVar, met;
  int outVar_int, rewVar_int;

  tree->SetBranchAddress("pfMET", &met);

  if(outVariable != "pfMET") {
    if(outIsAFloat) tree->SetBranchAddress(outVariable, &outVar);
    else tree->SetBranchAddress(outVariable, &outVar_int);
  }

  if(reweightVariable != outVariable) {
    if(rewByAFloat) tree->SetBranchAddress(reweightVariable, &rewVar);
    else tree->SetBranchAddress(reweightVariable, &rewVar_int);
  }

  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if(met >= 50.) continue;

    if(outVariable == "pfMET") outVar = met;

    if(reweightVariable == outVariable) {
      rewVar = outVar;
      rewVar_int = outVar_int;
    }

    Double_t ratioWeight = (rewByAFloat) ? ratio->GetBinContent(ratio->FindBin(rewVar)) : ratio->GetBinContent(ratio->FindBin(rewVar_int));
    Double_t ratioWeightErr = (rewByAFloat) ? ratio->GetBinError(ratio->FindBin(rewVar)) : ratio->GetBinError(ratio->FindBin(rewVar_int));

    if(ratioWeight == 0.) continue;

    Double_t oldError = (outIsAFloat) ? rew->GetBinError(rew->FindBin(outVar)) : rew->GetBinError(rew->FindBin(outVar_int));

    if(outIsAFloat) rew->Fill(outVar, 1.0*ratioWeight);
    else rew->Fill(outVar_int, 1.0*ratioWeight);

    Double_t newError = ratioWeightErr*ratioWeightErr / (ratioWeight*ratioWeight);

    if(outIsAFloat) rew->SetBinError(rew->FindBin(outVar), sqrt(oldError*oldError + newError*newError));
    else rew->SetBinError(rew->FindBin(outVar_int), sqrt(oldError*oldError + newError*newError));
    
  }  // loop over entries

  if(save) {
    TCanvas * canv = new TCanvas("canv", "Plot", 10, 10, 2000, 2000);

    if(ratio->Integral() > 0.) {
      ratio->Draw();
      canv->SaveAs("ratio_"+name+gifOrPdf);
    }

    delete canv;
  }

  tree->ResetBranchAddresses();

  return rew;
}

// Get number of events for fake rate calculations
pair<Float_t, Float_t> FitYield(TH1D* h_invmass, bool save, TString name) {

  TH1D * invmass = (TH1D*)h_invmass->Clone();
  invmass->Sumw2();

  RooRealVar x("x", name+" Invariant mass", 71, 111);
  RooDataHist data("data", "dataset", x, invmass);

  // Linear background
  //RooRealVar bkgd_poly_c1("slope of bkgd", "slope of background", 7.5, 0, 1500);
  //RooPolynomial bkgd_linear("bkgd_linear", "linear background", x, RooArgList(bkgd_poly_c1));
  //RooRealVar bkgd_yield("bkgd events", "background events", 1000, 0, 1000000);

  // RooCMSShape background
  RooRealVar aBkg("aBkg", "aBkg", 90.0, 75.0, 105.0, "GeV");
  RooRealVar bBkg("bBkg", "bBkg", 0.05 /*0.05, 0., 0.15*/);
  RooRealVar c("c", "c", 0.05 /*0.05, 0.0, 0.2*/);
  RooRealVar ZMass("bkgd_peak", "ZMass", 80.0, 70.0, 90.0, "GeV");  
  RooRealVar backgroundYield("bkgd events", "backgroundYield", 5000, 0, 12000);
  RooCMSShape background("background", "background", x, aBkg, bBkg, c, ZMass);

  RooRealVar cbmean("mass", "cbmean", 90, 85, 95);
  RooRealVar cbsigma("cb_sigma", "cbsigma", 4, 1, 10);
  RooRealVar cbsig("signal events", "cbsignal", 10000, 0, 1000000);
  RooRealVar n("n", "n", 50, 0, 1000);
  RooRealVar alpha("alpha", "alpha", 1, 0, 2);
  RooCBShape cball("cball", "crystal ball", x, cbmean, cbsigma, alpha, n);

  //RooAddPdf totalPDF("totalPDF", "totalPDF", RooArgList(bkgd_linear, cball), RooArgList(bkgd_yield, cbsig));
  RooAddPdf totalPDF("totalPDF", "totalPDF", RooArgList(background, cball), RooArgList(backgroundYield, cbsig));

  totalPDF.fitTo(data, Extended(), PrintLevel(-1));

  if(save) {
    TCanvas * canv = new TCanvas("canv", "Plot", 10, 10, 2000, 2000);

    RooPlot * xframe = x.frame(Title("Crstal Ball plus RooCMSShape background"));
    totalPDF.paramOn(xframe, Format("NE", AutoPrecision(2)), Layout(0.6, 0.98, 0.9));
    
    data.plotOn(xframe, LineColor(kBlack));
    totalPDF.plotOn(xframe, Components(background), LineColor(kRed));
    totalPDF.plotOn(xframe, LineColor(kBlue));
    xframe->Draw();
    canv->SaveAs(name+"_Mass_fitted"+gifOrPdf);
    
    delete canv;
    delete xframe;
  }

  pair<Float_t, Float_t> yield = make_pair(cbsig.getVal(), cbsig.getError());

  return yield;
}

TH1D * HistoFromTree(bool isAFloat, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t xlo, Double_t xhi, double restrictMet = false) {

  TH1D * h = new TH1D(name, title, nBins, xlo, xhi);
  h->Sumw2();

  Float_t met;
  tree->SetBranchAddress("pfMET", &met);

  Float_t var;
  Int_t var_int;
  if(variable != "pfMET") {
    if(isAFloat) tree->SetBranchAddress(variable, &var);
    else tree->SetBranchAddress(variable, &var_int);
  }
  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if(restrictMet && met >= 50.) continue;

    if(variable != "pfMET") {
      if(isAFloat) h->Fill(var);
      else h->Fill(var_int);
    }
    else h->Fill(met);

  }

  tree->ResetBranchAddresses();

  return h;
}

TH1D * SignalHistoFromTree(Float_t scale, bool isAFloat, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t xlo, Double_t xhi) {

  TH1D * h = new TH1D(name, title, nBins, xlo, xhi);
  h->Sumw2();

  Float_t var;
  Int_t var_int;
  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr, btagWeightUp, btagWeightDown;
  if(isAFloat) tree->SetBranchAddress(variable, &var);
  else tree->SetBranchAddress(variable, &var_int);
  tree->SetBranchAddress("pileupWeight", &puWeight);
  tree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  tree->SetBranchAddress("btagWeight", &btagWeight);
  tree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  tree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  tree->SetBranchAddress("btagWeightDown", &btagWeightDown);
  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    Float_t olderror = 0.;

    if(isAFloat) {
      olderror = h->GetBinError(h->FindBin(var));
      h->Fill(var, puWeight * btagWeight);
    }
    else {
      olderror = h->GetBinError(h->FindBin(var_int));
      h->Fill(var_int, puWeight * btagWeight);
    }
    
    // protection from weird 1200 weight errors...
    if(btagWeightErr > 20.) btagWeightErr = btagWeight;

    Float_t btagSFsys = (fabs(btagWeight - btagWeightUp) + fabs(btagWeight - btagWeightDown))/2.;
    Float_t btag_toterr = sqrt(btagWeightErr*btagWeightErr + btagSFsys*btagSFsys);

    Float_t addError2 = puWeight*puWeight*btag_toterr*btag_toterr + btagWeight*btagWeight*puWeightErr*puWeightErr;

    Float_t newerror = sqrt(olderror*olderror + addError2);

    if(isAFloat) h->SetBinError(h->FindBin(var), newerror);
    else h->SetBinError(h->FindBin(var_int), newerror);
      
  }

  h->Scale(scale);

  tree->ResetBranchAddresses();

  return h;
}

void calculateROC(TH1D * sig_a, TH1D * sig_b, TH1D * bkg, TString req, TString title) {

  TH1D * srootb_a = (TH1D*)sig_a->Clone("roc_a"); srootb_a->Reset();
  TH1D * srootb_b = (TH1D*)sig_b->Clone("roc_b"); srootb_b->Reset();

  int nbins = sig_a->GetNbinsX();
  Double_t x_a[nbins], x_b[nbins], y_a[nbins], y_b[nbins];

  Double_t s, b, serr, berr;

  for(int i = 0; i < nbins; i++) {
    s = sig_a->IntegralAndError(i+1, -1, serr, "");
    b = bkg->IntegralAndError(i+1, -1, berr, "");
    
    if(b == 0.) continue;

    x_a[i] = s / sig_a->Integral();
    y_a[i] = b / bkg->Integral();
    srootb_a->SetBinContent(i+1, x_a[i] / sqrt(y_a[i]));

    s = sig_b->IntegralAndError(i+1, -1, serr, "");
    x_b[i] = s / sig_b->Integral();
    y_b[i] = b / bkg->Integral();
    srootb_b->SetBinContent(i+1, x_b[i] / sqrt(y_b[i]));

  }

  TGraph * roc_a = new TGraph(nbins, x_a, y_a);
  TGraph * roc_b = new TGraph(nbins, x_b, y_b);

  roc_a->SetLineColor(kMagenta);
  srootb_a->SetLineColor(kMagenta);
  roc_b->SetLineColor(kBlue);
  srootb_b->SetLineColor(kBlue);

  TCanvas * canv = new TCanvas("roc_can_"+title+"_"+req, "ROC Plot", 10, 10, 2000, 2000);

  TPad * padhi = new TPad("padhi", "padhi", 0, 0.5, 1, 1);
  TPad * padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.5);

  padhi->Draw();
  padlo->Draw();
  padhi->cd();

  padhi->SetTickx(true);
  padhi->SetTicky(true);
  padhi->SetGridx(true);
  padhi->SetGridy(true);

  TH2D * blank = new TH2D("blank_"+title+"_"+req, "blank;#epsilon_{S};#epsilon_{B}", 1, 0, 1, 1, 0, 1);
  blank->Draw();
  roc_a->Draw("same L");
  roc_b->Draw("same L");

  padlo->cd();
  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);

  Int_t lastBin = 0;
  Int_t bestCutBin_a = 0;
  Int_t bestCutBin_b = 0;
  float bestDiscrim_a = 0;
  float bestDiscrim_b = 0;
  for(int i = 0; i < srootb_a->GetNbinsX(); i++) {
    if(srootb_a->GetBinContent(i+1) > 0 || srootb_b->GetBinContent(i+1) > 0) lastBin = i+1;
    if(srootb_a->GetBinContent(i+1) > bestDiscrim_a) {
      bestDiscrim_a = srootb_a->GetBinContent(i+1);
      bestCutBin_a = i+1;
    }
    if(srootb_b->GetBinContent(i+1) > bestDiscrim_b) {
      bestDiscrim_b = srootb_b->GetBinContent(i+1);
      bestCutBin_b = i+1;
    }
  }
  
  if(lastBin < srootb_a->GetNbinsX()) srootb_a->GetXaxis()->SetRangeUser(srootb_a->GetBinLowEdge(1), srootb_a->GetBinLowEdge(lastBin+1) * 1.1);
  if(srootb_b->GetMaximum() > srootb_a->GetMaximum()) srootb_a->GetYaxis()->SetRangeUser(0, 1.1 * srootb_b->GetMaximum());
  else srootb_a->GetYaxis()->SetRangeUser(0, 1.1 * srootb_a->GetMaximum());
  srootb_a->GetXaxis()->SetTitle("Lower cut on "+title);
  srootb_a->GetYaxis()->SetTitle("S / #sqrt{B}");
  srootb_a->Draw("hist");
  srootb_b->Draw("hist same");

  float lineMaxY = (srootb_b->GetMaximum() > srootb_a->GetMaximum()) ? 1.1 * srootb_b->GetMaximum() : 1.1 * srootb_a->GetMaximum();

  TLine * bestCutLine_a = new TLine(srootb_a->GetBinLowEdge(bestCutBin_a), 0, srootb_a->GetBinLowEdge(bestCutBin_a), lineMaxY);
  bestCutLine_a->SetLineColor(kMagenta);
  bestCutLine_a->SetLineWidth(2);
  bestCutLine_a->Draw("same");
  TLine * bestCutLine_b = new TLine(srootb_b->GetBinLowEdge(bestCutBin_b), 0, srootb_b->GetBinLowEdge(bestCutBin_b), lineMaxY);
  bestCutLine_b->SetLineColor(kBlue);
  bestCutLine_b->SetLineWidth(2);
  bestCutLine_b->Draw("same");

  canv->SaveAs("roc_"+title+"_"+req+".pdf");

}

vector<TH1D*> preparePlots(TTree * ggTree, TTree * egTree, TTree * ffTree, TTree * gfTree, TTree * ttgjetsTree,
			   TH1D * ratio_ff_0, TH1D * ratio_ff_1, TH1D * ratio_ff_2,
			   TH1D * ratio_gf_0, TH1D * ratio_gf_1, TH1D * ratio_gf_2,
			   Float_t egScale, Float_t egScaleErr,
			   Float_t ffScale, Float_t ffScaleErr,
			   Float_t gfScale, Float_t gfScaleErr, Float_t ttgjetsScale,
			   TString variable, bool isAFloat,
			   Int_t nBinsX, Float_t bin_lo, Float_t bin_hi,
			   TString req,
			   Int_t rebin,
			   bool compareReweighting) {

  TH1D * gg = HistoFromTree(isAFloat, variable, ggTree, variable+"_gg_"+req, variable, nBinsX, bin_lo, bin_hi);
  TH1D * ewk = HistoFromTree(isAFloat, variable, egTree, variable+"_eg_"+req, variable, nBinsX, bin_lo, bin_hi);
  TH1D * ttg = SignalHistoFromTree(ttgjetsScale, isAFloat, variable, ttgjetsTree, variable+"_ttgjets_"+req, variable, nBinsX, bin_lo, bin_hi);
  TH1D * qcd = Reweight(isAFloat, gg, variable, ratio_ff_0, ratio_ff_1, ratio_ff_2, ffTree, 0, false, "ff_"+req, false);
  TH1D * qcd_gf = Reweight(isAFloat, gg, variable, ratio_gf_0, ratio_gf_1, ratio_gf_2, gfTree, 0, false, "gf_"+req, false);
 
  TH1D * ewk_noNorm = (TH1D*)ewk->Clone();
  ewk->Scale(egScale);
  for(int i = 0; i < ewk->GetNbinsX(); i++) {
    Float_t normerr = egScaleErr*(ewk_noNorm->GetBinContent(i+1));
    Float_t staterr = ewk->GetBinError(i+1);
    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);
    ewk->SetBinError(i+1, new_err);
  }

  TH1D * qcd_noNorm = (TH1D*)qcd->Clone();
  qcd->Scale(ffScale);
  for(int i = 0; i < qcd->GetNbinsX(); i++) {
    Float_t normerr = ffScaleErr*(qcd_noNorm->GetBinContent(i+1));
    Float_t staterr = qcd->GetBinError(i+1);
    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);
    qcd->SetBinError(i+1, new_err);
  }

  TH1D * qcd_gf_noNorm = (TH1D*)qcd_gf->Clone();
  qcd_gf->Scale(gfScale);
  for(int i = 0; i < qcd_gf->GetNbinsX(); i++) {
    Float_t normerr = gfScaleErr*(qcd_gf_noNorm->GetBinContent(i+1));
    Float_t staterr = qcd_gf->GetBinError(i+1);
    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);
    qcd_gf->SetBinError(i+1, new_err);
  }

  if(rebin >= 2) {
    gg->Rebin(rebin);
    ewk->Rebin(rebin);
    ttg->Rebin(rebin);
    qcd->Rebin(rebin);
    qcd_gf->Rebin(rebin);
  }

  if(compareReweighting) {
    TH1D * qcd_unrew = HistoFromTree(isAFloat, variable, ffTree, variable+"_ff_unrew_"+req, variable, nBinsX, bin_lo, bin_hi);
    if(rebin >= 2) qcd_unrew->Rebin(rebin);
    qcd_unrew->Scale(1./qcd_unrew->Integral());

    TH1D * ratio_compare = (TH1D*)qcd->Clone("ratio_compare_"+variable+"_ff_"+req);
    ratio_compare->Scale(1./ratio_compare->Integral());
    ratio_compare->Divide(qcd_unrew);

    float var;
    int var_int;
    if(isAFloat) ffTree->SetBranchAddress(variable, &var);
    else ffTree->SetBranchAddress(variable, &var_int);

    float varmax = -999.;
    for(int i = 0; i < ffTree->GetEntries(); i++) {
      ffTree->GetEntry(i);
      if(isAFloat) {
	if(var > varmax) varmax = var;
      }
      else {
	if(var_int > varmax) varmax = var_int;
      }
    }

    ffTree->ResetBranchAddresses();

    TCanvas * canv = new TCanvas("canv", "Plot", 10, 10, 2000, 2000);
    ratio_compare->GetXaxis()->SetRangeUser(0, varmax * 1.1);
    ratio_compare->Draw();
    canv->SaveAs("avg_weight_ff_"+variable+"_"+req+gifOrPdf);
    delete canv;
  }

  vector<TH1D*> vec;
  vec.push_back(gg);
  vec.push_back(ewk);
  vec.push_back(qcd);
  vec.push_back(qcd_gf);
  vec.push_back(ttg);

  return vec;

}

void prettyPlot(TH1D* gg_, TH1D* eg_, TH1D* ff_, TH1D* gf_, TH1D* ttgjets,
		bool useQCDSystematic, bool useTTGJets, bool useFF, bool useMCforQCD,
                TH1D* sig_a_, TH1D* sig_b_,
		TH1D * mc_qcd30to40, TH1D * mc_qcd40,
                TString intLumi, bool ff_failure,
                TString req, TString title, TString xaxisTitle, TString yaxisTitle,
                Float_t xmin, Float_t xmax, Float_t ymin, Float_t ymax,
                Float_t ratiomin, Float_t ratiomax,
                bool drawMC, bool drawLegend, bool drawPrelim) {

  TH1D * gg = (TH1D*)gg_->Clone("gg");
  TH1D * ewk = (TH1D*)eg_->Clone("ewk");
  TH1D * ttg = (TH1D*)ttgjets->Clone("ttg");

  TH1D * qcd;
  if(!useMCforQCD) {
    // use FF as data-driven QCD
    if(useFF) {
      qcd = (TH1D*)ff_->Clone("qcd");
      if(useQCDSystematic) {
	for(int i = 0; i < qcd->GetNbinsX(); i++) {
	  Double_t syserror = qcd->GetBinContent(i+1) - gf_->GetBinContent(i+1);
	  Double_t totalerr = qcd->GetBinError(i+1);
	  qcd->SetBinError(i+1, sqrt(totalerr*totalerr + syserror*syserror));
	}
      }
    }
    // use GF as data-driven QCD
    else {
      qcd = (TH1D*)gf_->Clone("qcd");
      if(useQCDSystematic) {
	for(int i = 0; i < qcd->GetNbinsX(); i++) {
	  Double_t syserror = qcd->GetBinContent(i+1) - ff_->GetBinContent(i+1);
	  Double_t totalerr = qcd->GetBinError(i+1);
	  qcd->SetBinError(i+1, sqrt(totalerr*totalerr + syserror*syserror));
	}
      }
    }
  }
  // use monte-carlo QCD
  else {
    qcd = (TH1D*)mc_qcd30to40->Clone("qcd");
    qcd->Add(mc_qcd40);
  } 

  TH1D * bkg = (TH1D*)ewk->Clone("bkg");
  if(useTTGJets) bkg->Add(ttg);
  bkg->Add(qcd);  
  TH1D * errors = (TH1D*)bkg->Clone("errors");

  ttg->Add(ewk);

  TH1D * sig_a;
  TH1D * sig_b;
  if(drawMC) {
    sig_a = (TH1D*)sig_a_->Clone("sig_a");
    sig_b= (TH1D*)sig_b_->Clone("sig_b");
  }

  if(drawMC) calculateROC(sig_a, sig_b, bkg, req, title);

  TLegend * leg = new TLegend(0.50, 0.65, 0.85, 0.85, NULL, "brNDC");
  leg->AddEntry(gg, "#gamma#gamma Candidate Sample", "LP");
  leg->AddEntry(errors, "Total Background Uncertainty", "F");
  if(!ff_failure) leg->AddEntry(bkg, "QCD", "F");
  else leg->AddEntry(bkg, "QCD (UN-NORMALIZED)", "F");
  if(useTTGJets) leg->AddEntry(ttg, "t#bar{t}#gamma+jets MC", "F");
  leg->AddEntry(ewk, "Electroweak", "F");
  leg->SetFillColor(0);
  leg->SetTextSize(0.028);

  TPaveText * prelim = new TPaveText(0.50, 0.42, 0.85, 0.62, "NDC");
  prelim->SetFillColor(0);
  prelim->SetFillStyle(0);
  prelim->SetLineColor(0);
  prelim->AddText("CMS Preliminary 2013");
  prelim->AddText(" ");
  prelim->AddText("#sqrt{s} = 8 TeV, #intL = "+intLumi+" fb^{-1}");
  prelim->AddText(req+" Requirement");

  gg->SetMarkerStyle(20); 
  gg->SetMarkerSize(1.5);

  errors->SetFillColor(kOrange+10);
  errors->SetFillStyle(3154);
  errors->SetMarkerSize(0);

  ewk->SetFillColor(8);
  ewk->SetMarkerSize(0);
  ewk->SetLineColor(1);
  
  ttg->SetFillColor(kOrange-3);
  ttg->SetMarkerSize(0);
  ttg->SetLineColor(1);

  bkg->SetFillColor(kGray);
  bkg->SetMarkerSize(0);
  bkg->SetLineColor(1);

  TCanvas * can = new TCanvas("can", "Plot", 10, 10, 2000, 2000);

  TPad * padhi = new TPad("padhi", "padhi", 0, 0.3, 1, 1);
  TPad * padlo = new TPad("padlo", "padlo", 0, 0, 1, 0.3);

  padhi->Draw();
  padlo->Draw();
  padhi->cd();

  padhi->SetLogy(false);
  padhi->SetTickx(true);
  padhi->SetTicky(true);
  //padhi->SetGridx(true);
  //padhi->SetGridy(true);
  padhi->SetBottomMargin(0);

  bkg->SetTitle(title);
  bkg->GetXaxis()->SetTitle(xaxisTitle);
  bkg->GetYaxis()->SetTitle(yaxisTitle);

  if(xmax > xmin) bkg->GetXaxis()->SetRangeUser(xmin, xmax);
  bkg->GetYaxis()->SetRangeUser(ymin, ymax);

  bkg->Draw("hist");
  if(useTTGJets) ttg->Draw("same hist");
  ewk->Draw("same hist");
  errors->Draw("same e2");
  gg->Draw("same e1");
  ewk->Draw("same axis");

  if(drawMC) {
    sig_a->SetLineColor(kMagenta);
    sig_a->SetLineWidth(3);
    leg->AddEntry(sig_a, "GGM #gamma#gamma (460_175)", "L");
    sig_a->Draw("same hist");
    
    sig_b->SetLineColor(kBlue);
    sig_b->SetLineWidth(3);
    leg->AddEntry(sig_b, "GGM #gamma#gamma (560_325)", "L");
    sig_b->Draw("same hist");
  }

  if(drawLegend) leg->Draw("same");
  if(drawPrelim && drawLegend) prelim->Draw("same");

  padlo->cd();
  padlo->SetTopMargin(0);
  padlo->SetBottomMargin(0.2);

  TH1D * ratio = (TH1D*)gg->Clone("ratio");
  ratio->Reset();
  ratio->SetTitle("Data / Background");
  for(int i = 0; i < ratio->GetNbinsX(); i++) {
    if(bkg->GetBinContent(i+1) == 0.) continue;
    ratio->SetBinContent(i+1, gg->GetBinContent(i+1) / bkg->GetBinContent(i+1));
    ratio->SetBinError(i+1, gg->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  TH1D * ratio_sys;
  ratio_sys = (TH1D*)bkg->Clone("ratio_sys");
  for(int i = 0; i < ratio_sys->GetNbinsX(); i++) {
    ratio_sys->SetBinContent(i+1, 1.);
    if(bkg->GetBinContent(i+1) == 0.) ratio_sys->SetBinError(i+1, 0.);
    else ratio_sys->SetBinError(i+1, ratio_sys->GetBinError(i+1) / bkg->GetBinContent(i+1));
  }

  if(xmax > xmin) ratio->GetXaxis()->SetRangeUser(xmin, xmax);

  ratio_sys->SetFillStyle(1001);
  ratio_sys->SetFillColor(kGray);
  ratio_sys->SetLineColor(kGray);
  ratio_sys->SetMarkerColor(kGray);

  ratio->GetXaxis()->SetTitle(xaxisTitle);
  ratio->GetXaxis()->SetLabelFont(63);
  ratio->GetXaxis()->SetLabelSize(48);
  ratio->GetXaxis()->SetTitleSize(0.12);
  ratio->GetXaxis()->SetTitleOffset(0.6);
  ratio->GetYaxis()->SetTitle("Data / Background");
  ratio->GetYaxis()->SetLabelFont(63);
  ratio->GetYaxis()->SetLabelSize(48);
  ratio->GetYaxis()->SetTitleSize(0.08);
  ratio->GetYaxis()->SetTitleOffset(0.5);
  ratio->GetYaxis()->SetRangeUser(ratiomin, ratiomax);
  ratio->GetYaxis()->SetNdivisions(508);

  ratio->Draw("e1");
  ratio_sys->Draw("e2 same");
  ratio->Draw("e1 same");
  ratio->Draw("axis same");

  TLine * oneLine = new TLine(xmin, 1, xmax, 1);
  oneLine->SetLineStyle(2);
  oneLine->Draw();  

  padhi->cd();
  padhi->SetLogy(true);
  can->SaveAs(title+"_"+req+".pdf");

  delete can;

}

void prep_signal(TString req) {

  Double_t mst[29] = {110, 160, 185, 210, 235, 260, 285, 310, 335, 360, 385, 410, 460, 510, 560, 610, 660, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 1710, 2010, 5010};
  Double_t mBino[31] = {25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 375, 425, 475, 525, 575, 625, 675, 725, 825, 925, 1025, 1125, 1225, 1325, 1425, 1525, 1725, 2025};

  Double_t xbins[31];
  xbins[0] = 0;
  xbins[1] = 55;
  for(int i = 1; i < 29; i++) xbins[i+1] = (mst[i] + mst[i-1])/2.;
  xbins[30] = 6510;

  Double_t ybins[33];
  ybins[0] = 0;
  ybins[1] = 12.5;
  for(int i = 1; i < 31; i++) ybins[i+1] = (mBino[i] + mBino[i-1])/2.;
  ybins[32] = 2175;

  char code[100];
  int index1, index2;

  TH2D * h_acc = new TH2D("acc_"+req, "acc_"+req, 30, xbins, 32, ybins);
  TH2D * h_contam = new TH2D("contam_"+req, "contam_"+req, 30, xbins, 32, ybins);

  TFile * out = new TFile("signal_"+req+".root", "RECREATE");

  for(int i = 0; i < 899; i++) {

    index1 = mst[int(i)/31];
    index2 = mBino[int(i)%31];
    sprintf(code, "_mst_%d_m1_%d", index1, index2);
    TString code_t = code;

// DURP
if(index1 != 460 || index2 != 175) continue;

    TFile * f = new TFile("../acceptance/signal_contamination"+code_t+".root", "READ");
    if(f->IsZombie()) {
      f->Close();
      continue;
    }
    
    TTree * ggTree = (TTree*)f->Get("gg_"+req+"_EvtTree"+code_t);
    TTree * ffTree = (TTree*)f->Get("ff_"+req+"_EvtTree"+code_t);
    if(ggTree->GetEntries() > 0) {
      TH1D * gg = (TH1D*)SignalHistoFromTree(1.0, true, "pfMET", ggTree, "met_gg_"+req+code_t, "met_gg_"+req+code_t, 400, 0., 2000.);
      TH1D * ff = (TH1D*)SignalHistoFromTree(1.0, true, "pfMET", ffTree, "met_ff_"+req+code_t, "met_ff_"+req+code_t, 400, 0., 2000.);

      out->cd();
      gg->Write();
      ff->Write();
    }
    else {
      f->Close();
      continue;
    }

    //TH1D * h = (TH1D*)f->Get("ngen"+code_t);
    //double n = h->GetBinContent(1);
    double n = 15000.;

    TH1D * met_gg = (TH1D*)f->Get("met"+code_t+"_gg_"+req);
    double acceptance = met_gg->Integral();
    if(n > 0) h_acc->Fill(index1, index2, acceptance / n);

    TH1D * met_ff = (TH1D*)f->Get("met"+code_t+"_ff_"+req);
    double contamination = met_ff->Integral();
    if(n > 0) h_contam->Fill(index1, index2, contamination / acceptance);

    f->Close();
  }

  TCanvas * can = new TCanvas("canvas", "Plot", 10, 10, 2000, 2000);

  h_acc->GetXaxis()->SetTitle("#tilde{t} mass (GeV/c^{2})");
  h_acc->GetXaxis()->SetRangeUser(0, 1600);
  h_acc->GetXaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetTitle("Bino} mass (GeV/c^{2})");
  h_acc->GetYaxis()->SetTitleOffset(1.3);
  h_acc->GetYaxis()->SetLabelSize(0.03);
  h_acc->GetYaxis()->SetRangeUser(0, 1600);
  h_acc->GetZaxis()->SetLabelSize(0.02);
  h_acc->Draw("colz");
  can->SaveAs("acceptance_"+req+".pdf");
  
  h_contam->GetXaxis()->SetTitle("#tilde{t} mass (GeV/c^{2})");
  h_contam->GetXaxis()->SetRangeUser(0, 1600);
  h_contam->GetXaxis()->SetLabelSize(0.03);
  h_contam->GetYaxis()->SetTitle("Bino mass (GeV/c^{2})");
  h_contam->GetYaxis()->SetTitleOffset(1.3);
  h_contam->GetYaxis()->SetLabelSize(0.03);
  h_contam->GetYaxis()->SetRangeUser(0, 1600);
  h_contam->GetZaxis()->SetLabelSize(0.02);
  h_contam->Draw("colz");
  can->SaveAs("contamination_"+req+".pdf");

  out->Write();
  out->Close();

}

void plotReducedChi2(vector<TH1D*> gg, vector<TH1D*> gf, vector<TH1D*> ff,
		     TH2D*& gf_chi2,
		     TH2D*& ff_chi2,
		     Int_t binx) {

  for(unsigned int i = 0; i < gg.size(); i++) {
    if(gg[i]->Integral() >= 1.) gg[i]->Scale(1./gg[i]->Integral());
  }
  for(unsigned int i = 0; i < gf.size(); i++) {
    if(gf[i]->Integral() >= 1.) gf[i]->Scale(1./gf[i]->Integral());
  }
  for(unsigned int i = 0; i < ff.size(); i++) {
    if(ff[i]->Integral() >= 1.) ff[i]->Scale(1./ff[i]->Integral());
  }

  for(int i = 0; i < 10; i++) {

    Float_t chi2 = 0.;
    Int_t nBins = 0;
    
    for(int j = 0; j < gg[i]->GetNbinsX(); j++) {
      Float_t val_num = gg[i]->GetBinContent(j+1) - gf[i]->GetBinContent(j+1);
      Float_t val_den = gg[i]->GetBinError(j+1)*gg[i]->GetBinError(j+1) + gf[i]->GetBinError(j+1)*gf[i]->GetBinError(j+1);

      if(val_den == 0.) continue;

      chi2 += val_num * val_num / val_den;
      nBins++;
    }
    chi2 /= nBins;
    gf_chi2->SetBinContent(gf_chi2->FindBin(binx, i), chi2);

    chi2 = 0.;
    nBins = 0;
    
    for(int j = 0; j < gg[i]->GetNbinsX(); j++) {
      Float_t val_num = gg[i]->GetBinContent(j+1) - ff[i]->GetBinContent(j+1);
      Float_t val_den = gg[i]->GetBinError(j+1)*gg[i]->GetBinError(j+1) + ff[i]->GetBinError(j+1)*ff[i]->GetBinError(j+1);

      if(val_den == 0.) continue;

      chi2 += val_num * val_num / val_den;
      nBins++;
    }
    chi2 /= nBins;
    ff_chi2->SetBinContent(gf_chi2->FindBin(binx, i), chi2);
  }

  // now fill bin 10, chi2 across all variables
  Float_t chi2_all = 0.;
  Float_t nBins_all = 0;
  for(unsigned int i = 0; i < gg.size(); i++) {
    for(int j = 0; j < gg[i]->GetNbinsX(); j++) {
      Float_t val_num = gg[i]->GetBinContent(j+1) - gf[i]->GetBinContent(j+1);
      Float_t val_den = gg[i]->GetBinError(j+1)*gg[i]->GetBinError(j+1) + gf[i]->GetBinError(j+1)*gf[i]->GetBinError(j+1);

      if(val_den == 0.) continue;

      chi2_all += val_num * val_num / val_den;
      nBins_all++;
    }
  }
  chi2_all /= nBins_all;
  gf_chi2->SetBinContent(gf_chi2->FindBin(binx, 10), chi2_all);

  chi2_all = 0.;
  nBins_all = 0;
  for(unsigned int i = 0; i < gg.size(); i++) {
    for(int j = 0; j < gg[i]->GetNbinsX(); j++) {
      Float_t val_num = gg[i]->GetBinContent(j+1) - ff[i]->GetBinContent(j+1);
      Float_t val_den = gg[i]->GetBinError(j+1)*gg[i]->GetBinError(j+1) + ff[i]->GetBinError(j+1)*ff[i]->GetBinError(j+1);
    
      if(val_den == 0.) continue;
  
      chi2_all += val_num * val_num / val_den;
      nBins_all++;
    }
  }
  chi2_all /= nBins_all;
  ff_chi2->SetBinContent(ff_chi2->FindBin(binx, 10), chi2_all);

}

void analyzeReweighting(TTree * ggtree, TTree * fftree, TTree * gftree) {

  TString fVars[7] = {"leadPhotonEt", "trailPhotonEt", "invmass", "HT", "HT_jets", "diEMpT", "diJetPt"};
  TString iVars[2] = {"Njets", "Nbtags"};

  TH2D * gf_chi2 = new TH2D("gf_chi2", "gf_chi2;Reweighting Variable;#chi^{2} Evaluatation Variable", 10, 0, 10, 11, 0, 11);
  TH2D * ff_chi2 = new TH2D("ff_chi2", "ff_chi2;Reweighting Variable;#chi^{2} Evaluatation Variable", 10, 0, 10, 11, 0, 11);

  vector<TH1D*> gg;
  for(int i = 0; i < 7; i++) gg.push_back((TH1D*)HistoFromTree(true, fVars[i], ggtree, "gg_"+fVars[i]+"_anaRew", "gg_"+fVars[i]+"_anaRew", 200, 0., 2000., true));
  gg.push_back((TH1D*)HistoFromTree(true, "pfMET", ggtree, "gg_pfMET_anaRew", "gg_pfMET_anaRew", 10, 0., 50., true));
  for(int i = 0; i < 2; i++) gg.push_back((TH1D*)HistoFromTree(false, iVars[i], ggtree, "gg_"+iVars[i]+"_anaRew", "gg_"+iVars[i]+"_anaRew", 20, 0., 20., true));

  vector<TH1D*> gf;
  for(int i = 0; i < 7; i++) gf.push_back((TH1D*)HistoFromTree(true, fVars[i], gftree, "gf_"+fVars[i]+"_anaRew", "gf_"+fVars[i]+"_anaRew", 200, 0., 2000., true));
  gf.push_back((TH1D*)HistoFromTree(true, "pfMET", gftree, "gf_pfMET_anaRew", "gf_pfMET_anaRew", 10, 0., 50., true));
  for(int i = 0; i < 2; i++) gf.push_back((TH1D*)HistoFromTree(false, iVars[i], gftree, "gf_"+iVars[i]+"_anaRew", "gf_"+iVars[i]+"_anaRew", 20, 0., 20., true));

  vector<TH1D*> ff;
  for(int i = 0; i < 7; i++) ff.push_back((TH1D*)HistoFromTree(true, fVars[i], fftree, "ff_"+fVars[i]+"_anaRew", "ff_"+fVars[i]+"_anaRew", 200, 0., 2000., true));
  ff.push_back((TH1D*)HistoFromTree(true, "pfMET", fftree, "ff_pfMET_anaRew", "ff_pfMET_anaRew", 10, 0., 50., true));
  for(int i = 0; i < 2; i++) ff.push_back((TH1D*)HistoFromTree(false, iVars[i], fftree, "ff_"+iVars[i]+"_anaRew", "ff_"+iVars[i]+"_anaRew", 20, 0., 20., true));

  plotReducedChi2(gg, gf, ff,
		  gf_chi2,
		  ff_chi2,
		  0);

  const Int_t nBins_photonEt = 16;
  Double_t bins_photonEt[nBins_photonEt+1] = {25, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 250, 300, 400};
  TH1D * gf_weight_leadPhotonEt = GetAlternativeWeights(ggtree, gftree, "leadPhotonEt", true, nBins_photonEt, bins_photonEt, "gf");
  TH1D * ff_weight_leadPhotonEt = GetAlternativeWeights(ggtree, fftree, "leadPhotonEt", true, nBins_photonEt, bins_photonEt, "ff");
  
  vector<TH1D*> gf_leadPhotonEt;
  for(int i = 0; i < 7; i++) gf_leadPhotonEt.push_back((TH1D*)AlternativeReweight(fVars[i], true, "leadPhotonEt", true, gg[i], gf_weight_leadPhotonEt, gftree, 0, false, "gf_"+fVars[i]+"_rewLeadPhotonEt", false));
  gf_leadPhotonEt.push_back((TH1D*)AlternativeReweight("pfMET", true, "leadPhotonEt", true, gg[7], gf_weight_leadPhotonEt, gftree, 0, false, "gf_pfMET_rewLeadPhotonEt", false));
  for(int i = 0; i < 2; i++) gf_leadPhotonEt.push_back((TH1D*)AlternativeReweight(iVars[i], false, "leadPhotonEt", true, gg[i+8], gf_weight_leadPhotonEt, gftree, 0, false, "gf_"+iVars[i]+"_rewLeadPhotonEt", false));

  vector<TH1D*> ff_leadPhotonEt;
  for(int i = 0; i < 7; i++) ff_leadPhotonEt.push_back((TH1D*)AlternativeReweight(fVars[i], true, "leadPhotonEt", true, gg[i], ff_weight_leadPhotonEt, fftree, 0, false, "ff_"+fVars[i]+"_rewLeadPhotonEt", false));
  ff_leadPhotonEt.push_back((TH1D*)AlternativeReweight("pfMET", true, "leadPhotonEt", true, gg[7], ff_weight_leadPhotonEt, fftree, 0, false, "ff_pfMET_rewLeadPhotonEt", false));
  for(int i = 0; i < 2; i++) ff_leadPhotonEt.push_back((TH1D*)AlternativeReweight(iVars[i], false, "leadPhotonEt", true, gg[i+8], ff_weight_leadPhotonEt, fftree, 0, false, "ff_"+iVars[i]+"_rewLeadPhotonEt", false));

  plotReducedChi2(gg, gf_leadPhotonEt, ff_leadPhotonEt,
		  gf_chi2,
		  ff_chi2,
		  1);

  TH1D * gf_weight_trailPhotonEt = GetAlternativeWeights(ggtree, gftree, "trailPhotonEt", true, nBins_photonEt, bins_photonEt, "gf");
  TH1D * ff_weight_trailPhotonEt = GetAlternativeWeights(ggtree, fftree, "trailPhotonEt", true, nBins_photonEt, bins_photonEt, "ff");

  vector<TH1D*> gf_trailPhotonEt;
  for(int i = 0; i < 7; i++) gf_trailPhotonEt.push_back((TH1D*)AlternativeReweight(fVars[i], true, "trailPhotonEt", true, gg[i], gf_weight_trailPhotonEt, gftree, 0, false, "gf_"+fVars[i]+"_rewTrailPhotonEt", false));
  gf_trailPhotonEt.push_back((TH1D*)AlternativeReweight("pfMET", true, "trailPhotonEt", true, gg[7], gf_weight_trailPhotonEt, gftree, 0, false, "gf_pfMET_rewTrailPhotonEt", false));
  for(int i = 0; i < 2; i++) gf_trailPhotonEt.push_back((TH1D*)AlternativeReweight(iVars[i], false, "trailPhotonEt", true, gg[i+8], gf_weight_trailPhotonEt, gftree, 0, false, "gf_"+iVars[i]+"_rewTrailPhotonEt", false));

  vector<TH1D*> ff_trailPhotonEt;
  for(int i = 0; i < 7; i++) ff_trailPhotonEt.push_back((TH1D*)AlternativeReweight(fVars[i], true, "trailPhotonEt", true, gg[i], ff_weight_trailPhotonEt, fftree, 0, false, "ff_"+fVars[i]+"_rewTrailPhotonEt", false));
  ff_trailPhotonEt.push_back((TH1D*)AlternativeReweight("pfMET", true, "trailPhotonEt", true, gg[7], ff_weight_trailPhotonEt, fftree, 0, false, "ff_pfMET_rewTrailPhotonEt", false));
  for(int i = 0; i < 2; i++) ff_trailPhotonEt.push_back((TH1D*)AlternativeReweight(iVars[i], false, "trailPhotonEt", true, gg[i+8], ff_weight_trailPhotonEt, fftree, 0, false, "ff_"+iVars[i]+"_rewTrailPhotonEt", false));

    plotReducedChi2(gg, gf_trailPhotonEt, ff_trailPhotonEt,
		    gf_chi2,
		    ff_chi2,
		    2);

  const Int_t nBins_invmass = 26;
  Double_t bins_invmass[nBins_invmass+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 500, 600, 800};
  TH1D * gf_weight_invmass = GetAlternativeWeights(ggtree, gftree, "invmass", true, nBins_invmass, bins_invmass, "gf");
  TH1D * ff_weight_invmass = GetAlternativeWeights(ggtree, fftree, "invmass", true, nBins_invmass, bins_invmass, "ff");

  vector<TH1D*> gf_invmass;
  for(int i = 0; i < 7; i++) gf_invmass.push_back((TH1D*)AlternativeReweight(fVars[i], true, "invmass", true, gg[i], gf_weight_invmass, gftree, 0, false, "gf_"+fVars[i]+"_rewInvmass", false));
  gf_invmass.push_back((TH1D*)AlternativeReweight("pfMET", true, "invmass", true, gg[7], gf_weight_invmass, gftree, 0, false, "gf_pfMET_rewInvmass", false));
  for(int i = 0; i < 2; i++) gf_invmass.push_back((TH1D*)AlternativeReweight(iVars[i], false, "invmass", true, gg[i+8], gf_weight_invmass, gftree, 0, false, "gf_"+iVars[i]+"_rewInvmass", false));

  vector<TH1D*> ff_invmass;
  for(int i = 0; i < 7; i++) ff_invmass.push_back((TH1D*)AlternativeReweight(fVars[i], true, "invmass", true, gg[i], ff_weight_invmass, fftree, 0, false, "ff_"+fVars[i]+"_rewInvmass", false));
  ff_invmass.push_back((TH1D*)AlternativeReweight("pfMET", true, "invmass", true, gg[7], ff_weight_invmass, fftree, 0, false, "ff_pfMET_rewInvmass", false));
  for(int i = 0; i < 2; i++) ff_invmass.push_back((TH1D*)AlternativeReweight(iVars[i], false, "invmass", true, gg[i+8], ff_weight_invmass, fftree, 0, false, "ff_"+iVars[i]+"_rewInvmass", false));

  plotReducedChi2(gg, gf_invmass, ff_invmass,
		  gf_chi2,
		  ff_chi2,
		  3);

  const Int_t nBins_HT = 32;
  Double_t bins_HT[nBins_HT+1] = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 160, 180, 200, 220, 240, 260, 280, 300, 325, 350, 375, 400, 450, 500, 550, 600, 800, 1100};
  TH1D * gf_weight_HT = GetAlternativeWeights(ggtree, gftree, "HT", true, nBins_HT, bins_HT, "gf");
  TH1D * ff_weight_HT = GetAlternativeWeights(ggtree, fftree, "HT", true, nBins_HT, bins_HT, "ff");

  vector<TH1D*> gf_HT;
  for(int i = 0; i < 7; i++) gf_HT.push_back((TH1D*)AlternativeReweight(fVars[i], true, "HT", true, gg[i], gf_weight_HT, gftree, 0, false, "gf_"+fVars[i]+"_rewHT", false));
  gf_HT.push_back((TH1D*)AlternativeReweight("pfMET", true, "HT", true, gg[7], gf_weight_HT, gftree, 0, false, "gf_pfMET_rewHT", false));
  for(int i = 0; i < 2; i++) gf_HT.push_back((TH1D*)AlternativeReweight(iVars[i], false, "HT", true, gg[i+8], gf_weight_HT, gftree, 0, false, "gf_"+iVars[i]+"_rewHT", false));

  vector<TH1D*> ff_HT;
  for(int i = 0; i < 7; i++) ff_HT.push_back((TH1D*)AlternativeReweight(fVars[i], true, "HT", true, gg[i], ff_weight_HT, fftree, 0, false, "ff_"+fVars[i]+"_rewHT", false));
  ff_HT.push_back((TH1D*)AlternativeReweight("pfMET", true, "HT", true, gg[7], ff_weight_HT, fftree, 0, false, "ff_pfMET_rewHT", false));
  for(int i = 0; i < 2; i++) ff_HT.push_back((TH1D*)AlternativeReweight(iVars[i], false, "HT", true, gg[i+8], ff_weight_HT, fftree, 0, false, "ff_"+iVars[i]+"_rewHT", false));

  plotReducedChi2(gg, gf_HT, ff_HT,
		  gf_chi2,
		  ff_chi2,
		  4);

  TH1D * gf_weight_HT_jets = GetAlternativeWeights(ggtree, gftree, "HT_jets", true, nBins_HT, bins_HT, "gf");
  TH1D * ff_weight_HT_jets = GetAlternativeWeights(ggtree, fftree, "HT_jets", true, nBins_HT, bins_HT, "ff");
  
  vector<TH1D*> gf_HT_jets;
  for(int i = 0; i < 7; i++) gf_HT_jets.push_back((TH1D*)AlternativeReweight(fVars[i], true, "HT_jets", true, gg[i], gf_weight_HT_jets, gftree, 0, false, "gf_"+fVars[i]+"_rewHT_jets", false));
  gf_HT_jets.push_back((TH1D*)AlternativeReweight("pfMET", true, "HT_jets", true, gg[7], gf_weight_HT_jets, gftree, 0, false, "gf_pfMET_rewHT_jets", false));
  for(int i = 0; i < 2; i++) gf_HT_jets.push_back((TH1D*)AlternativeReweight(iVars[i], false, "HT_jets", true, gg[i+8], gf_weight_HT_jets, gftree, 0, false, "gf_"+iVars[i]+"_rewHT_jets", false));

  vector<TH1D*> ff_HT_jets;
  for(int i = 0; i < 7; i++) ff_HT_jets.push_back((TH1D*)AlternativeReweight(fVars[i], true, "HT_jets", true, gg[i], ff_weight_HT_jets, fftree, 0, false, "ff_"+fVars[i]+"_rewHT_jets", false));
  ff_HT_jets.push_back((TH1D*)AlternativeReweight("pfMET", true, "HT_jets", true, gg[7], ff_weight_HT_jets, fftree, 0, false, "ff_pfMET_rewHT_jets", false));
  for(int i = 0; i < 2; i++) ff_HT_jets.push_back((TH1D*)AlternativeReweight(iVars[i], false, "HT_jets", true, gg[i+8], ff_weight_HT_jets, fftree, 0, false, "ff_"+iVars[i]+"_rewHT_jets", false));

  plotReducedChi2(gg, gf_HT_jets, ff_HT_jets,
		  gf_chi2,
		  ff_chi2,
		  5);

  const Int_t nBins_diEMpT = 30;
  Double_t bins_diEMpT[nBins_diEMpT+1] = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 110, 120, 130, 140, 150, 200, 300, 400, 600, 1000};
  TH1D * gf_weight_diEMpT = GetAlternativeWeights(ggtree, gftree, "diEMpT", true, nBins_diEMpT, bins_diEMpT, "gf");
  TH1D * ff_weight_diEMpT = GetAlternativeWeights(ggtree, fftree, "diEMpT", true, nBins_diEMpT, bins_diEMpT, "ff");

  vector<TH1D*> gf_diEMpT;
  for(int i = 0; i < 7; i++) gf_diEMpT.push_back((TH1D*)AlternativeReweight(fVars[i], true, "diEMpT", true, gg[i], gf_weight_diEMpT, gftree, 0, false, "gf_"+fVars[i]+"_rewdiEMpT", false));
  gf_diEMpT.push_back((TH1D*)AlternativeReweight("pfMET", true, "diEMpT", true, gg[7], gf_weight_diEMpT, gftree, 0, false, "gf_pfMET_rewdiEMpT", false));
  for(int i = 0; i < 2; i++) gf_diEMpT.push_back((TH1D*)AlternativeReweight(iVars[i], false, "diEMpT", true, gg[i+8], gf_weight_diEMpT, gftree, 0, false, "gf_"+iVars[i]+"_rewdiEMpT", false));

  vector<TH1D*> ff_diEMpT;
  for(int i = 0; i < 7; i++) ff_diEMpT.push_back((TH1D*)AlternativeReweight(fVars[i], true, "diEMpT", true, gg[i], ff_weight_diEMpT, fftree, 0, false, "ff_"+fVars[i]+"_rewdiEMpT", false));
  ff_diEMpT.push_back((TH1D*)AlternativeReweight("pfMET", true, "diEMpT", true, gg[7], ff_weight_diEMpT, fftree, 0, false, "ff_pfMET_rewdiEMpT", false));
  for(int i = 0; i < 2; i++) ff_diEMpT.push_back((TH1D*)AlternativeReweight(iVars[i], false, "diEMpT", true, gg[i+8], ff_weight_diEMpT, fftree, 0, false, "ff_"+iVars[i]+"_rewdiEMpT", false));

  plotReducedChi2(gg, gf_diEMpT, ff_diEMpT,
		  gf_chi2,
		  ff_chi2,
		  6);

  TH1D * gf_weight_diJetPt = GetAlternativeWeights(ggtree, gftree, "diJetPt", true, nBins_diEMpT, bins_diEMpT, "gf");
  TH1D * ff_weight_diJetPt = GetAlternativeWeights(ggtree, fftree, "diJetPt", true, nBins_diEMpT, bins_diEMpT, "ff");

  vector<TH1D*> gf_diJetPt;
  for(int i = 0; i < 7; i++) gf_diJetPt.push_back((TH1D*)AlternativeReweight(fVars[i], true, "diJetPt", true, gg[i], gf_weight_diJetPt, gftree, 0, false, "gf_"+fVars[i]+"_rewdiJetPt", false));
  gf_diJetPt.push_back((TH1D*)AlternativeReweight("pfMET", true, "diJetPt", true, gg[7], gf_weight_diJetPt, gftree, 0, false, "gf_pfMET_rewdiJetPt", false));
  for(int i = 0; i < 2; i++) gf_diJetPt.push_back((TH1D*)AlternativeReweight(iVars[i], false, "diJetPt", true, gg[i+8], gf_weight_diJetPt, gftree, 0, false, "gf_"+iVars[i]+"_rewdiJetPt", false));

  vector<TH1D*> ff_diJetPt;
  for(int i = 0; i < 7; i++) ff_diJetPt.push_back((TH1D*)AlternativeReweight(fVars[i], true, "diJetPt", true, gg[i], ff_weight_diJetPt, fftree, 0, false, "ff_"+fVars[i]+"_rewdiJetPt", false));
  ff_diJetPt.push_back((TH1D*)AlternativeReweight("pfMET", true, "diJetPt", true, gg[7], ff_weight_diJetPt, fftree, 0, false, "ff_pfMET_rewdiJetPt", false));
  for(int i = 0; i < 2; i++) ff_diJetPt.push_back((TH1D*)AlternativeReweight(iVars[i], false, "diJetPt", true, gg[i+8], ff_weight_diJetPt, fftree, 0, false, "ff_"+iVars[i]+"_rewdiJetPt", false));

  plotReducedChi2(gg, gf_diJetPt, ff_diJetPt,
		  gf_chi2,
		  ff_chi2,
		  7);

  const Int_t nBins_Njets = 9;
  Double_t bins_Njets[nBins_Njets+1] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  TH1D * gf_weight_Njets = GetAlternativeWeights(ggtree, gftree, "Njets", false, nBins_Njets, bins_Njets, "gf");
  TH1D * ff_weight_Njets = GetAlternativeWeights(ggtree, fftree, "Njets", false, nBins_Njets, bins_Njets, "ff");

  vector<TH1D*> gf_Njets;
  for(int i = 0; i < 7; i++) gf_Njets.push_back((TH1D*)AlternativeReweight(fVars[i], true, "Njets", false, gg[i], gf_weight_Njets, gftree, 0, false, "gf_"+fVars[i]+"_rewNjets", false));
  gf_Njets.push_back((TH1D*)AlternativeReweight("pfMET", true, "Njets", false, gg[7], gf_weight_Njets, gftree, 0, false, "gf_pfMET_rewNjets", false));
  for(int i = 0; i < 2; i++) gf_Njets.push_back((TH1D*)AlternativeReweight(iVars[i], false, "Njets", false, gg[i+8], gf_weight_Njets, gftree, 0, false, "gf_"+iVars[i]+"_rewNjets", false));

  vector<TH1D*> ff_Njets;
  for(int i = 0; i < 7; i++) ff_Njets.push_back((TH1D*)AlternativeReweight(fVars[i], true, "Njets", false, gg[i], ff_weight_Njets, fftree, 0, false, "ff_"+fVars[i]+"_rewNjets", false));
  ff_Njets.push_back((TH1D*)AlternativeReweight("pfMET", true, "Njets", false, gg[7], ff_weight_Njets, fftree, 0, false, "ff_pfMET_rewNjets", false));
  for(int i = 0; i < 2; i++) ff_Njets.push_back((TH1D*)AlternativeReweight(iVars[i], false, "Njets", false, gg[i+8], ff_weight_Njets, fftree, 0, false, "ff_"+iVars[i]+"_rewNjets", false));

  plotReducedChi2(gg, gf_Njets, ff_Njets,
		  gf_chi2,
		  ff_chi2,
		  8);

  const Int_t nBins_Nbtags = 4;
  Double_t bins_Nbtags[nBins_Nbtags+1] = {0, 1, 2, 3, 4};
  TH1D * gf_weight_Nbtags = GetAlternativeWeights(ggtree, gftree, "Nbtags", false, nBins_Nbtags, bins_Nbtags, "gf");
  TH1D * ff_weight_Nbtags = GetAlternativeWeights(ggtree, fftree, "Nbtags", false, nBins_Nbtags, bins_Nbtags, "ff");

  vector<TH1D*> gf_Nbtags;
  for(int i = 0; i < 7; i++) gf_Nbtags.push_back((TH1D*)AlternativeReweight(fVars[i], true, "Nbtags", false, gg[i], gf_weight_Nbtags, gftree, 0, false, "gf_"+fVars[i]+"_rewNbtags", false));
  gf_Nbtags.push_back((TH1D*)AlternativeReweight("pfMET", true, "Nbtags", false, gg[7], gf_weight_Nbtags, gftree, 0, false, "gf_pfMET_rewNbtags", false));
  for(int i = 0; i < 2; i++) gf_Nbtags.push_back((TH1D*)AlternativeReweight(iVars[i], false, "Nbtags", false, gg[i+8], gf_weight_Nbtags, gftree, 0, false, "gf_"+iVars[i]+"_rewNbtags", false));

  vector<TH1D*> ff_Nbtags;
  for(int i = 0; i < 7; i++) ff_Nbtags.push_back((TH1D*)AlternativeReweight(fVars[i], true, "Nbtags", false, gg[i], ff_weight_Nbtags, fftree, 0, false, "ff_"+fVars[i]+"_rewNbtags", false));
  ff_Nbtags.push_back((TH1D*)AlternativeReweight("pfMET", true, "Nbtags", false, gg[7], ff_weight_Nbtags, fftree, 0, false, "ff_pfMET_rewNbtags", false));
  for(int i = 0; i < 2; i++) ff_Nbtags.push_back((TH1D*)AlternativeReweight(iVars[i], false, "Nbtags", false, gg[i+8], ff_weight_Nbtags, fftree, 0, false, "ff_"+iVars[i]+"_rewNbtags", false));

  plotReducedChi2(gg, gf_Nbtags, ff_Nbtags,
		  gf_chi2,
		  ff_chi2,
		  9);

  ggtree->ResetBranchAddresses();
  gftree->ResetBranchAddresses();
  fftree->ResetBranchAddresses();

  TString xlabels[10] = {"No reweighting", "leadPhotonEt", "trailPhotonEt", "invmass", "HT", "HT_jets", "di-EM Pt", "di-Jet Pt", "nJets", "nBtags"};
  TString ylabels[11] = {"leadPhotonEt", "trailPhotonEt", "invmass", "HT", "HT_jets", "di-EM Pt", "di-Jet Pt", "pfMET", "nJets", "nBtags", "all"};

  for(int i = 0; i < 10; i++) {
    gf_chi2->GetXaxis()->SetBinLabel(i+1, xlabels[i]);
    ff_chi2->GetXaxis()->SetBinLabel(i+1, xlabels[i]);
  }

  for(int i = 0; i < 11; i++) {
    gf_chi2->GetYaxis()->SetBinLabel(i+1, ylabels[i]);
    ff_chi2->GetYaxis()->SetBinLabel(i+1, ylabels[i]);
  }

  gf_chi2->GetXaxis()->SetLabelSize(0.035);
  gf_chi2->GetXaxis()->SetTitleSize(0.02);
  gf_chi2->GetXaxis()->SetTitleOffset(1.5);
  gf_chi2->GetYaxis()->SetLabelSize(0.022);
  gf_chi2->GetYaxis()->SetTitleSize(0.02);
  gf_chi2->GetYaxis()->SetTitleOffset(1.5);
  gf_chi2->GetZaxis()->SetLabelSize(0.02);
  
  ff_chi2->GetXaxis()->SetLabelSize(0.035);
  ff_chi2->GetXaxis()->SetTitleSize(0.02);
  ff_chi2->GetXaxis()->SetTitleOffset(1.5);
  ff_chi2->GetYaxis()->SetLabelSize(0.022);
  ff_chi2->GetYaxis()->SetTitleSize(0.02);
  ff_chi2->GetYaxis()->SetTitleOffset(1.5);
  ff_chi2->GetZaxis()->SetLabelSize(0.02);

  TCanvas * can = new TCanvas("can_reweight_gf", "can_reweight_gf", 10, 10, 2000, 2000);
  //can->SetLogz(true);
  gf_chi2->Draw("colz");
  gf_chi2->SetMarkerColor(kWhite);
  gf_chi2->Draw("text same");
  can->SaveAs("analyze_reweighting_gf.pdf");

  TCanvas * can2 = new TCanvas("can_reweight_ff", "can_reweight_ff", 10, 10, 2000, 2000);
  //can2->SetLogz(true);
  ff_chi2->Draw("colz");
  ff_chi2->SetMarkerColor(kWhite);
  ff_chi2->Draw("text same");
  can2->SaveAs("analyze_reweighting_ff.pdf");

  return;

}
