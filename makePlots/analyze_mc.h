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

#include "rootRoutines.h"

using namespace std;

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

TH1D * HistoFromTree(bool isAFloat, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t xlo, Double_t xhi, double metCut = -1.) {

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

    if(metCut > 0. && met >= metCut) continue;

    if(variable != "pfMET") {
      if(isAFloat) h->Fill(var);
      else h->Fill(var_int);
    }
    else h->Fill(met);

  }

  tree->ResetBranchAddresses();

  return h;
}

TH1D * HistoFromTree(bool isAFloat, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t* customBins, double metCut = -1.) {

  TH1D * h = new TH1D(name, title, nBins, customBins);
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

    if(metCut > 0. && met >= metCut) continue;

    if(variable != "pfMET") {
      if(isAFloat) h->Fill(var);
      else h->Fill(var_int);
    }
    else h->Fill(met);

  }

  tree->ResetBranchAddresses();

  return h;
}

TH1D * SignalHistoFromTree(Float_t scale, bool isAFloat, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t xlo, Double_t xhi, double metCut = -1.) {

  TH1D * h = new TH1D(name, title, nBins, xlo, xhi);
  h->Sumw2();

  Float_t var, met;
  Int_t var_int;
  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr, btagWeightUp, btagWeightDown;

  tree->SetBranchAddress("pfMET", &met);
  if(variable != "pfMET") {
    if(isAFloat) tree->SetBranchAddress(variable, &var);
    else tree->SetBranchAddress(variable, &var_int);
  }

  tree->SetBranchAddress("pileupWeight", &puWeight);
  tree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  tree->SetBranchAddress("btagWeight", &btagWeight);
  tree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  tree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  tree->SetBranchAddress("btagWeightDown", &btagWeightDown);

  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if(metCut > 0. && met >= metCut) continue;

    Float_t olderror = 0.;

    if(variable == "pfMET") var = met;

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

TH1D * SignalHistoFromTree(Float_t scale, bool isAFloat, TString variable, TTree * tree, TString name, TString title, Int_t nBins, Double_t* customBins, double metCut = -1.) {

  TH1D * h = new TH1D(name, title, nBins, customBins);
  h->Sumw2();

  Float_t var, met;
  Int_t var_int;
  Float_t puWeight, btagWeight;
  Float_t puWeightErr, btagWeightErr, btagWeightUp, btagWeightDown;

  tree->SetBranchAddress("pfMET", &met);
  if(variable != "pfMET") {
    if(isAFloat) tree->SetBranchAddress(variable, &var);
    else tree->SetBranchAddress(variable, &var_int);
  }

  tree->SetBranchAddress("pileupWeight", &puWeight);
  tree->SetBranchAddress("pileupWeightErr", &puWeightErr);
  tree->SetBranchAddress("btagWeight", &btagWeight);
  tree->SetBranchAddress("btagWeightErr", &btagWeightErr);
  tree->SetBranchAddress("btagWeightUp", &btagWeightUp);
  tree->SetBranchAddress("btagWeightDown", &btagWeightDown);
  for(int i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    if(metCut > 0. && met >= metCut) continue;

    Float_t olderror = 0.;

    if(variable == "pfMET") var = met;

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

class PlotMaker : public TObject {
  
  ClassDef(PlotMaker, 1);

 public:
  PlotMaker(Int_t lumi,
	    Float_t fakeRate, Float_t fakeRateErr,
	    TString requirement);
  virtual ~PlotMaker() { 

    delete ggTree;
    delete egTree;
    delete qcd30to40Tree;
    delete qcd40Tree;
    delete gjet20to40Tree;
    delete gjet40Tree;
    delete ttHadronicTree;
    delete ttSemiLepTree;
    delete ttgjetsTree;
    delete sigaTree;
    delete sigbTree;
    
  }

  void SetTrees(TTree * gg, TTree * eg,
		TTree * qcd30to40, TTree * qcd40,
		TTree * gjet20to40, TTree * gjet40,
		TTree * ttHadronic, TTree * ttSemiLep,
		TTree * ttgjets,
		TTree * sig_a, TTree * sig_b);

  void CreatePlot(TString variable, bool isAFloat,
		  Int_t nBinsX, Float_t bin_lo, Float_t bin_hi,
		  TString xaxisTitle, TString yaxisTitle,
		  Float_t xmin, Float_t xmax,
		  Float_t ymin, Float_t ymax,
		  Float_t ratiomin, Float_t ratiomax,
		  bool drawSignal, bool drawLegend, bool drawPrelim,
		  TFile*& out, double metCut);

  void CreatePlot(TString variable, bool isAFloat,
		  Int_t nBinsX, Double_t* customBins,
		  TString xaxisTitle,
		  Float_t xmin, Float_t xmax,
		  Float_t ymin, Float_t ymax,
		  Float_t ratiomin, Float_t ratiomax,
		  bool drawSignal, bool drawLegend, bool drawPrelim,
		  TFile*& out, double metCut);

 private:
  TTree * ggTree;
  TTree * egTree;
  
  TTree * qcd30to40Tree;
  TTree * qcd40Tree;
  TTree * gjet20to40Tree;
  TTree * gjet40Tree;
  TTree * ttHadronicTree;
  TTree * ttSemiLepTree;
  TTree * ttgjetsTree;
  TTree * sigaTree;
  TTree * sigbTree;

  Int_t intLumi_int;
  TString intLumi;
  Float_t egScale, egScaleErr;
  TString req;

};

PlotMaker::PlotMaker(Int_t lumi, Float_t fakeRate, Float_t fakeRateErr, TString requirement) :
  intLumi_int(lumi),
  egScale(fakeRate),
  egScaleErr(fakeRateErr),
  req(requirement)
{
  char buffer[50];
  sprintf(buffer, "%.3f", (float)intLumi_int / 1000.);
  intLumi = buffer;
}

void PlotMaker::SetTrees(TTree * gg, TTree * eg,
			 TTree * qcd30to40, TTree * qcd40,
			 TTree * gjet20to40, TTree * gjet40,
			 TTree * ttHadronic, TTree * ttSemiLep,
			 TTree * ttgjets,
			 TTree * sig_a, TTree * sig_b) {

  ggTree = gg;
  egTree = eg;
  
  qcd30to40Tree = qcd30to40;
  qcd40Tree = qcd40;
  gjet20to40Tree = gjet20to40;
  gjet40Tree = gjet40;
  ttHadronicTree = ttHadronic;
  ttSemiLepTree = ttSemiLep;
  ttgjetsTree = ttgjets;

  sigaTree = sig_a;
  sigbTree = sig_b;

}

void PlotMaker::CreatePlot(TString variable, bool isAFloat,
			   Int_t nBinsX, Float_t bin_lo, Float_t bin_hi,
			   TString xaxisTitle, TString yaxisTitle,
			   Float_t xmin, Float_t xmax,
			   Float_t ymin, Float_t ymax,
			   Float_t ratiomin, Float_t ratiomax,
			   bool drawSignal, bool drawLegend, bool drawPrelim,
			   TFile*& out, double metCut) {

  TH1D * gg = HistoFromTree(isAFloat, variable, ggTree, variable+"_gg_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * ewk = HistoFromTree(isAFloat, variable, egTree, variable+"_eg_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);

  TH1D * ewk_noNorm = (TH1D*)ewk->Clone();
  ewk->Scale(egScale);
  for(int i = 0; i < ewk->GetNbinsX(); i++) {
    Float_t normerr = egScaleErr*(ewk_noNorm->GetBinContent(i+1));
    Float_t staterr = ewk->GetBinError(i+1);
    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);
    ewk->SetBinError(i+1, new_err);
  }

  TH1D * qcd30to40 = SignalHistoFromTree(intLumi_int * 5.195E7 * 2.35E-4 * 1.019 * 1.019 / 6061407., isAFloat, variable, qcd30to40Tree, variable+"_qcd30to40_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * qcd40 = SignalHistoFromTree(intLumi_int * 5.195E7 * 0.002175 * 1.019 * 1.019 / 9782735., isAFloat, variable, qcd40Tree, variable+"_qcd40_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * qcd = (TH1D*)qcd30to40->Clone("qcd");
  qcd->Add(qcd40);

  TH1D * gjet20to40 = SignalHistoFromTree(intLumi_int * 81930.0 * 0.001835 * 1.019 * 1.019 / 5907942., isAFloat, variable, gjet20to40Tree, variable+"_gjet20to40_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * gjet40 = SignalHistoFromTree(intLumi_int * 8884.0 * 0.05387 * 1.019 * 1.019 / 5956149., isAFloat, variable, gjet40Tree, variable+"_gjet40_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * gjet = (TH1D*)gjet20to40->Clone("gjet");
  gjet->Add(gjet40);

  TH1D * ttHadronic = SignalHistoFromTree(intLumi_int * 53.4 * 1.019 * 1.019 / 10537444., isAFloat, variable, ttHadronicTree, variable+"_ttHadronic_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * ttSemiLep = SignalHistoFromTree(intLumi_int * 53.2 * 1.019 * 1.019 / 25424818., isAFloat, variable, ttSemiLepTree, variable+"_ttSemiLep_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  TH1D * ttbar = (TH1D*)ttHadronic->Clone("ttbar");
  ttbar->Add(ttSemiLep);

  TH1D * ttg = SignalHistoFromTree(intLumi_int * 1.019 * 1.019 * 14.0 / 1719954., isAFloat, variable, ttgjetsTree, variable+"_ttgjets_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);

  out->cd();
  gg->Write();
  ewk->Write();
  qcd->Write();
  gjet->Write();
  ttg->Write();
  ttbar->Write();

  TH1D * bkg = (TH1D*)qcd->Clone("bkg");

  bkg->Add(gjet);
  bkg->Add(ewk);
  bkg->Add(ttg);
  bkg->Add(ttbar);

  gjet->Add(ewk);
  gjet->Add(ttg);
  gjet->Add(ttbar);

  ewk->Add(ttg);
  ewk->Add(ttbar);
  
  ttg->Add(ttbar);

  TH1D * errors = (TH1D*)bkg->Clone("errors");

  TH1D * sig_a;
  TH1D * sig_b;
  if(drawSignal) {
    sig_a = SignalHistoFromTree(0.147492 * intLumi_int * 1.019 * 1.019 / 15000., isAFloat, variable, sigaTree, variable+"_a_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
    sig_b = SignalHistoFromTree(0.0399591 * intLumi_int * 1.019 * 1.019 / 15000., isAFloat, variable, sigbTree, variable+"_b_"+req, variable, nBinsX, bin_lo, bin_hi, metCut);
  }

  if(drawSignal) calculateROC(sig_a, sig_b, bkg, req, variable);

  TLegend * leg = new TLegend(0.50, 0.65, 0.85, 0.85, NULL, "brNDC");
  leg->AddEntry(gg, "#gamma#gamma Candidate Sample", "LP");
  leg->AddEntry(errors, "Total Background Uncertainty", "F");
  leg->AddEntry(bkg, "QCD", "F");
  leg->AddEntry(gjet, "#gamma + jets", "F");
  leg->AddEntry(ewk, "Electroweak", "F");
  leg->AddEntry(ttg, "t#bar{t}#gamma + jets", "F");
  leg->AddEntry(ttbar, "t#bar{t} + jets", "F");
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

  // new stack: qcd, gjet, ewk, ttg, ttbar
  bkg->SetFillColor(kGray);
  bkg->SetMarkerSize(0);
  bkg->SetLineColor(1);

  gjet->SetFillColor(kOrange-3);
  gjet->SetMarkerSize(0);
  gjet->SetLineColor(1);

  ewk->SetFillColor(kYellow);
  ewk->SetMarkerSize(0);
  ewk->SetLineColor(1);

  ttg->SetFillColor(8);
  ttg->SetMarkerSize(0);
  ttg->SetLineColor(1);

  ttbar->SetFillColor(kRed-3);
  ttbar->SetMarkerSize(0);
  ttbar->SetLineColor(1);

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

  bkg->SetTitle(variable);
  bkg->GetXaxis()->SetTitle(xaxisTitle);
  bkg->GetYaxis()->SetTitle(yaxisTitle);

  if(xmax > xmin) bkg->GetXaxis()->SetRangeUser(xmin, xmax);
  bkg->GetYaxis()->SetRangeUser(ymin, ymax);

  // new stack: qcd, gjet, ewk, ttg, ttbar
  bkg->Draw("hist");
  gjet->Draw("same hist");
  ewk->Draw("same hist");
  ttg->Draw("same hist");
  ttbar->Draw("same hist");
  errors->Draw("same e2");
  gg->Draw("same e1");
  bkg->Draw("same axis");

  if(drawSignal) {
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
  can->SaveAs(variable+"_"+req+".pdf");

  delete can;

}

void PlotMaker::CreatePlot(TString variable, bool isAFloat,
			   Int_t nBinsX, Double_t* customBins,
			   TString xaxisTitle,
			   Float_t xmin, Float_t xmax,
			   Float_t ymin, Float_t ymax,
			   Float_t ratiomin, Float_t ratiomax,
			   bool drawSignal, bool drawLegend, bool drawPrelim,
			   TFile*& out, double metCut) {

  TString yaxisTitle = "Number of Events / GeV";

  TH1D * gg = HistoFromTree(isAFloat, variable, ggTree, variable+"_gg_"+req, variable, nBinsX, customBins, metCut);
  gg = (TH1D*)DivideByBinWidth(gg);
  TH1D * ewk = HistoFromTree(isAFloat, variable, egTree, variable+"_eg_"+req, variable, nBinsX, customBins, metCut);
  
  TH1D * ewk_noNorm = (TH1D*)ewk->Clone();
  ewk->Scale(egScale);
  for(int i = 0; i < ewk->GetNbinsX(); i++) {
    Float_t normerr = egScaleErr*(ewk_noNorm->GetBinContent(i+1));
    Float_t staterr = ewk->GetBinError(i+1);
    Float_t new_err = sqrt(normerr*normerr + staterr*staterr);
    ewk->SetBinError(i+1, new_err);
  }
  ewk = (TH1D*)DivideByBinWidth(ewk);

  TH1D * qcd30to40 = SignalHistoFromTree(intLumi_int * 5.195E7 * 2.35E-4 * 1.019 * 1.019 / 6061407., isAFloat, variable, qcd30to40Tree, variable+"_qcd30to40_"+req, variable, nBinsX, customBins, metCut);
  TH1D * qcd40 = SignalHistoFromTree(intLumi_int * 5.195E7 * 0.002175 * 1.019 * 1.019 / 9782735., isAFloat, variable, qcd40Tree, variable+"_qcd40_"+req, variable, nBinsX, customBins, metCut);
  TH1D * qcd = (TH1D*)qcd30to40->Clone(variable+"_qcd_"+req);
  qcd->Add(qcd40);
  qcd = (TH1D*)DivideByBinWidth(qcd);

  TH1D * gjet20to40 = SignalHistoFromTree(intLumi_int * 81930.0 * 0.001835 * 1.019 * 1.019 / 5907942., isAFloat, variable, gjet20to40Tree, variable+"_gjet20to40_"+req, variable, nBinsX, customBins, metCut);
  TH1D * gjet40 = SignalHistoFromTree(intLumi_int * 8884.0 * 0.05387 * 1.019 * 1.019 / 5956149., isAFloat, variable, gjet40Tree, variable+"_gjet40_"+req, variable, nBinsX, customBins, metCut);
  TH1D * gjet = (TH1D*)gjet20to40->Clone(variable+"_gjet_"+req);
  gjet->Add(gjet40);
  gjet = (TH1D*)DivideByBinWidth(gjet);

  TH1D * ttHadronic = SignalHistoFromTree(intLumi_int * 53.4 * 1.019 * 1.019 / 10537444., isAFloat, variable, ttHadronicTree, variable+"_ttHadronic_"+req, variable, nBinsX, customBins, metCut);
  TH1D * ttSemiLep = SignalHistoFromTree(intLumi_int * 53.2 * 1.019 * 1.019 / 25424818., isAFloat, variable, ttSemiLepTree, variable+"_ttSemiLep_"+req, variable, nBinsX, customBins, metCut);
  TH1D * ttbar = (TH1D*)ttHadronic->Clone(variable+"_ttbar_"+req);
  ttbar->Add(ttSemiLep);
  ttbar = (TH1D*)DivideByBinWidth(ttbar);

  TH1D * ttg = SignalHistoFromTree(intLumi_int * 1.019 * 1.019 * 14.0 / 1719954., isAFloat, variable, ttgjetsTree, variable+"_ttgjets_"+req, variable, nBinsX, customBins, metCut);
  ttg = (TH1D*)DivideByBinWidth(ttg);

  out->cd();
  gg->Write();
  ewk->Write();
  qcd->Write();
  gjet->Write();
  ttg->Write();
  ttbar->Write();

  TH1D * bkg = (TH1D*)qcd->Clone("bkg");

  bkg->Add(gjet);
  bkg->Add(ewk);
  bkg->Add(ttg);
  bkg->Add(ttbar);

  gjet->Add(ewk);
  gjet->Add(ttg);
  gjet->Add(ttbar);

  ewk->Add(ttg);
  ewk->Add(ttbar);
  
  ttg->Add(ttbar);

  TH1D * errors = (TH1D*)bkg->Clone("errors");

  TH1D * sig_a;
  TH1D * sig_b;
  if(drawSignal) {
    sig_a = SignalHistoFromTree(0.147492 * intLumi_int * 1.019 * 1.019 / 15000., true, variable, sigaTree, variable+"_a_"+req, variable, nBinsX, customBins, metCut);
    sig_a = (TH1D*)DivideByBinWidth(sig_a);
    sig_b = SignalHistoFromTree(0.0399591 * intLumi_int * 1.019 * 1.019 / 15000., true, variable, sigbTree, variable+"_b_"+req, variable, nBinsX, customBins, metCut);
    sig_b = (TH1D*)DivideByBinWidth(sig_b);
  }

  if(drawSignal) calculateROC(sig_a, sig_b, bkg, req, variable);

  TLegend * leg = new TLegend(0.50, 0.65, 0.85, 0.85, NULL, "brNDC");
  leg->AddEntry(gg, "#gamma#gamma Candidate Sample", "LP");
  leg->AddEntry(errors, "Total Background Uncertainty", "F");
  leg->AddEntry(bkg, "QCD", "F");
  leg->AddEntry(gjet, "#gamma + jets", "F");
  leg->AddEntry(ewk, "Electroweak", "F");
  leg->AddEntry(ttg, "t#bar{t}#gamma + jets", "F");
  leg->AddEntry(ttbar, "t#bar{t} + jets", "F");
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

  // new stack: qcd, gjet, ewk, ttg, ttbar
  bkg->SetFillColor(kGray);
  bkg->SetMarkerSize(0);
  bkg->SetLineColor(1);

  gjet->SetFillColor(kOrange-3);
  gjet->SetMarkerSize(0);
  gjet->SetLineColor(1);

  ewk->SetFillColor(kYellow);
  ewk->SetMarkerSize(0);
  ewk->SetLineColor(1);

  ttg->SetFillColor(8);
  ttg->SetMarkerSize(0);
  ttg->SetLineColor(1);

  ttbar->SetFillColor(kRed-3);
  ttbar->SetMarkerSize(0);
  ttbar->SetLineColor(1);

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

  bkg->SetTitle(variable);
  bkg->GetXaxis()->SetTitle(xaxisTitle);
  bkg->GetYaxis()->SetTitle(yaxisTitle);

  if(xmax > xmin) bkg->GetXaxis()->SetRangeUser(xmin, xmax);
  bkg->GetYaxis()->SetRangeUser(ymin, ymax);

  // new stack: qcd, gjet, ewk, ttg, ttbar
  bkg->Draw("hist");
  gjet->Draw("same hist");
  ewk->Draw("same hist");
  ttg->Draw("same hist");
  ttbar->Draw("same hist");
  errors->Draw("same e2");
  gg->Draw("same e1");
  bkg->Draw("same axis");

  if(drawSignal) {
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
  can->SaveAs(variable+"_"+req+".pdf");

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
