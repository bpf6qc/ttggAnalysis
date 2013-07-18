#ifndef BtagInfo_h
#define BtagInfo_h

#include <math.h>

#include "ScaleFactorInfo.h"
#include "../src/SusyEvent.h"

using namespace std;

// Container for btagging-related values
// These are calculated from a ScaleFactorInfo object for each jet
class BtagInfo {
 public:
  BtagInfo(susy::PFJet& iJet, TLorentzVector& iCorrP4, TString& iTagger, double iScale, bool& iIsMC, bool& iIsFastSim, ScaleFactorInfo& iSfInfo, TString DataPeriod, double fractionAB, double fractionC, double fractionD);
  virtual ~BtagInfo() {};

  bool isTagged() { return (discr > discr_cut); };
  //bool isRandomizedTagged(Float_t coinflip);

  Float_t GetDiscriminant() { return discr; };

  Float_t GetSFlight() { return SFl; };
  Float_t GetSFlight_up() { return SFl + scale*(SFl_max - SFl); };
  Float_t GetSFlight_down() { return SFl - scale*(SFl - SFl_min); };

  Float_t GetSFbottom() { return SFb; };
  Float_t GetSFbottom_up() { return SFb + scale*SFb_error; };
  Float_t GetSFbottom_down() { return SFb - scale*SFb_error; };

  Float_t GetSFcharm() { return SFc; };
  Float_t GetSFcharm_up() { return SFc + scale*SFc_error; };
  Float_t GetSFcharm_down() { return SFc - scale*SFc_error; };

  Float_t GetScaleFactor(double iScale, bool doScale);

  Float_t GetTaggingEfficiency() { return eff; };
  Float_t GetTaggingEfficiencyError() { return eff_error; };

 private:
  TString tagger;
  Float_t discr_cut;

  Float_t discr;  
  Float_t pdgId;
  Float_t pt;
  Float_t eta;

  Float_t SFl;
  Float_t SFl_max;
  Float_t SFl_min;
  Float_t SFb;
  Float_t SFb_error;
  Float_t SFc;
  Float_t SFc_error;

  Float_t eff;
  Float_t eff_error;

  bool isMC;
  bool isFastSim;

  double scale;

};
  
BtagInfo::BtagInfo(susy::PFJet& iJet, TLorentzVector& iCorrP4, TString& iTagger, double iScale, bool& iIsMC, bool& iIsFastSim, ScaleFactorInfo& iSfInfo, TString DataPeriod, double fractionAB, double fractionC, double fractionD) {

  isMC = iIsMC;
  isFastSim = iIsFastSim;

  scale = iScale;

  Float_t ptmin[16] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
  Float_t ptmax[16] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};

  tagger = iTagger;

  if(tagger == "CSVL") discr_cut = 0.244;
  else if(tagger == "CSVM") discr_cut = 0.679;
  else if(tagger == "CSVT") discr_cut = 0.898;
  else discr_cut = 1.e6;

  if(tagger.Contains("CSV")) discr = iJet.bTagDiscriminators[susy::kCSV];
  else discr = -1.e6;

  pdgId = fabs(iJet.algDefFlavour);
  pt = iCorrP4.Pt();
  eta = fabs(iCorrP4.Eta());

  if(DataPeriod == "2012") {

    // for pt > 800 GeV and 0<abs(eta)<1.5: use the SFlight value at 800 GeV with twice the quoted uncertainty 
    if(pt > 800. && eta < 1.5) {
      SFl = iSfInfo.GetSFlightMean(eta, "AB", 800.) * fractionAB +
	iSfInfo.GetSFlightMean(eta, "C", 800.) * fractionC +
	iSfInfo.GetSFlightMean(eta, "D", 800.) * fractionD;
      
      SFl_max = iSfInfo.GetSFlightMax(eta, "AB", 800.) * fractionAB +
	iSfInfo.GetSFlightMax(eta, "C", 800.) * fractionC +
	iSfInfo.GetSFlightMax(eta, "D", 800.) * fractionD;
      SFl_max = 2.*SFl_max - SFl;

      SFl_min = iSfInfo.GetSFlightMin(eta, "AB", 800.) * fractionAB +
	iSfInfo.GetSFlightMin(eta, "C", 800.) * fractionC +
	iSfInfo.GetSFlightMin(eta, "D", 800.) * fractionD;
      SFl_min = SFl - 2.*SFl_min;
    }

    // for pt > 700 GeV and 1.5<abs(eta): use the SFlight value at 700 GeV with twice the quoted uncertainty 
    else if(pt > 700. && eta >= 1.5) {
      SFl = iSfInfo.GetSFlightMean(eta, "AB", 700.) * fractionAB +
	iSfInfo.GetSFlightMean(eta, "C", 700.) * fractionC +
	iSfInfo.GetSFlightMean(eta, "D", 700.) * fractionD;
      
      SFl_max = iSfInfo.GetSFlightMax(eta, "AB", 700.) * fractionAB +
	iSfInfo.GetSFlightMax(eta, "C", 700.) * fractionC +
	iSfInfo.GetSFlightMax(eta, "D", 700.) * fractionD;
      SFl_max = 2.*SFl_max - SFl;

      SFl_min = iSfInfo.GetSFlightMin(eta, "AB", 700.) * fractionAB +
	iSfInfo.GetSFlightMin(eta, "C", 700.) * fractionC +
	iSfInfo.GetSFlightMin(eta, "D", 700.) * fractionD;
      SFl_min = SFl - 2.*SFl_min;
    }
    
    else {
      SFl = iSfInfo.GetSFlightMean(eta, "AB", pt) * fractionAB +
	iSfInfo.GetSFlightMean(eta, "C", pt) * fractionC +
	iSfInfo.GetSFlightMean(eta, "D", pt) * fractionD;
      
      SFl_max = iSfInfo.GetSFlightMax(eta, "AB", pt) * fractionAB +
	iSfInfo.GetSFlightMax(eta, "C", pt) * fractionC +
	iSfInfo.GetSFlightMax(eta, "D", pt) * fractionD;
      
      SFl_min = iSfInfo.GetSFlightMin(eta, "AB", pt) * fractionAB +
	iSfInfo.GetSFlightMin(eta, "C", pt) * fractionC +
	iSfInfo.GetSFlightMin(eta, "D", pt) * fractionD;
    }

  } // if using all of 2012: AB + C + D and have fractions for those periods

  else {

    // for pt > 800 GeV and 0<abs(eta)<1.5: use the SFlight value at 800 GeV with twice the quoted uncertainty 
    if(pt > 800. && eta < 1.5) {
      SFl = iSfInfo.GetSFlightMean(eta, DataPeriod, 800.);
      
      SFl_max = iSfInfo.GetSFlightMax(eta, DataPeriod, 800.);
      SFl_max = 2.*SFl_max - SFl;

      SFl_min = iSfInfo.GetSFlightMin(eta, DataPeriod, 800.);
      SFl_min = SFl - 2.*SFl_min;
    }

    // for pt > 700 GeV and 1.5<abs(eta): use the SFlight value at 700 GeV with twice the quoted uncertainty 
    else if(pt > 700. && eta >= 1.5) {
      SFl = iSfInfo.GetSFlightMean(eta, DataPeriod, 700.);
      
      SFl_max = iSfInfo.GetSFlightMax(eta, DataPeriod, 700.);
      SFl_max = 2.*SFl_max - SFl;

      SFl_min = iSfInfo.GetSFlightMin(eta, DataPeriod, 700.);
      SFl_min = SFl - 2.*SFl_min;
    }
    
    else {
      SFl = iSfInfo.GetSFlightMean(eta, DataPeriod, pt);
      
      SFl_max = iSfInfo.GetSFlightMax(eta, DataPeriod, pt);
      
      SFl_min = iSfInfo.GetSFlightMin(eta, DataPeriod, pt);
    }

  } // if you're using only one of AB, C, or D

  // for pt > 800 GeV: use the SFb value at 800 GeV with twice the quoted uncertainty 
  if(pt > 800.) {
    SFb = iSfInfo.GetSFbottom(800.);
    SFb_error = (iSfInfo.GetSFbottomErrors())[15] * 2.;
  }
  
  // for pt < 20 GeV: use the SFb value at 20 GeV with twice the quoted uncertainty
  else if(pt < 20.) {
    SFb = iSfInfo.GetSFbottom(20.);
    SFb_error = (iSfInfo.GetSFbottomErrors())[0] * 2.;
  }
  
  else {
    SFb = iSfInfo.GetSFbottom(pt);
    for(int i = 0; i < 16; i++) {
      if(pt >= ptmin[i] && pt < ptmax[i]) SFb_error = (iSfInfo.GetSFbottomErrors())[i];
    }
  }
  
  // SFc = SFb with twice the quoted uncertainty
  SFc = SFb;
  SFc_error = 2. * SFb_error;
  
  if(pdgId == 5) {
    eff = iSfInfo.GetEffBottom(pt);
    eff_error = iSfInfo.GetEffErrorBottom(pt);
  }
  else if(pdgId == 4) {
    eff = iSfInfo.GetEffCharm(pt);
    eff_error = iSfInfo.GetEffErrorCharm(pt);
  }
  else if(pdgId == 1 || pdgId == 2 || pdgId == 3 || pdgId == 21) {
    eff = iSfInfo.GetEffLight(pt);
    eff_error = iSfInfo.GetEffErrorLight(pt);
  }

  // Additional correction to scale factors if this is FastSim
  if(isFastSim && isMC) {

    float ptmin_cf[14] = {30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500};
    float ptmax_cf[14] = {40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 670};

    Double_t CFb = 0., CFc = 0., CFl = 0.;
    Double_t CFb_err = 0., CFc_err = 0., CFl_err = 0.;
    // For jets with pT > 670 GeV, the correction factors from the highest jet pT bin (500 < pT< 670) are used, with the uncertainties doubled.
    if(pt > 670.) {
      //13
      CFb = (iSfInfo.GetCFbottom())[13];
      CFb_err = (iSfInfo.GetCFbottomErrors())[13] * 2.;
      CFc = (iSfInfo.GetCFcharm())[13];
      CFc_err = (iSfInfo.GetCFcharmErrors())[13] * 2.;
      CFl = (iSfInfo.GetCFlight(eta))[13];
      CFl_err = (iSfInfo.GetCFlightErrors(eta))[13] * 2.;
    }
    else {
      for(int i = 0; i < 14; i++) {
	if(pt >= ptmin_cf[i] && pt < ptmax_cf[i]) {
	  CFb = (iSfInfo.GetCFbottom())[i];
	  CFb_err = (iSfInfo.GetCFbottomErrors())[i];
	  CFc = (iSfInfo.GetCFcharm())[i];
	  CFc_err = (iSfInfo.GetCFcharmErrors())[i];
	  CFl = (iSfInfo.GetCFlight(eta))[i];
	  CFl_err = (iSfInfo.GetCFlightErrors(eta))[i];
	}
      }
    }

    // Add errors in quadrature
    SFb_error = CFb * SFb * sqrt(CFb_err*CFb_err/CFb/CFb + SFb_error*SFb_error/SFb/SFb);
    SFc_error = CFc * SFc * sqrt(CFc_err*CFc_err/CFc/CFc + SFc_error*SFc_error/SFc/SFc);
    SFl_max = SFl * CFl * (1. + sqrt(CFl_err*CFl_err/CFl/CFl + (SFl_max/SFl + 1.)*(SFl_max/SFl + 1.)));
    SFl_min = SFl * CFl * (1. - sqrt(CFl_err*CFl_err/CFl/CFl + (SFl_min/SFl + 1.)*(SFl_min/SFl + 1.)));
    
    // Apply the correction
    SFb = CFb * SFb;
    SFc = CFc * SFc;
    SFl = CFl * SFl;

  } // if FastSim


} // Constructor

Float_t BtagInfo::GetScaleFactor(double iScale, bool doScale) { 
  if(!doScale) {
    if(pdgId == 5) return GetSFbottom();
    else if(pdgId == 4) return GetSFcharm();
    else if(pdgId == 1 || pdgId == 2 || pdgId == 3 || pdgId == 21) return GetSFlight();
  }
  else {
    if(iScale > 0.) {
      scale = iScale;
      if(pdgId == 5) return GetSFbottom_up();
      else if(pdgId == 4) return GetSFcharm_up();
      else if(pdgId == 1 || pdgId == 2 || pdgId == 3 || pdgId == 21) return GetSFlight_up();
    }
    else if(iScale < 0.) {
      scale = -1. * iScale;
      if(pdgId == 5) return GetSFbottom_down();
      else if(pdgId == 4) return GetSFcharm_down();
      else if(pdgId == 1 || pdgId == 2 || pdgId == 3 || pdgId == 21) return GetSFlight_down();
    }
  }

  return 0.;
}

#endif
