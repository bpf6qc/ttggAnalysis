#ifndef ScaleFactorInfo_h
#define ScaleFactorInfo_h

#include <TString.h>
#include <TF1.h>
#include <TH1.h>

#include <iostream>
#include <vector>

using namespace std;

// Basic container for all scale factors, errors, and efficiencies
class ScaleFactorInfo {
 public:
  ScaleFactorInfo(TString tag);
  virtual ~ScaleFactorInfo() {};

  Float_t GetSFbottom(Float_t pt);
  Float_t GetSFcharm(Float_t pt) { return GetSFbottom(pt); }; // SFc = SFb with twice the quoted uncertainty

  // SFlightMin and SFlightMax correspond to SFlight +- (stat+syst error)
  Float_t GetSFlight(TString meanminmax, Float_t eta, TString DataPeriod, Float_t pt);
  Float_t GetSFlightMin(Float_t eta, TString DataPeriod, Float_t pt)  { return GetSFlight("min", eta, DataPeriod, pt); };
  Float_t GetSFlightMean(Float_t eta, TString DataPeriod, Float_t pt) { return GetSFlight("mean", eta, DataPeriod, pt); };
  Float_t GetSFlightMax(Float_t eta, TString DataPeriod, Float_t pt)  { return GetSFlight("max", eta, DataPeriod, pt); };
  
  // Tagging efficiencies. MC reweighting requires knowledge of all tagging efficiencies for that sample beforehand
  void SetTaggingEfficiencies(TH1F * eff_light, TH1F * eff_charm, TH1F * eff_bottom) {
    hasEfficiencies = true;
    lEff = (TH1F*)eff_light->Clone();
    cEff = (TH1F*)eff_charm->Clone();
    bEff = (TH1F*)eff_bottom->Clone();
  }

  Float_t GetEffLight(Float_t pt) { return (hasEfficiencies) ? lEff->GetBinContent(lEff->GetXaxis()->FindBin(pt)) : 0.; };
  Float_t GetEffCharm(Float_t pt) { return (hasEfficiencies) ? cEff->GetBinContent(cEff->GetXaxis()->FindBin(pt)) : 0.; };
  Float_t GetEffBottom(Float_t pt) { return (hasEfficiencies) ? bEff->GetBinContent(bEff->GetXaxis()->FindBin(pt)) : 0.; };

  Float_t GetEffErrorLight(Float_t pt) { return (hasEfficiencies) ? lEff->GetBinError(lEff->GetXaxis()->FindBin(pt)) : 0.; };
  Float_t GetEffErrorCharm(Float_t pt) { return (hasEfficiencies) ? cEff->GetBinError(cEff->GetXaxis()->FindBin(pt)) : 0.; };
  Float_t GetEffErrorBottom(Float_t pt) { return (hasEfficiencies) ? bEff->GetBinError(bEff->GetXaxis()->FindBin(pt)) : 0.; };

  // FastSim to FullSim correction factors
  vector<Float_t> GetCFbottom() { return CFb; };
  vector<Float_t> GetCFcharm() { return CFc; };
  vector<Float_t> GetCFlight(Float_t eta) {
    if(eta < 1.2) return CFudsg_12;
    return CFudsg_24;
  };

  vector<Float_t> GetCFbottomErrors() { return CFb_errors; };
  vector<Float_t> GetCFcharmErrors() { return CFc_errors; };
  vector<Float_t> GetCFlightErrors(Float_t eta) {
    if(eta < 1.2) return CFudsg_errors_12;
    return CFudsg_errors_24;
  };

  vector<Float_t> GetSFbottomErrors() { return SFb_error; };
  vector<Float_t> GetSFcharmErrors() { // SFc = SFb with twice the quoted uncertainty
    vector<Float_t> err = GetSFbottomErrors();
    for(unsigned int i = 0; i < err.size(); i++) err[i] *= 2.;
    return err;
  }
 
 private:
  TString tagger;

  bool hasEfficiencies;
  TH1F * lEff;
  TH1F * cEff;
  TH1F * bEff;

  vector<Float_t> CFb;
  vector<Float_t> CFc;
  vector<Float_t> CFudsg_12;
  vector<Float_t> CFudsg_24;

  vector<Float_t> CFb_errors;
  vector<Float_t> CFc_errors;
  vector<Float_t> CFudsg_errors_12;
  vector<Float_t> CFudsg_errors_24;

  vector<Float_t> SFb_error;

};

ScaleFactorInfo::ScaleFactorInfo(TString tag) {

  tagger = tag;
  hasEfficiencies = false;

  if(tagger == "CSVM") {
    CFb.push_back(0.982194); 
    CFb.push_back(0.980998); 
    CFb.push_back(0.992014); 
    CFb.push_back(0.994472); 
    CFb.push_back(0.996825); 
    CFb.push_back(0.999822); 
    CFb.push_back(1.00105); 
    CFb.push_back(1.00023); 
    CFb.push_back(0.991994); 
    CFb.push_back(0.979123); 
    CFb.push_back(0.947207); 
    CFb.push_back(0.928006); 
    CFb.push_back(0.874260); 
    CFb.push_back(0.839610);
  }
  else if(tagger == "CSVT") {
    CFb.push_back(0.968117); 
    CFb.push_back(0.967822); 
    CFb.push_back(0.978278); 
    CFb.push_back(0.981281); 
    CFb.push_back(0.987679); 
    CFb.push_back(0.986590); 
    CFb.push_back(0.990246); 
    CFb.push_back(0.984504); 
    CFb.push_back(0.967024); 
    CFb.push_back(0.940042); 
    CFb.push_back(0.873019); 
    CFb.push_back(0.850847); 
    CFb.push_back(0.769561); 
    CFb.push_back(0.650192);
  }
  else CFb.clear();

  if(tagger == "CSVM") {
    CFc.push_back(0.988545); 
    CFc.push_back(0.981714); 
    CFc.push_back(1.00946); 
    CFc.push_back(1.01591); 
    CFc.push_back(1.02810); 
    CFc.push_back(1.02195); 
    CFc.push_back(1.02590); 
    CFc.push_back(1.01936); 
    CFc.push_back(0.991228); 
    CFc.push_back(0.955343); 
    CFc.push_back(0.944433); 
    CFc.push_back(0.917282); 
    CFc.push_back(0.935018); 
    CFc.push_back(1.06375);
  }
  else if(tagger == "CSVT") {
    CFc.push_back(0.960959); 
    CFc.push_back(0.973876); 
    CFc.push_back(0.984323); 
    CFc.push_back(0.996344); 
    CFc.push_back(1.02418); 
    CFc.push_back(0.985580); 
    CFc.push_back(0.994745); 
    CFc.push_back(0.970489); 
    CFc.push_back(0.914155); 
    CFc.push_back(0.872072); 
    CFc.push_back(0.945289); 
    CFc.push_back(0.783816); 
    CFc.push_back(0.942773); 
    CFc.push_back(0.527354);
  }
  else CFc.clear();

  if(tagger == "CSVM") {
    CFudsg_12.push_back(1.21878); 
    CFudsg_12.push_back(1.28615); 
    CFudsg_12.push_back(1.37535); 
    CFudsg_12.push_back(1.38966); 
    CFudsg_12.push_back(1.40320); 
    CFudsg_12.push_back(1.49835); 
    CFudsg_12.push_back(1.44308); 
    CFudsg_12.push_back(1.58198); 
    CFudsg_12.push_back(1.55687); 
    CFudsg_12.push_back(1.65790); 
    CFudsg_12.push_back(1.90233); 
    CFudsg_12.push_back(1.92259); 
    CFudsg_12.push_back(2.66174); 
    CFudsg_12.push_back(3.08688);

    CFudsg_24.push_back(1.46970); 
    CFudsg_24.push_back(1.48732); 
    CFudsg_24.push_back(1.69024); 
    CFudsg_24.push_back(1.64494); 
    CFudsg_24.push_back(1.79297); 
    CFudsg_24.push_back(1.90760); 
    CFudsg_24.push_back(1.99867); 
    CFudsg_24.push_back(2.21659); 
    CFudsg_24.push_back(2.20103); 
    CFudsg_24.push_back(2.42645); 
    CFudsg_24.push_back(2.67594); 
    CFudsg_24.push_back(4.24735); 
    CFudsg_24.push_back(3.98979); 
    CFudsg_24.push_back(15.0457);
  }
  else if(tagger == "CSVT") {
    CFudsg_12.push_back(1.24890); 
    CFudsg_12.push_back(1.35145); 
    CFudsg_12.push_back(1.37205); 
    CFudsg_12.push_back(1.32472); 
    CFudsg_12.push_back(1.39976); 
    CFudsg_12.push_back(1.45884); 
    CFudsg_12.push_back(1.59912); 
    CFudsg_12.push_back(1.58971); 
    CFudsg_12.push_back(1.30841); 
    CFudsg_12.push_back(1.55936); 
    CFudsg_12.push_back(1.28346); 
    CFudsg_12.push_back(2.21265); 
    CFudsg_12.push_back(2.06927); 
    CFudsg_12.push_back(2.88109);

    CFudsg_24.push_back(1.67634); 
    CFudsg_24.push_back(1.70105); 
    CFudsg_24.push_back(1.75999); 
    CFudsg_24.push_back(1.78459); 
    CFudsg_24.push_back(2.19343); 
    CFudsg_24.push_back(2.73199); 
    CFudsg_24.push_back(3.49277); 
    CFudsg_24.push_back(2.58863); 
    CFudsg_24.push_back(2.48824); 
    CFudsg_24.push_back(4.01723); 
    CFudsg_24.push_back(3.86956); 
    CFudsg_24.push_back(0.000456049); 
    CFudsg_24.push_back(2.30988); 
    CFudsg_24.push_back(0.000855693);
  }
  else {
    CFudsg_12.clear();
    CFudsg_24.clear();
  }

  if(tagger == "CSVM") {
    CFb_errors.push_back(0.00253112); 
    CFb_errors.push_back(0.00296453); 
    CFb_errors.push_back(0.00113963); 
    CFb_errors.push_back(0.00128363); 
    CFb_errors.push_back(0.00232566); 
    CFb_errors.push_back(0.00232353); 
    CFb_errors.push_back(0.00219086); 
    CFb_errors.push_back(0.00156856); 
    CFb_errors.push_back(0.00322279); 
    CFb_errors.push_back(0.00400414); 
    CFb_errors.push_back(0.00737465); 
    CFb_errors.push_back(0.0105033); 
    CFb_errors.push_back(0.0171706); 
    CFb_errors.push_back(0.0344172);
  }
  else if(tagger == "CSVT") {
    CFb_errors.push_back(0.00223422); 
    CFb_errors.push_back(0.00367427); 
    CFb_errors.push_back(0.00145554); 
    CFb_errors.push_back(0.00337572); 
    CFb_errors.push_back(0.00344106); 
    CFb_errors.push_back(0.00591257); 
    CFb_errors.push_back(0.00218050); 
    CFb_errors.push_back(0.00472939); 
    CFb_errors.push_back(0.00353119); 
    CFb_errors.push_back(0.00739502); 
    CFb_errors.push_back(0.0193330); 
    CFb_errors.push_back(0.0158257); 
    CFb_errors.push_back(0.0306048); 
    CFb_errors.push_back(0.0603701);
  }
  else CFb_errors.clear();

  if(tagger == "CSVM") {
    CFc_errors.push_back(0.00746259); 
    CFc_errors.push_back(0.00661831); 
    CFc_errors.push_back(0.00968682); 
    CFc_errors.push_back(0.00751322); 
    CFc_errors.push_back(0.00675507); 
    CFc_errors.push_back(0.00562821); 
    CFc_errors.push_back(0.00862890); 
    CFc_errors.push_back(0.00768003); 
    CFc_errors.push_back(0.0188981); 
    CFc_errors.push_back(0.0261163); 
    CFc_errors.push_back(0.0450601); 
    CFc_errors.push_back(0.0448453); 
    CFc_errors.push_back(0.148805); 
    CFc_errors.push_back(0.177157);
  }
  else if(tagger == "CSVT") {
    CFc_errors.push_back(0.0155733); 
    CFc_errors.push_back(0.0121900); 
    CFc_errors.push_back(0.0131678); 
    CFc_errors.push_back(0.0113739); 
    CFc_errors.push_back(0.0213937); 
    CFc_errors.push_back(0.0123294); 
    CFc_errors.push_back(0.0153230); 
    CFc_errors.push_back(0.0156350); 
    CFc_errors.push_back(0.0409568); 
    CFc_errors.push_back(0.0654966); 
    CFc_errors.push_back(0.112785); 
    CFc_errors.push_back(0.187795); 
    CFc_errors.push_back(0.331301); 
    CFc_errors.push_back(0.162462);
  }
  else CFc_errors.clear();

  if(tagger == "CSVM") {
    CFudsg_errors_12.push_back(0.0182686); 
    CFudsg_errors_12.push_back(0.0373732); 
    CFudsg_errors_12.push_back(0.0461870); 
    CFudsg_errors_12.push_back(0.0288973); 
    CFudsg_errors_12.push_back(0.0333528); 
    CFudsg_errors_12.push_back(0.0513836); 
    CFudsg_errors_12.push_back(0.0420353); 
    CFudsg_errors_12.push_back(0.106627); 
    CFudsg_errors_12.push_back(0.0658359); 
    CFudsg_errors_12.push_back(0.117285); 
    CFudsg_errors_12.push_back(0.185533); 
    CFudsg_errors_12.push_back(0.214071); 
    CFudsg_errors_12.push_back(0.487274); 
    CFudsg_errors_12.push_back(0.871502);
    
    CFudsg_errors_24.push_back(0.104716); 
    CFudsg_errors_24.push_back(0.0392025); 
    CFudsg_errors_24.push_back(0.106315); 
    CFudsg_errors_24.push_back(0.115751); 
    CFudsg_errors_24.push_back(0.106807); 
    CFudsg_errors_24.push_back(0.0642086); 
    CFudsg_errors_24.push_back(0.138742); 
    CFudsg_errors_24.push_back(0.182345); 
    CFudsg_errors_24.push_back(0.169922); 
    CFudsg_errors_24.push_back(0.297889); 
    CFudsg_errors_24.push_back(0.320088); 
    CFudsg_errors_24.push_back(0.927736); 
    CFudsg_errors_24.push_back(1.24666); 
    CFudsg_errors_24.push_back(15.1860);
  }
  else if(tagger == "CSVT") {
    CFudsg_errors_12.push_back(0.0751438); 
    CFudsg_errors_12.push_back(0.0651619); 
    CFudsg_errors_12.push_back(0.0604241); 
    CFudsg_errors_12.push_back(0.0726285); 
    CFudsg_errors_12.push_back(0.0968158); 
    CFudsg_errors_12.push_back(0.0931768); 
    CFudsg_errors_12.push_back(0.163039); 
    CFudsg_errors_12.push_back(0.187749); 
    CFudsg_errors_12.push_back(0.198200); 
    CFudsg_errors_12.push_back(0.465354); 
    CFudsg_errors_12.push_back(0.339473); 
    CFudsg_errors_12.push_back(1.07079); 
    CFudsg_errors_12.push_back(1.07723); 
    CFudsg_errors_12.push_back(2.53188);
  
    CFudsg_errors_24.push_back(0.222165); 
    CFudsg_errors_24.push_back(0.161403); 
    CFudsg_errors_24.push_back(0.112342); 
    CFudsg_errors_24.push_back(0.275101); 
    CFudsg_errors_24.push_back(0.364229); 
    CFudsg_errors_24.push_back(0.330588); 
    CFudsg_errors_24.push_back(1.00953); 
    CFudsg_errors_24.push_back(0.404417); 
    CFudsg_errors_24.push_back(1.07731); 
    CFudsg_errors_24.push_back(2.65686); 
    CFudsg_errors_24.push_back(3.18286); 
    CFudsg_errors_24.push_back(5.25051e-05); 
    CFudsg_errors_24.push_back(2.38652); 
    CFudsg_errors_24.push_back(0.000438728);
  }
  else {
    CFudsg_errors_12.clear();
    CFudsg_errors_24.clear();
  }

  if(tagger == "CSVL") {
    SFb_error.push_back(0.0484285);
    SFb_error.push_back(0.0126178);
    SFb_error.push_back(0.0120027);
    SFb_error.push_back(0.0141137);
    SFb_error.push_back(0.0145441);
    SFb_error.push_back(0.0131145);
    SFb_error.push_back(0.0168479);
    SFb_error.push_back(0.0160836);
    SFb_error.push_back(0.0126209);
    SFb_error.push_back(0.0136017);
    SFb_error.push_back(0.019182);
    SFb_error.push_back(0.0198805);
    SFb_error.push_back(0.0386531);
    SFb_error.push_back(0.0392831);
    SFb_error.push_back(0.0481008);
    SFb_error.push_back(0.0474291);
  }
  else if(tagger == "CSVM") {
    SFb_error.push_back(0.0554504);
    SFb_error.push_back(0.0209663);
    SFb_error.push_back(0.0207019);
    SFb_error.push_back(0.0230073);
    SFb_error.push_back(0.0208719);
    SFb_error.push_back(0.0200453);
    SFb_error.push_back(0.0264232);
    SFb_error.push_back(0.0240102);
    SFb_error.push_back(0.0229375);
    SFb_error.push_back(0.0184615);
    SFb_error.push_back(0.0216242);
    SFb_error.push_back(0.0248119);
    SFb_error.push_back(0.0465748);
    SFb_error.push_back(0.0474666);
    SFb_error.push_back(0.0718173);
    SFb_error.push_back(0.0717567);
  }
  else if(tagger == "CSVT") {
    SFb_error.push_back(0.0567059);
    SFb_error.push_back(0.0266907);
    SFb_error.push_back(0.0263491);
    SFb_error.push_back(0.0342831);
    SFb_error.push_back(0.0303327);
    SFb_error.push_back(0.024608);
    SFb_error.push_back(0.0333786);
    SFb_error.push_back(0.0317642);
    SFb_error.push_back(0.031102);
    SFb_error.push_back(0.0295603);
    SFb_error.push_back(0.0474663);
    SFb_error.push_back(0.0503182);
    SFb_error.push_back(0.0580424);
    SFb_error.push_back(0.0575776);
    SFb_error.push_back(0.0769779);
    SFb_error.push_back(0.0898199);
  }
  else SFb_error.clear();

}



Float_t ScaleFactorInfo::GetSFbottom(Float_t pt) {

  Float_t x = pt;

  if(tagger == "CSVL") return 0.981149*((1.+(-0.000713295*x))/(1.+(-0.000703264*x)));
  else if(tagger == "CSVM") return 0.726981*((1.+(0.253238*x))/(1.+(0.188389*x)));
  else if(tagger == "CSVT") return 0.869965*((1.+(0.0335062*x))/(1.+(0.0304598*x)));

  return 0.;
}

Float_t ScaleFactorInfo::GetSFlight(TString meanminmax, Float_t eta, TString DataPeriod, Float_t pt) {
  
  Float_t eta_ = fabs(eta);
  Float_t x = pt;

  /*
  Double_t ptmax;
  if((tagger == "CSVL" && eta_ < 1.5) ||
     ((tagger == "CSVM" || tagger == "CSVT") && eta_ < 1.6)
     ) ptmax = 800.;
  else ptmax = 700.;
  */

  if(DataPeriod == "ABCD") {

    if(tagger == "CSVL") {
      if(eta_ < 0.5) {
	if( meanminmax == "mean" ) return ((1.04901+(0.00152181*x))+(-3.43568e-06*(x*x)))+(2.17219e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.973773+(0.00103049*x))+(-2.2277e-06*(x*x)))+(1.37208e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.12424+(0.00201136*x))+(-4.64021e-06*(x*x)))+(2.97219e-09*(x*(x*x)));
      }
      else if(eta_ >= 0.5 && eta_ < 1.0) {
	if( meanminmax == "mean" ) return ((0.991915+(0.00172552*x))+(-3.92652e-06*(x*x)))+(2.56816e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.921518+(0.00129098*x))+(-2.86488e-06*(x*x)))+(1.86022e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.06231+(0.00215815*x))+(-4.9844e-06*(x*x)))+(3.27623e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.0 && eta_ < 1.5) {
	if( meanminmax == "mean" ) return ((0.962127+(0.00192796*x))+(-4.53385e-06*(x*x)))+(3.0605e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.895419+(0.00153387*x))+(-3.48409e-06*(x*x)))+(2.30899e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.02883+(0.00231985*x))+(-5.57924e-06*(x*x)))+(3.81235e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.5 && eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.06121+(0.000332747*x))+(-8.81201e-07*(x*x)))+(7.43896e-10*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.983607+(0.000196747*x))+(-3.98327e-07*(x*x)))+(2.95764e-10*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.1388+(0.000468418*x))+(-1.36341e-06*(x*x)))+(1.19256e-09*(x*(x*x)));
      }
    } // if CSVL
    
    else if(tagger == "CSVM") {
      if(eta_ < 0.8) {
	if( meanminmax == "mean" ) return ((1.06238+(0.00198635*x))+(-4.89082e-06*(x*x)))+(3.29312e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.972746+(0.00104424*x))+(-2.36081e-06*(x*x)))+(1.53438e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.15201+(0.00292575*x))+(-7.41497e-06*(x*x)))+(5.0512e-09*(x*(x*x)));
      }
      else if(eta_ >= 0.8 && eta_ < 1.6) {
	if( meanminmax == "mean" ) return ((1.08048+(0.00110831*x))+(-2.96189e-06*(x*x)))+(2.16266e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.9836+(0.000649761*x))+(-1.59773e-06*(x*x)))+(1.14324e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.17735+(0.00156533*x))+(-4.32257e-06*(x*x)))+(3.18197e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.6 && eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.09145+(0.000687171*x))+(-2.45054e-06*(x*x)))+(1.7844e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((1.00616+(0.000358884*x))+(-1.23768e-06*(x*x)))+(6.86678e-10*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.17671+(0.0010147*x))+(-3.66269e-06*(x*x)))+(2.88425e-09*(x*(x*x)));
      }
    } // if CSVM

    else if(tagger == "CSVT") {
      if(eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.01739+(0.00283619*x))+(-7.93013e-06*(x*x)))+(5.97491e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.953587+(0.00124872*x))+(-3.97277e-06*(x*x)))+(3.23466e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.08119+(0.00441909*x))+(-1.18764e-05*(x*x)))+(8.71372e-09*(x*(x*x)));
      }
    } // if CSVT

  } // if SF_12 ABCD

  else if(DataPeriod == "AB") {

    if(tagger == "CSVL") {
      if(eta_ < 0.5) {
	if( meanminmax == "mean" ) return ((1.00989+(0.00155686*x))+(-3.72647e-06*(x*x)))+(2.47025e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.947488+(0.00105091*x))+(-2.43972e-06*(x*x)))+(1.58902e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.07229+(0.00206098*x))+(-5.00971e-06*(x*x)))+(3.35179e-09*(x*(x*x)));
      }
      else if(eta_ >= 0.5 && eta_ < 1.0) {
	if( meanminmax == "mean" ) return ((0.958598+(0.00173458*x))+(-4.12744e-06*(x*x)))+(2.83257e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.900024+(0.00129392*x))+(-3.01708e-06*(x*x)))+(2.06723e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.01716+(0.00217335*x))+(-5.23419e-06*(x*x)))+(3.5986e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.0 && eta_ < 1.5) {
	if( meanminmax == "mean" ) return ((0.963113+(0.00163674*x))+(-3.84776e-06*(x*x)))+(2.56918e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.90596+(0.00125465*x))+(-2.78863e-06*(x*x)))+(1.78602e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.02025+(0.0020171*x))+(-4.90389e-06*(x*x)))+(3.35329e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.5 && eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.04996+(0.00031979*x))+(-8.43322e-07*(x*x)))+(6.9451e-10*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.983472+(0.000169396*x))+(-2.82848e-07*(x*x)))+(1.52744e-10*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.11645+(0.000469873*x))+(-1.40321e-06*(x*x)))+(1.23681e-09*(x*(x*x)));
      }
    } // if CSVL

    else if(tagger == "CSVM") {
      if(eta_ < 0.8) {
	if( meanminmax == "mean" ) return ((1.02213+(0.00189078*x))+(-4.59419e-06*(x*x)))+(3.0366e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.946+(0.000940317*x))+(-1.99048e-06*(x*x)))+(1.18343e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.09827+(0.00283897*x))+(-7.19354e-06*(x*x)))+(4.89013e-09*(x*(x*x)));
      }
      else if(eta_ >= 0.8 && eta_ < 1.6) {
	if( meanminmax == "mean" ) return ((1.0596+(0.00102926*x))+(-2.70312e-06*(x*x)))+(1.82871e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.974966+(0.000545735*x))+(-1.23123e-06*(x*x)))+(7.05661e-10*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.14423+(0.00151156*x))+(-4.17277e-06*(x*x)))+(2.95233e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.6 && eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.04976+(0.000897158*x))+(-3.22829e-06*(x*x)))+(2.71316e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.977166+(0.000550586*x))+(-1.91114e-06*(x*x)))+(1.44817e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.12232+(0.00124269*x))+(-4.54368e-06*(x*x)))+(3.98079e-09*(x*(x*x)));
      }
    } // if CSVM

    else if(tagger == "CSVT") {
      if(eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((0.985589+(0.00302526*x))+(-8.73861e-06*(x*x)))+(6.65051e-09*(x*(x*x)));
        if( meanminmax == "min" ) return ((0.93612+(0.00131596*x))+(-4.30052e-06*(x*x)))+(3.45957e-09*(x*(x*x)));
        if( meanminmax == "max" ) return ((1.03505+(0.00472994*x))+(-1.31661e-05*(x*x)))+(9.84151e-09*(x*(x*x)));
      }
    } // if CSVT

  } // if SF_12 AB

  else if(DataPeriod == "C") {
    
    if(tagger == "CSVL") {
      if(eta_ < 0.5) {
	if( meanminmax == "mean" ) return ((1.03512+(0.00172098*x))+(-4.10286e-06*(x*x)))+(2.72413e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.971321+(0.00117532*x))+(-2.71334e-06*(x*x)))+(1.77294e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.0989+(0.0022646*x))+(-5.48834e-06*(x*x)))+(3.67551e-09*(x*(x*x)));
      }
      else if(eta_ >= 0.5 && eta_ < 1.0) {
	if( meanminmax == "mean" ) return ((0.977454+(0.00186222*x))+(-4.30874e-06*(x*x)))+(2.82227e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.917942+(0.00139264*x))+(-3.13422e-06*(x*x)))+(2.02475e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.03695+(0.00232982*x))+(-5.47968e-06*(x*x)))+(3.62048e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.0 && eta_ < 1.5) {
	if( meanminmax == "mean" ) return ((0.940154+(0.00214045*x))+(-5.30206e-06*(x*x)))+(3.75872e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.885078+(0.00170468*x))+(-4.08896e-06*(x*x)))+(2.85628e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((0.995215+(0.00257376*x))+(-6.5103e-06*(x*x)))+(4.66211e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.5 && eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.04882+(0.000373418*x))+(-1.00316e-06*(x*x)))+(8.52325e-10*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.982642+(0.000211816*x))+(-4.11471e-07*(x*x)))+(2.88443e-10*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.11499+(0.000534645*x))+(-1.59409e-06*(x*x)))+(1.41682e-09*(x*(x*x)));
      }
    } // if CSVL

    else if(tagger == "CSVM") {
      if(eta_ < 0.8) {
	if( meanminmax == "mean" ) return ((1.0444+(0.00216756*x))+(-5.4224e-06*(x*x)))+(3.69351e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.966203+(0.00112979*x))+(-2.56147e-06*(x*x)))+(1.65183e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.1226+(0.00320252*x))+(-8.27754e-06*(x*x)))+(5.73519e-09*(x*(x*x)));
      }
      else if(eta_ >= 0.8 && eta_ < 1.6) {
	if( meanminmax == "mean" ) return ((1.05203+(0.00138588*x))+(-3.97677e-06*(x*x)))+(3.13655e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.966774+(0.000855535*x))+(-2.33883e-06*(x*x)))+(1.86063e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.13729+(0.00191432*x))+(-5.61018e-06*(x*x)))+(4.41282e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.6 && eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.06547+(0.000850114*x))+(-2.76694e-06*(x*x)))+(1.75015e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.992673+(0.000455214*x))+(-1.29572e-06*(x*x)))+(3.89704e-10*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.13823+(0.00124422*x))+(-4.23813e-06*(x*x)))+(3.11339e-09*(x*(x*x)));
      }
    } // if CSVM

    else if(tagger == "CSVT") {
      if(eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.00197+(0.00266395*x))+(-6.95018e-06*(x*x)))+(4.91042e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.948887+(0.00103466*x))+(-2.88118e-06*(x*x)))+(2.07782e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.05505+(0.00428961*x))+(-1.10115e-05*(x*x)))+(7.74319e-09*(x*(x*x)));
      }
    } // if CSVT

  } // if SF_12 C

  else if(DataPeriod == "D") {
    
    if(tagger == "CSVL") {
      if(eta_ < 0.5) {
	if( meanminmax == "mean" ) return ((1.1121+(0.00156291*x))+(-3.72267e-06*(x*x)))+(2.54276e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((1.04345+(0.00100049*x))+(-2.27285e-06*(x*x)))+(1.53238e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.18074+(0.00212352*x))+(-5.16888e-06*(x*x)))+(3.55347e-09*(x*(x*x)));
      }
      else if(eta_ >= 0.5 && eta_ < 1.0) {
	if( meanminmax == "mean" ) return ((1.05107+(0.0018085*x))+(-4.42378e-06*(x*x)))+(3.12722e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.986937+(0.00132072*x))+(-3.17261e-06*(x*x)))+(2.25152e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.11519+(0.00229425*x))+(-5.67093e-06*(x*x)))+(4.00366e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.0 && eta_ < 1.5) {
	if( meanminmax == "mean" ) return ((0.984747+(0.00233796*x))+(-5.84283e-06*(x*x)))+(4.21798e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((0.927306+(0.00186598*x))+(-4.5141e-06*(x*x)))+(3.21483e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.04217+(0.0028073*x))+(-7.16639e-06*(x*x)))+(5.2225e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.5 && eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.0944+(0.000394694*x))+(-1.31095e-06*(x*x)))+(1.29556e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((1.0255+(0.000220197*x))+(-6.45505e-07*(x*x)))+(6.40579e-10*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.1633+(0.000568652*x))+(-1.97487e-06*(x*x)))+(1.95111e-09*(x*(x*x)));
      }
    } // if CSVL

    else if(tagger == "CSVM") {
      if(eta_ < 0.8) {
	if( meanminmax == "mean" ) return ((1.13626+(0.00209868*x))+(-5.54303e-06*(x*x)))+(3.9911e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((1.05089+(0.00100001*x))+(-2.44384e-06*(x*x)))+(1.72918e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.22164+(0.00319447*x))+(-8.63596e-06*(x*x)))+(6.25306e-09*(x*(x*x)));
      }
      else if(eta_ >= 0.8 && eta_ < 1.6) {
	if( meanminmax == "mean" ) return ((1.13291+(0.00128329*x))+(-3.79952e-06*(x*x)))+(3.03032e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((1.04112+(0.000728221*x))+(-2.04996e-06*(x*x)))+(1.64537e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.22469+(0.0018366*x))+(-5.54498e-06*(x*x)))+(4.4159e-09*(x*(x*x)));
      }
      else if(eta_ >= 1.6 && eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.18705+(0.000305854*x))+(-1.86925e-06*(x*x)))+(1.79183e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((1.10587+(-8.23503e-05*x))+(-3.06139e-07*(x*x)))+(2.38667e-10*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.26821+(0.000693404*x))+(-3.43071e-06*(x*x)))+(3.34622e-09*(x*(x*x)));
      }
    } // if CSVM

    else if(tagger == "CSVT") {
      if(eta_ < 2.4) {
	if( meanminmax == "mean" ) return ((1.08603+(0.0027994*x))+(-8.44182e-06*(x*x)))+(6.847e-09*(x*(x*x)));
	if( meanminmax == "min" ) return ((1.02837+(0.00104078*x))+(-3.81136e-06*(x*x)))+(3.43028e-09*(x*(x*x)));
	if( meanminmax == "max" ) return ((1.14368+(0.00455363*x))+(-1.30615e-05*(x*x)))+(1.0264e-08*(x*(x*x)));
      }
    } // if CSVT

  } // if SF_12 D

  return 0.;

}

#endif
