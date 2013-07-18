#include "math.h"
#include "TMath.h"

#include "../src/SusyEvent.h"

using namespace std;

float chargedHadronIso_corrected(susy::Photon gamma, float rho) {
  float eta = fabs(gamma.caloPosition.Eta());
  float ea;

  if(eta < 1.0) ea = 0.012;
  else if(eta < 1.479) ea = 0.010;
  else if(eta < 2.0) ea = 0.014;
  else if(eta < 2.2) ea = 0.012;
  else if(eta < 2.3) ea = 0.016;
  else if(eta < 2.4) ea = 0.020;
  else ea = 0.012;

  float iso = gamma.chargedHadronIso;
  iso = max(iso - rho*ea, (float)0.);

  return iso;
}
  
float neutralHadronIso_corrected(susy::Photon gamma, float rho) {
  float eta = fabs(gamma.caloPosition.Eta());
  float ea;

  if(eta < 1.0) ea = 0.030;
  else if(eta < 1.479) ea = 0.057;
  else if(eta < 2.0) ea = 0.039;
  else if(eta < 2.2) ea = 0.015;
  else if(eta < 2.3) ea = 0.024;
  else if(eta < 2.4) ea = 0.039;
  else ea = 0.072;

  float iso = gamma.neutralHadronIso;
  iso = max(iso - rho*ea, (float)0.);

  return iso;
}

float photonIso_corrected(susy::Photon gamma, float rho) {
  float eta = fabs(gamma.caloPosition.Eta());
  float ea;

  if(eta < 1.0) ea = 0.148;
  else if(eta < 1.479) ea = 0.130;
  else if(eta < 2.0) ea = 0.112;
  else if(eta < 2.2) ea = 0.216;
  else if(eta < 2.3) ea = 0.262;
  else if(eta < 2.4) ea = 0.260;
  else ea = 0.266;

  float iso = gamma.photonIso;
  iso = max(iso - rho*ea, (float)0.);

  return iso;
}

bool isLoosePhoton(susy::Photon gamma, float rho) {

  bool loose_barrel = gamma.momentum.Et() > 25.0 &&
    gamma.hadTowOverEm < 0.05 &&
    neutralHadronIso_corrected(gamma, rho) < 3.5 + 0.04*gamma.momentum.Pt() &&
    photonIso_corrected(gamma, rho) < 1.3 + 0.005*gamma.momentum.Pt() &&
    chargedHadronIso_corrected(gamma, rho) < 2.6 &&
    gamma.sigmaIetaIeta < 0.012 &&
    gamma.passelectronveto;

  bool loose_endcap = gamma.momentum.Et() > 25.0 &&
    gamma.hadTowOverEm < 0.05 &&
    neutralHadronIso_corrected(gamma, rho) < 2.9 + 0.04*gamma.momentum.Pt() &&
    chargedHadronIso_corrected(gamma, rho) < 2.3 &&
    gamma.sigmaIetaIeta < 0.034 &&
    gamma.passelectronveto;

  float pho_eta = fabs(gamma.caloPosition.Eta());

  if(pho_eta < 1.4442) return loose_barrel;
  else if(pho_eta > 1.566 && pho_eta < 2.5) return loose_endcap;

  return false;
}

bool is_egf(susy::Photon gamma, float rho) {

  if(fabs(gamma.caloPosition.Eta()) < 1.4442 &&
     gamma.momentum.Et() > 25.0 &&
     gamma.hadTowOverEm < 0.05 &&
     neutralHadronIso_corrected(gamma, rho) < 3.5 + 0.04*gamma.momentum.Pt() &&
     photonIso_corrected(gamma, rho) < 1.3 + 0.005*gamma.momentum.Pt() &&
     chargedHadronIso_corrected(gamma, rho) < 15.0 &&
     gamma.r9 < 1.0 &&
     gamma.sigmaIetaIeta > 0.001 &&
     gamma.sigmaIphiIphi > 0.001) {
    
    return true;

  }
  
  return false;
}

bool is_eg(susy::Photon gamma, float rho) {
  if(is_egf(gamma, rho) &&
     chargedHadronIso_corrected(gamma, rho) < 2.6 &&
     gamma.sigmaIetaIeta < 0.012) {
    return true;
  }
  
  return false;
}

bool is_f(susy::Photon gamma, float rho) {

  if(is_egf(gamma, rho) &&
     gamma.nPixelSeeds == 0 &&
     (gamma.sigmaIetaIeta >= 0.012 || chargedHadronIso_corrected(gamma, rho) >= 2.6 )) {
    return true;
  }
  return false;
}
