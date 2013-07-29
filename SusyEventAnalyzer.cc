#define SusyEventAnalyzer_cxx

#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TObject.h>

#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility>

#include "SusyEventAnalyzer.h"
#include "BtagWeight.h"
#include "EventQuality.h"
#include "Metadata.h"

using namespace std;

bool sortTriggers(pair<TString, int> i, pair<TString, int> j) { return (i.second > j.second); }

void SusyEventAnalyzer::PileupWeights(TString puFile) {

  TFile * in = new TFile(puFile, "READ");
  TH1F * _data = (TH1F*)in->Get("pileup");
  
  TString output_code_t = FormatName(scan);

  TH1F * data = (TH1F*)_data->Clone("pu_data"+output_code_t); data->Sumw2();
  TH1F * mc = new TH1F("pu_mc"+output_code_t, "pu_mc"+output_code_t, 70, 0, 70); mc->Sumw2();
  TH1F * mc_nPVertex = new TH1F("mc_nPVertex"+output_code_t, "mc_nPVertex"+output_code_t, 70, 0, 70);

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    int nPV = -1;
    susy::PUSummaryInfoCollection::const_iterator iBX = event.pu.begin();
    bool foundInTimeBX = false;
    while((iBX != event.pu.end()) && !foundInTimeBX) {
      if(iBX->BX == 0) {
	nPV = iBX->trueNumInteractions;
	foundInTimeBX = true;
      }
      ++iBX;
    }
    
    if(foundInTimeBX) mc->Fill(nPV);

    // Now find the nPV from reconstruction
    int nPV_reco = GetNumberPV(event);
    mc_nPVertex->Fill(nPV_reco);

  } // end event loop

  TH1D * data_nonorm = (TH1D*)data->Clone("pu_data_nonorm"+output_code_t);
  TH1D * mc_nonorm = (TH1D*)mc->Clone("pu_mc_nonorm"+output_code_t);

  Double_t intData = data->Integral();
  Double_t intMC = mc->Integral();

  data->Scale(1./intData);
  mc->Scale(1./intMC);

  TH1F * weights = (TH1F*)data->Clone("puWeights"+output_code_t);
  weights->Divide(mc);

  TFile * out = new TFile("pileupReweighting"+output_code_t+".root", "RECREATE");
  out->cd();

  data->Write();
  data_nonorm->Write();
  mc->Write();
  mc_nonorm->Write();
  weights->Write();
  out->Write();
  out->Close();

  in->Close();

  return;
}

void SusyEventAnalyzer::CalculateBtagEfficiency() {

  const int nChannels = 1;

  const int NCNT = 50;
  int nCnt[NCNT][nChannels];
  for(int i = 0; i < NCNT; i++) {
    for(int j = 0; j < nChannels; j++) {
    nCnt[i][j] = 0;
    }
  }
  
  TString output_code_t = FormatName(scan);

  // open histogram file and define histograms
  TFile * out = new TFile("btagEfficiency"+output_code_t+".root", "RECREATE");
  out->cd();

  TH1F * num_bjets = new TH1F("bjets"+output_code_t, "bjets"+output_code_t, 200, 0, 1000); num_bjets->Sumw2();
  TH1F * num_btags = new TH1F("btags"+output_code_t, "btags"+output_code_t, 200, 0, 1000); num_btags->Sumw2();
  TH1F * num_cjets = new TH1F("cjets"+output_code_t, "cjets"+output_code_t, 200, 0, 1000); num_cjets->Sumw2();
  TH1F * num_ctags = new TH1F("ctags"+output_code_t, "ctags"+output_code_t, 200, 0, 1000); num_ctags->Sumw2();
  TH1F * num_ljets = new TH1F("ljets"+output_code_t, "ljets"+output_code_t, 200, 0, 1000); num_ljets->Sumw2();
  TH1F * num_ltags = new TH1F("ltags"+output_code_t, "ltags"+output_code_t, 200, 0, 1000); num_ltags->Sumw2();

  ScaleFactorInfo sf(btagger);

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = " << event.runNumber << ", event = " << event.eventNumber << endl;
    }

    nCnt[0][0]++; // events

    vector<susy::Photon*> candidate_pair;
    vector<susy::PFJet*> pfJets, btags;
    vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
    vector<float> csvValues;
    vector<susy::Muon*> isoMuons, looseMuons;
    vector<susy::Electron*> isoEles, looseEles;
    vector<BtagInfo> tagInfos;

    int event_type = 0;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    findPhotons_prioritizeCount(event, candidate_pair, event_type);
    //findPhotons_prioritizeEt(event, candidate_pair, event_type);

    if(event_type != 1) {
      nCnt[28][0]++;
      continue;
    }

    float HT = 0.;
    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    findMuons(event, candidate_pair, isoMuons, looseMuons, HT);
    findElectrons(event, candidate_pair, isoEles, looseEles, HT);
    findJets(event, candidate_pair,
	     isoMuons, looseMuons,
	     isoEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT, hadronicSystem);

    HT += candidate_pair[0]->momentum.Pt();
    HT += candidate_pair[1]->momentum.Pt();

    ////////////////////

    bool passHLT = useTrigger ? PassTriggers(4) : true;
    if(!passHLT) {nCnt[36][0]++;continue;}
	
    for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
      map<TString, Float_t>::iterator s_it = pfJets[iJet]->jecScaleFactors.find("L1FastL2L3");
      if(s_it == pfJets[iJet]->jecScaleFactors.end()) {
	continue;
      }
      float scale = s_it->second;
      TLorentzVector corrP4 = scale * pfJets[iJet]->momentum;
      if(fabs(corrP4.Eta()) >= 2.4) continue;
	  
      if(fabs(pfJets[iJet]->algDefFlavour) == 5) {
	num_bjets->Fill(corrP4.Pt());
	if((btagger == "CSVL" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.244) ||
	   (btagger == "CSVM" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.679) ||
	   (btagger == "CSVT" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.898)) 
	  num_btags->Fill(corrP4.Pt());
      }
	  
      if(fabs(pfJets[iJet]->algDefFlavour) == 4) {
	num_cjets->Fill(corrP4.Pt());
	if((btagger == "CSVL" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.244) ||
	   (btagger == "CSVM" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.679) ||
	   (btagger == "CSVT" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.898)) 
	  num_ctags->Fill(corrP4.Pt());
      }
	  
      if(fabs(pfJets[iJet]->algDefFlavour) == 1 ||
	 fabs(pfJets[iJet]->algDefFlavour) == 2 ||
	 fabs(pfJets[iJet]->algDefFlavour) == 3 ||
	 fabs(pfJets[iJet]->algDefFlavour) == 21) {
	num_ljets->Fill(corrP4.Pt());
	if((btagger == "CSVL" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.244) ||
	   (btagger == "CSVM" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.679) ||
	   (btagger == "CSVT" && pfJets[iJet]->bTagDiscriminators[susy::kCSV] > 0.898)) 
	  num_ltags->Fill(corrP4.Pt());
      }
	  
    } // for jets
	
  } // for entries

  TH1F * bEff = (TH1F*)num_btags->Clone("bEff"+output_code_t);
  bEff->Divide(num_bjets);

  TH1F * cEff = (TH1F*)num_ctags->Clone("cEff"+output_code_t);
  cEff->Divide(num_cjets);

  TH1F * lEff = (TH1F*)num_ltags->Clone("lEff"+output_code_t);
  lEff->Divide(num_ljets);

  out->Write();
  out->Close();

}

void SusyEventAnalyzer::Data() {

  TFile* out = new TFile("hist_"+outputName+"_"+btagger+".root", "RECREATE");
  out->cd();

  const int nCategories = 8;
  TString categories[nCategories] = {"gg", "eg", "ff", "ee", "ee_loMass", "ee_onMass", "ee_hiMass", "gf"};
  const int nChannels = 7;

  TString channels[nChannels] = {
    "nojet",
    "j", "b",
    "bj",
    "muJets", // gg+mu+bj + X (dilep veto)
    "eleJets",
    "hadronic" // gg+5j1b + X (lep veto)
  };

  unsigned int nJetReq[nChannels] = {
    0,
    1, 1,
    2,
    2,
    2,
    4};
  
  unsigned int nBtagReq[nChannels] = {
    0,
    0, 1,
    1,
    1,
    1,
    1};
  
  unsigned int nEleReq[nChannels] = {
    0,
    0, 0,
    0,
    0,
    1,
    0};

  unsigned int nMuonReq[nChannels] = {
    0,
    0, 0,
    0,
    1,
    0,
    0};

  bool looseLeptonVeto[nChannels] = {
    false,
    false, false,
    false,
    true,
    true,
    true};
  
  const int NCNT = 50;
  int nCnt[NCNT][nChannels];
  for(int i = 0; i < NCNT; i++) {
    for(int j = 0; j < nChannels; j++) {
      nCnt[i][j] = 0;
    }
  }

  ///////////////////////////////////////////////////
  // Define histograms to be filled for all events
  ///////////////////////////////////////////////////

  TString metFilterNames[susy::nMetFilters] = {
    "CSCBeamHalo",
    "HcalNoise",
    "EcalDeadCellTP",
    "EcalDeadCellBE",
    "TrackingFailure",
    "EEBadSC",
    "HcalLaserOccupancy",
    "HcalLaserEventList",
    "HcalLaserRECOUserStep",
    "EcalLaserCorr",
    "ManyStripClus53X",
    "TooManyStripClus53X",
    "LogErrorTooManyClusters",
    "LogErrorTooManyTripletsPairs",
    "LogErrorTooManySeeds",
    "EERingOfFire",
    "InconsistentMuon",
    "GreedyMuon"};

  TH2F* h_metFilter = new TH2F("metFilter", "MET Filter Failures", susy::nMetFilters, 0, susy::nMetFilters, susy::nMetFilters, 0, susy::nMetFilters);
  for(int i = 0; i < susy::nMetFilters; i++) {
    h_metFilter->GetXaxis()->SetBinLabel(i+1, metFilterNames[i]);
    h_metFilter->GetYaxis()->SetBinLabel(i+1, metFilterNames[i]);
  }

  ///////////////////////////////////////////////////
  // Define histograms
  ///////////////////////////////////////////////////

  // MET
  VTH1F h_met = BookTH1FVector("met", "MET;#slash{E}_{T} (GeV);Events", 400, 0., 2000., nCategories, categories, nChannels, channels);

  // di-em pt
  VTH2F h_diempt = BookTH2FVector("diempt", "di-EM pt vs nJets;diEMPt (GeV/c);nJets", 1000, 0., 2000., 200, 0., 200., nCategories, categories, nChannels, channels);

  // di-jet pt
  VTH2F h_dijetpt = BookTH2FVector("dijetpt", "di-Jet pt vs nJets;diJetPt (GeV/c);nJets", 1000, 0., 2000., 200, 0., 200., nCategories, categories, nChannels, channels);
  
  /////////////////////////////////
  // Miscellaneous histograms
  /////////////////////////////////
  TH2F * h_DR_jet_gg = new TH2F("DR_jet_gg", "#DeltaR between jets and lead/trailing #gamma#gamma candidates;#DeltaR_{lead #gamma, jet};#DeltaR_{trail #gamma, jet}", 50, 0, 5, 50, 0, 5);

  const int nDivisions_chi2 = 50;

  TH1F * h_met_varyCSVcut_ff_j[nDivisions_chi2];
  TH1F * h_met_varyCSVcut_gg_j[nDivisions_chi2];
  TH1F * h_met_squareCSVcut_ff_jj[nDivisions_chi2];
  TH1F * h_met_squareCSVcut_gg_jj[nDivisions_chi2];

  for(int i = 0; i < nDivisions_chi2; i++) {
    char tmp[10];
    double tmp_val = (double)i / (double)nDivisions_chi2;
    sprintf(tmp, "%f", tmp_val);
    TString tmp_t = tmp;

    h_met_varyCSVcut_ff_j[i] = new TH1F("met_varyCSVcut_ff_j_"+tmp_t, "MET for ff+j (maxCSV > "+tmp_t+")", 400, 0, 2000);
    h_met_varyCSVcut_gg_j[i] = new TH1F("met_varyCSVcut_gg_j_"+tmp_t, "MET for gg+j (maxCSV > "+tmp_t+")", 400, 0, 2000);
    h_met_squareCSVcut_ff_jj[i] = new TH1F("met_squareCSVcut_ff_jj_"+tmp_t, "MET for ff+jj (maxCSV & submaxCSV > "+tmp_t+")", 400, 0, 2000);
    h_met_squareCSVcut_gg_jj[i] = new TH1F("met_squareCSVcut_gg_jj_"+tmp_t, "MET for gg+jj (maxCSV & submaxCSV > "+tmp_t+")", 400, 0, 2000);
  }

  /////////////////////////////////
  // Reweighting trees
  /////////////////////////////////

  float pfMET_ = 0.;
  float diEMpT_ = 0.;
  float diJetPt_ = 0.;
  int Njets_ = 0;
  int Nbtags_ = 0;
  int Nelectrons_ = 0;
  int Nmuons_ = 0;
  float invmass_ = 0.;
  float HT_ = 0.;
  float HT_jets_ = 0.;
  float hadronic_pt_ = 0.;
  float minDPhi_gMET_ = 0;
  float minDPhi_jMET_ = 0;
  float lead_Et_ = 0;
  float trail_Et_ = 0;
  float lead_Eta_ = 0;
  float trail_Eta_ = 0;
  float lead_Phi_ = 0;
  float trail_Phi_ = 0;

  float lead_matched_jetpt_ = 0;
  float trail_matched_jetpt_ = 0;

  float leadptOverInvmass_ = 0;
  float trailptOverInvmass_ = 0;

  float lead_chHadIso_ = 0;
  float trail_chHadIso_ = 0;
  float lead_sIetaIeta_ = 0;
  float trail_sIetaIeta_ = 0;

  int nPV_ = 0;
  float photon_dR_ = 0;
  float photon_dPhi_ = 0;

  float jet1_pt_ = -1;
  float jet2_pt_ = -1;
  float jet3_pt_ = -1;
  float jet4_pt_ = -1;

  float btag1_pt_ = -1;
  float btag2_pt_ = -1;

  float max_csv_ = -1;
  float submax_csv_ = -1;
  float min_csv_ = 10;

  float minDR_leadPhoton_jets_ = -1;
  float minDR_trailPhoton_jets_ = -1;

  int runNumber_ = 0;
  ULong_t eventNumber_ = 0;
  int lumiBlock_ = 0;
  Long64_t jentry_ = 0;

  float dimuon_invmass_ = 0.;
  float diphodimu_invmass_ = 0;

  Int_t metFilterBit_ = 0;

  vector<TTree*> ffTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree("ff_"+channels[i]+"_EvtTree", "An event tree for final analysis");
    tree->Branch("pfMET", &pfMET_, "pfMET_/F");
    tree->Branch("diEMpT", &diEMpT_, "diEMpT_/F");
    tree->Branch("diJetPt", &diJetPt_, "diJetPt_/F");
    tree->Branch("Njets", &Njets_, "Njets_/I");
    tree->Branch("Nbtags", &Nbtags_, "Nbtags_/I");
    tree->Branch("Nelectrons", &Nelectrons_, "Nelectrons_/I");
    tree->Branch("Nmuons", &Nmuons_, "Nmuons_/I");
    tree->Branch("invmass", &invmass_, "invmass_/F");
    tree->Branch("HT", &HT_, "HT_/F");
    tree->Branch("HT_jets", &HT_jets_, "HT_jets_/F");
    tree->Branch("hadronic_pt", &hadronic_pt_, "hadronic_pt_/F");
    tree->Branch("minDPhi_gMET", &minDPhi_gMET_, "minDPhi_gMET_/F");
    tree->Branch("minDPhi_jMET", &minDPhi_jMET_, "minDPhi_jMET_/F");
    tree->Branch("minDR_leadPhoton_jets", &minDR_leadPhoton_jets_, "minDR_leadPhoton_jets_/F");
    tree->Branch("minDR_trailPhoton_jets", &minDR_trailPhoton_jets_, "minDR_trailPhoton_jets_/F");
    tree->Branch("leadPhotonEt", &lead_Et_, "lead_Et_/F");
    tree->Branch("trailPhotonEt", &trail_Et_, "trail_Et_/F");
    tree->Branch("leadPhotonEta", &lead_Eta_, "lead_Eta_/F");
    tree->Branch("leadMatchedJetPt", &lead_matched_jetpt_, "lead_matched_jetpt_/F");
    tree->Branch("trailMatchedJetPt", &trail_matched_jetpt_, "trail_matched_jetpt_/F");
    tree->Branch("trailPhotonEta", &trail_Eta_, "trail_Eta_/F");
    tree->Branch("leadPhotonPhi", &lead_Phi_, "lead_Phi_/F");
    tree->Branch("trailPhotonPhi", &trail_Phi_, "trail_Phi_/F");
    tree->Branch("leadptOverInvmass", &leadptOverInvmass_, "leadptOverInvmass_/F");
    tree->Branch("trailptOverInvmass", &trailptOverInvmass_, "trailptOverInvmass_/F");
    tree->Branch("leadChargedHadronIso", &lead_chHadIso_, "lead_chHadIso_/F");
    tree->Branch("trailChargedHadronIso", &trail_chHadIso_, "trail_chHadIso_/F");
    tree->Branch("leadSigmaIetaIeta", &lead_sIetaIeta_, "lead_sIetaIeta_/F");
    tree->Branch("trailSigmaIetaIeta", &trail_sIetaIeta_, "trail_sIetaIeta_/F");
    tree->Branch("jet1_pt", &jet1_pt_, "jet1_pt_/F");
    tree->Branch("jet2_pt", &jet2_pt_, "jet2_pt_/F");
    tree->Branch("jet3_pt", &jet3_pt_, "jet3_pt_/F");
    tree->Branch("jet4_pt", &jet4_pt_, "jet4_pt_/F");
    tree->Branch("max_csv", &max_csv_, "max_csv_/F");
    tree->Branch("btag1_pt", &btag1_pt_, "btag1_pt_/F");
    tree->Branch("btag2_pt", &btag2_pt_, "btag2_pt_/F");
    tree->Branch("submax_csv", &submax_csv_, "submax_csv_/F");
    tree->Branch("min_csv", &min_csv_, "min_csv_/F");
    tree->Branch("nPV", &nPV_, "nPV_/I");
    tree->Branch("photon_dR", &photon_dR_, "photon_dR_/F");
    tree->Branch("photon_dPhi", &photon_dPhi_, "photon_dPhi_/F");
    tree->Branch("runNumber", &runNumber_, "runNumber_/I");
    tree->Branch("eventNumber", &eventNumber_, "eventNumber_/l");
    tree->Branch("luminosityBlockNumber", &lumiBlock_, "lumiBlock_/I");
    tree->Branch("jentry", &jentry_, "jentry_/L");
    tree->Branch("dimuon_invmass", &dimuon_invmass_, "dimuon_invmass_/F");
    tree->Branch("diphodimu_invmass", &diphodimu_invmass_, "diphodimu_invmass_/F");
    tree->Branch("metFilterBit", &metFilterBit_, "metFilterBit_/I");

    ffTrees.push_back(tree);
  }

  vector<TTree*> eeTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree("ee_"+channels[i]+"_EvtTree", "An event tree for final analysis");
    tree->Branch("pfMET", &pfMET_, "pfMET_/F");
    tree->Branch("diEMpT", &diEMpT_, "diEMpT_/F");
    tree->Branch("diJetPt", &diJetPt_, "diJetPt_/F");
    tree->Branch("Njets", &Njets_, "Njets_/I");
    tree->Branch("Nbtags", &Nbtags_, "Nbtags_/I");
    tree->Branch("Nelectrons", &Nelectrons_, "Nelectrons_/I");
    tree->Branch("Nmuons", &Nmuons_, "Nmuons_/I");
    tree->Branch("invmass", &invmass_, "invmass_/F");
    tree->Branch("HT", &HT_, "HT_/F");
    tree->Branch("HT_jets", &HT_jets_, "HT_jets_/F");
    tree->Branch("hadronic_pt", &hadronic_pt_, "hadronic_pt_/F");
    tree->Branch("minDPhi_gMET", &minDPhi_gMET_, "minDPhi_gMET_/F");
    tree->Branch("minDPhi_jMET", &minDPhi_jMET_, "minDPhi_jMET_/F");
    tree->Branch("minDR_leadPhoton_jets", &minDR_leadPhoton_jets_, "minDR_leadPhoton_jets_/F");
    tree->Branch("minDR_trailPhoton_jets", &minDR_trailPhoton_jets_, "minDR_trailPhoton_jets_/F");
    tree->Branch("leadPhotonEt", &lead_Et_, "lead_Et_/F");
    tree->Branch("trailPhotonEt", &trail_Et_, "trail_Et_/F");
    tree->Branch("leadMatchedJetPt", &lead_matched_jetpt_, "lead_matched_jetpt_/F");
    tree->Branch("trailMatchedJetPt", &trail_matched_jetpt_, "trail_matched_jetpt_/F");
    tree->Branch("leadPhotonEta", &lead_Eta_, "lead_Eta_/F");
    tree->Branch("trailPhotonEta", &trail_Eta_, "trail_Eta_/F");
    tree->Branch("leadPhotonPhi", &lead_Phi_, "lead_Phi_/F");
    tree->Branch("trailPhotonPhi", &trail_Phi_, "trail_Phi_/F");
    tree->Branch("leadptOverInvmass", &leadptOverInvmass_, "leadptOverInvmass_/F");
    tree->Branch("trailptOverInvmass", &trailptOverInvmass_, "trailptOverInvmass_/F");
    tree->Branch("leadChargedHadronIso", &lead_chHadIso_, "lead_chHadIso_/F");
    tree->Branch("trailChargedHadronIso", &trail_chHadIso_, "trail_chHadIso_/F");
    tree->Branch("leadSigmaIetaIeta", &lead_sIetaIeta_, "lead_sIetaIeta_/F");
    tree->Branch("trailSigmaIetaIeta", &trail_sIetaIeta_, "trail_sIetaIeta_/F");
    tree->Branch("jet1_pt", &jet1_pt_, "jet1_pt_/F");
    tree->Branch("jet2_pt", &jet2_pt_, "jet2_pt_/F");
    tree->Branch("jet3_pt", &jet3_pt_, "jet3_pt_/F");
    tree->Branch("jet4_pt", &jet4_pt_, "jet4_pt_/F");
    tree->Branch("btag1_pt", &btag1_pt_, "btag1_pt_/F");
    tree->Branch("btag2_pt", &btag2_pt_, "btag2_pt_/F");
    tree->Branch("max_csv", &max_csv_, "max_csv_/F");
    tree->Branch("submax_csv", &submax_csv_, "submax_csv_/F");
    tree->Branch("min_csv", &min_csv_, "min_csv_/F");
    tree->Branch("nPV", &nPV_, "nPV_/I");
    tree->Branch("photon_dR", &photon_dR_, "photon_dR_/F");
    tree->Branch("photon_dPhi", &photon_dPhi_, "photon_dPhi_/F");
    tree->Branch("runNumber", &runNumber_, "runNumber_/I");
    tree->Branch("eventNumber", &eventNumber_, "eventNumber_/l");
    tree->Branch("luminosityBlockNumber", &lumiBlock_, "lumiBlock_/I");
    tree->Branch("jentry", &jentry_, "jentry_/L");
    tree->Branch("dimuon_invmass", &dimuon_invmass_, "dimuon_invmass_/F");
    tree->Branch("diphodimu_invmass", &diphodimu_invmass_, "diphodimu_invmass_/F");
    tree->Branch("metFilterBit", &metFilterBit_, "metFilterBit_/I");

    eeTrees.push_back(tree);
  }

  vector<TTree*> gfTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree("gf_"+channels[i]+"_EvtTree", "An event tree for final analysis");
    tree->Branch("pfMET", &pfMET_, "pfMET_/F");
    tree->Branch("diEMpT", &diEMpT_, "diEMpT_/F");
    tree->Branch("diJetPt", &diJetPt_, "diJetPt_/F");
    tree->Branch("Njets", &Njets_, "Njets_/I");
    tree->Branch("Nbtags", &Nbtags_, "Nbtags_/I");
    tree->Branch("Nelectrons", &Nelectrons_, "Nelectrons_/I");
    tree->Branch("Nmuons", &Nmuons_, "Nmuons_/I");
    tree->Branch("invmass", &invmass_, "invmass_/F");
    tree->Branch("HT", &HT_, "HT_/F");
    tree->Branch("HT_jets", &HT_jets_, "HT_jets_/F");
    tree->Branch("hadronic_pt", &hadronic_pt_, "hadronic_pt_/F");
    tree->Branch("minDPhi_gMET", &minDPhi_gMET_, "minDPhi_gMET_/F");
    tree->Branch("minDPhi_jMET", &minDPhi_jMET_, "minDPhi_jMET_/F");
    tree->Branch("minDR_leadPhoton_jets", &minDR_leadPhoton_jets_, "minDR_leadPhoton_jets_/F");
    tree->Branch("minDR_trailPhoton_jets", &minDR_trailPhoton_jets_, "minDR_trailPhoton_jets_/F");
    tree->Branch("leadPhotonEt", &lead_Et_, "lead_Et_/F");
    tree->Branch("trailPhotonEt", &trail_Et_, "trail_Et_/F");
    tree->Branch("leadMatchedJetPt", &lead_matched_jetpt_, "lead_matched_jetpt_/F");
    tree->Branch("trailMatchedJetPt", &trail_matched_jetpt_, "trail_matched_jetpt_/F");
    tree->Branch("leadPhotonEta", &lead_Eta_, "lead_Eta_/F");
    tree->Branch("trailPhotonEta", &trail_Eta_, "trail_Eta_/F");
    tree->Branch("leadPhotonPhi", &lead_Phi_, "lead_Phi_/F");
    tree->Branch("trailPhotonPhi", &trail_Phi_, "trail_Phi_/F");
    tree->Branch("leadptOverInvmass", &leadptOverInvmass_, "leadptOverInvmass_/F");
    tree->Branch("trailptOverInvmass", &trailptOverInvmass_, "trailptOverInvmass_/F");
    tree->Branch("leadChargedHadronIso", &lead_chHadIso_, "lead_chHadIso_/F");
    tree->Branch("trailChargedHadronIso", &trail_chHadIso_, "trail_chHadIso_/F");
    tree->Branch("leadSigmaIetaIeta", &lead_sIetaIeta_, "lead_sIetaIeta_/F");
    tree->Branch("trailSigmaIetaIeta", &trail_sIetaIeta_, "trail_sIetaIeta_/F");
    tree->Branch("jet1_pt", &jet1_pt_, "jet1_pt_/F");
    tree->Branch("jet2_pt", &jet2_pt_, "jet2_pt_/F");
    tree->Branch("jet3_pt", &jet3_pt_, "jet3_pt_/F");
    tree->Branch("jet4_pt", &jet4_pt_, "jet4_pt_/F");
    tree->Branch("btag1_pt", &btag1_pt_, "btag1_pt_/F");
    tree->Branch("btag2_pt", &btag2_pt_, "btag2_pt_/F");
    tree->Branch("max_csv", &max_csv_, "max_csv_/F");
    tree->Branch("submax_csv", &submax_csv_, "submax_csv_/F");
    tree->Branch("min_csv", &min_csv_, "min_csv_/F");
    tree->Branch("nPV", &nPV_, "nPV_/I");
    tree->Branch("photon_dR", &photon_dR_, "photon_dR_/F");
    tree->Branch("photon_dPhi", &photon_dPhi_, "photon_dPhi_/F");
    tree->Branch("runNumber", &runNumber_, "runNumber_/I");
    tree->Branch("eventNumber", &eventNumber_, "eventNumber_/l");
    tree->Branch("luminosityBlockNumber", &lumiBlock_, "lumiBlock_/I");
    tree->Branch("jentry", &jentry_, "jentry_/L");
    tree->Branch("dimuon_invmass", &dimuon_invmass_, "dimuon_invmass_/F");
    tree->Branch("diphodimu_invmass", &diphodimu_invmass_, "diphodimu_invmass_/F");
    tree->Branch("metFilterBit", &metFilterBit_, "metFilterBit_/I");

    gfTrees.push_back(tree);
  }

  vector<TTree*> ggTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree("gg_"+channels[i]+"_EvtTree", "An event tree for final analysis");
    tree->Branch("pfMET", &pfMET_, "pfMET_/F");
    tree->Branch("diEMpT", &diEMpT_, "diEMpT_/F");
    tree->Branch("diJetPt", &diJetPt_, "diJetPt_/F");
    tree->Branch("Njets", &Njets_, "Njets_/I");
    tree->Branch("Nbtags", &Nbtags_, "Nbtags_/I");
    tree->Branch("Nelectrons", &Nelectrons_, "Nelectrons_/I");
    tree->Branch("Nmuons", &Nmuons_, "Nmuons_/I");
    tree->Branch("invmass", &invmass_, "invmass_/F");
    tree->Branch("HT", &HT_, "HT_/F");
    tree->Branch("HT_jets", &HT_jets_, "HT_jets_/F");
    tree->Branch("hadronic_pt", &hadronic_pt_, "hadronic_pt_/F");
    tree->Branch("minDPhi_gMET", &minDPhi_gMET_, "minDPhi_gMET_/F");
    tree->Branch("minDPhi_jMET", &minDPhi_jMET_, "minDPhi_jMET_/F");
    tree->Branch("minDR_leadPhoton_jets", &minDR_leadPhoton_jets_, "minDR_leadPhoton_jets_/F");
    tree->Branch("minDR_trailPhoton_jets", &minDR_trailPhoton_jets_, "minDR_trailPhoton_jets_/F");
    tree->Branch("leadPhotonEt", &lead_Et_, "lead_Et_/F");
    tree->Branch("trailPhotonEt", &trail_Et_, "trail_Et_/F");
    tree->Branch("leadMatchedJetPt", &lead_matched_jetpt_, "lead_matched_jetpt_/F");
    tree->Branch("trailMatchedJetPt", &trail_matched_jetpt_, "trail_matched_jetpt_/F");
    tree->Branch("leadPhotonEta", &lead_Eta_, "lead_Eta_/F");
    tree->Branch("trailPhotonEta", &trail_Eta_, "trail_Eta_/F");
    tree->Branch("leadPhotonPhi", &lead_Phi_, "lead_Phi_/F");
    tree->Branch("trailPhotonPhi", &trail_Phi_, "trail_Phi_/F");
    tree->Branch("leadptOverInvmass", &leadptOverInvmass_, "leadptOverInvmass_/F");
    tree->Branch("trailptOverInvmass", &trailptOverInvmass_, "trailptOverInvmass_/F");
    tree->Branch("leadChargedHadronIso", &lead_chHadIso_, "lead_chHadIso_/F");
    tree->Branch("trailChargedHadronIso", &trail_chHadIso_, "trail_chHadIso_/F");
    tree->Branch("leadSigmaIetaIeta", &lead_sIetaIeta_, "lead_sIetaIeta_/F");
    tree->Branch("trailSigmaIetaIeta", &trail_sIetaIeta_, "trail_sIetaIeta_/F");
    tree->Branch("jet1_pt", &jet1_pt_, "jet1_pt_/F");
    tree->Branch("jet2_pt", &jet2_pt_, "jet2_pt_/F");
    tree->Branch("jet3_pt", &jet3_pt_, "jet3_pt_/F");
    tree->Branch("jet4_pt", &jet4_pt_, "jet4_pt_/F");
    tree->Branch("btag1_pt", &btag1_pt_, "btag1_pt_/F");
    tree->Branch("btag2_pt", &btag2_pt_, "btag2_pt_/F");
    tree->Branch("max_csv", &max_csv_, "max_csv_/F");
    tree->Branch("submax_csv", &submax_csv_, "submax_csv_/F");
    tree->Branch("min_csv", &min_csv_, "min_csv_/F");
    tree->Branch("nPV", &nPV_, "nPV_/I");
    tree->Branch("photon_dR", &photon_dR_, "photon_dR_/F");
    tree->Branch("photon_dPhi", &photon_dPhi_, "photon_dPhi_/F");
    tree->Branch("runNumber", &runNumber_, "runNumber_/I");
    tree->Branch("eventNumber", &eventNumber_, "eventNumber_/l");
    tree->Branch("luminosityBlockNumber", &lumiBlock_, "lumiBlock_/I");
    tree->Branch("jentry", &jentry_, "jentry_/L");
    tree->Branch("dimuon_invmass", &dimuon_invmass_, "dimuon_invmass_/F");
    tree->Branch("diphodimu_invmass", &diphodimu_invmass_, "diphodimu_invmass_/F");
    tree->Branch("metFilterBit", &metFilterBit_, "metFilterBit_/I");

    ggTrees.push_back(tree);
  }

  vector<TTree*> egTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree("eg_"+channels[i]+"_EvtTree", "An event tree for final analysis");
    tree->Branch("pfMET", &pfMET_, "pfMET_/F");
    tree->Branch("diEMpT", &diEMpT_, "diEMpT_/F");
    tree->Branch("diJetPt", &diJetPt_, "diJetPt_/F");
    tree->Branch("Njets", &Njets_, "Njets_/I");
    tree->Branch("Nbtags", &Nbtags_, "Nbtags_/I");
    tree->Branch("Nelectrons", &Nelectrons_, "Nelectrons_/I");
    tree->Branch("Nmuons", &Nmuons_, "Nmuons_/I");
    tree->Branch("invmass", &invmass_, "invmass_/F");
    tree->Branch("HT", &HT_, "HT_/F");
    tree->Branch("HT_jets", &HT_jets_, "HT_jets_/F");
    tree->Branch("hadronic_pt", &hadronic_pt_, "hadronic_pt_/F");
    tree->Branch("minDPhi_gMET", &minDPhi_gMET_, "minDPhi_gMET_/F");
    tree->Branch("minDPhi_jMET", &minDPhi_jMET_, "minDPhi_jMET_/F");
    tree->Branch("minDR_leadPhoton_jets", &minDR_leadPhoton_jets_, "minDR_leadPhoton_jets_/F");
    tree->Branch("minDR_trailPhoton_jets", &minDR_trailPhoton_jets_, "minDR_trailPhoton_jets_/F");
    tree->Branch("leadPhotonEt", &lead_Et_, "lead_Et_/F");
    tree->Branch("trailPhotonEt", &trail_Et_, "trail_Et_/F");
    tree->Branch("leadMatchedJetPt", &lead_matched_jetpt_, "lead_matched_jetpt_/F");
    tree->Branch("trailMatchedJetPt", &trail_matched_jetpt_, "trail_matched_jetpt_/F");
    tree->Branch("leadPhotonEta", &lead_Eta_, "lead_Eta_/F");
    tree->Branch("trailPhotonEta", &trail_Eta_, "trail_Eta_/F");
    tree->Branch("leadPhotonPhi", &lead_Phi_, "lead_Phi_/F");
    tree->Branch("trailPhotonPhi", &trail_Phi_, "trail_Phi_/F");
    tree->Branch("leadptOverInvmass", &leadptOverInvmass_, "leadptOverInvmass_/F");
    tree->Branch("trailptOverInvmass", &trailptOverInvmass_, "trailptOverInvmass_/F");
    tree->Branch("leadChargedHadronIso", &lead_chHadIso_, "lead_chHadIso_/F");
    tree->Branch("trailChargedHadronIso", &trail_chHadIso_, "trail_chHadIso_/F");
    tree->Branch("leadSigmaIetaIeta", &lead_sIetaIeta_, "lead_sIetaIeta_/F");
    tree->Branch("trailSigmaIetaIeta", &trail_sIetaIeta_, "trail_sIetaIeta_/F");
    tree->Branch("jet1_pt", &jet1_pt_, "jet1_pt_/F");
    tree->Branch("jet2_pt", &jet2_pt_, "jet2_pt_/F");
    tree->Branch("jet3_pt", &jet3_pt_, "jet3_pt_/F");
    tree->Branch("jet4_pt", &jet4_pt_, "jet4_pt_/F");
    tree->Branch("btag1_pt", &btag1_pt_, "btag1_pt_/F");
    tree->Branch("btag2_pt", &btag2_pt_, "btag2_pt_/F");
    tree->Branch("max_csv", &max_csv_, "max_csv_/F");
    tree->Branch("submax_csv", &submax_csv_, "submax_csv_/F");
    tree->Branch("min_csv", &min_csv_, "min_csv_/F");
    tree->Branch("nPV", &nPV_, "nPV_/I");
    tree->Branch("photon_dR", &photon_dR_, "photon_dR_/F");
    tree->Branch("photon_dPhi", &photon_dPhi_, "photon_dPhi_/F");
    tree->Branch("runNumber", &runNumber_, "runNumber_/I");
    tree->Branch("eventNumber", &eventNumber_, "eventNumber_/l");
    tree->Branch("luminosityBlockNumber", &lumiBlock_, "lumiBlock_/I");
    tree->Branch("jentry", &jentry_, "jentry_/L");
    tree->Branch("dimuon_invmass", &dimuon_invmass_, "dimuon_invmass_/F");
    tree->Branch("diphodimu_invmass", &diphodimu_invmass_, "diphodimu_invmass_/F");
    tree->Branch("metFilterBit", &metFilterBit_, "metFilterBit_/I");

    egTrees.push_back(tree);
  }

  ScaleFactorInfo sf(btagger);

  // to check duplicate events
  map<int, set<int> > allEvents;

  bool quitAfterProcessing = false;

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    jentry_ = jentry;

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = " << event.runNumber << ", event = " << event.eventNumber << endl;
    }
    
    if(useSyncFile) {
      bool sync = false;
      for(unsigned int i = 0; i < syncRuns.size(); i++) {
	//if(event.runNumber == syncRuns[i] && event.luminosityBlockNumber == syncLumi[i] && event.eventNumber == syncEvents[i]) {
	if(event.runNumber == syncRuns[i] && event.eventNumber == syncEvents[i]) {
	  sync = true;
	  //Print(*event);
	  break;
	}
      }
      if(!sync) continue;

      //if(nCnt[0][0] == (syncRuns.size() - 1)) quitAfterProcessing = true;
    }

    if(singleEvent) {
      if(event.runNumber != single_run || event.luminosityBlockNumber != single_lumi || event.eventNumber != single_event) continue;
      //Print(event);
      quitAfterProcessing = true;
    }

    FillMetFilter2D(event, h_metFilter);

    nCnt[0][0]++; // events

    if(useJson && event.isRealData && !IsGoodLumi(event.runNumber, event.luminosityBlockNumber)) continue;
    nCnt[1][0]++;

    if(event.isRealData) {
      if(event.passMetFilters() != 1 ||
	 event.passMetFilter(susy::kEcalLaserCorr) != 1 ||
	 event.passMetFilter(susy::kManyStripClus53X) != 1 ||
	 event.passMetFilter(susy::kTooManyStripClus53X) != 1) {
	nCnt[21][0]++;
	continue;
      }
    }

    vector<susy::Photon*> candidate_pair;
    vector<susy::PFJet*> pfJets, btags;
    vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
    vector<float> csvValues;
    vector<susy::Muon*> isoMuons, looseMuons;
    vector<susy::Electron*> isoEles, looseEles;
    vector<BtagInfo> tagInfos;

    int event_type = 0;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    map<TString, susy::MET>::iterator met_it = event.metMap.find("pfMet");
    susy::MET* pfMet = &(met_it->second);

    findPhotons_prioritizeCount(event, candidate_pair, event_type);
    //findPhotons_prioritizeEt(event, candidate_pair, event_type);

    if(event_type == 0) {
      nCnt[28][0]++;
      continue;
    }

    bool duplicateEvent = ! (allEvents[event.runNumber].insert(event.eventNumber)).second;
    if(event.isRealData && duplicateEvent) continue;

    lead_Et_ = candidate_pair[0]->momentum.Et();
    lead_Eta_ = candidate_pair[0]->caloPosition.Eta();
    lead_Phi_ = candidate_pair[0]->caloPosition.Phi();
    trail_Et_ = candidate_pair[1]->momentum.Et();
    trail_Eta_ = candidate_pair[1]->caloPosition.Eta();
    trail_Phi_ = candidate_pair[1]->caloPosition.Phi();
    
    metFilterBit_ = event.metFilterBit;

    nPV_ = nPVertex;
    float dEta_ = candidate_pair[0]->caloPosition.Eta() - candidate_pair[1]->caloPosition.Eta();
    photon_dPhi_ = TVector2::Phi_mpi_pi(candidate_pair[0]->caloPosition.Phi() - candidate_pair[1]->caloPosition.Phi());
    photon_dR_ = sqrt(dEta_*dEta_ + photon_dPhi_*photon_dPhi_);
    runNumber_ = event.runNumber;
    eventNumber_ = event.eventNumber;
    lumiBlock_ = event.luminosityBlockNumber;

    leadptOverInvmass_ = candidate_pair[0]->momentum.Pt() / ((candidate_pair[0]->momentum + candidate_pair[1]->momentum).M());
    trailptOverInvmass_ = candidate_pair[1]->momentum.Pt() / ((candidate_pair[0]->momentum + candidate_pair[1]->momentum).M());

    lead_chHadIso_ = chargedHadronIso_corrected(*candidate_pair[0], event.rho25);
    trail_chHadIso_ = chargedHadronIso_corrected(*candidate_pair[1], event.rho25);

    lead_sIetaIeta_ = candidate_pair[0]->sigmaIetaIeta;
    trail_sIetaIeta_ = candidate_pair[1]->sigmaIetaIeta;

    float HT = 0.;
    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    // need to change these to also find veto leptons
    findMuons(event, candidate_pair, isoMuons, looseMuons, HT);
    findElectrons(event, candidate_pair, isoEles, looseEles, HT);

    findJets(event, candidate_pair, 
	     isoMuons, looseMuons,
	     isoEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT, hadronicSystem);

    for(unsigned int i = 0; i < pfJets_corrP4.size(); i++) h_DR_jet_gg->Fill(deltaR(pfJets_corrP4[i], candidate_pair[0]->caloPosition), deltaR(pfJets_corrP4[i], candidate_pair[1]->caloPosition));

    HT_jets_ = HT;
    hadronic_pt_ = hadronicSystem.Pt();

    max_csv_ = (csvValues.size() >= 1) ? csvValues[0] : -1.;
    submax_csv_ = (csvValues.size() >= 2) ? csvValues[1] : -1.;
    min_csv_ = (csvValues.size() >= 1) ? csvValues.back() : -1.;

    HT += candidate_pair[0]->momentum.Pt();
    HT += candidate_pair[1]->momentum.Pt();

    jet1_pt_ = (pfJets_corrP4.size() >= 1) ? pfJets_corrP4[0].Pt() : -1.;
    jet2_pt_ = (pfJets_corrP4.size() >= 2) ? pfJets_corrP4[1].Pt() : -1.;
    jet3_pt_ = (pfJets_corrP4.size() >= 3) ? pfJets_corrP4[2].Pt() : -1.;
    jet4_pt_ = (pfJets_corrP4.size() >= 4) ? pfJets_corrP4[3].Pt() : -1.;

    btag1_pt_ = (btags_corrP4.size() >= 1) ? btags_corrP4[0].Pt() : -1.;
    btag2_pt_ = (btags_corrP4.size() >= 2) ? btags_corrP4[1].Pt() : -1.;
    
    // Calculate dPhi_min(g, MET)
    float dPhi_gMET_lead = TVector2::Phi_mpi_pi(candidate_pair[0]->caloPosition.Phi() - pfMet->mEt.Phi());
    float dPhi_gMET_trail = TVector2::Phi_mpi_pi(candidate_pair[1]->caloPosition.Phi() - pfMet->mEt.Phi());
    minDPhi_gMET_ = min(fabs(dPhi_gMET_lead), fabs(dPhi_gMET_trail));

    float min_deltaPhi_jetsMET = 999.;
    for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
      float dPhi_jetMET = fabs(TVector2::Phi_mpi_pi(pfJets[iJet]->momentum.Phi() - pfMet->mEt.Phi()));
      if(dPhi_jetMET < min_deltaPhi_jetsMET) min_deltaPhi_jetsMET = dPhi_jetMET;
    }
    minDPhi_jMET_ = (pfJets.size() > 0) ? min_deltaPhi_jetsMET : -1.;

    float min_deltaR_lead_jet = 999.;
    for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
      float dEta_x = candidate_pair[0]->caloPosition.Eta() - pfJets[iJet]->momentum.Eta();
      float dPhi_x = TVector2::Phi_mpi_pi(candidate_pair[0]->caloPosition.Phi() - pfJets[iJet]->momentum.Phi());
      float dR_x = sqrt(dEta_x*dEta_x + dPhi_x*dPhi_x);

      if(dR_x < min_deltaR_lead_jet) min_deltaR_lead_jet = dR_x;
    }
    minDR_leadPhoton_jets_ = (pfJets.size() > 0) ? min_deltaR_lead_jet : -1.;

    float min_deltaR_trail_jet = 999.;
    for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
      float dEta_x = candidate_pair[1]->caloPosition.Eta() - pfJets[iJet]->momentum.Eta();
      float dPhi_x = TVector2::Phi_mpi_pi(candidate_pair[1]->caloPosition.Phi() - pfJets[iJet]->momentum.Phi());
      float dR_x = sqrt(dEta_x*dEta_x + dPhi_x*dPhi_x);

      if(dR_x < min_deltaR_trail_jet) min_deltaR_trail_jet = dR_x;
    }
    minDR_trailPhoton_jets_ = (pfJets.size() > 0) ? min_deltaR_trail_jet : -1.;
      
    ////////////////////

    if(event_type != 0) {

      for(int chan = 0; chan < nChannels; chan++) {

	if(pfJets.size() < nJetReq[chan]) continue;
	if(btags.size() < nBtagReq[chan]) continue;
	if(looseLeptonVeto[chan]) {
	  if(isoEles.size() != nEleReq[chan]) continue;
	  if(looseEles.size() != 0) continue;
	  if(isoMuons.size() != nMuonReq[chan]) continue;
	  if(looseMuons.size() != 0) continue;
	}

	if(event_type == 1) {
	
	  bool passHLT = useTrigger ? PassTriggers(1) : true;
	  if(!passHLT) {nCnt[31][chan]++;continue;}
	
	  h_diempt[0][chan]->Fill((candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt(), pfJets.size());
	  h_met[0][chan]->Fill(pfMet->met());

	  float diJetPt = 0.;
	  bool matchingWorked = GetDiJetPt(event, candidate_pair, diJetPt, lead_matched_jetpt_, trail_matched_jetpt_);
	  if(!matchingWorked) nCnt[46][chan]++;

	  h_dijetpt[0][chan]->Fill(diJetPt, pfJets.size());
		  
	  diJetPt_ = diJetPt;
	  
	  pfMET_ = pfMet->met();
	  diEMpT_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt();
	  Njets_ = pfJets.size();
	  Nbtags_ = btags.size();
	  Nelectrons_ = isoEles.size();
	  Nmuons_ = isoMuons.size();
	  invmass_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).M();
	  HT_ = HT;

	  ggTrees[chan]->Fill();

	  if(chan == 1) {
	    for(int iCut = 0; iCut < nDivisions_chi2; iCut++) {
	      double cut_val = (double)iCut / (double)nDivisions_chi2;
	      bool wouldPass = false;
    
	      for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
		if(pfJets[iJet]->bTagDiscriminators[susy::kCSV] > cut_val) {
		  wouldPass = true;
		  break;
		}
	      }
	      if(wouldPass) h_met_varyCSVcut_gg_j[iCut]->Fill(pfMet->met());

	    }
	  }

	  if(chan == 4) {
	    for(int iCut = 0; iCut < nDivisions_chi2; iCut++) {
	      double cut_val = (double)iCut / (double)nDivisions_chi2;
	      int nPassing = 0;
    
	      for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
		if(pfJets[iJet]->bTagDiscriminators[susy::kCSV] > cut_val) nPassing++;
	      }
	      if(nPassing >= 2) h_met_squareCSVcut_gg_jj[iCut]->Fill(pfMet->met());

	    }
	  }
	      
	  nCnt[2][chan]++;
	  
	} // if gg event

	else if(event_type == 2) {

	  bool passHLT = useTrigger ? PassTriggers(2) : true;
	  if(!passHLT) {nCnt[32][chan]++;continue;}

	  h_met[1][chan]->Fill(pfMet->met());

	  float diJetPt = 0.;
	  bool matchingWorked = GetDiJetPt(event, candidate_pair, diJetPt, lead_matched_jetpt_, trail_matched_jetpt_);
	  if(!matchingWorked) nCnt[46][chan]++;

	  diJetPt_ = diJetPt;
	  pfMET_ = pfMet->met();
	  diEMpT_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt();
	  Njets_ = pfJets.size();
	  Nbtags_ = btags.size();
	  Nelectrons_ = isoEles.size();
	  Nmuons_ = isoMuons.size();
	  invmass_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).M();
	  HT_ = HT;

	  egTrees[chan]->Fill();

	  nCnt[3][chan]++;
	  
	} // if ge event

	else if(event_type == -2) {
	  
	  bool passHLT = useTrigger ? PassTriggers(2) : true;
	  if(!passHLT) {nCnt[33][chan]++;continue;}
	  
	  h_met[1][chan]->Fill(pfMet->met());

	  float diJetPt = 0.;
	  bool matchingWorked = GetDiJetPt(event, candidate_pair, diJetPt, lead_matched_jetpt_, trail_matched_jetpt_);
	  if(!matchingWorked) nCnt[46][chan]++;

	  diJetPt_ = diJetPt;
	  pfMET_ = pfMet->met();
	  diEMpT_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt();
	  Njets_ = pfJets.size();
	  Nbtags_ = btags.size();
	  Nelectrons_ = isoEles.size();
	  Nmuons_ = isoMuons.size();
	  invmass_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).M();
	  HT_ = HT;

	  egTrees[chan]->Fill();

	  nCnt[3][chan]++;
	  
	} // if eg event

	else if(event_type == 3) {
	  
	  bool passHLT = useTrigger ? PassTriggers(2) : true;
	  if(!passHLT) {nCnt[34][chan]++;continue;}
	  
	  float zmass_ee = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).M();
	  h_diempt[3][chan]->Fill((candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt(), pfJets.size());
	  h_met[3][chan]->Fill(pfMet->met());
	  
	  float diJetPt = 0.;
	  bool matchingWorked = GetDiJetPt(event, candidate_pair, diJetPt, lead_matched_jetpt_, trail_matched_jetpt_);
	  if(!matchingWorked) nCnt[46][chan]++;

	  h_dijetpt[3][chan]->Fill(diJetPt, pfJets.size());

	  if(zmass_ee > 81.0 && zmass_ee < 101.0) {
	    h_met[5][chan]->Fill(pfMet->met());
	    h_diempt[5][chan]->Fill((candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt(), pfJets.size());
	    h_dijetpt[5][chan]->Fill(diJetPt, pfJets.size());
	  }
	  else if(zmass_ee > 71.0 && zmass_ee < 81.0) {
	    h_met[4][chan]->Fill(pfMet->met());
	    h_diempt[4][chan]->Fill((candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt(), pfJets.size());
	    h_dijetpt[4][chan]->Fill(diJetPt, pfJets.size());
	  }
	  else if(zmass_ee > 101.0 && zmass_ee < 111.0) {
	    h_met[6][chan]->Fill(pfMet->met());
	    h_diempt[6][chan]->Fill((candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt(), pfJets.size());
	    h_dijetpt[6][chan]->Fill(diJetPt, pfJets.size());
	  }
	  
	  diJetPt_ = diJetPt;
	  
	  pfMET_ = pfMet->met();
	  diEMpT_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt();
	  Njets_ = pfJets.size();
	  Nbtags_ = btags.size();
	  Nelectrons_ = isoEles.size();
	  Nmuons_ = isoMuons.size();
	  invmass_ = zmass_ee;
	  HT_ = HT;

	  eeTrees[chan]->Fill();

	  nCnt[4][chan]++;
	  if(zmass_ee > 81.0 && zmass_ee < 101.0) {
	    nCnt[14][chan]++;
	  }
	  
	} // if ee event
	
	else if(event_type == 4) {
	  
	  bool passHLT = useTrigger ? PassTriggers(4) : true;
	  if(!passHLT) {nCnt[35][chan]++;continue;}
	  
	  h_diempt[2][chan]->Fill((candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt(), pfJets.size());
	  h_met[2][chan]->Fill(pfMet->met());

	  float diJetPt = 0.;
	  bool matchingWorked = GetDiJetPt(event, candidate_pair, diJetPt, lead_matched_jetpt_, trail_matched_jetpt_);
	  if(!matchingWorked) nCnt[46][chan]++;

	  h_dijetpt[2][chan]->Fill(diJetPt, pfJets.size());

	  diJetPt_ = diJetPt;
	  pfMET_ = pfMet->met();
	  diEMpT_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt();
	  Njets_ = pfJets.size();
	  Nbtags_ = btags.size();
	  Nelectrons_ = isoEles.size();
	  Nmuons_ = isoMuons.size();
	  invmass_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).M();
	  HT_ = HT;

	  ffTrees[chan]->Fill();
	  
	  if(chan == 1) {
	    for(int iCut = 0; iCut < nDivisions_chi2; iCut++) {
	      double cut_val = (double)iCut / (double)nDivisions_chi2;
	      bool wouldPass = false;
	      
	      for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
		if(pfJets[iJet]->bTagDiscriminators[susy::kCSV] > cut_val) {
		  wouldPass = true;
		  break;
		}
	      }
	      if(wouldPass) h_met_varyCSVcut_ff_j[iCut]->Fill(pfMet->met());

	    }
	  }

	  if(chan == 4) {
	    for(int iCut = 0; iCut < nDivisions_chi2; iCut++) {
	      double cut_val = (double)iCut / (double)nDivisions_chi2;
	      int nPassing = 0;
    
	      for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
		if(pfJets[iJet]->bTagDiscriminators[susy::kCSV] > cut_val) nPassing++;
	      }
	      if(nPassing >= 2) h_met_squareCSVcut_ff_jj[iCut]->Fill(pfMet->met());

	    }
	  }

	  nCnt[5][chan]++;

	} // if chan jet event
	
	else if(event_type == 5 || event_type == -5) {
	  
	  bool passHLT = useTrigger ? PassTriggers(4) : true;
	  if(!passHLT) {/*nCnt[35][chan]++;*/continue;}
	  
	  h_diempt[7][chan]->Fill((candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt(), pfJets.size());
	  h_met[7][chan]->Fill(pfMet->met());

	  float diJetPt = 0.;
	  bool matchingWorked = GetDiJetPt(event, candidate_pair, diJetPt, lead_matched_jetpt_, trail_matched_jetpt_);
	  if(!matchingWorked) nCnt[46][chan]++;

	  h_dijetpt[7][chan]->Fill(diJetPt, pfJets.size());

	  diJetPt_ = diJetPt;
	  pfMET_ = pfMet->met();
	  diEMpT_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt();
	  Njets_ = pfJets.size();
	  Nbtags_ = btags.size();
	  Nelectrons_ = isoEles.size();
	  Nmuons_ = isoMuons.size();
	  invmass_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).M();
	  HT_ = HT;

	  gfTrees[chan]->Fill();
	  
	  nCnt[7][chan]++;
	  
	} // if chan jet event

      } // loop over jet/btag req channels

    } // if good event without dPhi cut

    ///////////////////////////////////
    
    if(quitAfterProcessing) break;
  } // for entries
  
  cout << "-------------------Job Summary-----------------" << endl;
  cout << "Total_events         : " << nCnt[0][0] << endl;
  cout << "in_JSON              : " << nCnt[1][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  for(int i = 0; i < nChannels; i++) {
    cout << "----------------" << channels[i] << " Requirement-------------" << endl;
    cout << "gg+" << channels[i] << " events              : " << nCnt[2][i] << endl;
    cout << "eg+" << channels[i] << " events              : " << nCnt[3][i] << endl;
    cout << "ee+" << channels[i] << " events              : " << nCnt[4][i] << endl;
    cout << "ee+" << channels[i] << " events (onMass)     : " << nCnt[14][i] << endl;
    cout << "ff+" << channels[i] << " events              : " << nCnt[5][i] << endl;
    cout << "gf+" << channels[i] << " events              : " << nCnt[7][i] << endl;
  }
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  cout << "----------------Continues, info----------------" << endl;
  cout << "duplicate events          : " << nCnt[20][0] << endl;
  cout << "fail MET filters          : " << nCnt[21][0] << endl;
  cout << "no good PV                : " << nCnt[22][0] << endl;
  cout << "found a good photon       : " << nCnt[23][0] << endl;
  cout << "pfMet not available       : " << nCnt[24][0] << endl;
  cout << ">=2 good photons          : " << nCnt[25][0] << endl;
  cout << "lead et is good (no dphi) : " << nCnt[26][0] << endl;
  cout << "lead et is good (dphi)    : " << nCnt[27][0] << endl;
  cout << "no passing candidates     : " << nCnt[28][0] << endl;
  cout << "JEC not available         : " << nCnt[29][0] << endl;
  cout << "bad jet                   : " << nCnt[30][0] << endl;
  /*
  for(int i = 0; i < nChannels; i++) {
    cout << "----------------" << channels[i] << " b-tagged matched jets-------" << endl;
    cout << "gg+" << channels[i] << " b-tagged matches       : " << nCnt[36][i] << endl;
    cout << "ee+" << channels[i] << " b-tagged matches       : " << nCnt[37][i] << endl;
    cout << "ff+" << channels[i] << " b-tagged matches       : " << nCnt[38][i] << endl;
  }
  */
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  cout << "events with no dijetpt    : " << nCnt[46][0] << endl;
  cout << "ee+j events with different e vertices : " << nCnt[47][1] << endl;

  out->cd();
  out->Write();
  out->Close();

}


void SusyEventAnalyzer::Filter() {

  int nFiltered = 0;
  TTree* filterTree = 0;
  
  TFile* filterFile = new TFile(outputName+".root", "RECREATE");
  filterTree = (TTree*)fTree->GetTree()->CloneTree(0);
  filterTree->SetAutoSave();
  
  std::map<int, std::set<int> > allEvents;
  
  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0)) ) {
      cout << int(jentry) << " events processed, " << nFiltered << " events skimmed" << endl;
    }

    //if(printLevel > 0) cout << "Check duplicated events for data only." << endl;
    //bool duplicateEvent = ! (allEvents[event.runNumber].insert(event.eventNumber)).second;
    //if(event.isRealData && duplicateEvent) continue;

    if(!IsGoodLumi(event.runNumber, event.luminosityBlockNumber)) continue;
    
    bool passHLT = useTrigger ? PassTriggers(4) : true;
    if(!passHLT) continue;

    int nVtx = GetNumberPV(event);
    if(nVtx == 0) continue;
    
    int nTrail = 0;
    int nLead = 0;
    
    map< TString, vector<susy::Photon> >::iterator phoMap = event.photons.find("photons");
    if(phoMap != event.photons.end()) {
      for(vector<susy::Photon>::iterator p_it = phoMap->second.begin(); p_it != phoMap->second.end(); p_it++) {
	if(p_it->momentum.Et() > 25.0 &&
	   fabs(p_it->caloPosition.Eta()) < 1.4442 &&
	   p_it->r9 <= 1.0 &&
	   p_it->hadTowOverEm < 0.05
	   ) 
	  {
	    nTrail++;
	    if(p_it->momentum.Et() > 40.0) nLead++;
	  }
      }
    }
    
    bool filterThis = (nTrail > 1 && nLead > 0 && nVtx > 0);
    if(filterThis) {
      nFiltered++;
      filterTree->Fill();
    }
    
  }

  cout << "All events      : " << processNEvents << endl;
  cout << "Filtered events : " << nFiltered << " (" << (float)nFiltered/processNEvents << "%)" << endl;

  filterTree->GetCurrentFile()->cd();
  filterTree->GetCurrentFile()->Write();
  filterTree->GetCurrentFile()->Close();

}

void SusyEventAnalyzer::Acceptance() {

  const int nCategories = 4;
  TString categories[nCategories] = {"gg", "eg", "ff", "gf"};
  const int nChannels = 7;

  TString channels[nChannels] = {
    "nojet",
    "j", "b",
    "bj",
    "muJets", // gg+mu+bj + X (dilep veto)
    "eleJets",
    "hadronic" // gg+5j1b + X (lep veto)
  };

  unsigned int nJetReq[nChannels] = {
    0,
    1, 1,
    2,
    2,
    2,
    4};
  
  unsigned int nBtagReq[nChannels] = {
    0,
    0, 1,
    1,
    1,
    1,
    1};
  
  unsigned int nEleReq[nChannels] = {
    0,
    0, 0,
    0,
    0,
    1,
    0};

  unsigned int nMuonReq[nChannels] = {
    0,
    0, 0,
    0,
    1,
    0,
    0};

  bool looseLeptonVeto[nChannels] = {
    false,
    false, false,
    false,
    true,
    true,
    true};

  const int NCNT = 50;
  int nCnt[NCNT][nChannels];
  for(int i = 0; i < NCNT; i++) {
    for(int j = 0; j < nChannels; j++) {
    nCnt[i][j] = 0;
    }
  }
  
  TString output_code_t = FormatName(scan);

  // open histogram file and define histograms
  TFile * out = new TFile("signal_contamination"+output_code_t+".root", "RECREATE");
  out->cd();

  VTH1F h_met = BookTH1FVector("met"+output_code_t, "MET;#slash{E}_{T} (GeV);Events", 400, 0., 2000., nCategories, categories, nChannels, channels);
  
  TH2F * h_dR_ele_gamma = new TH2F("dR_ele_gamma", "#DeltaR between gsf electrons and lead/trailing #gamma#gamma candidates;#DeltaR_{lead #gamma, ele};#DeltaR_{trail #gamma, ele}", 50, 0, 5, 50, 0, 5);
  TH2F * h_dR_mu_gamma = new TH2F("dR_mu_gamma", "#DeltaR between muons and lead/trailing #gamma#gamma candidates;#DeltaR_{lead #gamma, #mu};#DeltaR_{trail #gamma, #mu}", 50, 0, 5, 50, 0, 5);

  TH1F * h_jet_lF = new TH1F("jet_lF", "jet loose Full", 2, 0, 2);
  TH1F * h_jet_lS = new TH1F("jet_lS", "jet loose Simple", 2, 0, 2);
  TH1F * h_jet_lC = new TH1F("jet_lC", "jet loose CutBased", 2, 0, 2);
  TH1F * h_jet_mF = new TH1F("jet_mF", "jet medium Full", 2, 0, 2);
  TH1F * h_jet_mS = new TH1F("jet_mS", "jet medium Simple", 2, 0, 2);
  TH1F * h_jet_mC = new TH1F("jet_mC", "jet medium CutBased", 2, 0, 2);
  TH1F * h_jet_tF = new TH1F("jet_tF", "jet tight Full", 2, 0, 2);
  TH1F * h_jet_tS = new TH1F("jet_tS", "jet tight Simple", 2, 0, 2);
  TH1F * h_jet_tC = new TH1F("jet_tC", "jet tight CutBased", 2, 0, 2);

  TH1F * h_btag_lF = new TH1F("btag_lF", "btag loose Full", 2, 0, 2);
  TH1F * h_btag_lS = new TH1F("btag_lS", "btag loose Simple", 2, 0, 2);
  TH1F * h_btag_lC = new TH1F("btag_lC", "btag loose CutBased", 2, 0, 2);
  TH1F * h_btag_mF = new TH1F("btag_mF", "btag medium Full", 2, 0, 2);
  TH1F * h_btag_mS = new TH1F("btag_mS", "btag medium Simple", 2, 0, 2);
  TH1F * h_btag_mC = new TH1F("btag_mC", "btag medium CutBased", 2, 0, 2);
  TH1F * h_btag_tF = new TH1F("btag_tF", "btag tight Full", 2, 0, 2);
  TH1F * h_btag_tS = new TH1F("btag_tS", "btag tight Simple", 2, 0, 2);
  TH1F * h_btag_tC = new TH1F("btag_tC", "btag tight CutBased", 2, 0, 2);

  float pfMET_ = 0.;
  float diEMpT_ = 0.;
  float diJetPt_ = 0.;
  int Njets_ = 0;
  int Nbtags_ = 0;
  int Nelectrons_ = 0;
  int Nmuons_ = 0;
  float invmass_ = 0.;
  float HT_ = 0.;
  float HT_jets_ = 0.;
  float hadronic_pt_ = 0.;
  float minDPhi_gMET_ = 0;
  float minDPhi_jMET_ = 0;
  float lead_Et_ = 0;
  float trail_Et_ = 0;
  float lead_Eta_ = 0;
  float trail_Eta_ = 0;
  float lead_Phi_ = 0;
  float trail_Phi_ = 0;

  float lead_matched_jetpt_ = 0;
  float trail_matched_jetpt_ = 0;

  float leadptOverInvmass_ = 0;
  float trailptOverInvmass_ = 0;

  float lead_chHadIso_ = 0;
  float trail_chHadIso_ = 0;
  float lead_sIetaIeta_ = 0;
  float trail_sIetaIeta_ = 0;

  int nPV_ = 0;
  float photon_dR_ = 0;
  float photon_dPhi_ = 0;

  float jet1_pt_ = -1;
  float jet2_pt_ = -1;
  float jet3_pt_ = -1;
  float jet4_pt_ = -1;

  float btag1_pt_ = -1;
  float btag2_pt_ = -1;

  float max_csv_ = -1;
  float submax_csv_ = -1;
  float min_csv_ = 10;

  float minDR_leadPhoton_jets_ = -1;
  float minDR_trailPhoton_jets_ = -1;

  float pileupWeight_ = 0.;
  float pileupWeightErr_ = 0.;
  float btagWeight_ = 0.;
  float btagWeightUp_ = 0.;
  float btagWeightDown_ = 0.;
  float btagWeightErr_ = 0.;

  Int_t metFilterBit_ = 0;

  vector<TTree*> ggTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree("gg_"+channels[i]+"_EvtTree"+output_code_t, "An event tree for final analysis");
    tree->Branch("pfMET", &pfMET_, "pfMET_/F");
    tree->Branch("diEMpT", &diEMpT_, "diEMpT_/F");
    tree->Branch("diJetPt", &diJetPt_, "diJetPt_/F");
    tree->Branch("Njets", &Njets_, "Njets_/I");
    tree->Branch("Nbtags", &Nbtags_, "Nbtags_/I");
    tree->Branch("Nelectrons", &Nelectrons_, "Nelectrons_/I");
    tree->Branch("Nmuons", &Nmuons_, "Nmuons_/I");
    tree->Branch("invmass", &invmass_, "invmass_/F");
    tree->Branch("HT", &HT_, "HT_/F");
    tree->Branch("HT_jets", &HT_jets_, "HT_jets_/F");
    tree->Branch("hadronic_pt", &hadronic_pt_, "hadronic_pt_/F");
    tree->Branch("minDPhi_gMET", &minDPhi_gMET_, "minDPhi_gMET_/F");
    tree->Branch("minDPhi_jMET", &minDPhi_jMET_, "minDPhi_jMET_/F");
    tree->Branch("minDR_leadPhoton_jets", &minDR_leadPhoton_jets_, "minDR_leadPhoton_jets_/F");
    tree->Branch("minDR_trailPhoton_jets", &minDR_trailPhoton_jets_, "minDR_trailPhoton_jets_/F");
    tree->Branch("leadPhotonEt", &lead_Et_, "lead_Et_/F");
    tree->Branch("trailPhotonEt", &trail_Et_, "trail_Et_/F");
    tree->Branch("leadMatchedJetPt", &lead_matched_jetpt_, "lead_matched_jetpt_/F");
    tree->Branch("trailMatchedJetPt", &trail_matched_jetpt_, "trail_matched_jetpt_/F");
    tree->Branch("leadPhotonEta", &lead_Eta_, "lead_Eta_/F");
    tree->Branch("trailPhotonEta", &trail_Eta_, "trail_Eta_/F");
    tree->Branch("leadPhotonPhi", &lead_Phi_, "lead_Phi_/F");
    tree->Branch("trailPhotonPhi", &trail_Phi_, "trail_Phi_/F");
    tree->Branch("leadptOverInvmass", &leadptOverInvmass_, "leadptOverInvmass_/F");
    tree->Branch("trailptOverInvmass", &trailptOverInvmass_, "trailptOverInvmass_/F");
    tree->Branch("leadChargedHadronIso", &lead_chHadIso_, "lead_chHadIso_/F");
    tree->Branch("trailChargedHadronIso", &trail_chHadIso_, "trail_chHadIso_/F");
    tree->Branch("leadSigmaIetaIeta", &lead_sIetaIeta_, "lead_sIetaIeta_/F");
    tree->Branch("trailSigmaIetaIeta", &trail_sIetaIeta_, "trail_sIetaIeta_/F");
    tree->Branch("jet1_pt", &jet1_pt_, "jet1_pt_/F");
    tree->Branch("jet2_pt", &jet2_pt_, "jet2_pt_/F");
    tree->Branch("jet3_pt", &jet3_pt_, "jet3_pt_/F");
    tree->Branch("jet4_pt", &jet4_pt_, "jet4_pt_/F");
    tree->Branch("btag1_pt", &btag1_pt_, "btag1_pt_/F");
    tree->Branch("btag2_pt", &btag2_pt_, "btag2_pt_/F");
    tree->Branch("max_csv", &max_csv_, "max_csv_/F");
    tree->Branch("submax_csv", &submax_csv_, "submax_csv_/F");
    tree->Branch("min_csv", &min_csv_, "min_csv_/F");
    tree->Branch("nPV", &nPV_, "nPV_/I");
    tree->Branch("photon_dR", &photon_dR_, "photon_dR_/F");
    tree->Branch("photon_dPhi", &photon_dPhi_, "photon_dPhi_/F");
    tree->Branch("pileupWeight", &pileupWeight_, "pileupWeight_/F");
    tree->Branch("pileupWeightErr", &pileupWeightErr_, "pileupWeightErr_/F");
    tree->Branch("btagWeight", &btagWeight_, "btagWeight_/F");
    tree->Branch("btagWeightUp", &btagWeightUp_, "btagWeightUp_/F");
    tree->Branch("btagWeightDown", &btagWeightDown_, "btagWeightDown_/F");
    tree->Branch("btagWeightErr", &btagWeightErr_, "btagWeightErr_/F");
    tree->Branch("metFilterBit", &metFilterBit_, "metFilterBit_/I");

    ggTrees.push_back(tree);
  }

  vector<TTree*> ffTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree("ff_"+channels[i]+"_EvtTree"+output_code_t, "An event tree for final analysis");
    tree->Branch("pfMET", &pfMET_, "pfMET_/F");
    tree->Branch("diEMpT", &diEMpT_, "diEMpT_/F");
    tree->Branch("diJetPt", &diJetPt_, "diJetPt_/F");
    tree->Branch("Njets", &Njets_, "Njets_/I");
    tree->Branch("Nbtags", &Nbtags_, "Nbtags_/I");
    tree->Branch("Nelectrons", &Nelectrons_, "Nelectrons_/I");
    tree->Branch("Nmuons", &Nmuons_, "Nmuons_/I");
    tree->Branch("invmass", &invmass_, "invmass_/F");
    tree->Branch("HT", &HT_, "HT_/F");
    tree->Branch("HT_jets", &HT_jets_, "HT_jets_/F");
    tree->Branch("hadronic_pt", &hadronic_pt_, "hadronic_pt_/F");
    tree->Branch("minDPhi_gMET", &minDPhi_gMET_, "minDPhi_gMET_/F");
    tree->Branch("minDPhi_jMET", &minDPhi_jMET_, "minDPhi_jMET_/F");
    tree->Branch("minDR_leadPhoton_jets", &minDR_leadPhoton_jets_, "minDR_leadPhoton_jets_/F");
    tree->Branch("minDR_trailPhoton_jets", &minDR_trailPhoton_jets_, "minDR_trailPhoton_jets_/F");
    tree->Branch("leadPhotonEt", &lead_Et_, "lead_Et_/F");
    tree->Branch("trailPhotonEt", &trail_Et_, "trail_Et_/F");
    tree->Branch("leadMatchedJetPt", &lead_matched_jetpt_, "lead_matched_jetpt_/F");
    tree->Branch("trailMatchedJetPt", &trail_matched_jetpt_, "trail_matched_jetpt_/F");
    tree->Branch("leadPhotonEta", &lead_Eta_, "lead_Eta_/F");
    tree->Branch("trailPhotonEta", &trail_Eta_, "trail_Eta_/F");
    tree->Branch("leadPhotonPhi", &lead_Phi_, "lead_Phi_/F");
    tree->Branch("trailPhotonPhi", &trail_Phi_, "trail_Phi_/F");
    tree->Branch("leadptOverInvmass", &leadptOverInvmass_, "leadptOverInvmass_/F");
    tree->Branch("trailptOverInvmass", &trailptOverInvmass_, "trailptOverInvmass_/F");
    tree->Branch("leadChargedHadronIso", &lead_chHadIso_, "lead_chHadIso_/F");
    tree->Branch("trailChargedHadronIso", &trail_chHadIso_, "trail_chHadIso_/F");
    tree->Branch("leadSigmaIetaIeta", &lead_sIetaIeta_, "lead_sIetaIeta_/F");
    tree->Branch("trailSigmaIetaIeta", &trail_sIetaIeta_, "trail_sIetaIeta_/F");
    tree->Branch("jet1_pt", &jet1_pt_, "jet1_pt_/F");
    tree->Branch("jet2_pt", &jet2_pt_, "jet2_pt_/F");
    tree->Branch("jet3_pt", &jet3_pt_, "jet3_pt_/F");
    tree->Branch("jet4_pt", &jet4_pt_, "jet4_pt_/F");
    tree->Branch("btag1_pt", &btag1_pt_, "btag1_pt_/F");
    tree->Branch("btag2_pt", &btag2_pt_, "btag2_pt_/F");
    tree->Branch("max_csv", &max_csv_, "max_csv_/F");
    tree->Branch("submax_csv", &submax_csv_, "submax_csv_/F");
    tree->Branch("min_csv", &min_csv_, "min_csv_/F");
    tree->Branch("nPV", &nPV_, "nPV_/I");
    tree->Branch("photon_dR", &photon_dR_, "photon_dR_/F");
    tree->Branch("photon_dPhi", &photon_dPhi_, "photon_dPhi_/F");
    tree->Branch("pileupWeight", &pileupWeight_, "pileupWeight_/F");
    tree->Branch("pileupWeightErr", &pileupWeightErr_, "pileupWeightErr_/F");
    tree->Branch("btagWeight", &btagWeight_, "btagWeight_/F");
    tree->Branch("btagWeightUp", &btagWeightUp_, "btagWeightUp_/F");
    tree->Branch("btagWeightDown", &btagWeightDown_, "btagWeightDown_/F");
    tree->Branch("btagWeightErr", &btagWeightErr_, "btagWeightErr_/F");
    tree->Branch("metFilterBit", &metFilterBit_, "metFilterBit_/I");

    ffTrees.push_back(tree);
  }

  vector<TTree*> gfTrees;
  for(int i = 0; i < nChannels; i++) {
    TTree * tree = new TTree("gf_"+channels[i]+"_EvtTree"+output_code_t, "An event tree for final analysis");
    tree->Branch("pfMET", &pfMET_, "pfMET_/F");
    tree->Branch("diEMpT", &diEMpT_, "diEMpT_/F");
    tree->Branch("diJetPt", &diJetPt_, "diJetPt_/F");
    tree->Branch("Njets", &Njets_, "Njets_/I");
    tree->Branch("Nbtags", &Nbtags_, "Nbtags_/I");
    tree->Branch("Nelectrons", &Nelectrons_, "Nelectrons_/I");
    tree->Branch("Nmuons", &Nmuons_, "Nmuons_/I");
    tree->Branch("invmass", &invmass_, "invmass_/F");
    tree->Branch("HT", &HT_, "HT_/F");
    tree->Branch("HT_jets", &HT_jets_, "HT_jets_/F");
    tree->Branch("hadronic_pt", &hadronic_pt_, "hadronic_pt_/F");
    tree->Branch("minDPhi_gMET", &minDPhi_gMET_, "minDPhi_gMET_/F");
    tree->Branch("minDPhi_jMET", &minDPhi_jMET_, "minDPhi_jMET_/F");
    tree->Branch("minDR_leadPhoton_jets", &minDR_leadPhoton_jets_, "minDR_leadPhoton_jets_/F");
    tree->Branch("minDR_trailPhoton_jets", &minDR_trailPhoton_jets_, "minDR_trailPhoton_jets_/F");
    tree->Branch("leadPhotonEt", &lead_Et_, "lead_Et_/F");
    tree->Branch("trailPhotonEt", &trail_Et_, "trail_Et_/F");
    tree->Branch("leadMatchedJetPt", &lead_matched_jetpt_, "lead_matched_jetpt_/F");
    tree->Branch("trailMatchedJetPt", &trail_matched_jetpt_, "trail_matched_jetpt_/F");
    tree->Branch("leadPhotonEta", &lead_Eta_, "lead_Eta_/F");
    tree->Branch("trailPhotonEta", &trail_Eta_, "trail_Eta_/F");
    tree->Branch("leadPhotonPhi", &lead_Phi_, "lead_Phi_/F");
    tree->Branch("trailPhotonPhi", &trail_Phi_, "trail_Phi_/F");
    tree->Branch("leadptOverInvmass", &leadptOverInvmass_, "leadptOverInvmass_/F");
    tree->Branch("trailptOverInvmass", &trailptOverInvmass_, "trailptOverInvmass_/F");
    tree->Branch("leadChargedHadronIso", &lead_chHadIso_, "lead_chHadIso_/F");
    tree->Branch("trailChargedHadronIso", &trail_chHadIso_, "trail_chHadIso_/F");
    tree->Branch("leadSigmaIetaIeta", &lead_sIetaIeta_, "lead_sIetaIeta_/F");
    tree->Branch("trailSigmaIetaIeta", &trail_sIetaIeta_, "trail_sIetaIeta_/F");
    tree->Branch("jet1_pt", &jet1_pt_, "jet1_pt_/F");
    tree->Branch("jet2_pt", &jet2_pt_, "jet2_pt_/F");
    tree->Branch("jet3_pt", &jet3_pt_, "jet3_pt_/F");
    tree->Branch("jet4_pt", &jet4_pt_, "jet4_pt_/F");
    tree->Branch("btag1_pt", &btag1_pt_, "btag1_pt_/F");
    tree->Branch("btag2_pt", &btag2_pt_, "btag2_pt_/F");
    tree->Branch("max_csv", &max_csv_, "max_csv_/F");
    tree->Branch("submax_csv", &submax_csv_, "submax_csv_/F");
    tree->Branch("min_csv", &min_csv_, "min_csv_/F");
    tree->Branch("nPV", &nPV_, "nPV_/I");
    tree->Branch("photon_dR", &photon_dR_, "photon_dR_/F");
    tree->Branch("photon_dPhi", &photon_dPhi_, "photon_dPhi_/F");
    tree->Branch("pileupWeight", &pileupWeight_, "pileupWeight_/F");
    tree->Branch("pileupWeightErr", &pileupWeightErr_, "pileupWeightErr_/F");
    tree->Branch("btagWeight", &btagWeight_, "btagWeight_/F");
    tree->Branch("btagWeightUp", &btagWeightUp_, "btagWeightUp_/F");
    tree->Branch("btagWeightDown", &btagWeightDown_, "btagWeightDown_/F");
    tree->Branch("btagWeightErr", &btagWeightErr_, "btagWeightErr_/F");
    tree->Branch("metFilterBit", &metFilterBit_, "metFilterBit_/I");

    gfTrees.push_back(tree);
  }

  ScaleFactorInfo sf(btagger);
  TFile * btagEfficiency = new TFile("btagEfficiency"+output_code_t+".root", "READ");
  sf.SetTaggingEfficiencies((TH1F*)btagEfficiency->Get("lEff"+output_code_t), (TH1F*)btagEfficiency->Get("cEff"+output_code_t), (TH1F*)btagEfficiency->Get("bEff"+output_code_t));

  // get pileup weights
  TFile * puFile = new TFile("pileupReweighting"+output_code_t+".root", "READ");
  TH1F * puWeights = (TH1F*)puFile->Get("puWeights"+output_code_t);

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = " << event.runNumber << ", event = " << event.eventNumber << endl;
    }

    nCnt[0][0]++; // events

    float numTrueInt = -1.;
    susy::PUSummaryInfoCollection::const_iterator iBX = event.pu.begin();
    bool foundInTimeBX = false;
    while((iBX != event.pu.end()) && !foundInTimeBX) {
      if(iBX->BX == 0) {
	numTrueInt = iBX->trueNumInteractions;
	foundInTimeBX = true;
      }
      ++iBX;
    }

    float eventWeight = 0.;
    float eventWeightErr = 0.;
    if(numTrueInt >= 0.) {
      int binNum = puWeights->GetXaxis()->FindBin(numTrueInt);
      eventWeight = puWeights->GetBinContent(binNum);
      eventWeightErr = puWeights->GetBinError(binNum);
    }

    if(!doPileupReweighting) {
      eventWeight = 1.;
      eventWeightErr = 0.;
    }

    vector<susy::Photon*> candidate_pair;
    vector<susy::PFJet*> pfJets, btags;
    vector<TLorentzVector> pfJets_corrP4, btags_corrP4;
    vector<float> csvValues;
    vector<susy::Muon*> isoMuons, looseMuons;
    vector<susy::Electron*> isoEles, looseEles;
    vector<BtagInfo> tagInfos;

    int event_type = 0;

    int nPVertex = GetNumberPV(event);
    if(nPVertex == 0) continue;

    map<TString, susy::MET>::iterator met_it = event.metMap.find("pfMet");
    susy::MET* pfMet = &(met_it->second);

    findPhotons_prioritizeCount(event, candidate_pair, event_type);
    //findPhotons_prioritizeEt(event, candidate_pair, event_type);

    if(event_type == 0) {
      nCnt[28][0]++;
      continue;
    }

    lead_Et_ = candidate_pair[0]->momentum.Et();
    lead_Eta_ = candidate_pair[0]->caloPosition.Eta();
    lead_Phi_ = candidate_pair[0]->caloPosition.Phi();
    trail_Et_ = candidate_pair[1]->momentum.Et();
    trail_Eta_ = candidate_pair[1]->caloPosition.Eta();
    trail_Phi_ = candidate_pair[1]->caloPosition.Phi();
    
    nPV_ = nPVertex;
    float dEta_ = candidate_pair[0]->caloPosition.Eta() - candidate_pair[1]->caloPosition.Eta();
    photon_dPhi_ = TVector2::Phi_mpi_pi(candidate_pair[0]->caloPosition.Phi() - candidate_pair[1]->caloPosition.Phi());
    photon_dR_ = sqrt(dEta_*dEta_ + photon_dPhi_*photon_dPhi_);

    lead_chHadIso_ = chargedHadronIso_corrected(*candidate_pair[0], event.rho25);
    trail_chHadIso_ = chargedHadronIso_corrected(*candidate_pair[1], event.rho25);
    lead_sIetaIeta_ = candidate_pair[0]->sigmaIetaIeta;
    trail_sIetaIeta_ = candidate_pair[1]->sigmaIetaIeta;

    leadptOverInvmass_ = candidate_pair[0]->momentum.Pt() / ((candidate_pair[0]->momentum + candidate_pair[1]->momentum).M());
    trailptOverInvmass_ = candidate_pair[1]->momentum.Pt() / ((candidate_pair[0]->momentum + candidate_pair[1]->momentum).M());

    float HT = 0.;
    TLorentzVector hadronicSystem(0., 0., 0., 0.);

    findMuons(event, candidate_pair, isoMuons, looseMuons, HT);
    findElectrons(event, candidate_pair, isoEles, looseEles, HT);

    for(unsigned int i = 0; i < isoMuons.size(); i++) h_dR_mu_gamma->Fill(deltaR(isoMuons[i]->momentum, candidate_pair[0]->caloPosition), deltaR(isoMuons[i]->momentum, candidate_pair[1]->caloPosition));
    for(unsigned int i = 0; i < isoEles.size(); i++) h_dR_ele_gamma->Fill(deltaR(isoEles[i]->momentum, candidate_pair[0]->caloPosition), deltaR(isoEles[i]->momentum, candidate_pair[1]->caloPosition));

    findJets(event, candidate_pair, 
	     isoMuons, looseMuons,
	     isoEles, looseEles,
	     pfJets, btags,
	     sf,
	     tagInfos, csvValues, 
	     pfJets_corrP4, btags_corrP4, 
	     HT, hadronicSystem);

    HT_jets_ = HT;
    hadronic_pt_ = hadronicSystem.Pt();

    max_csv_ = (csvValues.size() >= 1) ? csvValues[0] : -1.;
    submax_csv_ = (csvValues.size() >= 2) ? csvValues[1] : -1.;
    min_csv_ = (csvValues.size() >= 1) ? csvValues.back() : -1.;


    HT += candidate_pair[0]->momentum.Pt();
    HT += candidate_pair[1]->momentum.Pt();

    jet1_pt_ = (pfJets_corrP4.size() >= 1) ? pfJets_corrP4[0].Pt() : -1.;
    jet2_pt_ = (pfJets_corrP4.size() >= 2) ? pfJets_corrP4[1].Pt() : -1.;
    jet3_pt_ = (pfJets_corrP4.size() >= 3) ? pfJets_corrP4[2].Pt() : -1.;
    jet4_pt_ = (pfJets_corrP4.size() >= 4) ? pfJets_corrP4[3].Pt() : -1.;

    btag1_pt_ = (btags_corrP4.size() >= 1) ? btags_corrP4[0].Pt() : -1.;
    btag2_pt_ = (btags_corrP4.size() >= 2) ? btags_corrP4[1].Pt() : -1.;

    for(unsigned int i = 0; i < pfJets.size(); i++) {
      h_jet_lF->Fill(pfJets[i]->passPuJetIdLoose(susy::kPUJetIdFull));
      h_jet_lS->Fill(pfJets[i]->passPuJetIdLoose(susy::kPUJetIdSimple));
      h_jet_lC->Fill(pfJets[i]->passPuJetIdLoose(susy::kPUJetIdCutBased));
      h_jet_mF->Fill(pfJets[i]->passPuJetIdMedium(susy::kPUJetIdFull));
      h_jet_mS->Fill(pfJets[i]->passPuJetIdMedium(susy::kPUJetIdSimple));
      h_jet_mC->Fill(pfJets[i]->passPuJetIdMedium(susy::kPUJetIdCutBased));
      h_jet_tF->Fill(pfJets[i]->passPuJetIdTight(susy::kPUJetIdFull));
      h_jet_tS->Fill(pfJets[i]->passPuJetIdTight(susy::kPUJetIdSimple));
      h_jet_tC->Fill(pfJets[i]->passPuJetIdTight(susy::kPUJetIdCutBased));
    }

    for(unsigned int i = 0; i < btags.size(); i++) {
      h_btag_lF->Fill(btags[i]->passPuJetIdLoose(susy::kPUJetIdFull));
      h_btag_lS->Fill(btags[i]->passPuJetIdLoose(susy::kPUJetIdSimple));
      h_btag_lC->Fill(btags[i]->passPuJetIdLoose(susy::kPUJetIdCutBased));
      h_btag_mF->Fill(btags[i]->passPuJetIdMedium(susy::kPUJetIdFull));
      h_btag_mS->Fill(btags[i]->passPuJetIdMedium(susy::kPUJetIdSimple));
      h_btag_mC->Fill(btags[i]->passPuJetIdMedium(susy::kPUJetIdCutBased));
      h_btag_tF->Fill(btags[i]->passPuJetIdTight(susy::kPUJetIdFull));
      h_btag_tS->Fill(btags[i]->passPuJetIdTight(susy::kPUJetIdSimple));
      h_btag_tC->Fill(btags[i]->passPuJetIdTight(susy::kPUJetIdCutBased));
    }

    // Calculate dPhi_min(g, MET)
    float dPhi_gMET_lead = TVector2::Phi_mpi_pi(candidate_pair[0]->caloPosition.Phi() - pfMet->mEt.Phi());
    float dPhi_gMET_trail = TVector2::Phi_mpi_pi(candidate_pair[1]->caloPosition.Phi() - pfMet->mEt.Phi());
    minDPhi_gMET_ = min(fabs(dPhi_gMET_lead), fabs(dPhi_gMET_trail));

    float min_deltaPhi_jetsMET = 999.;
    for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
      float dPhi_jetMET = fabs(TVector2::Phi_mpi_pi(pfJets[iJet]->momentum.Phi() - pfMet->mEt.Phi()));
      if(dPhi_jetMET < min_deltaPhi_jetsMET) min_deltaPhi_jetsMET = dPhi_jetMET;
    }
    minDPhi_jMET_ = (pfJets.size() > 0) ? min_deltaPhi_jetsMET : -1.;

    float min_deltaR_lead_jet = 999.;
    for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
      float dEta_x = candidate_pair[0]->caloPosition.Eta() - pfJets[iJet]->momentum.Eta();
      float dPhi_x = TVector2::Phi_mpi_pi(candidate_pair[0]->caloPosition.Phi() - pfJets[iJet]->momentum.Phi());
      float dR_x = sqrt(dEta_x*dEta_x + dPhi_x*dPhi_x);

      if(dR_x < min_deltaR_lead_jet) min_deltaR_lead_jet = dR_x;
    }
    minDR_leadPhoton_jets_ = (pfJets.size() > 0) ? min_deltaR_lead_jet : -1.;

    float min_deltaR_trail_jet = 999.;
    for(unsigned int iJet = 0; iJet < pfJets.size(); iJet++) {
      float dEta_x = candidate_pair[1]->caloPosition.Eta() - pfJets[iJet]->momentum.Eta();
      float dPhi_x = TVector2::Phi_mpi_pi(candidate_pair[1]->caloPosition.Phi() - pfJets[iJet]->momentum.Phi());
      float dR_x = sqrt(dEta_x*dEta_x + dPhi_x*dPhi_x);

      if(dR_x < min_deltaR_trail_jet) min_deltaR_trail_jet = dR_x;
    }
    minDR_trailPhoton_jets_ = (pfJets.size() > 0) ? min_deltaR_trail_jet : -1.;
      
    float btagWeight[nChannels];
    float btagWeightUp[nChannels];
    float btagWeightDown[nChannels];
    float btagWeightError[nChannels];
    for(int chan = 0; chan < nChannels; chan++) {
      BtagWeight * tagWeight = new BtagWeight(nBtagReq[chan]);
      pair<float, float> weightResult = tagWeight->weight(tagInfos, btags.size(), 0., false);
      btagWeight[chan] = weightResult.first;
      btagWeightError[chan] = weightResult.second;

      btagWeightUp[chan] = (tagWeight->weight(tagInfos, btags.size(), 1., true)).first;
      btagWeightDown[chan] = (tagWeight->weight(tagInfos, btags.size(), -1., true)).first;

      delete tagWeight;
    }
    tagInfos.clear();

    pileupWeight_ = eventWeight;
    pileupWeightErr_ = eventWeightErr;
        
    pfMET_ = pfMet->met();
    diEMpT_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).Pt();

    float diJetPt = 0.;
    bool matchingWorked = GetDiJetPt(event, candidate_pair, diJetPt, lead_matched_jetpt_, trail_matched_jetpt_);
    if(!matchingWorked) nCnt[31][0]++;
    diJetPt_ = diJetPt;

    Njets_ = pfJets.size();
    Nbtags_ = btags.size();
    Nelectrons_ = isoEles.size();
    Nmuons_ = isoMuons.size();
    invmass_ = (candidate_pair[0]->momentum + candidate_pair[1]->momentum).M();
    HT_ = HT;

    ////////////////////

    if(event_type != 0) {

      for(int chan = 0; chan < nChannels; chan++) {

	if(pfJets.size() < nJetReq[chan]) continue;
	if(btags.size() < nBtagReq[chan]) continue;
	if(looseLeptonVeto[chan]) {
	  if(isoEles.size() != nEleReq[chan]) continue;
	  if(looseEles.size() != 0) continue;
	  if(isoMuons.size() != nMuonReq[chan]) continue;
	  if(looseMuons.size() != 0) continue;
	}

	btagWeight_ = btagWeight[chan];
	btagWeightErr_ = btagWeightError[chan];
	btagWeightUp_ = btagWeightUp[chan];
	btagWeightDown_ = btagWeightDown[chan];

	if(event_type == 1) {
	
	  bool passHLT = useTrigger ? PassTriggers(1) : true;
	  if(!passHLT) {nCnt[36][0]++;continue;}
	
	  h_met[0][chan]->Fill(pfMet->met(), eventWeight * btagWeight[chan]);
	  nCnt[2][chan]++;

	  //diJetPt_ = diJetPt;

	  ggTrees[chan]->Fill();

	} // if gg jet event

	else if(event_type == 2 || event_type == -2) {
	  
	  bool passHLT = useTrigger ? PassTriggers(1) : true;
	  if(!passHLT) continue;
	  
	  h_met[1][chan]->Fill(pfMet->met(), eventWeight * btagWeight[chan]);
	  nCnt[4][chan]++;
	}

	else if(event_type == 4) {
	  
	  bool passHLT = useTrigger ? PassTriggers(4) : true;
	  if(!passHLT) {nCnt[40][0]++;continue;}
	
	  h_met[2][chan]->Fill(pfMet->met(), eventWeight * btagWeight[chan]);
	  nCnt[3][chan]++;

	  ffTrees[chan]->Fill();

	} // if ff jet event

	else if(event_type == 5 || event_type == -5) {
	  
	  bool passHLT = useTrigger ? PassTriggers(4) : true;
	  if(!passHLT) {/*nCnt[40][0]++;*/continue;}
	
	  h_met[3][chan]->Fill(pfMet->met(), eventWeight * btagWeight[chan]);
	  nCnt[5][chan]++;

	  gfTrees[chan]->Fill();

	} // if ff jet event

      } // for channels

    } // if chan event

  } // for entries

  cout << "-------------------Job Summary-----------------" << endl;
  cout << "Total_events         : " << nCnt[0][0] << endl;
  cout << "in_JSON              : " << nCnt[1][0] << endl;
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  for(int i = 0; i < nChannels; i++) {
    cout << "----------------" << channels[i] << " Requirement-------------" << endl;
    cout << "gg+" << channels[i] << " events              : " << nCnt[2][i] << endl;
    cout << "ff+" << channels[i] << " events              : " << nCnt[3][i] << endl;
    cout << "eg+" << channels[i] << " events              : " << nCnt[4][i] << endl;
    cout << "gf+" << channels[i] << " events              : " << nCnt[5][i] << endl;
  }
  cout << "-----------------------------------------------" << endl;
  cout << endl;
  cout << "----------------Continues, info----------------" << endl;
  cout << "fail MET filters          : " << nCnt[21][0] << endl;
  cout << "no good PV                : " << nCnt[22][0] << endl;
  cout << "found a good photon       : " << nCnt[23][0] << endl;
  cout << "pfMet not available       : " << nCnt[24][0] << endl;
  cout << ">=2 good photons          : " << nCnt[25][0] << endl;
  cout << "lead et is good (no dphi) : " << nCnt[26][0] << endl;
  cout << "lead et is good (dphi)    : " << nCnt[27][0] << endl;
  cout << "no passing candidates     : " << nCnt[28][0] << endl;
  cout << "JEC not available         : " << nCnt[29][0] << endl;
  cout << "bad jet                   : " << nCnt[30][0] << endl;
  cout << "diJet matching failed     : " << nCnt[31][0] << endl;

  puFile->Close();
  btagEfficiency->Close();

  out->Write();
  out->Close();

}

void SusyEventAnalyzer::SignalContent_gg() {

  char * tmp = getenv("CONDOR_SECTION");
  cout << "tmp = " << tmp << endl;
  int index = atoi(tmp);
  cout << "index = " << index << endl;
  
  // stop-bino scan
  int index1 = mst[int(index)/31];
  int index2 = mBino[int(index)%31];
  
  // bino/wino
  //int index1 = int(index)/17*100 + 400;
  //int index2 = int(index)%17*100 + 420;
  
  cout << "index1 = " << index1 << endl;
  cout << "index2 = " << index2 << endl;
  
  char output_file_name[100];
  
  ScaleFactorInfo sf(btagger);

  sprintf(output_file_name, "signalContent_mst_%d_m1_%d.root", index1, index2);

  TFile * out = new TFile(output_file_name, "RECREATE");
  out->cd();

  sprintf(output_file_name, "_mst_%d_m1_%d", index1, index2);
  TString code = output_file_name;

  TH2F * elePt_leadPt = new TH2F("elePt_leadPt"+code, "Electron pt vs leading photon pt", 400, 0, 2000, 400, 0, 2000);
  TH2F * elePt_subPt = new TH2F("elePt_subPt"+code, "Electron pt vs sub-leading photon pt", 400, 0, 2000, 400, 0, 2000);

  TH1F * h_w_mass = new TH1F("w_mass"+code, "jj mass", 1000, 0, 2000);
  TH1F * h_top_mass = new TH1F("top_mass"+code, "bjj mass", 1000, 0, 2000);

  // branching counts
  int nStops, tFromStop, bFromStop, wFromStop, nBinos;
  int nLeptonicW, nHadronicW, nLeptonicZ, nHadronicZ, nInvisibleZ;
  int gFromBino, zFromBino, wFromBino;

  // object counts
  int nPhotons_80, nPhotons_40, nPhotons_25, nPhotons;
  int nMuons, nElectrons;
  int nJets, nBjets;
  
  TTree * genTree = new TTree("genInfo"+code, "gen-level information tree");
  genTree->Branch("nStops", &nStops, "nStops/I");
  genTree->Branch("tFromStop", &tFromStop, "tFromStop/I");
  genTree->Branch("bFromStop", &bFromStop, "bFromStop/I");
  genTree->Branch("wFromStop", &wFromStop, "wFromStop/I");
  genTree->Branch("nBinos", &nBinos, "nBinos/I");
  genTree->Branch("nLeptonicW", &nLeptonicW, "nLeptonicW/I");
  genTree->Branch("nHadronicW", &nHadronicW, "nHadronicW/I");
  genTree->Branch("nLeptonicZ", &nLeptonicZ, "nLeptonicZ/I");
  genTree->Branch("nHadronicZ", &nHadronicZ, "nHadronicZ/I");
  genTree->Branch("nInvisibleZ", &nInvisibleZ, "nInvisibleZ/I");
  genTree->Branch("gFromBino", &gFromBino, "gFromBino/I");
  genTree->Branch("zFromBino", &zFromBino, "zFromBino/I");
  genTree->Branch("wFromBino", &wFromBino, "wFromBino/I");
  genTree->Branch("nPhotons_80", &nPhotons_80, "nPhotons_80/I");
  genTree->Branch("nPhotons_40", &nPhotons_40, "nPhotons_40/I");
  genTree->Branch("nPhotons_25", &nPhotons_25, "nPhotons_25/I");
  genTree->Branch("nPhotons", &nPhotons, "nPhotons/I");
  genTree->Branch("nMuons", &nMuons, "nMuons/I");
  genTree->Branch("nElectrons", &nElectrons, "nElectrons/I");
  genTree->Branch("nJets", &nJets, "nJets/I");
  genTree->Branch("nBjets", &nBjets, "nBjets/I");
  
  float ele_pt, ele_eta, ele_iso, ele_relIso;
  TTree * eleTree = new TTree("eleTree"+code, "gen-matched electron reco info");
  eleTree->Branch("pt", &ele_pt, "ele_pt/F");
  eleTree->Branch("eta", &ele_eta, "ele_eta/F");
  eleTree->Branch("iso", &ele_iso, "ele_iso/F");
  eleTree->Branch("relIso", &ele_relIso, "ele_relIso/F");

  float mu_pt, mu_eta, mu_iso, mu_relIso;
  TTree * muTree = new TTree("muTree"+code, "gen-matched muon reco info");
  muTree->Branch("pt", &mu_pt, "mu_pt/F");
  muTree->Branch("eta", &mu_eta, "mu_eta/F");
  muTree->Branch("iso", &mu_iso, "mu_iso/F");
  muTree->Branch("relIso", &mu_relIso, "mu_relIso/F");

  float g_pt, g_eta;
  TTree * photonTree = new TTree("photonTree"+code, "gen-matched photon reco info");
  photonTree->Branch("pt", &g_pt, "g_pt/F");
  photonTree->Branch("eta", &g_eta, "g_eta/F");

  float jet_pt, jet_eta, jet_csv, jet_jp, jet_tchp;
  int jet_flavor, jet_algDef, jet_phyDef, jet_mother;
  bool jet_pujid_lF, jet_pujid_lS, jet_pujid_lC, jet_pujid_mF, jet_pujid_mS, jet_pujid_mC, jet_pujid_tF, jet_pujid_tS, jet_pujid_tC;
  TTree * jetTree = new TTree("jetTree"+code, "gen-matched jet reco info");
  jetTree->Branch("pt", &jet_pt, "jet_pt/F");
  jetTree->Branch("eta", &jet_eta, "jet_eta/F");
  jetTree->Branch("csv", &jet_csv, "jet_csv/F");
  jetTree->Branch("jp", &jet_jp, "jet_jp/F");
  jetTree->Branch("tchp", &jet_tchp, "jet_tchp/F");
  jetTree->Branch("flavor", &jet_flavor, "jet_flavor/I");
  jetTree->Branch("algDef", &jet_algDef, "jet_algDef/I");
  jetTree->Branch("phyDef", &jet_phyDef, "jet_phyDef/I");
  jetTree->Branch("mother", &jet_mother, "jet_mother/I");
  jetTree->Branch("puJetId_lF", &jet_pujid_lF, "jet_pujid_lF/O");
  jetTree->Branch("puJetId_lS", &jet_pujid_lS, "jet_pujid_lS/O");
  jetTree->Branch("puJetId_lC", &jet_pujid_lC, "jet_pujid_lC/O");
  jetTree->Branch("puJetId_mF", &jet_pujid_mF, "jet_pujid_mF/O");
  jetTree->Branch("puJetId_mS", &jet_pujid_mS, "jet_pujid_mS/O");
  jetTree->Branch("puJetId_mC", &jet_pujid_mC, "jet_pujid_mC/O");
  jetTree->Branch("puJetId_tF", &jet_pujid_tF, "jet_pujid_tF/O");
  jetTree->Branch("puJetId_tS", &jet_pujid_tS, "jet_pujid_tS/O");
  jetTree->Branch("puJetId_tC", &jet_pujid_tC, "jet_pujid_tC/O");

  float btag_pt, btag_eta;
  int btag_flavor, btag_algDef, btag_phyDef, btag_mother;
  bool btag_pujid_lF, btag_pujid_lS, btag_pujid_lC, btag_pujid_mF, btag_pujid_mS, btag_pujid_mC, btag_pujid_tF, btag_pujid_tS, btag_pujid_tC;
  TTree * btagTree = new TTree("btagTree"+code, "gen-matched btag reco info");
  btagTree->Branch("pt", &btag_pt, "btag_pt/F");
  btagTree->Branch("eta", &btag_eta, "btag_eta/F");
  btagTree->Branch("flavor", &btag_flavor, "btag_flavor/I");
  btagTree->Branch("algDef", &btag_algDef, "btag_algDef/I");
  btagTree->Branch("phyDef", &btag_phyDef, "btag_phyDef/I");
  btagTree->Branch("mother", &btag_mother, "btag_mother/I");
  btagTree->Branch("puJetId_lF", &btag_pujid_lF, "btag_pujid_lF/O");
  btagTree->Branch("puJetId_lS", &btag_pujid_lS, "btag_pujid_lS/O");
  btagTree->Branch("puJetId_lC", &btag_pujid_lC, "btag_pujid_lC/O");
  btagTree->Branch("puJetId_mF", &btag_pujid_mF, "btag_pujid_mF/O");
  btagTree->Branch("puJetId_mS", &btag_pujid_mS, "btag_pujid_mS/O");
  btagTree->Branch("puJetId_mC", &btag_pujid_mC, "btag_pujid_mC/O");
  btagTree->Branch("puJetId_tF", &btag_pujid_tF, "btag_pujid_tF/O");
  btagTree->Branch("puJetId_tS", &btag_pujid_tS, "btag_pujid_tS/O");
  btagTree->Branch("puJetId_tC", &btag_pujid_tC, "btag_pujid_tC/O");
  
  float maxCSV, maxJP, maxTCHP;
  float submaxCSV, submaxJP, submaxTCHP;
  bool isDiPhoton, isSinglePhoton;
  bool isGGL, isGL;
  TTree * eventTree = new TTree("eventTree"+code, "event-wide info");
  eventTree->Branch("maxCSV", &maxCSV, "maxCSV/F");
  eventTree->Branch("submaxCSV", &submaxCSV, "submaxCSV/F");
  eventTree->Branch("maxJP", &maxJP, "maxJP/F");
  eventTree->Branch("submaxJP", &submaxJP, "submaxJP/F");
  eventTree->Branch("maxTCHP", &maxTCHP, "maxTCHP/F");
  eventTree->Branch("submaxTCHP", &submaxTCHP, "submaxTCHP/F");
  eventTree->Branch("isDiPhoton", &isDiPhoton, "isDiPhoton/O");
  eventTree->Branch("isSinglePhoton", &isSinglePhoton, "isSinglePhoton/O");
  eventTree->Branch("isGGL", &isGGL, "isGGL/O");
  eventTree->Branch("isGL", &isGL, "isGL/O");
  
  vector<TString> triggerNames_;
  vector<int> triggerCounts_;

  vector<pair<TString, int> > triggerFires;

  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;
  
  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = "
	   << event.runNumber << ", event = " << event.eventNumber << ", integ lumi = " 
	   << event.intgRecLumi << endl;
    }

    for(susy::TriggerMap::iterator it = event.hltMap.begin(); it != event.hltMap.end(); it++) {
      if((int(it->second.second)) && it->first.Contains("HLT_") && it->second.first == 1) {

	bool alreadyFound = false;
	for(unsigned int k = 0; k < triggerFires.size(); k++) {
	  if((triggerFires[k].first).Contains(it->first)) {
	    triggerFires[k].second++;
	    alreadyFound = true;
	    break;
	  }
	}

	if(!alreadyFound) {
	  triggerFires.push_back(make_pair(it->first, 1));
	}
	
	
      }
    }



    nStops = tFromStop = bFromStop = nBinos = 0;
    nLeptonicW = nHadronicW = nLeptonicZ = nHadronicZ = nInvisibleZ = 0;
    gFromBino = zFromBino = wFromBino = 0;
    nPhotons_80 = nPhotons_40 = nPhotons_25 = nPhotons = 0;
    nMuons = nElectrons = 0;
    nJets = nBjets = 0;
    ele_pt = ele_eta = ele_iso = ele_relIso = 0;
    mu_pt = mu_eta = mu_iso = mu_relIso = 0;
    g_pt = g_eta = 0;
    jet_pt = jet_eta = jet_csv = jet_jp = jet_tchp = 0;
    jet_flavor = jet_algDef = jet_phyDef = jet_mother = 0;
    jet_pujid_lF = jet_pujid_lS = jet_pujid_lC = jet_pujid_mF = jet_pujid_mS = jet_pujid_mC = jet_pujid_tF = jet_pujid_tS = jet_pujid_tC = false;
    btag_pujid_lF = btag_pujid_lS = btag_pujid_lC = btag_pujid_mF = btag_pujid_mS = btag_pujid_mC = btag_pujid_tF = btag_pujid_tS = btag_pujid_tC = false;
    btag_pt = btag_eta = 0;
    btag_flavor = btag_algDef = btag_phyDef = btag_mother = 0;
    maxCSV = maxJP = maxTCHP = 0;
    submaxCSV = submaxJP = submaxTCHP = 0;
    isDiPhoton = false;
    isSinglePhoton = false;
    isGGL = false;
    isGL = false;

    vector<double> pfJets_CSV, pfJets_JP, pfJets_TCHP;

    vector<TLorentzVector> bMomenta, jMomenta;

    map<TString, susy::PFJetCollection>::iterator pfJets_it = event.pfJets.find("ak5");
    if(pfJets_it == event.pfJets.end()) cout << "Wat @ find ak5" << endl;

    vector<float> v_photon_pt;

    for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {

      if(fabs(it->pdgId) == 1000006 && it->status == 3) nStops++;
      if(fabs(it->pdgId) == 1000022 && it->status == 3) nBinos++;
      if(fabs(it->pdgId) == 6 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000006 && it->status == 3) tFromStop++;
      if(fabs(it->pdgId) == 5 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000006 && it->status == 3) bFromStop++;
      if(fabs(it->pdgId) == 24 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000006 && it->status == 3) wFromStop++;
      if(fabs(it->pdgId) == 22 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000022 && it->status == 1) gFromBino++;
      if(fabs(it->pdgId) == 23 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000022 && it->status == 3) zFromBino++;
      if(fabs(it->pdgId) == 24 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000022 && it->status == 3) wFromBino++;

      // Gen-level photons from binos in acceptance
      if(fabs(it->pdgId) == 22 && it->status == 1 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000022 && fabs(it->momentum.Eta()) < 1.4442) {

	  if(it->momentum.Pt() > 25.0) nPhotons_25++;
	  if(it->momentum.Pt() > 40.0) nPhotons_40++;
	  if(it->momentum.Pt() > 80.0) nPhotons_80++;
	  nPhotons++;

	// find photons
	map<TString, vector<susy::Photon> >::iterator phoMap = event.photons.find("photons");
	if(phoMap != event.photons.end()) {
	  
	  for(vector<susy::Photon>::iterator pho_it = phoMap->second.begin();
	      pho_it != phoMap->second.end(); pho_it++) {
	    
	    if(deltaR(pho_it->momentum, it->momentum) >= 0.3) continue;

	    g_pt = pho_it->momentum.Pt();
	    g_eta = pho_it->caloPosition.Eta();

	    v_photon_pt.push_back(g_pt);

	    photonTree->Fill();
	    break;
	    
	  } // for photon
	} // if
    
      }

      float leading_photon_pt = -1.;
      float subleading_photon_pt = -1.;
      sort(v_photon_pt.begin(), v_photon_pt.end(), greater<float>());
      if(v_photon_pt.size() >= 1) leading_photon_pt = v_photon_pt[0];
      if(v_photon_pt.size() >= 2) subleading_photon_pt = v_photon_pt[1];

      bool isFromWfromStopOrTop = fabs(event.genParticles[it->motherIndex].pdgId) == 24 && 
	 (
	  fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 6 || 
	  fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000006
	  );

      bool isFromWorZfromBino = (fabs(event.genParticles[it->motherIndex].pdgId) == 23 || fabs(event.genParticles[it->motherIndex].pdgId) == 24) &&
	fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000022;

      if((isFromWfromStopOrTop || isFromWorZfromBino) && 
	 fabs(it->momentum.Eta()) < 2.6 &&
	 !(fabs(it->momentum.Eta()) >= 1.4442 && fabs(it->momentum.Eta()) <= 1.566)) {

	if(fabs(it->pdgId) == 11 ||
	   fabs(it->pdgId) == 13) {

	  if(fabs(event.genParticles[it->motherIndex].pdgId) == 24) nLeptonicW++;
	  if(fabs(event.genParticles[it->motherIndex].pdgId) == 23 && nLeptonicZ == 0) nLeptonicZ++;

	  if(fabs(it->pdgId) == 13) {

	    nMuons++;

	    // lepton from stop->top->W; look at kinematics of reconstructed lepton
	    map<TString, vector<susy::Muon> >::iterator muMap = event.muons.find("muons");
	    if(muMap != event.muons.end()) {
	      for(vector<susy::Muon>::iterator mu_it = muMap->second.begin(); mu_it != muMap->second.end(); mu_it++) {
	
		if(deltaR(mu_it->momentum, it->momentum) >= 0.3) continue;

		mu_pt = mu_it->momentum.Pt();
		mu_eta = mu_it->momentum.Eta();
		mu_iso = max(0., (mu_it->sumNeutralHadronEt04 + mu_it->sumPhotonEt04 - 0.5*(mu_it->sumPUPt04)));
		mu_iso += mu_it->sumChargedHadronPt04;
		mu_relIso = mu_iso / mu_pt;
		
		muTree->Fill();
		break;

	      }
	    }
	  } // if muon

	  if(fabs(it->pdgId) == 11) {

	    nElectrons++;

	    map<TString, vector<susy::Electron> >::iterator eleMap = event.electrons.find("gsfElectrons");
	    if(eleMap != event.electrons.end()) {
	      for(vector<susy::Electron>::iterator ele_it = eleMap->second.begin(); ele_it != eleMap->second.end(); ele_it++) {
		
		if(deltaR(ele_it->momentum, it->momentum) >= 0.3) continue;
		
		ele_iso = ele_it->chargedHadronIso + ele_it->photonIso + ele_it->neutralHadronIso;
		ele_pt = ele_it->momentum.Pt();
		ele_eta = ele_it->momentum.Eta();
		ele_relIso = ele_iso / ele_pt;

		elePt_leadPt->Fill(ele_pt, leading_photon_pt);
		elePt_subPt->Fill(ele_pt, subleading_photon_pt);

		eleTree->Fill();
		break;
	      }
	      
	    }
	  } // if electron
	    
	} // if lepton from stop->top->W
	else if(fabs(it->pdgId) != 12 &&
		fabs(it->pdgId) != 14 &&
		fabs(it->pdgId) != 16) {
	  if(fabs(event.genParticles[it->motherIndex].pdgId) == 24) nHadronicW++;
	  if(fabs(event.genParticles[it->motherIndex].pdgId) == 23 && nHadronicZ == 0) nHadronicZ++;
	}
	else if(fabs(it->pdgId) == 12 &&
		fabs(it->pdgId) == 14 &&
		fabs(it->pdgId) == 16 &&
		nInvisibleZ == 0) nInvisibleZ++;

      }

      bool jetFromCascade = fabs(event.genParticles[it->motherIndex].pdgId) >= 1000001 && 
	fabs(event.genParticles[it->motherIndex].pdgId) <= 1000038;

      bool jetFromWfromTop = fabs(event.genParticles[it->motherIndex].pdgId) == 24 && 
	fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 6;

      bool jetFromTop = fabs(event.genParticles[it->motherIndex].pdgId) == 6;

      bool jetFromWfromBino = fabs(event.genParticles[it->motherIndex].pdgId) == 24 && 
	fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000022;

      bool jetFromZfromBino = fabs(event.genParticles[it->motherIndex].pdgId) == 23 && 
	fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId) == 1000022;

      bool jetFromInteraction = jetFromCascade || jetFromWfromTop || jetFromTop || jetFromWfromBino || jetFromZfromBino;

      // Gen-level jets
      if((fabs(it->pdgId) == 1 ||
	  fabs(it->pdgId) == 2 ||
	  fabs(it->pdgId) == 3 ||
	  fabs(it->pdgId) == 4 ||
	  fabs(it->pdgId) == 5 ||
	  fabs(it->pdgId) == 21) &&
	 //it->status == 3 &&
	 jetFromInteraction &&
	 it->momentum.Pt() > 30.0 &&
	 fabs(it->momentum.Eta()) < 2.6) {

	nJets++;

	if(fabs(it->pdgId) == 5 && fabs(it->momentum.Eta()) < 2.4) {
	  nBjets++;
	}

	// Match this jet to a good reco pfjet
	if(pfJets_it != event.pfJets.end()) {
	  for(vector<susy::PFJet>::iterator jet = pfJets_it->second.begin(); jet != pfJets_it->second.end(); jet++) {
	    map<TString, Float_t>::iterator s_it = jet->jecScaleFactors.find("L1FastL2L3");
	    if(s_it == jet->jecScaleFactors.end()) { cout << "Wat @ jec" << endl; continue; }
	    TLorentzVector corrP4 = s_it->second * jet->momentum;
	    
	    if(deltaR(corrP4, it->momentum) >= 0.3) continue;

	    jet_pt = corrP4.Pt();
	    jet_eta = corrP4.Eta();
	    jet_csv = jet->bTagDiscriminators[susy::kCSV];
	    jet_jp = jet->bTagDiscriminators[susy::kJP];
	    jet_tchp = jet->bTagDiscriminators[susy::kTCHP];

	    jet_flavor = it->pdgId;
	    jet_algDef = jet->algDefFlavour;
	    jet_phyDef = jet->phyDefFlavour;
	    jet_mother = event.genParticles[it->motherIndex].pdgId;

	    jet_pujid_lF = jet->passPuJetIdLoose(susy::kPUJetIdFull);
	    jet_pujid_lS = jet->passPuJetIdLoose(susy::kPUJetIdSimple);
	    jet_pujid_lC = jet->passPuJetIdLoose(susy::kPUJetIdCutBased);

	    jet_pujid_mF = jet->passPuJetIdMedium(susy::kPUJetIdFull);
	    jet_pujid_mS = jet->passPuJetIdMedium(susy::kPUJetIdSimple);
	    jet_pujid_mC = jet->passPuJetIdMedium(susy::kPUJetIdCutBased);

	    jet_pujid_tF = jet->passPuJetIdTight(susy::kPUJetIdFull);
	    jet_pujid_tS = jet->passPuJetIdTight(susy::kPUJetIdSimple);
	    jet_pujid_tC = jet->passPuJetIdTight(susy::kPUJetIdCutBased);

	    jetTree->Fill();

	    bool is_btagged = false;
	    
	    if(fabs(corrP4.Eta()) < 2.4) {
	      pfJets_CSV.push_back(jet_csv);
	      pfJets_JP.push_back(jet_jp);
	      pfJets_TCHP.push_back(jet_tchp);
	   
	      if((btagger == "CSVL" && jet->bTagDiscriminators[susy::kCSV] > 0.244) ||
		 (btagger == "CSVM" && jet->bTagDiscriminators[susy::kCSV] > 0.679) ||
		 (btagger == "CSVT" && jet->bTagDiscriminators[susy::kCSV] > 0.898)) {

		btag_pt = jet_pt;
		btag_eta = jet_eta;

		btag_flavor = jet_flavor;
		btag_algDef = jet_algDef;
		btag_phyDef = jet_phyDef;
		btag_mother = jet_mother;

		btag_pujid_lF = jet->passPuJetIdLoose(susy::kPUJetIdFull);
		btag_pujid_lS = jet->passPuJetIdLoose(susy::kPUJetIdSimple);
		btag_pujid_lC = jet->passPuJetIdLoose(susy::kPUJetIdCutBased);

		btag_pujid_mF = jet->passPuJetIdMedium(susy::kPUJetIdFull);
		btag_pujid_mS = jet->passPuJetIdMedium(susy::kPUJetIdSimple);
		btag_pujid_mC = jet->passPuJetIdMedium(susy::kPUJetIdCutBased);

		btag_pujid_tF = jet->passPuJetIdTight(susy::kPUJetIdFull);
		btag_pujid_tS = jet->passPuJetIdTight(susy::kPUJetIdSimple);
		btag_pujid_tC = jet->passPuJetIdTight(susy::kPUJetIdCutBased);

		btagTree->Fill();

		bMomenta.push_back(corrP4);
		is_btagged = true;
	      }
	      
	    }

	    if(!is_btagged)  jMomenta.push_back(corrP4);

	    break;
	  }
	}

      }

    }

    sort(pfJets_CSV.begin(), pfJets_CSV.end(), greater<double>());
    if(pfJets_CSV.size() >= 1) maxCSV = pfJets_CSV[0];
    else maxCSV = -9999.;
    if(pfJets_CSV.size() >= 2) submaxCSV = pfJets_CSV[1];
    else submaxCSV = -9999.;

    sort(pfJets_JP.begin(), pfJets_JP.end(), greater<double>());
    if(pfJets_JP.size() >= 1) maxJP = pfJets_JP[0];
    else maxJP = -9999.;
    if(pfJets_JP.size() >= 2) submaxJP = pfJets_JP[1];
    else submaxJP = -9999.;

    sort(pfJets_TCHP.begin(), pfJets_TCHP.end(), greater<double>());
    if(pfJets_TCHP.size() >= 1) maxTCHP = pfJets_TCHP[0];
    else maxTCHP = -9999.;
    if(pfJets_TCHP.size() >= 2) submaxTCHP = pfJets_TCHP[1];
    else submaxTCHP = -9999.;

    isDiPhoton = nPhotons_25 > 1 && nPhotons_40 > 0;
    isSinglePhoton = nPhotons_80 > 0;

    isGL = isSinglePhoton && (nMuons > 0 || nElectrons > 0);
    isGGL = isDiPhoton && (nMuons > 0 || nElectrons > 0);

    genTree->Fill();
    eventTree->Fill();

    if(bMomenta.size() >= 1 && jMomenta.size() >= 2) {
      for(unsigned int ib = 0; ib < bMomenta.size(); ib++) {
	for(unsigned int j1 = 0; j1 < jMomenta.size() - 1; j1++) {
	  for(unsigned int j2 = j1; j2 < jMomenta.size(); j2++) {

	    double w_mass = (jMomenta[j1] + jMomenta[j2]).M();
	    h_w_mass->Fill(w_mass);

	    double top_mass = (bMomenta[ib] + jMomenta[j1] + jMomenta[j2]).M();
	    h_top_mass->Fill(top_mass);

	  }
	}
      }

    }
	 

  } // event loop

  cout << endl << endl << "Found trigger counts: " << endl << endl;

  sort(triggerFires.begin(), triggerFires.end(), sortTriggers);
  for(unsigned int i = 0; i < triggerFires.size(); i++) {
    cout << triggerFires[i].first << " -- " << triggerFires[i].second << " times" << endl;
  }

  out->cd();
  out->Write();
  out->Close();

  
}

void SusyEventAnalyzer::Step0() {

  char * tmp = getenv("CONDOR_SECTION");
  int index = atoi(tmp);
  
  bool stopScan = (scan == "stop-bino");
  bool stopScan_ext = (scan == "stop-bino extended");
  bool squarkGluinoScan = (scan == "squark-gluino");

  int index1, index2;

  if(stopScan) {
    index1 = mst[int(index)/31];
    index2 = mBino[int(index)%31];
  }
  else if(stopScan_ext) {
    index1 = mst_ext[int(index)/5];
    index2 = mBino_ext[int(index)%5];
  }
  else if(squarkGluinoScan) {
    index1 = int(index)/17*100 + 400;
    index2 = int(index)%17*100 + 420;
  }
  else {
    index1 = 0;
    index2 = 0;
  }

  const int NCNT = 50;
  int nCnt[NCNT];
  for(int i = 0; i < NCNT; i++) nCnt[i] = 0;

  TString output_code_t = FormatName(scan);

  vector<TString> hltList;

  // L1_DoubleEG_13_7 seed
  hltList.push_back("L1_DoubleEG_13_7");
  hltList.push_back("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v");
  hltList.push_back("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v");
  hltList.push_back("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v");
  hltList.push_back("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v");
  hltList.push_back("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v");
  hltList.push_back("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v");
  hltList.push_back("HLT_Photon26_Photon18_v"); // 7 (pre-scaled)
  hltList.push_back("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v"); // 8 (pre-scaled)

  // L1_SingleEG_22 seed
  hltList.push_back("L1_SingleEG_22");
  hltList.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
  hltList.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v");
  hltList.push_back("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");
  hltList.push_back("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v");
  hltList.push_back("HLT_Photon36_R9Id85_Photon22_R9Id85_v");
  hltList.push_back("HLT_Photon36_Photon22_v"); // 15 (pre-scaled)
  hltList.push_back("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_v"); // 16 (pre-scaled)

  // top pag stuff
  hltList.push_back("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFNoPUJet50_40_30_v");
  hltList.push_back("HLT_Ele25_CaloIdVT_CaloIsoVL_TrkIdVL_TrkIsoT_TriCentralPFNoPUJet50_40_30_v");
  hltList.push_back("HLT_IsoMu17_eta2p1_TriCentralPFNoPUJet50_40_30_v");
  hltList.push_back("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFNoPUJet30_BTagIPIter_v");
  hltList.push_back("HLT_IsoMu17_eta2p1_CentralPFNoPUJet30_BTagIPIter_v");
  hltList.push_back("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  hltList.push_back("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");
  hltList.push_back("HLT_Mu17_Mu8_v");
  hltList.push_back("HLT_Mu17_TkMu8_v");
  hltList.push_back("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v");

  hltList.push_back("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TriCentralPFJet30_v");
  hltList.push_back("HLT_IsoMu17_eta2p1_TriCentralPFJet30_v");
  hltList.push_back("HLT_Ele25_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralPFJet30_BTagIPIter_v");
  hltList.push_back("HLT_IsoMu20_eta2p1_CentralPFJet30_BTagIPIter_v");

  // joey's list
  hltList.push_back("HLT_IsoMu24_eta2p1_v");
  hltList.push_back("HLT_Ele27_WP80_v");

  vector<pair<TString, int> > triggerFires;
  for(unsigned int i = 0; i < hltList.size(); i++) triggerFires.push_back(make_pair(hltList[i], 0));

  // open histogram file and define histograms
  TFile * out = new TFile("step0"+output_code_t+".root", "RECREATE");
  out->cd();

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

  vector<TH2D*> h_gg, h_gg_muJets, h_gg_eleJets, h_gg_hadronic, h_gg_eleMu, h_gg_diEle, h_gg_diMu;
  for(unsigned int i = 0; i < triggerFires.size(); i++) {
    char buff[10];
    sprintf(buff, "_%d", i);
    TString buff_t = buff;
    h_gg.push_back(new TH2D("gg"+buff_t, triggerFires[i].first, 30, xbins, 32, ybins));
    h_gg_muJets.push_back(new TH2D("gg_muJets"+buff_t, triggerFires[i].first, 30, xbins, 32, ybins));
    h_gg_eleJets.push_back(new TH2D("gg_eleJets"+buff_t, triggerFires[i].first, 30, xbins, 32, ybins));
    h_gg_hadronic.push_back(new TH2D("gg_hadronic"+buff_t, triggerFires[i].first, 30, xbins, 32, ybins));
    h_gg_eleMu.push_back(new TH2D("gg_eleMu"+buff_t, triggerFires[i].first, 30, xbins, 32, ybins));
    h_gg_diEle.push_back(new TH2D("gg_diEle"+buff_t, triggerFires[i].first, 30, xbins, 32, ybins));
    h_gg_diMu.push_back(new TH2D("gg_diMu"+buff_t, triggerFires[i].first, 30, xbins, 32, ybins));
  }

  Double_t xbins_ext[24];
  xbins_ext[0] = 0;
  xbins_ext[1] = 50;
  for(int i = 1; i < 22; i++) xbins_ext[i+1] = (mst_ext[i] + mst_ext[i-1])/2.;
  xbins_ext[23] = 8000;
  
  Double_t ybins_ext[7] = {0, 12.5, 37.5, 62.5, 89.5, 112.5, 137.5};

  vector<TH2D*> h_ext_gg, h_ext_gg_muJets, h_ext_gg_eleJets, h_ext_gg_hadronic, h_ext_gg_eleMu, h_ext_gg_diEle, h_ext_gg_diMu;
  for(unsigned int i = 0; i < triggerFires.size(); i++) {
    char buff[10];
    sprintf(buff, "_%d", i);
    TString buff_t = buff;
    h_ext_gg.push_back(new TH2D("gg_ext"+buff_t, triggerFires[i].first, 23, xbins_ext, 6, ybins_ext));
    h_ext_gg_muJets.push_back(new TH2D("gg_ext_muJets"+buff_t, triggerFires[i].first, 23, xbins_ext, 6, ybins_ext));
    h_ext_gg_eleJets.push_back(new TH2D("gg_ext_eleJets"+buff_t, triggerFires[i].first, 23, xbins_ext, 6, ybins_ext));
    h_ext_gg_hadronic.push_back(new TH2D("gg_ext_hadronic"+buff_t, triggerFires[i].first, 23, xbins_ext, 6, ybins_ext));
    h_ext_gg_eleMu.push_back(new TH2D("gg_ext_eleMu"+buff_t, triggerFires[i].first, 23, xbins_ext, 6, ybins_ext));
    h_ext_gg_diEle.push_back(new TH2D("gg_ext_diEle"+buff_t, triggerFires[i].first, 23, xbins_ext, 6, ybins_ext));
    h_ext_gg_diMu.push_back(new TH2D("gg_ext_diMu"+buff_t, triggerFires[i].first, 23, xbins_ext, 6, ybins_ext));
  }

  float lead_hOe, lead_sIetaIeta, lead_chHadIso, lead_neutralHadIso, lead_photonIso, lead_r9, lead_Et;
  float trail_hOe, trail_sIetaIeta, trail_chHadIso, trail_neutralHadIso, trail_photonIso, trail_r9, trail_Et;
    
  bool pass_36iso_22iso,
    pass_36or_22or,
    pass_26iso_18iso_m60,
    pass_26or_18or_m60,
    pass_26or_18or_m70;

  TTree * ggTree = new TTree("gg"+output_code_t+"_EvtTree", "An event tree for final analysis");
  ggTree->Branch("lead_hOe", &lead_hOe, "lead_hOe/F");
  ggTree->Branch("lead_sIetaIeta", &lead_sIetaIeta, "lead_sIetaIeta/F");
  ggTree->Branch("lead_chHadIso", &lead_chHadIso, "lead_chHadIso/F");
  ggTree->Branch("lead_neutralHadIso", &lead_neutralHadIso, "lead_neutralHadIso/F");
  ggTree->Branch("lead_photonIso", &lead_photonIso, "lead_photonIso/F");
  ggTree->Branch("lead_r9", &lead_r9, "lead_r9/F");
  ggTree->Branch("lead_Et", &lead_Et, "lead_Et/F");
  ggTree->Branch("trail_hOe", &trail_hOe, "trail_hOe/F");
  ggTree->Branch("trail_sIetaIeta", &trail_sIetaIeta, "trail_sIetaIeta/F");
  ggTree->Branch("trail_chHadIso", &trail_chHadIso, "trail_chHadIso/F");
  ggTree->Branch("trail_neutralHadIso", &trail_neutralHadIso, "trail_neutralHadIso/F");
  ggTree->Branch("trail_photonIso", &trail_photonIso, "trail_photonIso/F");
  ggTree->Branch("trail_r9", &trail_r9, "trail_r9/F");
  ggTree->Branch("trail_Et", &trail_Et, "trail_Et/F");
  ggTree->Branch("pass_36iso_22iso", &pass_36iso_22iso, "pass_36iso_22iso/O");
  ggTree->Branch("pass_36or_22or", &pass_36or_22or, "pass_36or_22or/O");
  ggTree->Branch("pass_26iso_18iso_m60", &pass_26iso_18iso_m60, "pass_26iso_18iso_m60/O");
  ggTree->Branch("pass_26or_18or_m60", &pass_26or_18or_m60, "pass_26or_18or_m60/O");
  ggTree->Branch("pass_26or_18or_m70", &pass_26or_18or_m70, "pass_26or_18or_m670/O");
  
  Long64_t nEntries = fTree->GetEntries();
  cout << "Total events in files : " << nEntries << endl;
  cout << "Events to be processed : " << processNEvents << endl;

  // start event looping
  Long64_t jentry = 0;
  while(jentry != processNEvents && event.getEntry(jentry++) != 0) {

    if(printLevel > 0 || (printInterval > 0 && (jentry >= printInterval && jentry%printInterval == 0))) {
      cout << int(jentry) << " events processed with run = "
	   << event.runNumber << ", event = " << event.eventNumber << ", integ lumi = " 
	   << event.intgRecLumi << endl;
    }

    nCnt[0]++; // events

    // truth-match what the ttbar decay is
    int nPhotons_fromBino = 0;
    int nElectrons_fromW = 0;
    int nMuons_fromW = 0;
    int nQuarks_fromW = 0;

    vector<susy::Photon*> reco_photons;

    for(vector<susy::Particle>::iterator it = event.genParticles.begin(); it != event.genParticles.end(); it++) {
      
      // photon from bino
      if(fabs(it->pdgId) == 22 && it->status == 1 && fabs(event.genParticles[it->motherIndex].pdgId) == 1000022) {
	nPhotons_fromBino++;

	map<TString, vector<susy::Photon> >::iterator phoMap = event.photons.find("photons");
	if(phoMap != event.photons.end()) {
	  for(vector<susy::Photon>::iterator itp = phoMap->second.begin();
	      itp != phoMap->second.end(); itp++) {
	    
	    if(is_eg(*itp, event.rho25) && itp->nPixelSeeds == 0) {
	      if(deltaR(it->momentum, itp->caloPosition) < 0.3) {
		reco_photons.push_back(&*itp);
	      }
	    }
	    
	  }
	}

      } // found gen photon
      sort(reco_photons.begin(), reco_photons.end(), EtGreater<susy::Photon>);

      // if it comes from a W
      if(fabs(event.genParticles[it->motherIndex].pdgId) == 24 && event.genParticles[it->motherIndex].status == 3) {

	int grandma_id = fabs(event.genParticles[event.genParticles[it->motherIndex].motherIndex].pdgId);
	// and if that W came from a stop or top
	if(grandma_id == 1000006 || grandma_id == 6) {

	  if(fabs(it->pdgId) == 11) nElectrons_fromW++;
	  if(fabs(it->pdgId) == 13) nMuons_fromW++;
	  if(fabs(it->pdgId) <= 6) nQuarks_fromW++;
	  if(fabs(it->pdgId) > 18) cout << "I don't know what this is, found a W --> id(" << it->pdgId << ")" << endl;
	}

      }

    }

    bool isGG = (nPhotons_fromBino == 2);
    bool isGG_hadronic = (isGG && nMuons_fromW == 0 && nElectrons_fromW == 0);
    bool isGG_muJets = (isGG && nMuons_fromW == 1 && nElectrons_fromW == 0);
    bool isGG_eleJets = (isGG && nElectrons_fromW == 1 && nMuons_fromW == 0);
    bool isGG_eleMu = (isGG && nElectrons_fromW == 1 && nMuons_fromW == 1);
    bool isGG_diEle = (isGG && nElectrons_fromW == 2 && nMuons_fromW == 0);
    bool isGG_diMu = (isGG && nElectrons_fromW == 0 && nMuons_fromW == 2);

    if(isGG && reco_photons.size() == 2) {
      lead_hOe = reco_photons[0]->hadTowOverEm;
      lead_sIetaIeta = reco_photons[0]->sigmaIetaIeta;
      lead_chHadIso = chargedHadronIso_corrected(*reco_photons[0], event.rho25);
      lead_neutralHadIso = neutralHadronIso_corrected(*reco_photons[0], event.rho25);
      lead_photonIso = photonIso_corrected(*reco_photons[0], event.rho25);
      lead_r9 = reco_photons[0]->r9;
      lead_Et = reco_photons[0]->momentum.Et();

      trail_hOe = reco_photons[1]->hadronicOverEm;
      trail_sIetaIeta = reco_photons[1]->sigmaIetaIeta;
      trail_chHadIso = chargedHadronIso_corrected(*reco_photons[1], event.rho25);
      trail_neutralHadIso = neutralHadronIso_corrected(*reco_photons[1], event.rho25);
      trail_photonIso = photonIso_corrected(*reco_photons[1], event.rho25);
      trail_r9 = reco_photons[1]->r9;
      trail_Et = reco_photons[1]->momentum.Et();
      
      pass_36iso_22iso = PassTrigger("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
      pass_36or_22or = PassTrigger("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v");
      
      pass_26iso_18iso_m60 = PassTrigger("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v");
      pass_26or_18or_m60 = PassTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v");
      pass_26or_18or_m70 = PassTrigger("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass70_v");

      ggTree->Fill();
    }

    for(susy::TriggerMap::iterator it = event.hltMap.begin(); it != event.hltMap.end(); it++) {
      for(unsigned int i = 0; i < triggerFires.size(); i++) {

	if(it->first.Contains(triggerFires[i].first) && (int(it->second.second))) {
	  triggerFires[i].second++;

	  if(scan == "stop-bino") {
	    if(isGG) h_gg[i]->Fill(index1, index2);
	    if(isGG_hadronic) h_gg_hadronic[i]->Fill(index1, index2);
	    if(isGG_muJets) h_gg_muJets[i]->Fill(index1, index2);
	    if(isGG_eleJets) h_gg_eleJets[i]->Fill(index1, index2);
	    if(isGG_eleMu) h_gg_eleMu[i]->Fill(index1, index2);
	    if(isGG_diEle) h_gg_diEle[i]->Fill(index1, index2);
	    if(isGG_diMu) h_gg_diMu[i]->Fill(index1, index2);
	  }

	  else if(scan == "stop-bino extended") {
	    if(isGG) h_ext_gg[i]->Fill(index1, index2);
	    if(isGG_hadronic) h_ext_gg_hadronic[i]->Fill(index1, index2);
	    if(isGG_muJets) h_ext_gg_muJets[i]->Fill(index1, index2);
	    if(isGG_eleJets) h_ext_gg_eleJets[i]->Fill(index1, index2);
	    if(isGG_eleMu) h_ext_gg_eleMu[i]->Fill(index1, index2);
	    if(isGG_diEle) h_ext_gg_diEle[i]->Fill(index1, index2);
	    if(isGG_diMu) h_ext_gg_diMu[i]->Fill(index1, index2);
	  }

	}

      }
      
    }

  } // end event loop

  cout << endl << "Out of " << nCnt[0] << " events:" << endl;

  for(unsigned int i = 0; i < triggerFires.size(); i++) cout << triggerFires[i].first << "  ==  " << triggerFires[i].second << endl;

  out->Write();
  out->Close();

}
