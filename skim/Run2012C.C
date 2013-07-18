void Run2012C(TString discriminant = "CSVM", bool isMC = false) {

  gROOT->Reset();
  gSystem->Load("libSusyEvent.so");

  gROOT->LoadMacro("SusyEventAnalyzer.cc++");

  TChain chain("susyTree");
  //chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v0p1/Run2012A-22Jan2013-v1/Photon/*.root");
  //chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v0p1/Run2012B-22Jan2013-v1/DoublePhoton/susyEvents_SOURCE_NUMBER*.root");
  chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v0p1/Run2012C-22Jan2013-v2/DoublePhoton/susyEvents_SOURCE_NUMBER*.root");
  //chain.Add("/eos/uscms/store/user/lpcpjm/SusyNtuples/cms538v0p1/Run2012D-22Jan2013-v1/DoublePhoton/*.root");

  chain.SetBranchStatus("*", 1);

  if(chain.LoadTree(0) != 0) {
    cerr << "Error with input chain. Do the files exist?" << endl;
    return;
  }

  SusyEventAnalyzer* sea = new SusyEventAnalyzer(chain);

  // configuration parameters
  // any values given here will replace the default values
  sea->SetPrintInterval(1e5);             // print frequency
  sea->SetPrintLevel(0);                  // print level for event contents

  sea->SetOutput("Run2012C-22Jan2013-DoublePhoton_SOURCE_NUMBER");

  std::vector<TString> eg_names;
  eg_names.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
  std::vector<int> eg_types;
  eg_types.push_back(1);
  eg_types.push_back(2);
  eg_types.push_back(3);
  sea->AddHlt(eg_names, eg_types);

  std::vector<TString> f_names;
  f_names.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50_v");
  f_names.push_back("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85_v");
  f_names.push_back("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50_v");
  f_names.push_back("HLT_Photon36_R9Id85_Photon22_R9Id85_v");
  std::vector<int> f_types;
  f_types.push_back(4);
  sea->AddHlt(f_names, f_types);

  sea->SetProcessNEvents(-1);      	  // number of events to be processed
  
  sea->IncludeAJson("Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt");

  sea->SetIsMC(isMC);
  
  sea->SetUseSyncFile(false);
  //sea->IncludeSyncFile("synchro/dmorse_ff.txt_not_brian_ff_nojet.txt");
  sea->SetCheckSingleEvent(false);
  sea->AddCheckSingleEvent(196203, 33, 27883630);

  sea->SetBtagger(discriminant);

  sea->AddValidTagger("TCHPT");
  sea->AddValidTagger("JPL");
  sea->AddValidTagger("JPM");
  sea->AddValidTagger("JPT");
  sea->AddValidTagger("CSVL");
  sea->AddValidTagger("CSVM");
  sea->AddValidTagger("CSVT");
  
  TStopwatch ts;

  ts.Start();
  
  sea->Filter();
  
  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;
  
}
