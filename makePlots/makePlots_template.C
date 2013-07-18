void makePlots() {

  gSystem->Load("libRooFit.so");
  gSystem->Load("libRooFitCore.so");
  gSystem->SetIncludePath("-I../../../.. -I/uscmst1/prod/sw/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms6/include");
  gSystem->Load("/uscms_data/d2/bfrancis/btagRA3/CMSSW_5_3_3/src/PhysicsTools/lib/libTagAndProbe.so");
  gSystem->Load("../libSusyEvent.so");

  gROOT->LoadMacro("analyze.C+");

  TStopwatch ts;
  ts.Start();

  TString input = "FILE_TO_RUN";
  bool addMC = true;
  TString intLumi = "19.789";
  int intLumi_int = 19789;
  //4.041/fb ichep
  //19.499/fb 2012
  //19.789/fb 2012 rereco

  bool useFF = true;
  bool useDifferenceSystematic = false;
  bool useTTGJets = false;
  bool useMCforQCD = false;

  for(int i = 0; i < 7; i++) {
//  for(int i = 0; i < 1; i++) {
    analyze(input, addMC, i, intLumi, intLumi_int, useFF, useDifferenceSystematic, useTTGJets, useMCforQCD);
  }  

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
