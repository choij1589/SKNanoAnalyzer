#include "GetEffLumi.h"
//#include "correction.h"
//using correction::CorrectionSet;
GetEffLumi::GetEffLumi(){};
GetEffLumi::~GetEffLumi(){};

void GetEffLumi::initializeAnalyzer(){
  useTH1F = false;
  cout << "[GetEffLumi] Using TH1D for histograms." << endl;

  myCorr = new MyCorrection(DataEra, DataPeriod, IsDATA?DataStream:MCSample ,IsDATA);
  fChain->SetBranchStatus("*", 0);
  if(!IsDATA){
    fChain->SetBranchStatus("genWeight", 1);

    // Read generator filter efficiency from LuminosityBlocks tree
    long long totalPassed = 0, totalTotal = 0;
    TObjArray *fileElements = fChain->GetListOfFiles();
    for (int iFile = 0; iFile < fileElements->GetEntries(); iFile++) {
      TString fileName = fileElements->At(iFile)->GetTitle();
      TFile *file = TFile::Open(fileName);
      if (!file || file->IsZombie()) {
        if (file) file->Close();
        continue;
      }
      TTree *lumiTree = (TTree *)file->Get("LuminosityBlocks");
      if (!lumiTree || !lumiTree->GetBranch("GenFilter_numEventsPassed")) {
        file->Close();
        continue;
      }
      Int_t passed, total;
      lumiTree->SetBranchAddress("GenFilter_numEventsPassed", &passed);
      lumiTree->SetBranchAddress("GenFilter_numEventsTotal", &total);
      for (Long64_t j = 0; j < lumiTree->GetEntries(); j++) {
        lumiTree->GetEntry(j);
        totalPassed += passed;
        totalTotal += total;
      }
      file->Close();
    }
    genFilterPassed = totalPassed;
    genFilterTotal = totalTotal;
    if (totalTotal > 0) {
      double filterEff = (double)totalPassed / (double)totalTotal;
      cout << "[GetEffLumi] GenFilter efficiency: " << totalPassed << " / " << totalTotal << " = " << filterEff << endl;
    }
  }
  sumW = 1.;
  sumSign = 1.;
}
void GetEffLumi::executeEvent() {
  double weight = 1.;
  double weight_sign = 1.;
  if(!IsDATA)
    {
      weight = MCweight(false,false);
      weight_sign = MCweight(true,false);
      FillHist("NEvents", 0, 1, 1, 0., 1.);
      FillHist("sumW", 0, weight, 1, 0., 1.);
      FillHist("sumSign", 0, weight_sign, 1, 0., 1.);
    }
  else{
    FillHist("NEvents", 0, 1, 1, 0., 1.);
  }


  return;
}

void GetEffLumi::WriteHist() {
  cout << "[GetEffLumi::WriteHist] genFilterPassed=" << genFilterPassed
       << " genFilterTotal=" << genFilterTotal << endl;
  if (genFilterTotal > 0) {
    GetOutfile()->cd();
    TH1D *hPassed = new TH1D("GenFilter_numEventsPassed", "", 1, 0., 1.);
    hPassed->SetBinContent(1, (double)genFilterPassed);
    hPassed->Write();
    TH1D *hTotal = new TH1D("GenFilter_numEventsTotal", "", 1, 0., 1.);
    hTotal->SetBinContent(1, (double)genFilterTotal);
    hTotal->Write();
    delete hPassed;
    delete hTotal;
  }
  AnalyzerCore::WriteHist();
}

