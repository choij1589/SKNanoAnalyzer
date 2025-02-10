#include "MeasureJetTaggingEff.h"
#include "JetTaggingParameter.h"
//////////

MeasureJetTaggingEff::MeasureJetTaggingEff()
{
} // MeasureJetTaggingEff::MeasureJetTaggingEff

//////////

MeasureJetTaggingEff::~MeasureJetTaggingEff()
{
} // MeasureJetTaggingEff::~MeasureJetTaggingEff()

//////////

void MeasureJetTaggingEff::initializeAnalyzer()
{
    fChain->SetBranchStatus("*", 0);
    fChain->SetBranchStatus("Jet_*", 1);
    fChain->SetBranchStatus("Muon_*", 1);
    fChain->SetBranchStatus("Pileup_*", 1);
    fChain->SetBranchStatus("PV_*", 1);
    fChain->SetBranchStatus("genWeight", 1);

    myCorr = new MyCorrection(DataEra, IsDATA?DataStream:MCSample ,IsDATA);
    TString datapath = getenv("DATA_DIR");
    TString btagpath = datapath + "/" + DataEra + "/BTag/";

    Taggers.clear();
    WPs.clear();

    Taggers.push_back(JetTagging::JetFlavTagger::DeepJet);
    Taggers.push_back(JetTagging::JetFlavTagger::ParticleNet);
    Taggers.push_back(JetTagging::JetFlavTagger::ParT);
    WPs.push_back(JetTagging::JetFlavTaggerWP::Loose);
    WPs.push_back(JetTagging::JetFlavTaggerWP::Medium);
    WPs.push_back(JetTagging::JetFlavTaggerWP::Tight);
    WPs.push_back(JetTagging::JetFlavTaggerWP::VeryTight);
    WPs.push_back(JetTagging::JetFlavTaggerWP::SuperTight);





    vec_etabins = {0.0, 0.8, 1.6, 2., 2.5};
    vec_ptbins = {20., 25., 30., 50., 70., 100., 140., 200., 300., 600., 1000.}; // PT bins used in POG SF measurements
    // for average users, this binning will be sufficient.
    // but eta-dependence of efficiency can be larger for |eta|>~2, where track & muon detector information of jet constituents starts to get lost, which is critical in tagging.
    // precision analysis with high-eta b may use finer binnings there, but beware of small number of b-jets in high-eta, high-pt bins if you use ttbar sample; proper optimization of bin size should be studied.

    PtMax = vec_ptbins.at(vec_ptbins.size() - 1);
    NEtaBin = vec_etabins.size() - 1;
    NPtBin = vec_ptbins.size() - 1;

    etabins = new float[NEtaBin + 1];
    for (int i = 0; i < NEtaBin + 1; i++)
        etabins[i] = vec_etabins.at(i);
    ptbins = new float[NPtBin + 1];
    for (int i = 0; i < NPtBin + 1; i++)
        ptbins[i] = vec_ptbins.at(i);
} // void MeasureJetTaggingEff::initializeAnalyzer

//////////

void MeasureJetTaggingEff::executeEvent()
{

    Event ev = GetEvent();
    RVec<Jet> AllJets = GetAllJets();
    RVec<Muon> AllMuons = GetAllMuons();
    if(!PassJetVetoMap(AllJets,AllMuons)) return;

    RVec<Electron> AllElectrons =  GetAllElectrons();
    AllElectrons = SelectElectrons(AllElectrons, Electron::ElectronID::POG_TIGHT, 30., 2.5);
    myCorr->GetElectronRECOSF(AllElectrons);

    float JetPtCut = 25.;
    float JetEtaCut = 2.4;

    //define your analysis phase space here
    throw std::runtime_error("Please define your analysis phase space");
    RVec<Jet> jets = SelectJets(AllJets, Jet::JetID::TIGHT, 25., 2.5);
    if (jets.size() < 4)
        return;

    float weight = 1.;
    float w_Gen = MCweight();
    float w_Norm = ev.GetTriggerLumi("Full");
    float w_PU = myCorr->GetPUWeight(ev.nTrueInt()); 
    weight *= w_Gen * w_Norm * w_PU;
    // tagging performance depends on PU, so it is better reweight to proper PU profile

    auto isWPAvailable= [&](const std::string flav, const std::string tagger, const JetTagging::JetFlavTaggerWP wp, const int run){ 
        if(wp != JetTagging::JetFlavTaggerWP::VeryTight && wp != JetTagging::JetFlavTaggerWP::SuperTight) return true;
        if(flav == "c") return false;
        if(run == 2) return false;
        return true;
    };
    //==== code to measure btag efficiencies in TT MC
    //==== Reference : https://github.com/rappoccio/usercode/blob/Dev_53x/EDSHyFT/plugins/BTaggingEffAnalyzer.cc
    for (unsigned int ij = 0; ij < jets.size(); ij++)
    {
        TString flav = "0";
        if (fabs(jets.at(ij).hadronFlavour()) == 4)
            flav = "4";
        if (fabs(jets.at(ij).hadronFlavour()) == 5)
            flav = "5";

        float this_Eta = fabs(jets.at(ij).Eta());                                // POG recommendation is to use |eta|
        float this_Pt = jets.at(ij).Pt() < PtMax ? jets.at(ij).Pt() : PtMax - 1; // put overflows in the last bin

        //==== First, fill the denominator
        FillHist(string("tagging#b") + "##era#" + DataEra.Data() + "##flavor#" + string(flav) + "##systematic#central##den", this_Eta, this_Pt, weight, NEtaBin, etabins, NPtBin, ptbins);
        FillHist(string("tagging#c") + "##era#" + DataEra.Data() + "##flavor#" + string(flav) + "##systematic#central##den", this_Eta, this_Pt, weight, NEtaBin, etabins, NPtBin, ptbins);
        FillHist("DeepJetBTaggingScore"+flav, jets.at(ij).GetBTaggerResult(JetTagging::JetFlavTagger::DeepJet), weight, 100, 0, 1);
        FillHist("ParticleNetBTaggingScore"+flav, jets.at(ij).GetBTaggerResult(JetTagging::JetFlavTagger::ParticleNet), weight, 100, 0, 1);
        FillHist("ParTBTaggingScore"+flav, jets.at(ij).GetBTaggerResult(JetTagging::JetFlavTagger::ParT), weight, 100, 0, 1);
        for(unsigned int i_tag=0; i_tag < Taggers.size(); i_tag++){
            JetTagging::JetFlavTagger this_tagger = Taggers.at(i_tag);
            //============ b-tagging
            for(unsigned int i_wp=0; i_wp < WPs.size(); i_wp++){
                JetTagging::JetFlavTaggerWP this_wp = WPs.at(i_wp);
                if (!isWPAvailable("b", JetTagging::GetTaggerCorrectionLibStr(this_tagger).Data(), this_wp, Run)) continue;
                myCorr->SetTaggingParam(this_tagger, this_wp);
                float this_bTaggingCut = myCorr->GetBTaggingWP();
                if (jets.at(ij).GetBTaggerResult(this_tagger) > this_bTaggingCut)
                    FillHist(string("tagging#b") + "##era#" + DataEra.Data() + "##tagger#" + JetTagging::GetTaggerCorrectionLibStr(this_tagger).Data() + "##working_point#" + JetTagging::GetTaggerCorrectionWPStr(this_wp).Data() + "##flavor#" + string(flav) + "##systematic#central##num", this_Eta, this_Pt, weight, NEtaBin, etabins, NPtBin, ptbins);
            }
            //============ c-tagging
            for(unsigned int i_wp=0; i_wp < WPs.size(); i_wp++){
                JetTagging::JetFlavTaggerWP this_wp = WPs.at(i_wp);
                if (!isWPAvailable("c", JetTagging::GetTaggerCorrectionLibStr(this_tagger).Data(), this_wp, Run)) continue;
                myCorr->SetTaggingParam(this_tagger, this_wp);
                float this_CvBCut = myCorr->GetCTaggingWP().first;
                float this_CvLCut = myCorr->GetCTaggingWP().second;
                if (jets.at(ij).GetCTaggerResult(this_tagger).first > this_CvBCut && jets.at(ij).GetCTaggerResult(this_tagger).second > this_CvLCut)
                    FillHist(string("tagging#c") + "##era#" + DataEra.Data() + "##tagger#" + JetTagging::GetTaggerCorrectionLibStr(this_tagger).Data() + "##working_point#" + JetTagging::GetTaggerCorrectionWPStr(this_wp).Data() + "##flavor#" + string(flav) + "##systematic#central##num", this_Eta, this_Pt, weight, NEtaBin, etabins, NPtBin, ptbins);
            }


        }
    }
} // void MeasureJetTaggingEff::executeEvent
