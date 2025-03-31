#include "ExampleRun.h"

ExampleRun::ExampleRun() {}
ExampleRun::~ExampleRun() {}

void ExampleRun::initializeAnalyzer() {
    //==== Example 1
    //==== Dimuon Z-peak events with two muon IDs, with systematics
    
    // if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return true
    RunSyst = HasFlag("RunSyst");
    cout << "[ExampleRun::initializeAnalyzer] RunSyst = " << RunSyst << endl;

    // Dimuon Z-peak with two muon IDs
    // "RVec<TString> MuonIDs;" is defined in Analyzers/include/ExampleRun.h
    MuonIDs = {Muon::MuonID::POG_MEDIUM, Muon::MuonID::POG_TIGHT};
    MuonIDSFKeys = {"NUM_MediumID_DEN_TrackerMuons", "NUM_TightID_DEN_TrackerMuons"};

    // At this point, sample informations (e.g. IsDATA, DataStream, MCSample, or DataEra) are all set
    // You can define sample-dependent or era-dependent variables here
    // example: era-dependent variables
    // "TString IsoMuTriggerName;" and "float TriggerSafePtCut;" are defined in Analyzers/include/ExampleRun.h
    // IsoMuTriggerName is an era-dependent variable, which is used throhoughout the events(let's make it global variable)
    if (DataEra == "2016preVFP" || DataEra == "2016postVFP" || DataEra == "2018") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2017") {
        IsoMuTriggerName = "HLT_IsoMu27";
        TriggerSafePtCut = 29.;
    } else if (DataEra == "2022") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2022EE") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    } else if (DataEra == "2023") {
        IsoMuTriggerName = "";
        TriggerSafePtCut = 26.;
    } else {
        cerr << "[ExampleRun::initializeAnalyzer] DataEra is not set properly" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "[ExampleRun::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
    cout << "[ExampleRun::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

    // B-tagging
    // add taggers and WP that you wnat to use in analysis
    // Not implemented yet
    
    // Correction
    myCorr = new MyCorrection(DataEra, IsDATA?DataStream:MCSample ,IsDATA);
    // SystematicHelper
    string SKNANO_HOME = getenv("SKNANO_HOME");
    if (IsDATA) {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/noSyst.yaml", DataStream);
    } else {
        systHelper = std::make_unique<SystematicHelper>(SKNANO_HOME + "/docs/ExampleSystematic.yaml", MCSample);
    }

    //==== Example 2
    //==== Using new PDF
    //==== It consumes so much time, so only being activated with --userflags RunNewPDF
    RunNewPDF = HasFlag("RunNewPDF");
    cout << "[ExampleRun::initializeAnalyzer] RunNewPDF = " << RunNewPDF << endl;
    if (RunNewPDF && !IsDATA) {
        LHAPDFHandler LHAPDFHandler_Prod;
        LHAPDFHandler_Prod.CentralPDFName = "NNPDF31_nnlo_hessian_pdfas";
        LHAPDFHandler_Prod.init();

        LHAPDFHandler LHAPDFHandler_New;
        LHAPDFHandler_New.CentralPDFName = "NNPDF31_nlo_hessian_pdfas";
        LHAPDFHandler_New.ErrorSetMember_Start = 1;
        LHAPDFHandler_New.ErrorSetMember_End = 100;
        LHAPDFHandler_New.AlphaSMember_Down = 101;
        LHAPDFHandler_New.AlphaSMember_Up = 102;
        LHAPDFHandler_New.init();
        
        pdfReweight->SetProdPDF( LHAPDFHandler_Prod.PDFCentral );
        pdfReweight->SetNewPDF( LHAPDFHandler_New.PDFCentral );
        pdfReweight->SetNewPDFErrorSet( LHAPDFHandler_New.PDFErrorSet );
        pdfReweight->SetNewPDFAlphaS( LHAPDFHandler_New.PDFAlphaS_Down, LHAPDFHandler_New.PDFAlphaS_Up );
    }

    //==== Example 3
    //==== How to estimate xsec errors (PDF & Scale)
    //==== For example, MET
    //RunXSecSyst = HasFlag("RunXSecSyst");
    RunXSecSyst = false;
    cout << "[ExampleRun::initializeAnalyzer] RunXSecSyst = " << RunXSecSyst << endl;
}

void ExampleRun::executeEvent() {
    //==== Example 1
    //==== Dimuon Z-peak events with two muon IDs, with systematics
    
    // *IMPORTANT TO SAVE CPU TIME*
    // Every GetMuon() funtion first collect ALL NANOAOD muons with GetAllMuons(), and then check ID booleans.
    // GetAllMuons not only loops over all NANOAOD muons, but also actually CONSTRUCT muon objects for each muons.
    // We are now running systematics, and you don't want to do this for every systematic sources
    // So, I defined "RVec<Muon> AllMuons;" in Analyzers/include/ExampleRun.h,
    // and save muons objects at the very beginning of executeEvent().
    // Later, do "SelectMuons(AllMuons, ID, pt, eta)" to get muons with ID cuts
    AllMuons = GetAllMuons();
    //AllJets = GetAllJets();
    
    // No prefire weight for Run3?
    weight_Prefire = 1.;  //조사해보기


    // Loop over systematic sources that require separate event loop
    for (const auto &syst_dummy : *systHelper)
    {
        executeEventFromParameter();
    }

    //==== Example 2
    //==== Using new PDF
    if (RunNewPDF && !IsDATA) {
        // cout << "[ExampleRun::executeEvent] PDF reweight = " << GetPDFReweight() << endl;
        FillHist("NewPDF_PDFReweight", GetPDFReweight(), 1., 2000, 0.90, 1.10);
        //cout << "[ExampleRun::executeEvent] PDF reweight for error set (NErrorSet = "<<pdfReweight->NErrorSet<< ") :" << endl;
        for(int i=0; i<pdfReweight->NErrorSet; i++){
            //cout << "[ExampleRun::executeEvent]   " << GetPDFReweight(i) << endl;
            FillHist("NewPDF_PDFErrorSet/PDFReweight_Member_"+TString::Itoa(i,10), GetPDFReweight(i), 1., 2000, 0.90, 1.10);
        }
    }

    //==== Example 3
    //==== How to estimate xsec errors (PDF & Scale)
    //==== For example, MET
    //==== Need Update
}

void ExampleRun::executeEventFromParameter() {
    // Assign systematic sources, In this example, only muon ID scale factors are considered
    Muon::MuonID this_muon_id;
    TString this_muon_id_sf_key;

    const TString this_syst = systHelper->getCurrentSysName();
    if(this_syst == "Central"){
        this_muon_id = MuonIDs[0];
        this_muon_id_sf_key = MuonIDSFKeys[0];        
    } else if (this_syst == "Muon_ID_Tight_Up") {
        // In this case, change of Muon ID is considered as systematic variation.
        // Of course this only has one varation, not up and down
        // I turned on OneSided option in the yaml file, then it variates in "Up" variation and stays same in "Down" variation
        this_muon_id = MuonIDs[1];
        this_muon_id_sf_key = MuonIDSFKeys[1];
    } else{
        //cerr << "[ExampleRun::executeEventFromParameter] check syst.yaml to set the systematic source for " << this_syst << endl;
        //exit(EXIT_FAILURE);
        return;
    }

    std::unordered_map<std::string, std::variant<std::function<float(MyCorrection::variation, TString)>, std::function<float()>>> weight_function_map;

    // No cut
    FillHist(this_syst+"/NoCut", 0., 1., 1, 0., 1.);

    // No MET filter for NanoAODv12? 
    Event ev = GetEvent();
    //Particle METv = ev.GetMETVector()
    if (! (ev.PassTrigger(IsoMuTriggerName))) return;
    

    // Copy All objects
    RVec<Muon> this_AllMuons = AllMuons;
    //RVec<Jet> this_AllJets = AllJets;

    // apply ID selections using this_AllXXX
    RVec<Muon> muons = SelectMuons(this_AllMuons, this_muon_id, 20., 2.4);
    //RVec<Jet> jets = SelectJets(this_AllJets, param.Jet_ID, 30., 2.4);
    
    // sort in pt-order
    // 1) leptons : after scaling/smearing, pt ordring can differ from NANOAOD
    sort(muons.begin(), muons.end(), PtComparing);
    // jets : similar, but also when applying new JEC, ordering is changes. This is important if you use leading jets
    //sort(jets.begin(), jets.end(), PtComparing);

    // b-tagging
    // need update
    // int NBJets_NoSF(0), NBJets_WithSF_2a(0);
    // JetTagging::Parameters jtp_DeepCSV_Medium = JetTagging::Parameters(JetTagging::DeepCSV, JetTagging::Medium, JetTagging::incl, JetTagging::comb);
    

    // Event selection
    // dimuon
    if (muons.size() != 2) return;
    // leading muon has to pass trigger-safe cut
    if (muons.at(0).Pt() <= TriggerSafePtCut) return;
    // On-Z
    Particle ZCand = muons.at(0) + muons.at(1);
    if (! (fabs(ZCand.M() - 91.2) < 15.)) return;

    // example of applying muon scale factors
    auto mu_id_lambda = [&](MyCorrection::variation syst, TString source)
    { return myCorr->GetMuonIDSF(this_muon_id_sf_key, muons, syst, source); };
    weight_function_map["Muon_ID"] = mu_id_lambda;
    systHelper->assignWeightFunctionMap(weight_function_map);

    // Event weight
    float weight = 1.;
    if (!IsDATA) {
        weight *= MCweight();
        weight *= ev.GetTriggerLumi("Full");
        weight *= weight_Prefire;
        float default_weight = weight;
        // Below line will caculate weight for each systematic sources
        unordered_map<std::string, float> weight_map = systHelper->calculateWeight();
        for (const auto &w : weight_map)
        {
            // fill histograms
            FillHist(this_syst+"/ZCand_Mass_"+w.first, ZCand.M(), default_weight*w.second, 50, 70., 110.);
        }
    } else {
        // fill histograms
        FillHist(this_syst+"/ZCand_Mass_"+this_syst, ZCand.M(), weight, 50, 70., 110.);
    }
}
