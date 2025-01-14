#include "ExampleRun.h"

ExampleRun::ExampleRun() {}
ExampleRun::~ExampleRun() {}

void ExampleRun::initializeAnalyzer() {
    //==== Example 1
    //==== Dimuon Z-peak events with two muon IDs, with systematics
    
    // if you use "--userflags RunSyst" with SKFlat.py, HasFlag("RunSyst") will return true
    //RunSyst = HasFlag("RunSyst");
    RunSyst = false;
    cout << "[ExampleRun::initializeAnalyzer] RunSyst = " << RunSyst << endl;

    // Dimuon Z-peak with two muon IDs
    // "RVec<TString> MuonIDs;" is defined in Analyzers/include/ExampleRun.h
    MuonIDs = {"POGMedium", "POGTight"};
    MuonIDSFKeys = {"NUM_MediumID_DEN_TrackerMuons", "NUM_TightID_DEN_TrackerMuons"};

    // At this point, sample informations (e.g. IsDATA, DataStream, MCSample, or DataEra) are all set
    // You can define sample-dependent or era-dependent variables here
    // example: era-dependent variables
    // "TString IsoMuTriggerName;" and "float TriggerSafePtCut;" are defined in Analyzers/include/ExampleRun.h
    // IsoMuTriggerName is an era-dependent variable, which is used throhughout the events(let's make it global variable)
    if (DataEra == "2022") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    }
    else if (DataEra == "2022EE") {
        IsoMuTriggerName = "HLT_IsoMu24";
        TriggerSafePtCut = 26.;
    }
    else if (DataEra == "2023") {
        IsoMuTriggerName = "";
        TriggerSafePtCut = 26.;
    }
    else {
        cerr << "[ExampleRun::initializeAnalyzer] DataEra is not set properly" << endl;
        exit(EXIT_FAILURE);
    }

    cout << "[ExampleRun::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
    cout << "[ExampleRun::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

    // B-tagging
    // add taggers and WP that you wnat to use in analysis
    // Not implemented yet
    
    // Correction
    myCorr = new Correction(DataEra, IsDATA?DataStream:MCSample ,IsDATA);

    //==== Example 2
    //==== Using new PDF
    //==== It consumes so much time, so only being activated with --userflags RunNewPDF
    //RunNewPDF = HasFlag("RunNewPDF");
    RunNewPDF = false;
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

    // Declare parameter set
    ExampleParameter param;

    // Loop over muonIDs
    for (unsigned int it_MuonID=0; it_MuonID < MuonIDs.size(); it_MuonID++) {
        const TString MuonID = MuonIDs.at(it_MuonID);
        const TString MuonIDSFKey = MuonIDSFKeys.at(it_MuonID);

        // 1) First, let's run central values
        // clear parameter set
        param.Clear();

        // set which systematic sources you want to run this time
        param.syst = ExampleParameter::SYST::Central;

        // set name of the parameter set
        param.Name = MuonID+"_"+param.GetSystType();

        // you can define lepton ID string here
        param.Muon_Tight_ID = MuonID;
        param.Muon_ID_SF_Key = MuonIDSFKey;

        // and jet ID
        param.Jet_ID = "tight";

        // now all parameters are set, run executeEventFromParmeter(param)
        executeEventFromParameter(param);

        // 2) Now, loop over systematic sources
        if (RunSyst) {
            for (unsigned int syst=1; syst < ExampleParameter::SYST::NSyst; syst++) {
                param.syst = syst; 
                param.Name = MuonID+"_"+"Syst_"+param.GetSystType();
                executeEventFromParameter(param);
            }
        }
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

void ExampleRun::executeEventFromParameter(ExampleParameter param) {
    // No cut
    //FillHist(param.Name+"/NoCut_"+param.Name, 0., 1., 1, 0., 1.);

    // No MET filter for NanoAODv12? 
    Event ev = GetEvent();
    //Particle METv = ev.GetMETVector()
    //cout << "[ExampleRun::executeEventFromParameter] PassTrigger = " << ev.PassTrigger(IsoMuTriggerName) << endl;
    if (! (ev.PassTrigger(IsoMuTriggerName))) return;
    

    // Copy All objects
    RVec<Muon> this_AllMuons = AllMuons;
    //RVec<Jet> this_AllJets = AllJets;
    
    // For each systematic sources,
    // 1. Smear or scale them
    // 2. apply ID selections
    // This order should be explicitly followed
    // below are all varibales for available systematic sources
    // Update systematic sources!
    if (param.syst == ExampleParameter::SYST::Central) {}
    else {
        cerr << "[ExampleRun::executeEventFromParameter] syst = " << param.syst << " is not implemented" << endl;
        exit(EXIT_FAILURE);
    }

    // apply ID selections using this_AllXXX
    RVec<Muon> muons = SelectMuons(this_AllMuons, param.Muon_Tight_ID, 20., 2.4);
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

    // Event weight
    float weight = 1.;
    if (!IsDATA) {
        weight *= MCweight();
        weight *= ev.GetTriggerLumi("Full");
        weight *= weight_Prefire;

        // example of applying muon scale factors
        for (const auto &muon: muons) {
            float this_idsf = myCorr->GetMuonIDSF(param.Muon_ID_SF_Key, muon.Eta(), muon.Pt());
            float this_isosf = 1.;
            weight *= this_idsf * this_isosf;
        }
    }

    // fill histograms
    FillHist(param.Name+"/ZCand_Mass_"+param.Name, ZCand.M(), weight, 50, 70., 110.);
}
