#include "AnalyzerCore.h"

AnalyzerCore::AnalyzerCore() {
    myCorr = nullptr;
    outfile = nullptr;
    if(HasFlag("useTH1F")) {
        cout << "[AnalyzerCore::AnalyzerCore] Using TH1F" << endl;
        useTH1F = true;
    }
    else {
        cout << "[AnalyzerCore::AnalyzerCore] Using TH1D" << endl;
        useTH1F = false;
    }
    // pdfReweight = new PDFReweight();
}

AnalyzerCore::~AnalyzerCore() {
    for (const auto &pair: histmap1d) delete pair.second; histmap1d.clear();
    for (const auto &pair: histmap2d) delete pair.second; histmap2d.clear();
    for (const auto &pair: histmap3d) delete pair.second; histmap3d.clear();
    if (outfile) delete outfile;
    if (myCorr) delete myCorr;
    // if (pdfReweight) delete pdfReweight;
}

// https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
bool AnalyzerCore::PassNoiseFilter(const RVec<Jet> &Alljets, const Event &ev, const Event::MET_Type MET_type){
    bool passNoiseFilter = true;
    if (Run == 2) {
        passNoiseFilter = Flag_goodVertices
                       && Flag_globalSuperTightHalo2016Filter
                       && Flag_HBHENoiseFilter
                       && Flag_HBHENoiseIsoFilter
                       && Flag_EcalDeadCellTriggerPrimitiveFilter
                       && Flag_BadPFMuonFilter
                       && Flag_BadPFMuonDzFilter
                       && Flag_hfNoisyHitsFilter
                       && Flag_eeBadScFilter;
        if (DataEra.Contains("2017") || DataEra.Contains("2018")) {
            passNoiseFilter = passNoiseFilter && Flag_ecalBadCalibFilter;
        }
        return passNoiseFilter;
    } else { // Run3
        bool passNoiseFilter = Flag_goodVertices 
                            && Flag_globalSuperTightHalo2016Filter
                            && Flag_EcalDeadCellTriggerPrimitiveFilter
                            && Flag_BadPFMuonFilter 
                            && Flag_BadPFMuonDzFilter 
                            && Flag_hfNoisyHitsFilter 
                            && Flag_eeBadScFilter;
        if (!passNoiseFilter) return false;
    
        // Temporarily remove the Flag_ecalBadCalibFilter. Instead
        if (! IsDATA) return passNoiseFilter; // Not applying to MC
        if (! (362433 <= RunNumber && RunNumber <= 367144)) return passNoiseFilter;
    
        const Particle METv = ev.GetMETVector(MET_type); // Not re-corrected...
        if (! (METv.Pt() > 100.)) return passNoiseFilter;
        RVec<Jet> this_jet = SelectJets(Alljets, "NOCUT", 50., 0.5);
        for(const auto &jet: this_jet){
            bool badEcal = (jet.Pt() > 50);
            badEcal = badEcal && (-0.5 < jet.Eta() && jet.Eta() < -0.1);
            badEcal = badEcal && (-2.1 < jet.Phi() && jet.Phi() < -1.8);
            badEcal = badEcal && (jet.neEmEF() > 0.9 || jet.chEmEF() > 0.9) ;
            badEcal = badEcal && jet.DeltaPhi(METv) > 2.9;
            if(badEcal) return false;
        }
        return true;
    }
}

float AnalyzerCore::GetScaleVariation(const MyCorrection::variation &muF_syst, const MyCorrection::variation &muR_syst) {
    if(nLHEScaleWeight == 0) return 1.;
    if(muF_syst == MyCorrection::variation::down && muR_syst == MyCorrection::variation::down) return LHEScaleWeight[0];
    else if(muF_syst == MyCorrection::variation::nom && muR_syst == MyCorrection::variation::down) return LHEScaleWeight[1];
    else if(muF_syst == MyCorrection::variation::up && muR_syst == MyCorrection::variation::down) return LHEScaleWeight[2];
    else if(muF_syst == MyCorrection::variation::down && muR_syst == MyCorrection::variation::nom) return LHEScaleWeight[3];
    else if(muF_syst == MyCorrection::variation::up && muR_syst == MyCorrection::variation::nom) return LHEScaleWeight[4];
    else if(muF_syst == MyCorrection::variation::down && muR_syst == MyCorrection::variation::up) return LHEScaleWeight[5];
    else if(muF_syst == MyCorrection::variation::nom && muR_syst == MyCorrection::variation::up) return LHEScaleWeight[6];
    else if(muF_syst == MyCorrection::variation::up && muR_syst == MyCorrection::variation::up) return LHEScaleWeight[7];
    else if(muF_syst == MyCorrection::variation::nom && muR_syst == MyCorrection::variation::nom) return 1.f;
    else{
        cout << "[AnalyzerCore::GetScaleVariation] muF_syst = " << int(muF_syst) << " and muR_syst = " << int(muR_syst) << " is not implemented" << endl;
        exit(ENODATA);
    }
}

float AnalyzerCore::GetPSWeight(const MyCorrection::variation &ISR_syst, const MyCorrection::variation &FSR_syst){
    //[0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2[2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;
    if(nPSWeight == 1) return 1.;
    if(ISR_syst == MyCorrection::variation::up && FSR_syst == MyCorrection::variation::nom) return PSWeight[0];
    else if(ISR_syst == MyCorrection::variation::nom && FSR_syst == MyCorrection::variation::up) return PSWeight[1];
    else if(ISR_syst == MyCorrection::variation::down && FSR_syst == MyCorrection::variation::nom) return PSWeight[2];
    else if(ISR_syst == MyCorrection::variation::nom && FSR_syst == MyCorrection::variation::down) return PSWeight[3];
    else{
        cout << "[AnalyzerCore::GetPSWeight] ISR_syst = " << int(ISR_syst) << " and FSR_syst = " << int(FSR_syst) << " is not implemented" << endl;
        exit(ENODATA);
    }
}

float AnalyzerCore::GetHT(const RVec<Jet> &jets){
    float HT = 0.;
    for(const auto &jet: jets){
        HT += jet.Pt();
    }
    return HT;
}

void AnalyzerCore::METType1Propagation(Particle &MET, RVec<Particle> &original_objects, RVec<Particle> &corrected_objects){
   Particle Typle1Correction;
   Typle1Correction.SetPtEtaPhiM(0., 0., 0., 0.);
   for(const auto &obj: original_objects){
       Typle1Correction += obj;
   }
   for(const auto &obj: corrected_objects){
       Typle1Correction -= obj;
   }
    MET += Typle1Correction;
}

float AnalyzerCore::GetL1PrefireWeight(MyCorrection::variation syst){
    if(Run == 3) return 1.;
    if (IsDATA)
        return 1.;
    float weight = 1.;
    switch(syst){
        case MyCorrection::variation::nom:
            weight = L1PreFiringWeight_Nom;
            break;
        case MyCorrection::variation::up:
            weight = L1PreFiringWeight_Up;
            break;
        case MyCorrection::variation::down:
            weight = L1PreFiringWeight_Dn;
            break;
        default:
            exit(ENODATA);
    }
    return weight;
}

unordered_map<int, int> AnalyzerCore::GenJetMatching(const RVec<Jet> &jets, const RVec<GenJet> &genjets, const float &rho, const float dR, const float pTJerCut){
    unordered_map<int, int> matched_genjet_idx;
    RVec<tuple<int,int,float,float>> possible_matches;

    //All possible matches that pass the dR and pTJerCut
    for(size_t i = 0; i < jets.size(); i++){
        for(size_t j = 0; j < genjets.size(); j++){
            float this_DeltaR = jets[i].DeltaR(genjets[j]);
            float this_pt_diff = fabs(jets[i].Pt() - genjets[j].Pt());
            float this_jer = myCorr->GetJER(jets[i].Eta(), jets[i].Pt(), rho);
            this_jer *= jets[i].Pt();
            if(jets[i].DeltaR(genjets[j]) < dR && this_pt_diff < pTJerCut*this_jer){
                possible_matches.emplace_back(i, j, this_DeltaR, this_pt_diff);
            }
        }
    }
    //sort them by dR and then pt_diff
    sort(possible_matches.begin(), possible_matches.end(), [](const tuple<int,int,float,float> &a, const tuple<int,int,float,float> &b){
        if(get<2>(a) == get<2>(b)) return get<3>(a) < get<3>(b);
        return get<2>(a) < get<2>(b);
    });
    //greeedy matching
    RVec<bool> used_jet(jets.size(), false);
    RVec<bool> used_genjet(genjets.size(), false);
    for(const auto &match: possible_matches){
        int jet_idx = get<0>(match);
        int genjet_idx = get<1>(match);
        if(used_jet[jet_idx] || used_genjet[genjet_idx]) continue;
        matched_genjet_idx[jet_idx] = genjet_idx;
        used_jet[jet_idx] = true;
        used_genjet[genjet_idx] = true;
    }
    // assign -999 to unmatched jet
    for(size_t i = 0; i < jets.size(); i++){
        if(used_jet[i]) continue;
        else matched_genjet_idx[i] = -999;
    }
    return matched_genjet_idx;
}

unordered_map<int, int> AnalyzerCore::deltaRMatching(const RVec<Particle> &objs1, const RVec<Particle> &objs2, const float dR) {
    RVec<tuple<size_t, size_t, float>> possible_matches;
    RVec<bool> used_obj1(objs1.size(), false);
    RVec<bool> used_obj2(objs2.size(), false);
    for(size_t i = 0; i < objs1.size(); i++){
        for(size_t j = 0; j < objs2.size(); j++){
            float this_dR = objs1[i].DeltaR(objs2[j]);
            if(this_dR < dR){
                possible_matches.emplace_back(i, j, this_dR);
            }
        }
    }
    sort(possible_matches.begin(), possible_matches.end(), [](const tuple<size_t, size_t, float> &a, const tuple<size_t, size_t, float> &b){
        return get<2>(a) < get<2>(b);
    });
    unordered_map<int, int> matched_idx;
    //print all
    for(size_t i = 0; i < possible_matches.size(); i++){
        int obj1_idx = get<0>(possible_matches[i]);
        int obj2_idx = get<1>(possible_matches[i]);
        if(used_obj1[obj1_idx] || used_obj2[obj2_idx]) continue;
        matched_idx[obj1_idx] = obj2_idx;
        used_obj1[obj1_idx] = true;
        used_obj2[obj2_idx] = true;
    }

    for(size_t i = 0; i < objs1.size(); i++){
        if(used_obj1[i]) continue;
        matched_idx[i] = -999;
    }

    return matched_idx;
}

RVec<Jet> AnalyzerCore::SmearJets(const RVec<Jet> &jets, const RVec<GenJet> &genjets, const MyCorrection::variation &syst, const TString &source) {
    unordered_map<int, int> matched_idx = GenJetMatching(jets, genjets, fixedGridRhoFastjetAll);
    RVec<Jet> smeared_jets;
    gRandom->SetSeed(int(MET_pt*1e6));
    const float MIN_JET_ENERGY=1e-2;
    for(size_t i = 0; i < jets.size(); i++){
        Jet this_jet = jets.at(i);

        // backward smearing for systematic variation
        if (this_jet.GetUnsmearedP4().E() < 0) {
            // never smeared yet
            this_jet.SetUnsmearedP4(this_jet);
        } else {
            // already smeared
            const float backward_factor = this_jet.GetUnsmearedP4().Pt() / this_jet.Pt();
            this_jet *= backward_factor;
        }

        // Following the procedures in https://github.com/choij1589/SKFlatMaker/blob/master/SKFlatMaker/src/SKFlatMaker.cc#L3378-L3419
        float this_corr = 1.;
        const float this_jer = myCorr->GetJER(this_jet.Eta(), this_jet.Pt(), fixedGridRhoFastjetAll);
        const float this_sf = myCorr->GetJERSF(this_jet.Eta(), this_jet.Pt(), syst, source);
        if (matched_idx[i] > 0 && matched_idx[i] < genjets.size()) {
            // found matched jet
            const float matched_genjet_pt = genjets[matched_idx[i]].Pt();
            this_corr += (this_sf-1.) * (1.-matched_genjet_pt/this_jet.Pt());
        } else {
            this_corr += (gRandom->Gaus(0., this_jer))*sqrt(max(this_sf*this_sf-1., 0.));
        }
        // To avoid flipping direction (this_corr < 0)
        const float min_corr = MIN_JET_ENERGY/this_jet.E();
        this_corr = max(this_corr, min_corr);
        
        /*
        if(matched_idx[i] < 0) {
            // if the jet is not matched to any genjet, do stochastic smearing
            if(this_sf < 1.) {
                this_corr = MIN_JET_ENERGY / this_jet.E();
            } else{
                this_corr += (gRandom->Gaus(0, this_jer)) * sqrt(max(this_sf * this_sf - 1., 0.0));
                float new_corr = MIN_JET_ENERGY / this_jet.E();
                this_corr = max(this_corr, new_corr);
            }
        } else {
            float this_genjet_pt = genjets[matched_idx[i]].Pt();
            this_corr += (this_sf - 1.) * (1. - this_genjet_pt / this_jet.Pt());
            float new_corr = MIN_JET_ENERGY / this_jet.E();
            this_corr = max(this_corr, new_corr);
        }
        */
        this_jet *= this_corr;
        smeared_jets.push_back(this_jet);
    }
    return smeared_jets;
}

RVec<Jet> AnalyzerCore::SmearJets(const RVec<Jet> &jets, const RVec<GenJet> &genjets, const TString &syst, const TString &source){
    if (syst == "nominal") return jets;
    else if (syst == "up") return SmearJets(jets, genjets, MyCorrection::variation::up, source);
    else if (syst == "down") return SmearJets(jets, genjets, MyCorrection::variation::down, source);
    else throw runtime_error("[AnalyzerCore::SmearJets] Invalid syst value");
}

RVec<Jet> AnalyzerCore::ScaleJets(const RVec<Jet> &jets, const MyCorrection::variation &syst, const TString &source) {
    if(syst == MyCorrection::variation::nom) return jets;
    
    RVec<Jet> scaled_jets;
    RVec<TString> syst_sources = {"AbsoluteMPFBias",
                                  "AbsoluteScale",
                                  "AbsoluteStat",
                                  "FlavorQCD",
                                  "Fragmentation",
                                  //"PileUpDataMC",
                                  //"PileUpPtBB",
                                  //"PileUpPtEC1",
                                  //"PileUpPtEC2",
                                  //"PileUpPtHF",
                                  //"PileUpPtRef",
                                  "PileUpEnvelope",
                                  "RelativeJEREC1",
                                  "RelativeJEREC2",
                                  "RelativeJERHF",
                                  "RelativePtBB",
                                  "RelativePtEC1",
                                  "RelativePtEC2",
                                  "RelativePtHF",
                                  "RelativeBal",
                                  "RelativeSample",
                                  "RelativeFSR",
                                  "RelativeStatFSR",
                                  "RelativeStatEC",
                                  "RelativeStatHF",
                                  "SinglePionECAL",
                                  "SinglePionHCAL",
                                  "TimePtEta"};

    for(const auto &jet: jets){
        Jet this_jet = jet;

        TLorentzVector this_unsmearedP4 = jet.GetUnsmearedP4();
        if(this_unsmearedP4.E() < 0){
            cerr << "[AnalyzerCore::ScaleJets] It seems the attribute UnsmearedP4 has never been set for this Jet." << endl;
            cerr << "[AnalyzerCore::ScaleJets] Since all jets in MC are smeared when calling GetAllJets(), UnsmearedP4 should have been set by AnalyzerCore::SmearJets at the same time." << endl;
            cerr << "[AnalyzerCore::ScaleJets] This may indicate spurious behavior in your code. For now, Jet.pt() is used for scaling, but I recommend inspecting your code." << endl;

            this_unsmearedP4 = jet;
        }

        if(source == "total"){
            for(const auto &it: syst_sources) {
                this_jet *= myCorr->GetJESUncertainty(this_unsmearedP4.Eta(), this_unsmearedP4.Pt(), syst, it);
            }
            scaled_jets.push_back(this_jet);
        } else{
            this_jet *= myCorr->GetJESUncertainty(this_unsmearedP4.Eta(), this_unsmearedP4.Pt(), syst, source);
            scaled_jets.push_back(this_jet);
        }
    }
    return scaled_jets;
}

RVec<Jet> AnalyzerCore::ScaleJets(const RVec<Jet> &jets, const TString &syst, const TString &source){
    if(syst == "nom") return jets;
    else if (syst == "up") return ScaleJets(jets, MyCorrection::variation::up, source);
    else if (syst == "down") return ScaleJets(jets, MyCorrection::variation::down, source);
    else throw runtime_error("[AnalyzerCore::ScaleJets] Invalid syst value");
}

Particle AnalyzerCore::ApplyTypeICorrection(const Particle& MET,
                                            const RVec<Jet>& jets,
                                            const RVec<Electron>& electrons,
                                            const RVec<Muon>& muons,
                                            const MyCorrection::variation& unclustered_syst) {
    // Unclustered energy systematic is calculated separately by variation.
    // We keep this functionaltiy for the CHS MET usage, as there is no branch for UE_UP and UE_DOWN.
    // If you are using PUPPI MET, you can use the branches PuppiMET_ptUnclusteredUp and PuppiMET_ptUnclusteredDown
    // So only use up/down variations for unclustered energy if you are using CHS MET.
    // For the others (like JES, JER, etc.), we assume they are already applied to the input


    // Apply Type-I correction using TLorentzVector
    TLorentzVector metVector;
    metVector.SetPxPyPzE(MET.Px(), MET.Py(), 0.0, MET.Pt());

    // Sum of old corrected jets passing acceptance (pT > 15 GeV) - these were used in NanoAOD Type-I
    TLorentzVector sumOldCorrJets;
    for (const auto& jet : jets) {
        if (jet.GetOriginalPt() > 15.0 && fabs(jet.Eta()) < 5.0) {
            TLorentzVector oldCorrJet;
            // Original mass: Jet_mass[i] from NanoAOD (before our re-correction)
            float originalMass = jet.M();
            if (jet.Pt() > 0.0f) {
                // Reconstruct original mass: since both pt and mass scale by JESSF,
                // originalMass = currentMass * (originalPt / currentPt)
                originalMass = jet.M() * (jet.GetOriginalPt() / jet.Pt());
            }
            oldCorrJet.SetPtEtaPhiM(jet.GetOriginalPt(), jet.Eta(), jet.Phi(), originalMass);
            sumOldCorrJets += oldCorrJet;
        }
    }

    // Sum of corrected jets passing acceptance (pT > 15 GeV, |eta| < 5.0)
    TLorentzVector sumCorrJets;
    for (const auto& jet : jets) {
        if (jet.Pt() > 15.0 && fabs(jet.Eta()) < 5.0) {
            sumCorrJets += jet;
        }
    }
    metVector -= (sumCorrJets - sumOldCorrJets);

    // Sum of raw muons passing acceptance (pT > 5 GeV, |eta| < 2.4)
    TLorentzVector sumRawMuons;
    for (const auto& muon : muons) {
        if (muon.OriginalPt() > 5.0 && fabs(muon.Eta()) < 2.4) {
            TLorentzVector rawMuon;
            rawMuon.SetPtEtaPhiM(muon.OriginalPt(), muon.Eta(), muon.Phi(), muon.M());
            sumRawMuons += rawMuon;
        }
    }

    // Sum of corrected muons passing acceptance (pT > 5 GeV, |eta| < 2.4)
    TLorentzVector sumCorrMuons;
    for (const auto& muon : muons) {
        if (muon.Pt() > 5.0 && fabs(muon.Eta()) < 2.4) {
            sumCorrMuons += muon;
        }
    }
    metVector -= (sumCorrMuons - sumRawMuons);

    // Electrons: correlate MET with electron scale/resolution variations
    // Use nominal electrons from the event (NanoAOD nominal == raw reference)
    RVec<Electron> nominal_electrons = GetAllElectrons();
    TLorentzVector sumNominalElectrons;
    for (const auto& ele : nominal_electrons) {
        if (ele.Pt() > 7.0 && fabs(ele.Eta()) < 2.5) {
            sumNominalElectrons += ele;
        }
    }
    TLorentzVector sumCorrElectrons;
    for (const auto& ele : electrons) {
        if (ele.Pt() > 7.0 && fabs(ele.Eta()) < 2.5) {
            sumCorrElectrons += ele;
        }
    }
    metVector -= (sumCorrElectrons - sumNominalElectrons);

    // Unclustered energy systematic (±10%) applied on top of correlated object variations
    if (unclustered_syst != MyCorrection::variation::nom) {
        TLorentzVector totalObjects = sumCorrJets + sumCorrMuons + sumCorrElectrons;
        TLorentzVector unclusteredEnergy = metVector + totalObjects;
        float unc_factor = (unclustered_syst == MyCorrection::variation::up) ? 0.1f : -0.1f;
        TLorentzVector unclusteredVariation = unclusteredEnergy * unc_factor;
        metVector -= unclusteredVariation;
    }

    Particle correctedMET;
    correctedMET.SetPtEtaPhiM(metVector.Pt(), 0.0, metVector.Phi(), 0.0);
    return correctedMET;
}

// pdfs
float AnalyzerCore::GetPDFWeight(LHAPDF::PDF *pdf) {
    float pdf1 = pdf->xfxQ2(Generator_id1, Generator_x1, Generator_scalePDF);
    float pdf2 = pdf->xfxQ2(Generator_id2, Generator_x2, Generator_scalePDF);
    return pdf1 * pdf2;
}

float AnalyzerCore::GetPDFReweight() {
    return GetPDFWeight(pdfReweight->NewPDF) / GetPDFWeight(pdfReweight->ProdPDF);
}

float AnalyzerCore::GetPDFReweight(int member) {
    return GetPDFWeight(pdfReweight->PDFErrorSet.at(member)) / GetPDFWeight(pdfReweight->ProdPDF);
}

// MC weights
// TODO: Treat MiNNLO and Sherpa weights 
float AnalyzerCore::MCweight(bool usesign, bool norm_1invpb) const {
    if (IsDATA) return 1.;
    float weight = genWeight;
      if (usesign){
        if (weight > 0)
        weight = 1.0;
        else if (weight < 0)
        weight = -1.0;
        else
        weight = 0.0;
    }

    if (norm_1invpb)
    {
        if (usesign)
        weight *= xsec / sumSign;
        else
        weight *= xsec / sumW;
    }
    return weight;

}

// Objects
Event AnalyzerCore::GetEvent() {
    Event ev;
    ev.SetRunLumiEvent(RunNumber, LumiBlock, EventNumber);
    ev.SetnPileUp(Pileup_nPU);
    ev.SetnTrueInt(Pileup_nTrueInt);
    if (Run == 2) ev.SetnPVsGood(static_cast<int>(PV_npvsGood_RunII));
    if (Run == 3) ev.SetnPVsGood(PV_npvsGood);
    ev.SetGenMET(GenMET_pt, GenMET_phi);
    ev.SetMETVector(MET_pt, MET_phi, Event::MET_Type::CHS);
    ev.SetMETVector(PuppiMET_pt, PuppiMET_phi, Event::MET_Type::PUPPI);
    ev.SetMETVector(PuppiMET_ptUnclusteredUp, PuppiMET_phiUnclusteredUp, Event::MET_Type::PUPPI, Event::MET_Syst::UE_UP);
    ev.SetMETVector(PuppiMET_ptUnclusteredDown, PuppiMET_phiUnclusteredDown, Event::MET_Type::PUPPI, Event::MET_Syst::UE_DOWN);
    ev.SetTrigger(TriggerMap);
    ev.SetEra(GetEra());
    ev.setRho(fixedGridRhoFastjetAll);
    return ev;
}

RVec<Muon> AnalyzerCore::GetAllMuons() {
    RVec<Muon> muons;
    RVec<Gen> truth;
    if (!IsDATA) truth = GetAllGens();

    for (int i = 0; i < nMuon; i++) {
        Muon muon;
        muon.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_mass[i]);
        muon.SetCharge(Muon_charge[i]);
        muon.SetNTrackerLayers(Muon_nTrackerLayers[i]);
        float roccor = 1.;
        float roccor_err = 0.;
        if (IsDATA) {
            roccor = myCorr->GetMuonScaleSF(muon, MyCorrection::variation::nom);
            roccor_err = myCorr->GetMuonScaleSF(muon, MyCorrection::variation::up) - roccor;
        } else {
            Gen matched_gen = GetGenMatchedMuon(muon, truth);
            float matched_pt = matched_gen.Pt();
            roccor = myCorr->GetMuonScaleSF(muon, MyCorrection::variation::nom, matched_pt);
            roccor_err = myCorr->GetMuonScaleSF(muon, MyCorrection::variation::up, matched_pt) - roccor;
        }
        muon.SetOriginalPt(muon.Pt());
        muon.SetMomentumScaleUpDown(muon.Pt()*(roccor+roccor_err), muon.Pt()*(roccor-roccor_err)); 
        muon.SetPtEtaPhiM(muon.Pt()*roccor, muon.Eta(), muon.Phi(), muon.M());
        muon.SetTkRelIso(Muon_tkRelIso[i]);
        muon.SetPfRelIso03(Muon_pfRelIso03_all[i]);
        muon.SetPfRelIso04(Muon_pfRelIso04_all[i]);
        muon.SetMiniPFRelIso(Muon_miniPFRelIso_all[i]);
        muon.SetdXY(Muon_dxy[i], Muon_dxyErr[i]);
        muon.SetdZ(Muon_dz[i], Muon_dzErr[i]);
        muon.SetIP3D(Muon_ip3d[i], Muon_sip3d[i]);
        muon.SetNTrackerLayers(Muon_nTrackerLayers[i]);
        muon.SetGenPartFlav(Muon_genPartFlav[i]);
        muon.SetBIDBit(Muon::BooleanID::LOOSE, Muon_looseId[i]);
        muon.SetBIDBit(Muon::BooleanID::MEDIUM, Muon_mediumId[i]);
        muon.SetBIDBit(Muon::BooleanID::MEDIUMPROMPT, Muon_mediumPromptId[i]);
        muon.SetBIDBit(Muon::BooleanID::TIGHT, Muon_tightId[i]);
        muon.SetBIDBit(Muon::BooleanID::SOFT, Muon_softId[i]);
        muon.SetBIDBit(Muon::BooleanID::SOFTMVA, Muon_softMvaId[i]);
        muon.SetBIDBit(Muon::BooleanID::TRIGGERLOOSE, Muon_triggerIdLoose[i]);
        muon.SetWIDBit(Muon::WorkingPointID::HIGHPT, Muon_highPtId[i]);
        muon.SetWIDBit(Muon::WorkingPointID::MINIISO, Muon_miniIsoId[i]);
        muon.SetWIDBit(Muon::WorkingPointID::MULTIISO, Muon_multiIsoId[i]);
        if(Run == 3){
            muon.SetWIDBit(Muon::WorkingPointID::MVAMU, Muon_mvaMuID_WP[i]);
            muon.SetGenPartIdx(Muon_genPartIdx[i]);
            muon.SetJetIdx(Muon_jetIdx[i]);
        }
        else if(Run == 2){
            muon.SetWIDBit(Muon::WorkingPointID::MVAMU, Muon_mvaId[i]);
            muon.SetGenPartIdx(static_cast<short>(Muon_genPartIdx_RunII[i]));
            muon.SetJetIdx(Muon_jetIdx_RunII[i]);
        }
        //muon.SetWIDBit(Muon::WorkingPointID::MVALOWPT, Muon_mvaLowPtId[i]);
        muon.SetWIDBit(Muon::WorkingPointID::PFISO, Muon_pfIsoId[i]);
        muon.SetWIDBit(Muon::WorkingPointID::PUPPIISO, Muon_puppiIsoId[i]);
        muon.SetWIDBit(Muon::WorkingPointID::TKISO, Muon_tkIsoId[i]);
        muon.SetMVAID(Muon::MVAID::SOFTMVA, Muon_softMva[i]);
        muon.SetMVAID(Muon::MVAID::MVALOWPT, Muon_mvaLowPt[i]);
        muon.SetMVAID(Muon::MVAID::MVATTH, Muon_mvaTTH[i]);
        muon.SetIsTracker(Muon_isTracker[i]);
        muon.SetIsGlobal(Muon_isGlobal[i]);
        muon.SetIsStandalone(Muon_isStandalone[i]);

        muons.push_back(muon);
    }

    return muons;
}

RVec<Muon> AnalyzerCore::ScaleMuons(const RVec<Muon> &muons, const TString &syst) {
    RVec<Muon> scaled_muons;
    for (const auto &muon: muons) {
        Muon scaled_muon = muon;
        if (syst == "up") {
            scaled_muon.SetPtEtaPhiM(muon.MomentumScaleUp(), muon.Eta(), muon.Phi(), muon.M());
        } else if (syst == "down") {
            scaled_muon.SetPtEtaPhiM(muon.MomentumScaleDown(), muon.Eta(), muon.Phi(), muon.M());
        } else {
            throw runtime_error("[AnalyzerCore::ScaleMuons] Invalid variation");
        }
        scaled_muons.emplace_back(scaled_muon);
    }
    return scaled_muons;
}

RVec<Muon> AnalyzerCore::GetMuons(const TString ID, const float ptmin, const float fetamax) {
    RVec<Muon> muons = GetAllMuons();
    RVec<Muon> selected_muons;
    for (const auto &muon: muons) {
        if (! (muon.Pt() > ptmin)) continue;
        if (! (fabs(muon.Eta()) < fetamax)) continue;
        if (! muon.PassID(ID)) continue;
        selected_muons.push_back(muon);
    }
    return selected_muons;
}

RVec<Muon> AnalyzerCore::SelectMuons(const RVec<Muon> &muons, const TString ID, const float ptmin, const float fetamax) const {
    RVec<Muon> selected_muons;

    for (const auto &muon: muons) {
        if (! (muon.Pt() > ptmin)) continue;
        if (! (fabs(muon.Eta()) < fetamax)) continue;
        if (ID.Contains("&&"))
        {
            TObjArray *IDs = ID.Tokenize("&&");
            bool pass = true;
            for (int i = 0; i < IDs->GetEntries(); i++)
            {
                TString thisID = ((TObjString *)IDs->At(i))->GetString();
                if (!muon.PassID(thisID))
                {
                    pass = false;
                    break;
                }
            }
            delete IDs;
            if (!pass)
                continue;
        }
        else
        {
            if (!muon.PassID(ID))
                continue;
        }
        selected_muons.push_back(muon);
    }
    return selected_muons;
}

RVec<Muon> AnalyzerCore::SelectMuons(const RVec<Muon> &muons, const Muon::MuonID ID, const float ptmin, const float fetamax) const
{
    RVec<Muon> selected_muons;

    for (const auto &muon : muons)
    {
        if (!(muon.Pt() > ptmin))
            continue;
        if (!(fabs(muon.Eta()) < fetamax))
            continue;
        if (!muon.PassID(ID))
            continue;
        selected_muons.push_back(muon);
    }
    return selected_muons;
}

RVec<Electron> AnalyzerCore::GetAllElectrons(){
    RVec<Electron> electrons;
    for (int i = 0; i < nElectron; i++){
        // Reject GAP region electrons
        const float fscEta = fabs(Electron_scEta[i]);
        if (1.444 < fscEta && fscEta < 1.566) continue;

        Electron electron;
        electron.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_mass[i]);
        electron.SetCharge(Electron_charge[i]);
        electron.SetScEta(Electron_scEta[i]);
        electron.SetDeltaEtaInSC(Electron_deltaEtaInSC[i]);
        electron.SetDeltaEtaInSeed(Electron_deltaEtaInSeed[i]);
        electron.SetDeltaPhiInSC(Electron_deltaPhiInSC[i]);
        electron.SetDeltaPhiInSeed(Electron_deltaPhiInSeed[i]);
        electron.SetPFClusterIso(Electron_ecalPFClusterIso[i], Electron_hcalPFClusterIso[i]);
        electron.SetPfRelIso03(Electron_pfRelIso03_all[i]);
        electron.SetMiniPFRelIso(Electron_miniPFRelIso_all[i]);
        electron.SetdXY(Electron_dxy[i], Electron_dxyErr[i]);
        electron.SetdZ(Electron_dz[i], Electron_dzErr[i]);
        electron.SetIP3D(Electron_ip3d[i], Electron_sip3d[i]);
        electron.SetConvVeto(Electron_convVeto[i]);
        electron.SetLostHits(Electron_lostHits[i]);
        electron.SetSeedGain(Electron_seedGain[i]);
        electron.SetTightCharge(Electron_tightCharge[i]);
        electron.SetSieie(Electron_sieie[i]);
        electron.SetHoe(Electron_hoe[i]);
        electron.SetEInvMinusPInv(Electron_eInvMinusPInv[i]);
        electron.SetDr03EcalRecHitSumEt(Electron_dr03EcalRecHitSumEt[i]);
        electron.SetDr03HcalDepth1TowerSumEt(Electron_dr03HcalDepth1TowerSumEt[i]);
        electron.SetDr03TkSumPt(Electron_dr03TkSumPt[i]);
        electron.SetDr03TkSumPtHEEP(Electron_dr03TkSumPtHEEP[i]);
        electron.SetR9(Electron_r9[i]);
        electron.SetRho(fixedGridRhoFastjetAll);
        electron.SetEnergyErr(Electron_energyErr[i]);
        electron.SetBIDBit(Electron::BooleanID::CUTBASEDHEEP, Electron_cutBased_HEEP[i]);
        electron.SetMVA(Electron::MVATYPE::MVATTH, Electron_mvaTTH[i]);
        electron.SetGenPartFlav(Electron_genPartFlav[i]);

        if (Run == 2) {
            electron.SetEnergyResUnc(Electron_dEsigmaUp[i], Electron_dEsigmaDown[i]);
            electron.SetBIDBit(Electron::BooleanID::MVAISOWP80, Electron_mvaFall17V2Iso_WP80[i]);
            electron.SetBIDBit(Electron::BooleanID::MVAISOWP90, Electron_mvaFall17V2Iso_WP90[i]);
            electron.SetBIDBit(Electron::BooleanID::MVAISOWPL, Electron_mvaFall17V2Iso_WPL[i]);
            electron.SetBIDBit(Electron::BooleanID::MVANOISOWP80, Electron_mvaFall17V2noIso_WP80[i]);
            electron.SetBIDBit(Electron::BooleanID::MVANOISOWP90, Electron_mvaFall17V2noIso_WP90[i]);
            electron.SetBIDBit(Electron::BooleanID::MVANOISOWPL, Electron_mvaFall17V2noIso_WPL[i]);
            electron.SetMVA(Electron::MVATYPE::MVAISO, Electron_mvaFall17V2Iso[i]);
            electron.SetMVA(Electron::MVATYPE::MVANOISO, Electron_mvaFall17V2noIso[i]);
            electron.SetCBIDBit(Electron::CutBasedID::CUTBASED, Electron_cutBased_RunII[i]);
            electron.SetGenPartIdx(static_cast<short>(Electron_genPartIdx_RunII[i]));
            electron.SetJetIdx(Electron_jetIdx_RunII[i]);
        } else if (Run == 3) {
            electron.SetBIDBit(Electron::BooleanID::MVAISOWP80, Electron_mvaIso_WP80[i]);
            electron.SetBIDBit(Electron::BooleanID::MVAISOWP90, Electron_mvaIso_WP90[i]);
            electron.SetBIDBit(Electron::BooleanID::MVAISOWPL, Electron_mvaIso_WPL[i]);
            electron.SetBIDBit(Electron::BooleanID::MVANOISOWP80, Electron_mvaNoIso_WP80[i]);
            electron.SetBIDBit(Electron::BooleanID::MVANOISOWP90, Electron_mvaNoIso_WP90[i]);
            electron.SetBIDBit(Electron::BooleanID::MVANOISOWPL, Electron_mvaNoIso_WPL[i]);
            electron.SetMVA(Electron::MVATYPE::MVAISO, Electron_mvaIso[i]);
            electron.SetMVA(Electron::MVATYPE::MVANOISO, Electron_mvaNoIso[i]);
            electron.SetCBIDBit(Electron::CutBasedID::CUTBASED, Electron_cutBased[i]);
            electron.SetGenPartIdx(Electron_genPartIdx[i]);
            electron.SetJetIdx(Electron_jetIdx[i]); 
        } else {
            throw runtime_error("[AnalyzerCore::GetAllElectrons] Invalid run number");
        }


        electrons.push_back(electron);
    }

    return electrons;
}

RVec<Electron> AnalyzerCore::GetElectrons(const TString ID, const float ptmin, const float fetamax, bool vetoHEM) {
    RVec<Electron> electrons = GetAllElectrons();
    RVec<Electron> selected_electrons;
    for (const auto &electron: electrons) {
        if (! (electron.Pt() > ptmin)) continue;
        if (! (fabs(electron.Eta()) < fetamax)) continue;
        if (! electron.PassID(ID)) continue;
        if (vetoHEM && IsHEMElectron(electron)) continue;
        selected_electrons.push_back(electron);
    }
    return selected_electrons;
}

RVec<Electron> AnalyzerCore::SelectElectrons(const RVec<Electron> &electrons, const TString ID, const float ptmin, const float fetamax, bool vetoHEM) const {
    RVec<Electron> selected_electrons;
    for (const auto &electron: electrons) {
        if (! (electron.Pt() > ptmin)) continue;
        if (! (fabs(electron.Eta()) < fetamax)) continue;
        if (! electron.PassID(ID)) continue;
        if (vetoHEM && IsHEMElectron(electron)) continue;
        selected_electrons.push_back(electron);
    }
    return selected_electrons;
}

RVec<Electron> AnalyzerCore::SelectElectrons(const RVec<Electron> &electrons, const Electron::ElectronID ID, const float ptmin, const float fetamax, bool vetoHEM) const {
    RVec<Electron> selected_electrons;
    for (const auto &electron : electrons) {
        if (!(electron.Pt() > ptmin)) continue;
        if (!(fabs(electron.Eta()) < fetamax)) continue;
        if (!electron.PassID(ID)) continue;
        if (vetoHEM && IsHEMElectron(electron)) continue;
        selected_electrons.push_back(electron);
    }
    return selected_electrons;
}

RVec<Electron> AnalyzerCore::ScaleElectrons(const Event &ev, const RVec<Electron> &electrons, const TString &syst) {
    if (IsDATA || syst == "nom") return electrons;
    
    RVec<Electron> scaled_electrons;
    for (const auto &electron: electrons) {
        float scale_variation = 1.;
        if (syst == "up") {
            scale_variation = myCorr->GetElectronScaleUnc(electron.scEta(), electron.SeedGain(), ev.run(), electron.r9(), electron.Pt(), MyCorrection::variation::up);    
        } else if (syst == "down") {
            scale_variation = myCorr->GetElectronScaleUnc(electron.scEta(), electron.SeedGain(), ev.run(), electron.r9(), electron.Pt(), MyCorrection::variation::down);    
        } else {
            throw runtime_error("[AnalyzerCore::ScaleElectrons] Invalid variation");
        }
        Electron scaled_electron = electron;
        scaled_electron.SetPtEtaPhiM(electron.Pt()*scale_variation, electron.Eta(), electron.Phi(), electron.M());
        scaled_electrons.emplace_back(scaled_electron);
    }
    return scaled_electrons;
}

RVec<Electron> AnalyzerCore::SmearElectrons(const RVec<Electron> &electrons, const TString &syst) {
    RVec<Electron> smeared_electrons;
    switch(Run) {
        case 2:
            for (const auto &electron: electrons) {
                float smeared_pt = electron.Pt();
                if (syst == "up") {
                    smeared_pt *= (electron.E() - electron.dEsigmaUp())/electron.E();
                } else if (syst == "down") {
                    smeared_pt *= (electron.E() - electron.dEsigmaDown())/electron.E();
                } else {
                    throw runtime_error("[AnalyzerCore::SmearElectrons] Invalid variation");
                }
                Electron smeared_electron = electron;
                smeared_electron.SetPtEtaPhiM(smeared_pt, electron.Eta(), electron.Phi(), electron.M());
                smeared_electrons.emplace_back(smeared_electron);
            }
            break;
        case 3: 
            for (const auto &electron: electrons) {
                float smeared_pt = electron.Pt();
                if (syst == "nom") {
                    smeared_pt *= myCorr->GetElectronSmearUnc(electron, MyCorrection::variation::nom, int(electron.Rho()));
                } else if (syst == "up") {
                    smeared_pt *= myCorr->GetElectronSmearUnc(electron, MyCorrection::variation::up, int(electron.Rho()));
                } else if (syst == "down") {
                    smeared_pt *= myCorr->GetElectronSmearUnc(electron, MyCorrection::variation::down, int(electron.Rho()));
                } else {
                    throw runtime_error("[AnalyzerCore::SmearElectrons] Invalid variation");
                }
                Electron smeared_electron = electron;
                smeared_electron.SetPtEtaPhiM(smeared_pt, electron.Eta(), electron.Phi(), electron.M());
                smeared_electrons.emplace_back(smeared_electron);
            }
            break;
        default:
            throw runtime_error("[AnalyzerCore::SmearElectrons] Run " + to_string(Run) + " is not supported");
    }
    return smeared_electrons;
}

RVec<Gen> AnalyzerCore::GetAllGens(){
    RVec<Gen> Gens;
    if(IsDATA) return Gens;

    for (int i = 0; i < nGenPart; i++){

        Gen gen;

        gen.SetIsEmpty(false);
        gen.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
        gen.SetIndexPIDStatus(i, GenPart_pdgId[i], GenPart_status[i]);
        if(Run == 3){
            gen.SetMother(GenPart_genPartIdxMother[i]);
            gen.SetGenStatusFlags(GenPart_statusFlags[i]);
        }
        else if(Run == 2){
            gen.SetMother(static_cast<short>(GenPart_genPartIdxMother_RunII[i]));
            gen.SetGenStatusFlags(static_cast<unsigned short>(GenPart_statusFlags_RunII[i])); 
        }

        Gens.push_back(gen);
    }

    return Gens;
}

RVec<LHE> AnalyzerCore::GetAllLHEs() {
    RVec<LHE> lhes;
    if (IsDATA) return lhes;
    
    // Safety check for array sizes
    if (nLHEPart <= 0) return lhes;

    for (int i = 0; i < nLHEPart; i++) {
        LHE lhe;
        lhe.SetIsEmpty(false);
        lhe.SetPtEtaPhiM(LHEPart_pt[i], LHEPart_eta[i], LHEPart_phi[i], LHEPart_mass[i]);
        lhe.SetStatus(LHEPart_status[i]);
        lhe.SetSpin(LHEPart_spin[i]);
        lhe.SetIncomingPz(LHEPart_incomingpz[i]);
        lhe.SetPdgId(LHEPart_pdgId[i]);
        lhes.push_back(lhe);
    }
    return lhes;
}

RVec<Tau> AnalyzerCore::GetAllTaus(){

    RVec<Tau> taus;
    
    for (int i = 0; i < nTau; i++) {

        Tau tau;
        tau.SetPtEtaPhiM(Tau_pt[i], Tau_eta[i], Tau_phi[i], Tau_mass[i]);
        tau.SetdXY(Tau_dxy[i]);
        tau.SetdZ(Tau_dz[i]);
        tau.SetGenPartFlav(Tau_genPartFlav[i]);
        if(Run == 3){
            tau.SetCharge(Tau_charge[i]);
            tau.SetDecayMode(Tau_decayMode[i]);
            tau.SetIdDeepTau2018v2p5VSjet(Tau_idDeepTau2018v2p5VSjet[i]);
            tau.SetIdDeepTau2018v2p5VSmu(Tau_idDeepTau2018v2p5VSmu[i]);
            tau.SetIdDeepTau2018v2p5VSe(Tau_idDeepTau2018v2p5VSe[i]);
            if(!IsDATA)tau.SetGenPartIdx(Tau_genPartIdx[i]);
        }
        else if(Run == 2){
            tau.SetCharge(Tau_charge_RunII[i]);
            tau.SetDecayMode(Tau_decayMode_RunII[i]);
            tau.SetIdDeepTau2018v2p5VSjet(Tau_idDeepTau2017v2p1VSjet[i]);
            tau.SetIdDeepTau2018v2p5VSmu(Tau_idDeepTau2017v2p1VSmu[i]);
            tau.SetIdDeepTau2018v2p5VSe(Tau_idDeepTau2017v2p1VSe[i]);
            if(!IsDATA)tau.SetGenPartIdx(Tau_genPartIdx_RunII[i]);
        }
        tau.SetIdDecayModeNewDMs(Tau_idDecayModeNewDMs[i]);
        taus.push_back(tau);

    }

    return taus;

}

RVec<Tau> AnalyzerCore::SelectTaus(const RVec<Tau> &taus, const TString ID, const float ptmin, const float absetamax) const{
    
    RVec<Tau> selected_taus;
    
    for (const auto &tau: taus) {

        if (! (tau.Pt() > ptmin)) continue;
        if (! (fabs(tau.Eta()) < absetamax)) continue;
        if (! tau.PassID(ID)) continue;
        selected_taus.push_back(tau);

    }

    return selected_taus;

}

RVec<Jet> AnalyzerCore::GetAllJets() {
    RVec<Jet> Jets;
    
    for (int i = 0; i < nJet; i++) {
        Jet jet;
        const float rawPt = Jet_pt[i] * (1.-Jet_rawFactor[i]);
        const float rawMass = Jet_mass[i] * (1.-Jet_rawFactor[i]);
        const float JESSF = myCorr->GetJESSF(Jet_area[i], Jet_eta[i], rawPt, Jet_phi[i], fixedGridRhoFastjetAll, RunNumber);
        const float correctedPt = rawPt * JESSF;
        const float correctedMass = rawMass * JESSF;
        jet.SetPtEtaPhiM(correctedPt, Jet_eta[i], Jet_phi[i], correctedMass);
        jet.SetRawPt(rawPt);
        jet.SetOriginalPt(Jet_pt[i]);
        jet.SetArea(Jet_area[i]);
        jet.SetOriginalIndex(i);
        jet.SetEnergyFractions(Jet_chHEF[i], Jet_neHEF[i], Jet_neEmEF[i], Jet_chEmEF[i], Jet_muEF[i]);
        if(!IsDATA){
            if(Run == 3)
            jet.SetJetFlavours(Jet_partonFlavour[i] ,Jet_hadronFlavour[i]);
            else if(Run == 2)
            jet.SetJetFlavours(static_cast<short>(Jet_partonFlavour_RunII[i]) , static_cast<unsigned char>(Jet_hadronFlavour_RunII[i]));
        }
        RVec<float> tvs;
        RVec<float> tvs2;
        if(Run == 3){
            tvs = {Jet_btagDeepFlavB[i], 
                   Jet_btagDeepFlavCvB[i], 
                   Jet_btagDeepFlavCvL[i], 
                   Jet_btagDeepFlavQG[i],
                   Jet_btagPNetB[i], 
                   Jet_btagPNetCvB[i], 
                   Jet_btagPNetCvL[i], 
                   Jet_btagPNetQvG[i],
                   Jet_btagPNetTauVJet[i], 
                   Jet_btagRobustParTAK4B[i], 
                   Jet_btagRobustParTAK4CvB[i], 
                   Jet_btagRobustParTAK4CvL[i], 
                   Jet_btagRobustParTAK4QG[i]};
            jet.SetMultiplicities(Jet_nConstituents[i], 
                                  Jet_nElectrons[i], 
                                  Jet_nMuons[i], 
                                  Jet_nSVs[i]);
            jet.SetHadronMultiplicities(Jet_chMultiplicity[i], 
                                        Jet_neMultiplicity[i]);
            if(!IsDATA){
                jet.SetMatchingIndices(Jet_electronIdx1[i], 
                                       Jet_electronIdx2[i], 
                                       Jet_muonIdx1[i], 
                                       Jet_muonIdx2[i], 
                                       Jet_svIdx1[i], 
                                       Jet_svIdx2[i], 
                                       Jet_genJetIdx[i]);
            }
            else{
                jet.SetMatchingIndices(Jet_electronIdx1[i], 
                                       Jet_electronIdx2[i], 
                                       Jet_muonIdx1[i], 
                                       Jet_muonIdx2[i], 
                                       Jet_svIdx1[i], 
                                       Jet_svIdx2[i]);
            }
            jet.SetJetID(Jet_jetId[i], 3);
            jet.SetJetPuID(0b111);
            tvs2 = {Jet_PNetRegPtRawCorr[i], Jet_PNetRegPtRawCorrNeutrino[i], Jet_PNetRegPtRawRes[i], Jet_rawFactor[i], -999.0, -999.0, -999.0, -999.0};
        }
        //for Run 2 DeepJet is only option
        else if(Run == 2){
            tvs = {Jet_btagDeepFlavB[i], 
                   Jet_btagDeepFlavCvB[i], 
                   Jet_btagDeepFlavCvL[i], 
                   Jet_btagDeepFlavQG[i],
                   -999., -999., -999., -999.,
                   -999., -999., -999., -999., -999.};
            jet.SetMultiplicities(Jet_nConstituents[i], 
                                  Jet_nElectrons_RunII[i], 
                                  Jet_nMuons_RunII[i],
                                  0);
            if(!IsDATA){
                jet.SetMatchingIndices(Jet_electronIdx1_RunII[i], 
                                       Jet_electronIdx2_RunII[i], 
                                       Jet_muonIdx1_RunII[i], 
                                       Jet_muonIdx2_RunII[i], 
                                       -9, 
                                       -9, 
                                       Jet_genJetIdx_RunII[i]);
            }
            else{
                jet.SetMatchingIndices(Jet_electronIdx1_RunII[i], 
                                       Jet_electronIdx2_RunII[i], 
                                       Jet_muonIdx1_RunII[i], 
                                       Jet_muonIdx2_RunII[i], 
                                       -9, 
                                       -9);
            }
            jet.SetJetID(Jet_jetId_RunII[i], 2);
            if (DataYear == 2016) {
                // due to the bug in the NanoAODv9, the puId is stored in a wrong way
                int InterChanged = 0;
                InterChanged = Jet_puId[i] >> 2 | ((Jet_puId[i] & 0b001) << 2) | (Jet_puId[i] & 0b010);
                jet.SetJetPuID(InterChanged);
            } else {
                jet.SetJetPuID(Jet_puId[i]);
            }
            tvs2 = {-999.0, -999.0, -999.0, -999.0, Jet_rawFactor[i], Jet_bRegCorr[i], Jet_bRegRes[i], Jet_cRegCorr[i], Jet_cRegRes[i]};
        }
        jet.SetTaggerResults(tvs);
        jet.SetCorrections(tvs2);
        Jets.push_back(jet);
    }
    if(!IsDATA) Jets = SmearJets(Jets, GetAllGenJets());
    return Jets;
}

RVec<Photon> AnalyzerCore::GetAllPhotons() {
    RVec<Photon> Photons;
    for (int i = 0; i< nPhoton; i++) {
        Photon photon;
        photon.SetPtEtaPhiE(Photon_pt[i], Photon_eta[i], Photon_phi[i], 1.0);
        float photon_theta = photon.Theta();
        float photon_E = Photon_pt[i] / sin( photon_theta );
        photon.SetPtEtaPhiE( Photon_pt[i], Photon_eta[i], Photon_phi[i], photon_E);
        photon.SetEnergy(photon_E);
        photon.SetSigmaIetaIeta(Photon_sieie[i]);
        photon.SetHoe(Photon_hoe[i]);
        if(Run == 3){
            photon.SetEnergyRaw(Photon_energyRaw[i]);
        }
        photon.SetPixelSeed(Photon_pixelSeed[i]);
        photon.SetisScEtaEB(Photon_isScEtaEB[i]);
        photon.SetisScEtaEE(Photon_isScEtaEE[i]);
//        Need to Add PV Info
//        Photon.SetSCEta(Photon_eta[i], Photon_phi[i], PV_x, PV_y, PV_z, Photon_isScEtaEB[i], Photon_isScEtaEE[i]);  
        photon.SetSCEta(Photon_eta[i], Photon_phi[i], 0, 0, 0, false, false); 
        photon.SetBIDBit(Photon::BooleanID::MVAIDWP80, Photon_mvaID_WP80[i]);
        photon.SetBIDBit(Photon::BooleanID::MVAIDWP90, Photon_mvaID_WP90[i]);
        if(Run == 3){
            photon.SetCBIDBit(Photon::CutBasedID::CUTBASED, Photon_cutBased[i]);
        }
        else if(Run == 2){
            photon.SetCBIDBit(Photon::CutBasedID::CUTBASED, Photon_cutBased_RunII[i]);
        }


        Photons.push_back(photon);
    }
    return Photons;
}

RVec<Photon> AnalyzerCore::GetPhotons(TString id, double ptmin, double fetamax) {
    RVec<Photon> photons = GetAllPhotons();
    RVec<Photon> out;
    for(unsigned int i=0; i<photons.size(); i++) {
        if(!( photons.at(i).Pt()>ptmin )){
            continue;
        }
        if (!( fabs(photons.at(i).scEta())<fetamax)) {
            continue;
        }
        if(!( photons.at(i).PassID(id))) {
            continue;
        }
        out.push_back(photons.at(i));
    }
    return out;
}

RVec<Jet> AnalyzerCore::GetJets(const TString ID, const float ptmin, const float fetamax) {
    RVec<Jet> jets = GetAllJets();
    RVec<Jet> selected_jets;
    for (const auto &jet: jets) {
        if (jet.Pt() < ptmin) continue;
        if (fabs(jet.Eta()) > fetamax) continue;
        if (! jet.PassID(ID)) continue;
        selected_jets.push_back(jet);
    }
    return selected_jets;
}

RVec<Jet> AnalyzerCore::SelectJets(const RVec<Jet> &jets, const TString ID, const float ptmin, const float fetamax) const{
    RVec<Jet> selected_jets;
    for (const auto &jet: jets) {
        if (jet.Pt() < ptmin) continue;
        if (fabs(jet.Eta()) > fetamax) continue;
        if (! jet.PassID(ID)) continue;
        selected_jets.push_back(jet);
    }
    return selected_jets;
}

RVec<Jet> AnalyzerCore::SelectJets(const RVec<Jet> &jets, const Jet::JetID ID, const float ptmin, const float fetamax) const
{
    RVec<Jet> selected_jets;
    for (const auto &jet : jets)
    {
        if (jet.Pt() < ptmin)
            continue;
        if (fabs(jet.Eta()) > fetamax)
            continue;
        if (!jet.PassID(ID))
            continue;
        selected_jets.push_back(jet);
    }
    return selected_jets;
}

RVec<Jet> AnalyzerCore::JetsVetoLeptonInside(const RVec<Jet> &jets, const RVec<Electron> &electrons, const RVec<Muon> &muons, const float dR) const{
    RVec<Jet> selected_jets;
    for(const auto &jet: jets){
        bool isLeptonInside = false;
        for(const auto &electron: electrons){
            if(jet.DeltaR(electron) < dR){
                isLeptonInside = true;
                break;
            }
        }
        if(isLeptonInside) continue;
        for(const auto &muon: muons){
            if(jet.DeltaR(muon) < dR){
                isLeptonInside = true;
                break;
            }
        }
        if(isLeptonInside) continue;
        selected_jets.push_back(jet);
    }
    return selected_jets;
}

// For Run2, reject jets within the veto map
bool AnalyzerCore::PassVetoMap(const Jet &jet, const RVec<Muon> &AllMuons, const TString mapCategory){
    if (! (Run == 2)) return true;
    // Only apply to the jets with em energy fraction less than 0.9
    if (jet.chEmEF() + jet.neEmEF() > 0.9) return true;

    // Selections should be looser than Analysis jet selections
    bool pass_loose_selection = jet.Pt() > 15.;
    pass_loose_selection = pass_loose_selection && jet.PassID(Jet::JetID::TIGHT);
    pass_loose_selection = pass_loose_selection && (jet.Pt() > 50. || jet.PassID(Jet::JetID::PUID_LOOSE));
    for (const auto &muon: AllMuons){
        pass_loose_selection = pass_loose_selection && (jet.DeltaR(muon) > 0.2);
    }
    bool pass_veto_map = pass_loose_selection && (!myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory));
    return pass_veto_map;
}

// For Run3, reject events if any jet is within the veto map
bool AnalyzerCore::PassVetoMap(const RVec<Jet> &AllJets, const RVec<Muon> &AllMuons, const TString mapCategory) {
    if (! (Run == 3)) return true;
    RVec<Jet> selected_jets;
    RVec<Electron> empty_electrons;

    RVec<Jet> this_jet = SelectJets(AllJets, Jet::JetID::TIGHT, 15., 5.0);
    this_jet = JetsVetoLeptonInside(this_jet, empty_electrons, AllMuons, 0.2);
    for(const auto &jet: this_jet){
        if(jet.chEmEF() + jet.neEmEF() < 0.9) selected_jets.push_back(jet);
    }
    for(const auto &jet: selected_jets){
        if(myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory)) return false;
    }
    return true;
}

bool AnalyzerCore::PassJetVetoMap(const RVec<Jet> &AllJets, const RVec<Muon> &AllMuons, const TString mapCategory){
    RVec<Jet> this_jet = SelectJets(AllJets, Jet::JetID::TIGHT, 15., 5.0);
    this_jet = SelectJets(this_jet, Jet::JetID::PUID_TIGHT, 15., 5.0);
    RVec<Jet> selected_jets;
    RVec<Electron> empty_electrons;
    this_jet = JetsVetoLeptonInside(this_jet, empty_electrons, AllMuons, 0.2);
    for(const auto &jet: this_jet){
        if(jet.chEmEF() + jet.neEmEF() < 0.9) selected_jets.push_back(jet);
    }
    for(const auto &jet: selected_jets){
        if(myCorr->IsJetVetoZone(jet.Eta(), jet.Phi(), mapCategory)) return false;
    }
    return true;
}

RVec<FatJet> AnalyzerCore::GetAllFatJets() {
    
    RVec<FatJet> FatJets;

    for (int i = 0; i < nFatJet; i++) {

        FatJet fatjet;

        RVec<float> pnet_m;
        RVec<float> pnet;

        if(Run == 3){
            pnet_m = {FatJet_particleNetWithMass_H4qvsQCD[i], FatJet_particleNetWithMass_HccvsQCD[i],
               FatJet_particleNetWithMass_HbbvsQCD[i], FatJet_particleNetWithMass_QCD[i],
               FatJet_particleNetWithMass_TvsQCD[i], FatJet_particleNetWithMass_WvsQCD[i],
               FatJet_particleNetWithMass_ZvsQCD[i]};

            pnet = {FatJet_particleNet_QCD[i], FatJet_particleNet_QCD0HF[i],
                    FatJet_particleNet_QCD1HF[i], FatJet_particleNet_QCD2HF[i],
                    FatJet_particleNet_XbbVsQCD[i], FatJet_particleNet_XccVsQCD[i],
                    FatJet_particleNet_XqqVsQCD[i], FatJet_particleNet_XggVsQCD[i],
                    FatJet_particleNet_XteVsQCD[i], FatJet_particleNet_XtmVsQCD[i],
                    FatJet_particleNet_XttVsQCD[i], FatJet_particleNet_massCorr[i]};

            fatjet.SetJetID(FatJet_jetId[i]);
            if(!IsDATA) fatjet.SetGenMatchIDs(FatJet_genJetAK8Idx[i], FatJet_subJetIdx1[i], FatJet_subJetIdx2[i]);
        }
        else if(Run == 2){
            pnet_m = {FatJet_particleNet_H4qvsQCD[i], FatJet_particleNet_HccvsQCD[i],
               FatJet_particleNet_HbbvsQCD[i], FatJet_particleNet_QCD[i],
               FatJet_particleNet_TvsQCD[i], FatJet_particleNet_WvsQCD[i],
               FatJet_particleNet_ZvsQCD[i]};
            
            pnet = {FatJet_particleNetMD_QCD[i], -999., -999., -999.,
                    FatJet_particleNetMD_Xbb[i], FatJet_particleNetMD_Xcc[i],
                    FatJet_particleNetMD_Xqq[i], -999., -999., -999., -999., -999.};
            fatjet.SetJetID(static_cast<unsigned char>(FatJet_jetId_RunII[i]));
            if(!IsDATA) fatjet.SetGenMatchIDs(static_cast<short>(FatJet_genJetAK8Idx_RunII[i]), static_cast<short>(FatJet_subJetIdx1_RunII[i]), static_cast<short>(FatJet_subJetIdx2_RunII[i]));
        }

        fatjet.SetPtEtaPhiM(FatJet_pt[i], FatJet_eta[i], FatJet_phi[i], FatJet_mass[i]);
        fatjet.SetArea(FatJet_area[i]);
        fatjet.SetSDMass(FatJet_msoftdrop[i]);
        fatjet.SetLSF3(FatJet_lsf3[i]);
        fatjet.SetConstituents(FatJet_nBHadrons[i], FatJet_nCHadrons[i], FatJet_nConstituents[i]);
        fatjet.SetBTaggingInfo(FatJet_btagDDBvLV2[i], FatJet_btagDDCvBV2[i], FatJet_btagDDCvLV2[i], FatJet_btagDeepB[i], FatJet_btagHbb[i]);
        fatjet.SetPNetwithMassResults(pnet_m);
        fatjet.SetPNetResults(pnet);
        fatjet.SetSubjettiness(FatJet_tau1[i], FatJet_tau2[i], FatJet_tau3[i], FatJet_tau4[i]);

        FatJets.push_back(fatjet);
    }

    return FatJets;
}

RVec<GenJet> AnalyzerCore::GetAllGenJets() {
    
    RVec<GenJet> GenJets;
    if(IsDATA) return GenJets;

    for (int i = 0; i < nGenJet; i++) {

        GenJet genjet;

        genjet.SetPtEtaPhiM(GenJet_pt[i], GenJet_eta[i], GenJet_phi[i], GenJet_mass[i]);
        if(Run == 3)
        genjet.SetGenFlavours(GenJet_partonFlavour[i], GenJet_hadronFlavour[i]);
        else if(Run == 2)
        genjet.SetGenFlavours(static_cast<short>(GenJet_partonFlavour_RunII[i]), GenJet_hadronFlavour[i]); 
        GenJets.push_back(genjet);
    }

    return GenJets;
}

RVec<GenDressedLepton> AnalyzerCore::GetAllGenDressedLeptons() {
    RVec<GenDressedLepton> GenDressedLeptons;
    if(IsDATA) return GenDressedLeptons;

    for(int i = 0; i < nGenDressedLepton; i++) {
        GenDressedLepton lep;
        lep.SetPtEtaPhiM(GenDressedLepton_pt[i], GenDressedLepton_eta[i], GenDressedLepton_phi[i], GenDressedLepton_mass[i]);
        lep.SetPdgId(GenDressedLepton_pdgId[i]);
        lep.SetHasTauAnc(GenDressedLepton_hasTauAnc[i]);
        GenDressedLeptons.push_back(lep);
    }
    return GenDressedLeptons;
}

RVec<GenIsolatedPhoton> AnalyzerCore::GetAllGenIsolatedPhotons() {
    RVec<GenIsolatedPhoton> GenIsolatedPhotons;
    if(IsDATA) return GenIsolatedPhotons;

    for(int i = 0; i < nGenIsolatedPhoton; i++) {
        GenIsolatedPhoton photon;
        photon.SetPtEtaPhiM(GenIsolatedPhoton_pt[i], GenIsolatedPhoton_eta[i], GenIsolatedPhoton_phi[i], GenIsolatedPhoton_mass[i]);
        GenIsolatedPhotons.push_back(photon);
    }
    return GenIsolatedPhotons;
}

RVec<GenVisTau> AnalyzerCore::GetAllGenVisTaus() {
    RVec<GenVisTau> GenVisTaus;
    if(IsDATA) return GenVisTaus;

    for(int i = 0; i < nGenVisTau; i++) {
        GenVisTau tau;
        tau.SetPtEtaPhiM(GenVisTau_pt[i], GenVisTau_eta[i], GenVisTau_phi[i], GenVisTau_mass[i]);
        tau.SetCharge(GenVisTau_charge[i]);
        tau.SetGenPartIdxMother(GenVisTau_genPartIdxMother[i]);
        tau.SetStatus(GenVisTau_status[i]);
        GenVisTaus.push_back(tau);
    }
    return GenVisTaus;
}

RVec<TrigObj> AnalyzerCore::GetAllTrigObjs() {
    RVec<TrigObj> TrigObjs;
    
    for(int i = 0; i < nTrigObj; i++) {
        TrigObj trigObj;
        trigObj.SetRun(Run);
        trigObj.SetPtEtaPhiM(TrigObj_pt[i], TrigObj_eta[i], TrigObj_phi[i], 0.0); // TrigObj mass is typically 0
        if (Run == 3) {
            trigObj.SetId(static_cast<Int_t>(TrigObj_id[i]));
        } else {
            trigObj.SetId(TrigObj_id_RunII[i]);
        }
        trigObj.SetFilterBits(TrigObj_filterBits[i]);
        TrigObjs.push_back(trigObj);
    }
    return TrigObjs;
}

bool AnalyzerCore::IsHEMElectron(const Electron& electron) const {
    if (DataYear != 2018) return false;

    if (electron.Eta() < -1.25){
        if((electron.Phi() < -0.82) && (electron.Phi() > -1.62)) return true;
    }
    return false;
}

// Gen Matching
void AnalyzerCore::PrintGen(const RVec<Gen> &gens) {
    cout << "===========================================================" << endl;
    //cout << "RunNumber:EventNumber = " << Run << ":" << event << endl;
    cout << "index\tPID\tStatus\tMIdx\tMPID\tStart\tPt\tEta\tPhi\tM" << endl;
    for(unsigned int i=2; i<gens.size(); i++){
        Gen gen = gens.at(i);
        RVec<int> history = TrackGenSelfHistory(gen, gens);
        cout << i << "\t" 
             << gen.PID() << "\t" 
             << gen.Status() << "\t" 
             << gen.MotherIndex() << "\t" 
             << gens.at(gen.MotherIndex()).PID() << "\t" 
             << history[0] << "\t" 
             << fixed << setprecision(2) 
             << gen.Pt() << "\t" 
             << gen.Eta() << "\t" 
             << gen.Phi() << "\t" 
             << gen.M() << endl;
    }
}

RVec<int> AnalyzerCore::TrackGenSelfHistory(const Gen &me, const RVec<Gen> &gens) {
    //returns {index of the first history of the gen, 
    //         index of the last history of the gen's mother}
    int myindex = me.Index();
    if (myindex < 2) return {myindex, -1};

    int mypid = gens.at(myindex).PID();
    int currentidx = myindex;
    int motherindex = me.MotherIndex();
    if (motherindex < 0) return {myindex, -1};
    
    while(gens.at(motherindex).PID() == mypid){
        // Go one generation up
        currentidx = motherindex;
        motherindex = gens.at(motherindex).MotherIndex();
        if(motherindex<0) break;
    }
    RVec<int> out = {currentidx, motherindex};
    return out;
}

Gen AnalyzerCore::GetGenMatchedLepton(const Lepton& lep, const RVec<Gen>& gens){
    //==== find status 1 lepton
    int reco_PID = -999;
    if(lep.LeptonFlavour() == Lepton::ELECTRON) {
        reco_PID = 11;
    } else if(lep.LeptonFlavour() == Lepton::MUON) {
        reco_PID = 13;
    } else {
        cerr << "[AnalyzerCore::GetGenMatchedLepton] input lepton flavour not set" << endl;
        exit(EXIT_FAILURE);
    }

    float min_dR = 0.1;
    Gen gen_closest;
    bool found_match = false;
    for(unsigned int i=2; i<gens.size(); i++){
        const Gen &gen = gens.at(i);
        if( gen.Status() != 1 ) continue; 
        if( abs( gen.PID() ) != reco_PID ) continue; 
        if( gen.MotherIndex() < 0 ) continue; // reject ISR
        if( gen.DeltaR( lep ) < min_dR ){ // dR matching
            min_dR = gen.DeltaR( lep ) ;
            gen_closest = gen;
            found_match = true;
        }
    }
    if(!found_match){
        // Return a Gen with invalid index if no match found
        Gen dummy;
        dummy.SetIndexPIDStatus(-1, 0, -1); 
        return dummy;
    }
    return gen_closest;
}

Gen AnalyzerCore::GetGenMatchedMuon(const Muon& muon, const RVec<Gen>& gens){
    // Find status 1 muon
    int reco_PID = 13;

    float max_dR = 0.2;
    float min_dR = 0.1;
    Gen gen_closest;
    bool found_match = false;
    float distance = 100000;
    for (const auto &gen : gens) {
        if (gen.Status() != 1) continue;
        if (abs(gen.PID()) != reco_PID) continue;
        if (gen.MotherIndex() < 0) continue; // reject ISR
        if (gen.DeltaR(muon) > max_dR) continue;
        float this_distance = pow(gen.DeltaR(muon)/0.005,2) + pow((muon.Pt()/gen.Pt()-1)/0.02,2);
        if (this_distance < distance) {
            distance = this_distance;
            gen_closest = gen;
            found_match = true;
        }
    }
    if (!found_match) {
        Gen dummy;
        return dummy;
    }
    return gen_closest;
}

Gen AnalyzerCore::GetGenMatchedPhoton(const Lepton& lep, const RVec<Gen>& gens){
    //Find if there is a photon candidate for a source of external conversion (similar pt, direction of lepton)
    //1) Validity of cuts are checked for electron with PT>10. Note that the previous dPtRel<0.2 cut in the CatNtuple analysis was optimized for PT(e)>25.
    //   External conversion rate proportionate to M^{-2}, thus muon external conversion rate is negligible.
    //2) As the algorithm targets stable photon near lepton, only status-1 photon should be investigated theoretically. 
    //   However, in some events, last status in photon history is 23. This is believed to be due to skimming of GEN history between pythia & MiniAOD. 
    //   This case is also covered here. Note that IsFinalPhotonSt23 has to analyze whole history, i.e. computationally heavier part in this function.
    //   Therefore it is put at last part of criteria for efficient functioning.
    //For more details about the optimization, ask J.B..
    float min_dR = 0.2;//1)
    Gen gen_closest;
    bool found_match = false;
    float pt_min = 10., dPtRelmax=0.5, dRmax=0.2;//1)
    for(unsigned int i=2; i<gens.size(); i++){
        const Gen &gen = gens.at(i);
        if( gen.MotherIndex() < 0 ) continue;
        if( ! ( abs(gen.PID())==22 && (gen.Status()==1 || gen.Status()==23) ) ) continue;
        if( gen.Pt() < pt_min ) continue;
        if( !(lep.Pt()/gen.Pt()>(1.-dPtRelmax) && lep.Pt()/gen.Pt()<(1.+dPtRelmax)) ) continue;
        if( gen.DeltaR( lep ) > dRmax ) continue;
        if( gen.Status()==23 && !IsFinalPhotonSt23_Public(gens) ) continue;

        if( gen.DeltaR( lep ) < min_dR ){
            min_dR = gen.DeltaR( lep ) ;
            gen_closest = gen;
            found_match = true;
        }
    }
    if(!found_match){
        // Return a Gen with invalid index if no match found
        Gen dummy;
        dummy.SetIndexPIDStatus(-1, 0, -1); 
        return dummy;
    }
    return gen_closest;
}

bool AnalyzerCore::IsFinalPhotonSt23_Public(const RVec<Gen>& gens){
    //In Some XG proc events, it is observed that some of photons' last status is 23. Presumably due to skimming of generator history between pythia and MiniAOD.
    //The function returns if this is the case.
    //And this is designed only for 1 hard photon case as W+G or Z+G or TT+G
    bool IsFinalGammaStatus23 = false;
    bool HasStatus23Photon    = false;
    for(unsigned int i=2; i<gens.size(); i++){
        const Gen &gen = gens.at(i);
        int fpid  = fabs(gen.PID());
        int GenSt = gen.Status();
        int MPID_direct= gens.at(gen.MotherIndex()).PID();
        if( !((fpid!=22 && MPID_direct==22) || (fpid==22 && (GenSt==23||GenSt==1))) ) continue;

        RVec<int> my_history, mom_history;
        my_history  = TrackGenSelfHistory(gen, gens);
        int LastSelfIdx     = my_history[0];
        int MotherIdx       = my_history[1];
        int LastSelfSt      = gens.at(LastSelfIdx).Status();
        int LastSelfMIdx=-1, MStatus_orig=-1;
        if(MotherIdx!=-1){
            mom_history  = TrackGenSelfHistory(gens.at(MotherIdx), gens);
            LastSelfMIdx = mom_history[0];
            MStatus_orig = gens.at(LastSelfMIdx).Status();
        }

        if(fpid==22){
            if(GenSt==23) {
                HasStatus23Photon=true; 
                IsFinalGammaStatus23=true;
            } else if(GenSt==1 && LastSelfSt==23) {
                IsFinalGammaStatus23=false;
                break; //a
            }
        } else if( MPID_direct==22 && MStatus_orig==23 ){ 
            IsFinalGammaStatus23=false; 
            break; //b
        }
    }
    if(!HasStatus23Photon) return false;
    return IsFinalGammaStatus23;
    //**footnotes
    // a) Status-23 photon's last is 1. Thus status-23 photon is not the last history.
    // b) Daughter particle of status-23 photon is found. Thus status-23 photon is not the last history. 
}

bool AnalyzerCore::IsFromHadron(const Gen& me, const RVec<Gen>& gens) {
    bool out = false;
    int myindex = me.Index();
    if(myindex<2) return true;
    RVec<int> my_history = TrackGenSelfHistory(me, gens);
    if (my_history[1] < 0) return true;
    Gen          Start = gens.at( my_history[0] );
    Gen MotherOf_Start = gens.at( my_history[1] );

    //==== Status 21~29 are from hardprocess
    //==== Means it's lepton from hardprocess
    //==== e.g., leptons from Z start their lives with status 23
    if( 20 < Start.Status() && Start.Status() < 30 ) return false;

    Gen current_me = Start; // me will always be Start
    Gen current_mother = Start; // initializing
    while( current_mother.Index() >= 2 ){
        RVec<int> current_history = TrackGenSelfHistory(current_me, gens);
        if (current_history[1] == -1) { // current_me is the initial partons
            out = true;
            break;
        }
        //==== Go one generation up
        current_me = gens.at(current_history[1]);

       //==== Now look at mother of previous "me"
        current_mother = gens.at(current_history[1]);

        RVec<int> current_mother_history = TrackGenSelfHistory(current_mother, gens);
        Gen StartOf_current_mother = gens.at(current_mother_history[0]);
        int current_mother_PID = current_mother.PID();

        //==== From Z,W,H,t
        if( current_mother_PID==23 || current_mother_PID==24 || current_mother_PID==25 || current_mother_PID==6 ){
            out = false;
            break;
        }
    
        //==== From Signal
        if( IsSignalPID(current_mother_PID) ){
            out = false;
            break;
        }
        if( (current_mother_PID==11 || current_mother_PID==13 || current_mother_PID==15 || current_mother_PID==22) && (StartOf_current_mother.Status()>20 && StartOf_current_mother.Status()<30)){
            out = false;
            break;
        }
        if( current_mother_PID>50 ){
            out=true;
            break;
        }
        if( (current_mother_PID>=1 && current_mother_PID<=5) || current_mother_PID==21 ){
            out=true; 
            break;
        }
    }

    return out;
}

bool AnalyzerCore::IsSignalPID(const int &pid) {
    const unsigned int fpid = abs(pid);

    //==== HeavyNeutrino
    if(fpid>=9900000) return true;
    //==== CP-odd Higgs
    if(fpid==36) return true;
    return false;
}

//==== [Type]
//====  1 : EWPrompt
//====  2 : Signal Daughter
//====  3 : EWtau daughter
//====  4 : Internal Conversion daughter from t/EWV/EWlep(Implicit,Explicit)
//====  5 : Internal Conversion daughter from HardScatterPhoton
//====  6 : From Offshell W (mother == 37)
//==== -1 : Unmatched & not EW Conversion candidate
//==== -2 : Hadron daughter
//==== -3 : Daughter of tau from hadron or parton
//==== -4 : Internal conversion daughter(implicit,explicit) having hadronic origin 
//==== -5 : External conversion candidate(Hard scattered photon)
//==== -6 : External conversion from t/EWV/EWlep
//==== (-4: Daughter of Non-hard scattered photon & has parton or hadron ancestor OR implicit Conv from quark)
//==== -------------------------------------------------------------
//====  0 : Error
//==== >0 : Non-fake; Non-hadronic origin
//==== <0 : Fakes   ; Hadronic origin or external conversion

int AnalyzerCore::GetLeptonType(const Gen& gen, const RVec<Gen>& gens) {
    int LeptonType = GetLeptonType_Public(gen.Index(), gens);//4)
    if(LeptonType>=4 && LeptonType<=5){//5)
        int NearbyPrElType = GetPrElType_InSameSCRange_Public(gen.Index(), gens);
        if(NearbyPrElType>0) LeptonType = NearbyPrElType;
    }
    return LeptonType;
}

int AnalyzerCore::GetLeptonType(const Lepton& lep, const RVec<Gen>& gens) {
    int LeptonType=0, MatchedTruthIdx=-1; 
    Gen gen_closest = GetGenMatchedLepton(lep, gens);
    //cout << "[AnalyzerCore::GetLeptonType] [Reco] pt = " << lep.Pt() << "\t, eta = " << lep.Eta() << endl;
    //cout << "[AnalyzerCore::GetLeptonType] [Gen] Index = " << gen_closest.Index() << endl;

    //==== No matched gen lepton
    if( gen_closest.IsEmpty() ){
        Gen gen_photon_closest = GetGenMatchedPhoton(lep, gens);
        int NearPhotonType = GetGenPhotonType(gen_photon_closest,gens);
        if (NearPhotonType == 1)      LeptonType = -5;
        else if (NearPhotonType == 2) LeptonType = -6;
        else                          LeptonType = -1; // NearPhotonType <= 0
    } else { //==== Has macthed gen lepton
        MatchedTruthIdx = gen_closest.Index();
        LeptonType = GetLeptonType_Public(MatchedTruthIdx, gens);
        if(LeptonType>=4 && LeptonType<=5){
            int NearbyPrElType = GetPrElType_InSameSCRange_Public(MatchedTruthIdx, gens);
            if(NearbyPrElType>0) LeptonType = NearbyPrElType;
        }
    }
    return LeptonType;
    //**footnote
    //1) matched to no gen-lepton nor gen-photon -> mis-reco.
    //   matched to no gen-lepton, but to photon with hadronic origin -> mis-reco. (e.g. pions->e) or external conversion from photon in jets
    //2) matched to no gen-lepton, but to photon with non-hadronic origin (hard process) -> external conversion from photon with non-hadronic source (ME-level)
    //3) matched to no gen-lepton, but to photon with non-hadronic origin (soft QED radiation) -> external conversion from photon with non-hadronic source (PS-level)
    //4) matched to gen-lepton, categorize based on the truth categorization algo.: AnalyzerCore::GetLeptonType_Public(int TruthIdx, std::vector<Gen>& TruthColl)
    //5) collimated e/gm objects are merged in SC, hence if there is prompt electron within SC-merging range, reco-electron's properties actually represent pre-QED-FSR prompt-electron,
    //   rather than the closest internal conversion electron. Therefore shift the type to the prompt lepton's type.
    //- Note: distinction between type 4 vs. 5 and type -5 vs. -6 is unphysical. it is only for debugging.
}

int AnalyzerCore::GetLeptonType_Public(const int& genIdx, const RVec<Gen>& gens) {
    //Type : 1:EW-Prompt / 2: BSM-Prompt / 3:EW/BSM-Prompt-Tau Daughter / 6: From Offshell W (mother == 37)
    //       4:Internal Conversion from Soft QED Radiation (PS-level) / 5:Internal Conversion from Hard Process Photon (ME-level)
    //      -1:Unmatched & not EW Conversion candidate (mis-reco. or external conversion within jets)
    //      -2:Hadron Daughter / -3:Daughter of Tau from Hadron or Parton / -4:Internal Conversion Daughter having hadronic origin
    //      -5:External Conversion from Hard process photon / -6:External conversion from t/EWV/EWlep
    //      (-4:Daughter of Non-hard scattered photon & has parton or hadron ancestor OR implicit Conv from quark)
    //       0:Error / >0: Non-fake: Non-hadronic origin / <0 : Fakes: Hadronic origin or external conversion
    //Note: There is no physical distinction between type 4 vs. 5 and type -5 vs. -6. The distinction is only intended for straightforward debugging.

    //Only consider Status 1 lepton
    if (genIdx<2) return 0;
    if (gens.at(genIdx).Status()!=1) return 0;
    if( !(fabs(gens.at(genIdx).PID())==11 || fabs(gens.at(genIdx).PID())==13) ) return 0;

    int LeptonType=0;
    int MPID=0, GrMPID=0;
    RVec<int> my_history, mom_history, grmom_history;
    int LastSelfIdx=-1, MotherIdx=-1, LastSelfMIdx=-1, GrMotherIdx=-1, LastSelfGrMIdx=-1;
    int Status_orig=0, MStatus_orig=0, MStatus_last=0, GrMStatus_orig=0, GrMStatus_last=0;
    bool HadronicOrigin = false;

    my_history      = TrackGenSelfHistory(gens.at(genIdx), gens);
    LastSelfIdx     = my_history[0];
    MotherIdx       = my_history[1];
    Status_orig     = gens.at(LastSelfIdx).Status();
    HadronicOrigin  = IsFromHadron(gens.at(genIdx), gens);

    if(   MotherIdx!=-1   ){ 
        mom_history  = TrackGenSelfHistory(gens.at(MotherIdx), gens);
        LastSelfMIdx = mom_history[0];
        GrMotherIdx  = mom_history[1];
        MPID         = gens.at(MotherIdx).PID();
        MStatus_orig = gens.at(LastSelfMIdx).Status();
        MStatus_last = gens.at(MotherIdx).Status();
    }
    if(  GrMotherIdx!=-1  ){ 
        grmom_history  = TrackGenSelfHistory(gens.at(GrMotherIdx), gens);
        LastSelfGrMIdx = grmom_history[0];
        GrMPID         = gens.at(GrMotherIdx).PID();
        GrMStatus_orig = gens.at(LastSelfGrMIdx).Status();
        GrMStatus_last = gens.at(GrMotherIdx).Status();
    }

    if     ( genIdx==-1 )                                         LeptonType= 0;
    else if( fabs(MPID)==23 || fabs(MPID)==24 || fabs(MPID)==25 ) LeptonType= 1;
    else if( IsSignalPID(MPID) )                                  LeptonType= 2;
    else if( abs(MPID) == 37 )                                    LeptonType= 6; 
    else if( Status_orig>20 && Status_orig<30 )                   LeptonType= 1;//1)
    else if( fabs(MPID)>50 )                                      LeptonType=-2;
    else if( fabs(MPID)==15 && MStatus_last==2 ){
        if     ( fabs(GrMPID)==23 || fabs(GrMPID)==24 || fabs(GrMPID)==25 )    LeptonType= 3;
        else if( IsSignalPID(GrMPID) )                                         LeptonType= 3;
        else if( MStatus_orig>20  && MStatus_orig<30  )                        LeptonType= 3;//1)
        else if( HadronicOrigin )                                              LeptonType=-3;//2-a)
        else if( fabs(GrMPID)==22  && GrMStatus_orig>20 && GrMStatus_orig<30 ) LeptonType= 5;//2-b)
        else if( fabs(GrMPID)==22 )                                            LeptonType= 4;//2-c)
        else if( (fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15) && GrMStatus_last!=2 ) LeptonType= 4;//2-d)
        else                                                                   LeptonType= 0;
    }
    else if( fabs(MPID)==22 ){
        if( MStatus_orig>20 && MStatus_orig<30 )                            LeptonType= 5;//3-a)
        else if( HadronicOrigin )                                           LeptonType=-4;//3-b)
        else if( fabs(GrMPID)==24 || fabs(GrMPID)==23 || fabs(GrMPID)==6  ) LeptonType= 4;//3-c)
        else if( fabs(GrMPID)==11 || fabs(GrMPID)==13 || fabs(GrMPID)==15 ) LeptonType= 4;//3-d)
        else                                                                LeptonType= 0;
    }
    else if( (fabs(MPID)==11 || fabs(MPID)==13 || fabs(MPID)==15) && MStatus_last!=2 && !HadronicOrigin ) LeptonType= 4;//4-a)
    else if( ((fabs(MPID)>=1 && fabs(MPID)<=5) || fabs(MPID)==21) && MStatus_last!=2 )                    LeptonType=-4;//4-b)
    else if( fabs(MPID)==6 ) LeptonType=4;//4-c)
    else LeptonType=0;

    return LeptonType;
}

//==== [Type]
//====  0 : Invalid input or Error or HardScatter is input when hardscatter is not final state
//====  1 : HardScatter / 2: Else prompt daughter(l,V,t)
//==== -1 : Reserved for unmatched(Not used now) / -2: Hadronic origin
int AnalyzerCore::GetGenPhotonType(const Gen& genph, const RVec<Gen>& gens) {
    if(genph.IsEmpty()) return 0;
    int genph_index = genph.Index();
    if( genph_index<2 ) return 0;
    if( !(gens.at(genph_index).PID()==22 && (gens.at(genph_index).Status()==1 || gens.at(genph_index).Status()==23)) ) return 0;

    if(gens.at(genph_index).Status()==23){
        if(IsFinalPhotonSt23_Public(gens)) return 1;
        else                               return 0;
    }//From this pt, only St1 Photon is treated.

    RVec<int> phhist = TrackGenSelfHistory(genph, gens);
    const Gen&          Start = gens.at(phhist[0]);
    const Gen& MotherOf_Start = gens.at(phhist[1]);
    const int& MotherOf_Start_PID = abs(MotherOf_Start.PID()); // |PID|
    const bool fromhadron = IsFromHadron(genph, gens);

    if     ( 20 < Start.Status() && Start.Status() < 30 ) return 1;//1)
    else if( MotherOf_Start_PID==23 || MotherOf_Start_PID==25) return 1;//2)
    else if( MotherOf_Start_PID==24 || MotherOf_Start_PID==6 || IsSignalPID(MotherOf_Start_PID) ) return 2;//3)
    else if( fromhadron ) return -2;//4)
    else if( MotherOf_Start_PID==11 || MotherOf_Start_PID==13 || MotherOf_Start_PID==15 ) return 2;//5)
    else return 0;

    //**footnote
    //1) Hard process photon's original status is 20's, unless the history is skimmed in MiniAOD.
    //2) If part of intermediate history is skimmed at MinoAOD, then last photon status is 1 without any preceding history lines.
    //   In that case, intermediate mother is written in history, and this line catches the case.
    //3) top and charged bosons can radiate photons.
    // - Note that distinction between 1 and 2 is not physical. You should use both 1 & 2 for prompt photons.
    //4) this category does not include top, but photons from hadrons and quarks. Predominantly, in most of the cases they are daughter of pi0.
    //   Rarer, but other mesons as eta, B, or even some quarks can also radiate energetic photons.
    //5) Photons radiated from lepton FSR. Sometimes they are quite energetic.
}

int AnalyzerCore::GetPrElType_InSameSCRange_Public(int genIdx, const RVec<Gen>& gens){
    //Abbreviation: Get Prompt Electron Type In Same Supercluster Range (Public (shared) version of J.B.'s original gen-matching code)
    //Prompt e>eee (int. conv.) case, collimated electrons can be merged in one SC & track is selected among them, and reconstructed as single electron. 
    //In this case, still there will be a nearby prompt electron of LepType 1-3 within supercluster merging range.
    //Supercluster merging range: |dphi|<0.3/0.3 (EB/EE), |deta|<0/0.07 (EB/EE) (additionally put 0.03 for marginal difference from size of crystals)
    //Ref:JINST 10 (2015) P06005, arXiv:1502.02701
    //Return Value: 1/2/3:LeptonType-1/2/3 electron found in same SC range
    //             -1: No LeptonType-1/2/3 electron found in same SC range
    //Note: Hadronic origins are not included in this algo. as I am not sure of effect of many nearby particles in jets on the ele-reco performance.
    //      + it is not of my interest yet.

    if(genIdx<2) return -1;
    if(abs(gens.at(genIdx).PID())!=11) return -1;
    if(gens.at(genIdx).Status()!=1) return -1;

    float dPhiMax=0.3, dEtaMax=0.1;
    int NearbyElType=-1;

    for(unsigned int it_gen=2; it_gen<gens.size(); it_gen++){
        if(gens.at(it_gen).Status()!=1) continue;
        if(abs(gens.at(it_gen).PID())!=11) continue;
        if(fabs(gens.at(genIdx).Eta()-gens.at(it_gen).Eta())>dEtaMax) continue;
        if(gens.at(genIdx).DeltaPhi(gens.at(it_gen))>dPhiMax) continue;

        int LepType = GetLeptonType_Public(it_gen, gens);
        if(!(LepType>=1 && LepType<=3)) continue;
        NearbyElType  = LepType;
    }
    return NearbyElType;
}

// Histogram Handlers
void AnalyzerCore::FillHist(const TString &histname, float value, float weight, int n_bin, float x_min, float x_max) {
    auto histkey = string(histname);
    auto it = histmap1d.find(histkey);
    if (it == histmap1d.end()) {
        TH1 *this_hist;
        if(useTH1F) this_hist = new TH1F(histkey.c_str(), "", n_bin, x_min, x_max);
        else this_hist = new TH1D(histkey.c_str(), "", n_bin, x_min, x_max);
        this_hist->SetDirectory(nullptr);
        histmap1d[histkey] = this_hist;
        this_hist->Fill(value, weight);
    }
    else {
        it->second->Fill(value, weight);
    }
}

void AnalyzerCore::FillHist(const TString &histname, float value, float weight, int n_bin, float *xbins) {
    auto histkey = string(histname.Data());
    auto it = histmap1d.find(histkey);
    if (it == histmap1d.end()) {
        TH1 *this_hist;
        if(useTH1F) this_hist = new TH1F(histkey.c_str(), "", n_bin, xbins);
        else this_hist = new TH1D(histkey.c_str(), "", n_bin, xbins);
        this_hist->SetDirectory(nullptr);
        histmap1d[histkey] = this_hist;
        this_hist->Fill(value, weight);
    }
    else {
        it->second->Fill(value, weight);
    }
}

void AnalyzerCore::FillHist(const TString &histname, float value_x, float value_y, float weight,
                                                   int n_binx, float x_min, float x_max,
                                                   int n_biny, float y_min, float y_max) {
    auto histkey = string(histname);
    auto it = histmap2d.find(histkey);
    if (it == histmap2d.end()) {
        TH2 *this_hist;
        if(useTH1F) this_hist = new TH2F(histkey.c_str(), "", n_binx, x_min, x_max, n_biny, y_min, y_max);
        else this_hist = new TH2D(histkey.c_str(), "", n_binx, x_min, x_max, n_biny, y_min, y_max);
        this_hist->SetDirectory(nullptr);
        histmap2d[histkey] = this_hist;
        this_hist->Fill(value_x, value_y, weight);
    }
    else {
        it->second->Fill(value_x, value_y, weight);
    }
}

void AnalyzerCore::FillHist(const TString &histname, float value_x, float value_y, float weight,
                                                    int n_binx, float *xbins,
                                                    int n_biny, float *ybins) {
    auto histkey = string(histname);
    auto it = histmap2d.find(histkey);
    if (it == histmap2d.end()) {
        TH2 *this_hist;
        if(useTH1F) this_hist = new TH2F(histkey.c_str(), "", n_binx, xbins, n_biny, ybins);
        else this_hist = new TH2D(histkey.c_str(), "", n_binx, xbins, n_biny, ybins);
        this_hist->SetDirectory(nullptr);
        histmap2d[histkey] = this_hist;
        this_hist->Fill(value_x, value_y, weight);
    }
    else {
        it->second->Fill(value_x, value_y, weight);
    }
}

void AnalyzerCore::FillHist(const TString &histname, float value_x, float value_y, float value_z, 
                                      float weight, int n_binx, float x_min, float x_max,
                                                    int n_biny, float y_min, float y_max,
                                                    int n_binz, float z_min, float z_max) {
    auto histkey = string(histname);
    auto it = histmap3d.find(histkey);
    if (it == histmap3d.end()) {
        TH3 *this_hist;
        if(useTH1F) this_hist = new TH3F(histkey.c_str(), "", n_binx, x_min, x_max, n_biny, y_min, y_max, n_binz, z_min, z_max);
        else this_hist = new TH3D(histkey.c_str(), "", n_binx, x_min, x_max, n_biny, y_min, y_max, n_binz, z_min, z_max);
        this_hist->SetDirectory(nullptr);
        histmap3d[histkey] = this_hist;
        this_hist->Fill(value_x, value_y, value_z, weight);
    }
    else {
        it->second->Fill(value_x, value_y, value_z, weight);
    }
}

void AnalyzerCore::FillHist(const TString &histname, float value_x, float value_y, float value_z, 
                                      float weight, int n_binx, float *xbins,
                                                    int n_biny, float *ybins,
                                                    int n_binz, float *zbins) {
    auto histkey = string(histname);
    auto it = histmap3d.find(histkey);
    if (it == histmap3d.end()) {
        TH3 *this_hist;
        if(useTH1F) this_hist = new TH3F(histkey.c_str(), "", n_binx, xbins, n_biny, ybins, n_binz, zbins);
        else this_hist = new TH3D(histkey.c_str(), "", n_binx, xbins, n_biny, ybins, n_binz, zbins);
        this_hist->SetDirectory(nullptr);
        histmap3d[histkey] = this_hist;
        this_hist->Fill(value_x, value_y, value_z, weight);
    }
    else {
        it->second->Fill(value_x, value_y, value_z, weight);
    }
}
TTree* AnalyzerCore::NewTree(const TString &treename, const RVec<TString> &keeps, const RVec<TString> &drops){
    auto treekey = string(treename);
    auto it = treemap.find(treekey);
    if (it == treemap.end()){
        //if keeps and drops are empty, make new tree
        if(keeps.size() == 0 && drops.size() == 0){
            TTree *newtree = new TTree(treekey.c_str(), "");
            treemap[treekey] = newtree;
            return newtree;
        } else{
            //check tree is empty. 
            if(fChain->GetEntries() == 0){
                cout << "[AnalyzerCore::NewTree] fChain is empty." << endl;
                exit(0);
            }
            TTree *newtree = fChain->CloneTree(0);
            newtree->SetName(treekey.c_str());
            for (const auto &drop : drops) {
                newtree->SetBranchStatus(drop, 0);
            }
            for (const auto &keep : keeps) {
                newtree->SetBranchStatus(keep, 1);
            }
            treemap[treekey] = newtree;
            unordered_map<string, TBranch*> this_branchmap;
            branchmaps[newtree] = this_branchmap;
            return newtree;
        }
    } else{
        return it->second;
    }
}

TTree* AnalyzerCore::GetTree(const TString &treename){
    auto treekey = string(treename);
    auto it = treemap.find(treekey);
    if (it == treemap.end()){
        cout << "[AnalyzerCore::GetTree] Tree " << treename << " not found" << endl;
        exit(ENODATA);
    }
    return it->second;
}

void AnalyzerCore::SetBranch(const TString &treename, const TString &branchname, void *address, const TString &leaflist){
    try{
       void* this_address = address;
        TTree *tree = GetTree(treename);
        
        unordered_map<string, TBranch*>* this_branchmap = &branchmaps[tree];
        auto it = this_branchmap->find(string(branchname));

        if (it == this_branchmap->end()) {
            auto br = tree->Branch(branchname, this_address, leaflist);
            this_branchmap->insert({string(branchname), br});
        } else {
            it->second->SetAddress(this_address);
        }
    }
    catch(int e){
        cout << "[AnalyzerCore::SetBranch] Error get tree: " << treename.Data() << endl;
        exit(e);
    }
}

template void AnalyzerCore::SetBranch_Vector<int>(const TString &, const TString &, std::vector<int> &);
template void AnalyzerCore::SetBranch_Vector<float>(const TString &, const TString &, std::vector<float> &);
template void AnalyzerCore::SetBranch_Vector<double>(const TString &, const TString &, std::vector<double> &);
template void AnalyzerCore::SetBranch_Vector<bool>(const TString &, const TString &, std::vector<bool> &);


void AnalyzerCore::FillTrees(const TString &treename) {
    if (treename == "") {
        for (const auto &pair : treemap) {
            const string &treename = pair.first;
            TTree *tree = pair.second;
            tree->Fill();
        }
        this_floats.clear();
        this_ints.clear();
        this_bools.clear();
        this_floats.shrink_to_fit();
        this_ints.shrink_to_fit();
        this_bools.shrink_to_fit(); // Mandatory;
    } else {
        // Convert treeName to std::string for comparison
        std::string treeNameStr(treename.Data());

        auto it = treemap.find(treeNameStr);
        if (it != treemap.end()) {
            // Tree with the given name exists, fill it
            TTree *tree = it->second;
            tree->Fill();
            this_floats.clear();
            this_ints.clear();
            this_bools.clear();
            this_floats.shrink_to_fit();
            this_ints.shrink_to_fit();
            this_bools.shrink_to_fit(); // Mandatory;
        } else {
            // Handle the case where the treeName is not found in the map
            throw std::runtime_error("[AnalyzerCore::FillTrees] Tree with name '" + treeNameStr + "' not found in treemap.");
        }
    }
}


void AnalyzerCore::WriteHist() {
    cout << "[AnalyzerCore::WriteHist] Writing histograms to " << outfile->GetName() << endl;
    std::vector<std::pair<std::string, TH1 *>> sorted_histograms1d(histmap1d.begin(), histmap1d.end());
    std::vector<std::pair<std::string, TH2 *>> sorted_histograms2d(histmap2d.begin(), histmap2d.end());
    std::vector<std::pair<std::string, TH3 *>> sorted_histograms3d(histmap3d.begin(), histmap3d.end());
    std::sort(sorted_histograms1d.begin(), sorted_histograms1d.end(),
              [](const std::pair<std::string, TH1 *> &a, const std::pair<std::string, TH1 *> &b)
              {
                  return a.first < b.first;
              });
    std::sort(sorted_histograms2d.begin(), sorted_histograms2d.end(),
              [](const std::pair<std::string, TH2 *> &a, const std::pair<std::string, TH2 *> &b)
              {
                  return a.first < b.first;
              });
    std::sort(sorted_histograms3d.begin(), sorted_histograms3d.end(),
              [](const std::pair<std::string, TH3 *> &a, const std::pair<std::string, TH3 *> &b)
              {
                  return a.first < b.first;
              });
    for (const auto &pair: sorted_histograms1d) {
        const string &histname = pair.first;
        TH1 *hist = pair.second;
        cout << "[AnalyzerCore::WriteHist] Writing 1D histogram: " << histname << endl;
        // Split the directory and name
        // e.g. "dir1/dir2/histname" -> "dir1/dir2", "histname"
        // e.g. "histname" -> "", "histname"
        size_t last_slash = histname.find_last_of('/');
        string this_prefix, this_name;
        last_slash == string::npos ? this_prefix = "" : this_prefix = histname.substr(0, last_slash);
        last_slash == string::npos ? this_name = histname : this_name = histname.substr(last_slash + 1);

        TDirectory *this_dir = outfile->GetDirectory(this_prefix.c_str());
        if (!this_dir) outfile->mkdir(this_prefix.c_str());
        outfile->cd(this_prefix.c_str());
        hist->Write(this_name.c_str());
    }
    for (const auto &pair: sorted_histograms2d) {
        const string &histname = pair.first;
        cout << "[AnalyzerCore::WriteHist] Writing 2D histogram: " << histname << endl;
        TH2 *hist = pair.second;
        // Split the directory and name
        // e.g. "dir1/dir2/histname" -> "dir1/dir2", "histname"
        // e.g. "histname" -> "", "histname"
        size_t last_slash = histname.find_last_of('/');
        string this_prefix, this_name;
        last_slash == string::npos ? this_prefix = "" : this_prefix = histname.substr(0, last_slash);
        last_slash == string::npos ? this_name = histname : this_name = histname.substr(last_slash + 1);

        TDirectory *this_dir = outfile->GetDirectory(this_prefix.c_str());
        if (!this_dir) outfile->mkdir(this_prefix.c_str());
        outfile->cd(this_prefix.c_str());
        hist->Write(this_name.c_str());
    }
    for (const auto &pair: sorted_histograms3d) {
        const string &histname = pair.first;
        cout << "[AnalyzerCore::WriteHist] Writing 3D histogram: " << histname << endl;
        TH3 *hist = pair.second;
        // Split the directory and name
        // e.g. "dir1/dir2/histname" -> "dir1/dir2", "histname"
        // e.g. "histname" -> "", "histname"
        size_t last_slash = histname.find_last_of('/');
        string this_prefix, this_name;
        last_slash == string::npos ? this_prefix = "" : this_prefix = histname.substr(0, last_slash);
        last_slash == string::npos ? this_name = histname : this_name = histname.substr(last_slash + 1);

        TDirectory *this_dir = outfile->GetDirectory(this_prefix.c_str());
        if (!this_dir) outfile->mkdir(this_prefix.c_str());
        outfile->cd(this_prefix.c_str());
        hist->Write(this_name.c_str());
    }
    for (const auto &pair: treemap) {
        const string &treename = pair.first;
        cout << "[AnalyzerCore::WriteHist] Writing tree: " << treename << endl;
        TTree *tree = pair.second;

        size_t last_slash = treename.find_last_of('/');
        string this_prefix, this_name;
        last_slash == string::npos ? this_prefix = "" : this_prefix = treename.substr(0, last_slash);
        last_slash == string::npos ? this_name = treename : this_name = treename.substr(last_slash + 1);

        TDirectory *this_dir = outfile->GetDirectory(this_prefix.c_str());
        if (!this_dir) outfile->mkdir(this_prefix.c_str());
        outfile->cd(this_prefix.c_str());
        TTree* temptree = tree->CloneTree(-1);//this is because the tree contains lot of empty disabled branch. I don't know better way to handle this, so just keep memory-consuming way for now.
        temptree->Write();
        delete temptree;
        delete tree;
    }
    cout << "[AnalyzerCore::WriteHist] Writing histograms done" << endl;
    outfile->Close();
}
