#include "TrigObj.h"

ClassImp(TrigObj)

TrigObj::TrigObj() : Particle() {
    j_id = -999;
    j_filterBits = 0;
    j_run = 2; // Default to Run 2
}

TrigObj::~TrigObj() {}

// Electron filter bits
// Run2
// bit 0 = *CaloIdLTrackIdLIsoVL*TrackIso*Filter
// bit 1 = hltEle*WPTight*TrackIsoFilter*
// bit 2 = hltEle*WPLoose*TrackIsoFilter
// bit 3 = *OverlapFilter*IsoEle*PFTau*
// bit 4 = hltEle*Ele*CaloIdLTrackIdLIsoVL*Filter
// bit 5 = hltMu*TrkIsoVVL*Ele*CaloIdLTrackIdLIsoVL*Filter*
// bit 6 = hltOverlapFilterIsoEle*PFTau*
// bit 7 = hltEle*Ele*Ele*CaloIdLTrackIdLDphiLeg*Filter
// bit 8 = hltL3fL1Mu*DoubleEG*Filtered* || hltMu*DiEle*CaloIdLTrackIdLElectronleg*Filter
// bit 9 = hltL3fL1DoubleMu*EG*Filter* || hltDiMu*Ele*CaloIdLTrackIdLElectronleg*Filter
// bit 10 = hltEle32L1DoubleEGWPTightGsfTrackIsoFilter && hltEGL1SingleEGOrFilter
// bit 11 = hltEle*CaloIdVTGsfTrkIdTGsfDphiFilter
// bit 12 = HLT_Ele*PFJet*
// bit 13 = hltEG175HEFilter || hltEG200HEFilter
    
// Run3
// bit 0 = *CaloIdLTrackIdLIsoVL*TrackIso*Filter
// bit 1 = hltEle*WPTight*TrackIsoFilter*
// bit 2 = hltEle*WPLoose*TrackIsoFilter
// bit 3 = *OverlapFilter*IsoEle*PFTau*
// bit 4 = hltEle*Ele*CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter
// bit 5 = hltEle*Ele*CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter
// bit 6 = hltMu*TrkIsoVVL*Ele*CaloIdLTrackIdLIsoVL*Filter*
// bit 7 = hlt*OverlapFilterIsoEle*PFTau*
// bit 8 = hltEle*Ele*Ele*CaloIdLTrackIdLDphiLeg*Filter
// bit 9 = hltL3fL1Mu*DoubleEG*Filtered* || hltMu*DiEle*CaloIdLTrackIdLElectronleg*Filter
// bit 10 = hltL3fL1DoubleMu*EG*Filter* || hltDiMu*Ele*CaloIdLTrackIdLElectronleg*Filter
// bit 11 = hltEle32L1DoubleEGWPTightGsfTrackIsoFilter && hltEGL1SingleEGOrFilter
// bit 12 = hltEle*CaloIdVTGsfTrkIdTGsfDphiFilter
// bit 13 = HLT_Ele*PFJet*
// bit 14 = hltEG175HEFilter || hltEG200HEFilter
// bit 15 = hltEle*CaloIdLMWPMS2Filter
// bit 16 = hltDiEle*CaloIdLMWPMS2UnseededFilter

// Muon filter bits
// Run2
// bit 0 = *RelTrkIsoVVLFiltered0p4
// bit 1 = hltL3crIso*Filtered0p07
// bit 2 = *OverlapFilterIsoMu*PFTau*
// bit 3 = hltL3crIsoL1*SingleMu*Filtered0p07 || hltL3crIsoL1sMu*Filtered0p07
// bit 4 = hltDiMuon*Filtered*
// bit 5 = hltMu*TrkIsoVVL*Ele*CaloIdLTrackIdLIsoVL*Filter*
// bit 6 = hltOverlapFilterIsoMu*PFTau*
// bit 7 = hltL3fL1TripleMu*
// bit 8 = hltL3fL1DoubleMu*EG*Filtered* || hltDiMu*Ele*CaloIdLTrackIdLElectronleg*Filter
// bit 9 = hltL3fL1Mu*DoubleEG*Filtered* || hltMu*DiEle*CaloIdLTrackIdLElectronleg*Filter
// bit 10 = hltL3fL1sMu*L3Filtered50* || hltL3fL1sMu*TkFiltered50*
// bit 11 = hltL3fL1sMu*L3Filtered100* || hltL3fL1sMu*TkFiltered100*

// Run3
// bit 0 = *RelTrkIsoVVLFiltered0p4
// bit 1 = hltL3crIso*Filtered0p07 || hltL3crIso*IsoFiltered0p08 || hltL3crIso*IsoFiltered
// bit 2 = *OverlapFilterIsoMu*PFTau*
// bit 3 = hltL3crIsoL1*SingleMu*IsoFiltered0p07 || hltL3crIsoL1sMu*IsoFiltered0p07 || hltL3crIsoL1*SingleMu*IsoFiltered0p08 || hltL3crIsoL1sMu*IsoFiltered0p08 || hltL3crIsoL1*SingleMu*IsoFiltered || hltL3crIsoL1sMu*IsoFiltered
// bit 4 = hltDiMuon*Filtered*
// bit 5 = hltMu*TrkIsoVVL*Ele*CaloIdLTrackIdLIsoVL*Filter*
// bit 6 = hltOverlapFilterIsoMu*PFTau*
// bit 7 = hltL3fL1TripleMu*
// bit 8 = hltL3fL1DoubleMu*EG*Filtered* || hltDiMu*Ele*CaloIdLTrackIdLElectronleg*Filter
// bit 9 = hltL3fL1Mu*DoubleEG*Filtered* || hltMu*DiEle*CaloIdLTrackIdLElectronleg*Filter
// bit 10 = hltL3fL1sMu*L3Filtered50* || hltL3fL1sMu*TkFiltered50*
// bit 11 = hltL3fL1sMu*L3Filtered100* || hltL3fL1sMu*TkFiltered100*
// bit 12 = hltMu17Photon30IsoCaloIdMuonlegL3Filtered17Q 
