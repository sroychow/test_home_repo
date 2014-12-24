#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

#include "MuonEfficiency.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::abs;
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
MuonEfficiency::MuonEfficiency()
  : AnaBase(),
    _dumpEvent(false)
{
  _triggerPathTagList.clear();  
}
// ----------
// Destructor
// ----------
MuonEfficiency::~MuonEfficiency() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MuonEfficiency::beginJob() 
{ 
  AnaBase::beginJob();

  for (unsigned int i = 0; i < NEL(nProbe); i++){
    nProbe[i] = 0;
  }
  for (unsigned int i = 0; i < NEL(nSingleCut); i++){
    nSingleCut[i] = 0;
  }
  // Open the output ROOT file
  _histf->cd();
  bookHistograms();

  return true;
}
// ---------------
// Book histograms
// ---------------
void MuonEfficiency::bookHistograms() 
{
  new TH1D("evcounter", "Selected event counter", 9, -0.5, 8.5);

  //// FILLING HISTOGRAMS FOR FINDING OUT EFFICIENCY for abs(eta) < 1.44
  new TH1F("mass_ptle10_pass_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_pass_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_pass_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_pass_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_pass_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_pass_centralEta", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_ptle10_fail_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_fail_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_fail_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_fail_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_fail_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_fail_centralEta", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_ptle10_passIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_passIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_passIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_passIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_passIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_passIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_ptle10_failIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_failIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_failIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_failIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_failIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_failIso_centralEta", "Invariant mass of the two muons", 120, 20, 140);
  //// FILLING HISTOGRAMS FOR FINDING OUT EFFICIENCY for 1.44 < abs(eta) < 2.1
  //=========================================================================================
  new TH1F("mass_ptle10_pass_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_pass_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_pass_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_pass_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_pass_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_pass_endEta", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_ptle10_fail_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_fail_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_fail_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_fail_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_fail_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_fail_endEta", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_ptle10_passIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_passIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_passIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_passIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_passIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_passIso_endEta", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_ptle10_failIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_failIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_failIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_failIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_failIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_failIso_endEta", "Invariant mass of the two muons", 120, 20, 140);
  //========================================================================================
  new TH1F("mass_ptle10_pass", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_pass", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_pass", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_pass", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_pass", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_pass", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_ptle10_fail", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_fail", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_fail", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_fail", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_fail", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_fail", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_ptle10_passIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_passIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_passIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_passIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_passIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_passIso", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_ptle10_failIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle20_failIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle30_failIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle40_failIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle50_failIso", "Invariant mass of the two muons", 120, 20, 140);
  new TH1F("mass_ptle100_failIso", "Invariant mass of the two muons", 120, 20, 140);

  new TH1F("mass_outwindow", "Invariant mass of the two muons outside of the fixed window", 120, 0, 120);
  new TH1F("mass_inwindow", "Invariant mass of the two muons inside of the fixed window", 120, 0, 120);
  new TH1F("globalchi2_inwindow", "ch2 inside the fixed window", 90, 0, 30);
  new TH1F("globalchi2_all", "ch2 for all probe muon", 90, 0, 30);
  new TH1F("met_", "met befor puting total charge #neq 0 cut on tag muon", 100, 0, 100);
  new TH2F("TnP_correlation","pt correlation between tag and probe moun within acceptance cut", 100, 0, 100, 100, 0, 100);
  new TProfile("probe_trkD0_pt_before","profile plot of trkD0 as a func of pt just after the acceptance cut", 100, 0, 100, 0, 0.1);
  new TProfile("probe_dB_pt_before","profile plot of dB as a func of pt just after the acceptance cut", 100, 0, 100, 0, 0.1);
  new TProfile("probe_vtxDistZ_pt_before","profile plot of vtxDistZ as a func of pt just after the acceptance cut", 100, 0, 100, 0, 0.8);
  new TProfile("pt_pfRelIso", "pfRelIso of probe muon after ID and acceptance cuts as a function of pt", 100, 0, 100, 0,0.8);
  new TProfile("pt_pfRelIso_afterAllCuts", "pfRelIso of probe muon after ID, acceptance and Iso cuts as a function of pt", 100, 0, 100, 0,0.8);
  new TProfile("pt_pfRelIso_noIsoCut", "pfRelIso of probe muon after ID and acceptance cuts as a function of pt", 100, 0, 100, 0,0.8);
  new TProfile("eta_pfRelIso_noIsoCut", "pfRelIso of probe muon after ID and acceptance cuts as a function of pt", 50, -2.5, 2.5, 0,0.8);
  new TProfile("eta_pfRelIso_afterAllCuts", "pfRelIso of probe muon after ID, acceptance and Iso cuts as a function of pt", 50, -2.5, 2.5, 0, 0.8);

  new TH1F("pfRelIso","distribution of pfRelIso after id cuts",1000, 0, 10);
  new TH1F("mueta_probe_before","eta distribution of mu before id cuts",100,-3,3);
  new TH1F("mueta_probe_after","eta distribution of mu after id cuts",100,-3,3);

  new TH1F("tagpt","pt distribution of tag muon within acceptance cut) ", 100, 0, 100);

  new TH1F("mupt_probe_beforeEtaCut","pt distribution of mu before acceptance cut (ID cuts are applied)", 100, 0, 100);
  new TH1F("mupt_probe_before","pt distribution of mu before cuts (only #eta cut applied)", 100, 0, 100);
  new TH1F("mupt_probe_before_low","pt distribution of mu for #eta < 1.4 before cuts ", 100, 0, 100);
  new TH1F("mupt_probe_before_high","pt distribution of mu for 1.4 <= #eta < 2.1 before cuts", 100, 0, 100);
  new TH1F("mupt_probe_after","pt distribution of mu after cuts", 100, 0, 100);
  new TH1F("mupt_probe_after_low","pt distribution of mu for #eta < 1.4 after cuts", 100, 0, 100);
  new TH1F("mupt_probe_after_high","pt distribution of mu for 1.4 <= #eta < 2.1 after cuts", 100, 0, 100);

  new TH1F("mupt_probe_after_AllCuts","pt distribution of probe mu after all ID, acceptance and pfRelIso cuts", 100, 0, 100);
  new TH1F("mupt_probe_after_AllCuts_lowEta","pt distribution of probe mu after all ID, pfRelIso cuts for #eta < 1.4", 100, 0, 100);
  new TH1F("mupt_probe_after_AllCuts_highEta","pt distribution of probe mu after all Id and pfRelIso cuts for 1.4<= #eta < 2.1", 100, 0, 100);

  new TH1F("isTrackerMuon","pt distribution of probe muon after isTrackerMuon cut", 100, 0, 100);
  new TH1F("isGlobalMuonPromptTight","pt of probe muon after isGlobalMuonPromptTight cut", 100, 0, 100);
  new TH1F("isAllArbitrated","pt of probe muon after of isAllArbitrated cut", 100, 0, 100);
  new TH1F("vtxDistZ","pt of probe muon after vtxDistZ cut", 100, 0, 100);
  new TH1F("nChambers","pt of probe muon after nChambers cut", 100, 0, 100);
  new TH1F("nMatches","pt of probe muon after nMatches cut", 100, 0, 100);
  new TH1F("nMatchedStations","pt of probe muon after nMatchedStations cut", 100, 0, 100);
  new TH1F("pixHits","pt distribution of probe muon after trkHits cut", 100, 0, 100);
  new TH1F("trkHits","pt distribution of probe muon after trkHits cut", 100, 0, 100);
  new TH1F("globalChi2","pt distribution of probe muon after globalChi2 cut", 100, 0, 100);
  new TH1F("trkD0","pt distribution of probe muon after trkD0 cut", 100, 0, 100);
  new TH1F("dB","pt distribution of probe muon after dB cut", 100, 0, 100);
  new TH1F("iso","pt distribution of probe muon afterrelIso cut", 100, 0, 100);
}

// -------------------
// The main event loop
// -------------------
void MuonEfficiency::clearLists() {
  vtxList.clear();
  eleList.clear();
  tauList.clear();
  bjetList.clear();
  trigObjList.clear();
  probeMuonList.clear();
  tagMuonList.clear();
}
void MuonEfficiency::eventLoop() 
{
  // Initialize analysis
  if (!beginJob()) return;

  int nPrint = max(10000, nEvents/1000);

  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  for (int ev = 0; ev < nEvents; ev++) {
    clearEvent();
    int lflag = _chain->LoadTree(ev); 
    int nentries = getEntry(lflag);    // returns bytes read
    
    string currentFile(gSystem->BaseName(_chain->GetCurrentFile()->GetName())); 

    const Event* evt = dynamic_cast<Event*>(eventA->At(0));
    assert(evt);
    int run   = evt->run;
    int event = evt->event;

    // PileUP weight
    _puevWt = 1; // for data
    if (_isMC) {
      int npu = 0;
      _puevWt = wtPileUp(npu);
    }
    // Show status of the run
    if (currentFile != lastFile) 
      cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
           << " ==> " << currentFile 
           << " <<< Run# " << run
           << " Event# " << setw(8) << event << " >>> " 
           << " Events proc. " << setw(8) << ev << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
       cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
            << " ==> " << _chain->GetCurrentFile()->GetName() 
            << " <<< Run# " << run 
            << " Event# " << setw(8) << event << " >>> " 
            << " Events proc. " << setw(8) << ev << endl;

    AnaUtil::fillHist1D("evcounter", 0, _puevWt);

    // Trigger selection
    if (_useTrigger && !isTriggered()) continue;
    AnaUtil::fillHist1D("evcounter", 1, _puevWt);

    clearLists();
    op.verbose = (_logOption >> 1 & 0x1); 
    findVtxInfo(vtxList, op, _fLog);
    // Event Selection Starts here .....
    // presence of > 1 good vertex
    if (vtxList.size() < 1) continue;
    AnaUtil::fillHist1D("evcounter", 2, _puevWt);

    double vz = vtxList.at(0).z;
    op.verbose = (_logOption >> 4 & 0x1); 
    findElectronInfo(eleList, vz, op, _fLog);

    op.verbose = (_logOption >> 2 & 0x1); 
    findTauInfo(tauList, vz, op, _fLog);

    op.verbose = (_logOption >> 5 & 0x1); 
    findJetInfo(bjetList, op, _fLog);

    findTriggerObjectInfo(trigObjList);

    selectEvent();
  }  
  // Analysis is over
  endJob();
}

// ----------------------------------------------------------
// Perform event selection, For selection of Z -> e+e- events
// we need,
//   - > 0 Tight electron
//   - event within e+e- invariant mass window
// ----------------------------------------------------------
void MuonEfficiency::selectEvent() 
{
  if (eleList.size()) return;
  AnaUtil::fillHist1D("evcounter", 3, _puevWt);
  if (tauList.size()) return;
  AnaUtil::fillHist1D("evcounter", 4, _puevWt);
  if (bjetList.size()) return;
  AnaUtil::fillHist1D("evcounter", 5, _puevWt);

  findMuonIDInfo();
}
void MuonEfficiency::findMuonIDInfo() {
  // Selection of tag muon
  if (n_muon < 2) return;
  AnaUtil::fillHist1D("evcounter", 6, _puevWt);
  if (_dumpEvent) {
    dumpTriggerPaths(_fLog);
    dumpEvent("11111", _fLog, true);
    dumpTriggerObjectInfo(trigObjList, _fLog);
  }

  size_t tagindex;
  vector<Muon> preTagMuonList;
  for (int indx = 0; indx < n_muon; ++indx) {
    const Muon* muon = dynamic_cast<Muon*>(muonA->At(indx));
    if (muon) preTagMuonList.push_back(*muon);
  }
  random_shuffle (preTagMuonList.begin(), preTagMuonList.end());

  double maxPtDiff = AnaUtil::cutValue(_evselCutMap, "maxPtDiff");
  double vz = vtxList.at(0).z; // assuming events are selected using good vtx condition

  int indx = 0;
  for (vector<Muon>::const_iterator it  = preTagMuonList.begin(); 
                                    it != preTagMuonList.end(); ++it,++indx) {
    const Muon& muon = (*it);
    if (abs(muon.eta) >= AnaUtil::cutValue(_muonCutMap, "eta")) continue;

    int sbit = 0;
    if (!muon.isTrackerMuon)                                                    sbit |= (1 << 0);
    if (!muon.isGlobalMuonPromptTight)                                          sbit |= (1 << 1);
    if (!muon.isAllArbitrated)                                                  sbit |= (1 << 2);
    if (abs(muon.vtxDistZ) >= AnaUtil::cutValue(_muonCutMap, "vtxDistZ"))       sbit |= (1 << 3);
    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(_muonCutMap, "dz"))          sbit |= (1 << 4);
    if (abs(muon.nChambers) <= AnaUtil::cutValue(_muonCutMap, "nChambers"))     sbit |= (1 << 5);
    if (abs(muon.nMatches) <= AnaUtil::cutValue(_muonCutMap, "nMatches"))       sbit |= (1 << 6);
    if (abs(muon.nMatchedStations) <= AnaUtil::cutValue(_muonCutMap, "nMatchedStations"))     
                                                                                sbit |= (1 << 7);
    if (muon.pfRelIso >= AnaUtil::cutValue(_muonCutMap, "relIso"))              sbit |= (1 << 8);
    if (muon.pixHits <= AnaUtil::cutValue(_muonCutMap,"pixHits"))               sbit |= (1 << 9);
    if (muon.trkHits <= AnaUtil::cutValue(_muonCutMap,"trkHits"))               sbit |= (1 << 10);
    if (muon.globalChi2 >= AnaUtil::cutValue(_muonCutMap,"globalChi2"))         sbit |= (1 << 11);
    if (abs(muon.trkD0) >= AnaUtil::cutValue(_muonCutMap,"trkD0"))              sbit |= (1 << 12);
    if (abs(muon.dB) >= AnaUtil::cutValue(_muonCutMap,"dB"))                    sbit |= (1 << 13);

    // Now apply cuts
    if (sbit) continue;  

    // Trigger match
    TLorentzVector taglv;
    taglv.SetPtEtaPhiE(muon.pt, muon.eta, muon.phi, muon.energy);
    int ntobj = trigObjList.size();
    int tindx = -1;
    uint flag = 0;
    double drTag = matchTriggerObject(trigObjList, taglv, _triggerPathTagList, -1, maxPtDiff, tindx, flag);
    if (0) cout << "=> indx = " << indx
                << " tindx = " << tindx
                << " drTag = " << drTag
                << " flag = " << flag
                << endl;
    if (tindx < 0 || tindx >= ntobj) continue;
    if (drTag >= AnaUtil::cutValue(_evselCutMap, "maxDr") || flag != 1) continue;

    tagindex = indx;
    tagMuonList.push_back(muon);
    break;
  }
  if (tagMuonList.size() < 1) return;
  AnaUtil::fillHist1D("evcounter", 7, _puevWt);
  const Muon& tagmuon = tagMuonList.at(0);
    
  const MET* Met = dynamic_cast<MET*>(metA->At(0));
  if (!Met) return;
  AnaUtil::fillHist1D("met_", Met->met, _puevWt);
  if (Met->met > AnaUtil::cutValue(_evselCutMap, "maxMET")) return;
  AnaUtil::fillHist1D("evcounter", 8);
    
  TLorentzVector tag;
  tag.SetPtEtaPhiE(tagmuon.pt, tagmuon.eta, tagmuon.phi, tagmuon.energy);
  vector<Muon> probeAllList;
  size_t jndx = 0;
  for (vector<Muon>::const_iterator it  = preTagMuonList.begin(); 
                                    it != preTagMuonList.end(); ++it,++jndx) {
    if (jndx == tagindex) continue;
    const Muon& muon = (*it);
    TLorentzVector probe;
    probe.SetPtEtaPhiE(muon.pt, muon.eta, muon.phi, muon.energy); 
    TLorentzVector z = tag + probe;
    double mass = z.M();
    probeAllList.push_back(muon);
    if (tagmuon.charge + muon.charge != 0) continue;
    AnaUtil::fillHist1D("mass_outwindow", mass, _puevWt);
    if (mass > AnaUtil::cutValue(_evselCutMap, "massLow") && mass < AnaUtil::cutValue(_evselCutMap, "massHigh")) {
      AnaUtil::fillHist1D("tagpt", tagmuon.pt, _puevWt);
      if (tagmuon.pt > AnaUtil::cutValue(_evselCutMap, "tagPt")) {
        probeMuonList.push_back(muon);
	AnaUtil::fillHist1D("mass_inwindow", mass, _puevWt);
      }
    }
  }
  for (vector<Muon>::const_iterator it  = probeAllList.begin(); 
                                    it != probeAllList.end(); ++it) {
    const Muon& muon = (*it);
    AnaUtil::fillHist1D("globalchi2_all", muon.globalChi2, _puevWt);
  }
  if (probeMuonList.size() < 1) return;
  AnaUtil::fillHist1D("evcounter", 9);
  sort(probeMuonList.begin(), probeMuonList.end(), PtComparator<Muon>());

  double etaCut = AnaUtil::cutValue(_muonCutMap, "eta");
  if (abs(probeMuonList.at(0).eta) < etaCut && abs(tagmuon.eta) < etaCut && tagmuon.pt > 15)
    AnaUtil::fillHist2D("TnP_correlation", tagmuon.pt, probeMuonList.at(0).pt, _puevWt);

  computeMuonEff();
}
void MuonEfficiency::computeMuonEff() {
  TLorentzVector t, p;

  const Muon& tmuon = tagMuonList.at(0);
  t.SetPtEtaPhiE(tmuon.pt, tmuon.eta, tmuon.phi, tmuon.energy);

  const Muon& pmuon = probeMuonList.at(0);
  p.SetPtEtaPhiE(pmuon.pt, pmuon.eta, pmuon.phi, pmuon.energy);

  double prbeta = pmuon.eta;
  double prbpt  = pmuon.pt;

  double mass = (t+p).M();
  if (abs(prbeta) < 1.44)
    AnaUtil::fillHist1D("mupt_probe_before_low", prbpt, _puevWt);
  else if (1.44 <= abs(prbeta) && abs(prbeta) < 2.1)
    AnaUtil::fillHist1D("mupt_probe_before_high", prbpt, _puevWt);
   
  const string HistAfter[] = { // contains name of some histograms, used to calculate 
                               // the efficiency for each individual cuts
    "isTrackerMuon",
    "isGlobalMuonPromptTight",
    "isAllArbitrated",
    "vtxDistZ",
    "nChambers",
    "nMatches",
    "nMatchedStations",
    "pixHits",
    "trkHits",
    "globalChi2",
    "trkD0",
    "dB",
    "iso"
  }; 

  double vz = vtxList.at(0).z; // assuming events are selected using good vtx condition

  // Selection of probe muon
  int sbit = 0;
  if (!pmuon.isTrackerMuon)                                                    sbit |= (1 << 0);
  if (!pmuon.isGlobalMuonPromptTight)                                          sbit |= (1 << 1);
  if (!pmuon.isAllArbitrated)                                                  sbit |= (1 << 2);
  if (abs(pmuon.vtxDistZ) >= AnaUtil::cutValue(_muonCutMap, "vtxDistZ"))       sbit |= (1 << 3);
  bool isGoodVtx;
  TVector3 vmu = findLeptonVtx(pmuon.vtxIndex, isGoodVtx);
  double dz = vmu.z() - vz;
  if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(_muonCutMap, "dz"))           sbit |= (1 << 4);
  if (abs(pmuon.nChambers) <= AnaUtil::cutValue(_muonCutMap, "nChambers"))     sbit |= (1 << 5);
  if (abs(pmuon.nMatches) <= AnaUtil::cutValue(_muonCutMap, "nMatches"))       sbit |= (1 << 6);
  if (abs(pmuon.nMatchedStations) <= AnaUtil::cutValue(_muonCutMap, "nMatchedStations"))
                                                                               sbit |= (1 << 7);
  if (pmuon.pixHits <= AnaUtil::cutValue(_muonCutMap,"pixHits"))               sbit |= (1 << 8);
  if (pmuon.trkHits <= AnaUtil::cutValue(_muonCutMap,"trkHits"))               sbit |= (1 << 9);
  if (pmuon.globalChi2 >= AnaUtil::cutValue(_muonCutMap,"globalChi2"))         sbit |= (1 << 10);
  if (abs(pmuon.trkD0) >= AnaUtil::cutValue(_muonCutMap,"trkD0"))              sbit |= (1 << 11);
  if (abs(pmuon.dB) >= AnaUtil::cutValue(_muonCutMap, "dB"))                   sbit |= (1 << 12);
    
  if (pmuon.globalChi2 <= AnaUtil::cutValue(_muonCutMap, "globalChi2"))
    AnaUtil::fillHist1D("globalchi2_inwindow", pmuon.globalChi2, _puevWt);
    
  AnaUtil::fillHist1D("mueta_probe_before", prbeta, _puevWt);
    
  if (sbit == 0) {
    AnaUtil::fillHist1D("mueta_probe_after", prbeta, _puevWt);
    AnaUtil::fillHist1D("mupt_probe_beforeEtaCut", prbpt, _puevWt);
    if (abs(prbeta) < 1.44)
      AnaUtil::fillHist1D("mupt_probe_after_low", prbpt, _puevWt);
    else if (1.44 <= abs(prbeta) && abs(prbeta) < 2.1)
      AnaUtil::fillHist1D("mupt_probe_after_high", prbpt, _puevWt);
  }
  // Filling Histograms for finding out efficiency
  if (sbit == 0) {
    // Allowed pseudorapidity region
    if (8 < prbpt && prbpt <= 10)        AnaUtil::fillHist1D("mass_ptle10_pass", mass, _puevWt);
    else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_pass", mass, _puevWt);
    else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_pass", mass, _puevWt);
    else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_pass", mass, _puevWt);
    else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_pass", mass, _puevWt);
    else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_pass", mass, _puevWt);
    if (pmuon.pfRelIso < AnaUtil::cutValue(_muonCutMap, "relIso")) {
      if (8 < prbpt && prbpt <= 10)        AnaUtil::fillHist1D("mass_ptle10_passIso", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_passIso", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_passIso", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_passIso", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_passIso", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_passIso", mass, _puevWt);
    }
    else {
      if (8 < prbpt && prbpt <= 10)        AnaUtil::fillHist1D("mass_ptle10_fail", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_failIso", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_failIso", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_failIso", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_failIso", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_failIso", mass, _puevWt); 
    }
    // Central region
    if (abs(prbeta) < 1.44) {
      if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_pass_centralEta", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_pass_centralEta", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_pass_centralEta", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_pass_centralEta", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_pass_centralEta", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_pass_centralEta", mass, _puevWt);
      // study isolation only when sbit = 0
      if (pmuon.pfRelIso < AnaUtil::cutValue(_muonCutMap, "relIso")) {
        if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_passIso_centralEta", mass, _puevWt);
        else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_passIso_centralEta", mass, _puevWt);
        else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_passIso_centralEta", mass, _puevWt);
        else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_passIso_centralEta", mass, _puevWt);
        else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_passIso_centralEta", mass, _puevWt);
        else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_passIso_centralEta", mass, _puevWt);
      }
      else {
        if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_failIso_centralEta", mass, _puevWt);
        else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_failIso_centralEta", mass, _puevWt);
        else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_failIso_centralEta", mass, _puevWt);
        else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_failIso_centralEta", mass, _puevWt);
        else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_failIso_centralEta", mass, _puevWt);
        else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_failIso_centralEta", mass, _puevWt); 
      }
    }
    // Forward region
    else if (1.44 <= abs(prbeta) && abs(prbeta) < 2.1) {
      if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_pass_endEta", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_pass_endEta", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_pass_endEta", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_pass_endEta", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_pass_endEta", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_pass_endEta", mass, _puevWt);
     
      if (pmuon.pfRelIso < AnaUtil::cutValue(_muonCutMap, "relIso")) {
	if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_passIso_endEta", mass, _puevWt);
	else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_passIso_endEta", mass, _puevWt);
        else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_passIso_endEta", mass, _puevWt);
	else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_passIso_endEta", mass, _puevWt);
	else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_passIso_endEta", mass, _puevWt);
	else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_passIso_endEta", mass, _puevWt);
      }
      else {
	if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_fail_endEta", mass, _puevWt);
	else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_failIso_endEta", mass, _puevWt);
	else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_failIso_endEta", mass, _puevWt);
	else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_failIso_endEta", mass, _puevWt);
	else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_failIso_endEta", mass, _puevWt);
	else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_failIso_endEta", mass, _puevWt); 
      }
    }
  }
  else {
    // no need to study Isolation for sbit = 0
    // Allowed pseudorapidity region
    if (8 < prbpt && prbpt <= 10)        AnaUtil::fillHist1D("mass_ptle10_fail", mass, _puevWt);
    else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_fail", mass, _puevWt);
    else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_fail", mass, _puevWt);
    else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_fail", mass, _puevWt);
    else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_fail", mass, _puevWt);
    else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_fail", mass, _puevWt);
    // Central region
    if (abs(prbeta) < 1.44) {
      if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_fail_centralEta", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_fail_centralEta", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_fail_centralEta", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_fail_centralEta", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_fail_centralEta", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_fail_centralEta", mass, _puevWt);
    }
    // Forward region
    else if (1.44 <= abs(prbeta) && abs(prbeta) < 2.1) {
      if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_fail", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_fail_endEta", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_fail_endEta", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_fail_endEta", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_fail_endEta", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_fail_endEta", mass, _puevWt);
    }
  } 
  // =====================================================================================================
  if (abs(prbeta) >= AnaUtil::cutValue(_muonCutMap, "eta")) return;                //acceptance cut
  ++nProbe[0];
  ++nSingleCut[0];
  AnaUtil::fillHist1D("mupt_probe_before", prbpt, _puevWt);
  AnaUtil::fillProfile("probe_trkD0_pt_before", prbpt, pmuon.trkD0, _puevWt);
  AnaUtil::fillProfile("probe_dB_pt_before", prbpt, pmuon.dB, _puevWt);
  AnaUtil::fillProfile("probe_vtxDistZ_pt_before", prbpt, pmuon.vtxDistZ, _puevWt);

  if (sbit == 0) {
    AnaUtil::fillHist1D("mupt_probe_after",prbpt,_puevWt);
    AnaUtil::fillProfile("pt_pfRelIso_noIsoCut", prbpt, pmuon.pfRelIso, _puevWt);
    AnaUtil::fillProfile("eta_pfRelIso_noIsoCut", prbeta, pmuon.pfRelIso, _puevWt);
    if (pmuon.pfRelIso <= AnaUtil::cutValue(_muonCutMap, "relIso")) {
      AnaUtil::fillHist1D("mupt_probe_after_AllCuts", prbpt, _puevWt);
      AnaUtil::fillProfile("pt_pfRelIso_afterAllCuts", prbpt, pmuon.pfRelIso, _puevWt);
      AnaUtil::fillProfile("eta_pfRelIso_afterAllCuts", prbeta, pmuon.pfRelIso, _puevWt);
      if (abs(prbeta) < 1.44)
  	AnaUtil::fillHist1D("mupt_probe_after_AllCuts_lowEta", prbpt, _puevWt);
      else if (1.44 <= abs(prbeta) && abs(prbeta) < 2.1)
  	AnaUtil::fillHist1D("mupt_probe_after_AllCuts_highEta", prbpt, _puevWt);        
    }
  }
  AnaUtil::fillProfile("pt_pfRelIso", prbpt, pmuon.pfRelIso, _puevWt);
  AnaUtil::fillHist1D("pfRelIso", pmuon.pfRelIso, _puevWt);

  // Isolation Cut 
  if (pmuon.pfRelIso >= AnaUtil::cutValue(_muonCutMap, "relIso")) sbit |= (1 << 13);
  for (int i = 1; i < 15; ++i) {
    int j = pow(2,i) - 1;
    if ((sbit & j) == 0) 
      ++nProbe[i];
    if ((sbit & (1<< (i-1))) == 0) {
      ++nSingleCut[i];
      if (i > 1) AnaUtil::fillHist1D(HistAfter[i-2], prbpt, _puevWt);  
    }
  }  
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void MuonEfficiency::endJob() {  
  // Print Muon ID Efficiency Numbers
  const string tagnprobe[] =
  { 
    "abs(muon->eta) < 2.1 ",
    "probe muon with no quality cut",
    "isTrackerMuon",                       
    "isGlobalMuonPromptTight",
    "isAllArbitrated",
    "abs(vtxDistZ) <= 0.2 ",
    "nChambers > 0 ",
    "nMatches >0 ",                   
    "nMatchedStations > 1 ",
    "muon->pixHits >= 1 ",
    "muon->trkHits >= 10 ",
    "globalChi2 <= 10 ",
    "abs(muon->trkD0) <= 0.02 ",
    "abs(muon->dB) <= 0.02 ",
    "relIso > 0.1 "                     
  };
  _fLog << "Statistics for Muon Identification Efficiency" << endl;
  _fLog << setw(89) << "eff_total ->" 
        << setw(10) << "err"
        << setw(20) << "Single Cut" 
        << setw(20) << "singleCuteff ->"
        << setw(10) << "err"
        << endl;

  for (unsigned int i = 0; i < NEL(nProbe); i++) {
    if (!nProbe[1]) break; 
    float eff_total =  nProbe[i]*1.0/nProbe[1];
    float err1 = sqrt(nProbe[i]*(1-eff_total))/nProbe[1]; 
    if (i == 0)
      _fLog << setw(64) << tagnprobe[i] 
            << setw(10) << nProbe[i]
            << setw(15) << "-" 
            << setw(12) << "-"
            << setw(20) << nSingleCut[0]
            << setw(20) << "-"
            << setw(12) << "-"
            << endl;
    else if (i == 1)
      _fLog << setw(64) << tagnprobe[i] 
            << setw(10) << nProbe[i]
            << setw(15) << "-"
            << setw(12) << "-"
            << setw(20) << nSingleCut[0] 
            << setw(20) << "-"
            << setw(12) << "-"
            << endl;
    else {
      float eff_single = nSingleCut[i]*1.0/nSingleCut[1];
      float err2 = sqrt(nSingleCut[i]*(1-eff_single))/nSingleCut[1]; 
      _fLog << setw(64) << tagnprobe[i] 
            << setw(10) << nProbe[i]
            << setw(15) << eff_total
            << setw(12) << err1 
            << setw(20) << nSingleCut[i] 
            << setw(20) << eff_single 
            << setw(12) << err2 
            << endl;
    }
  }
  _fLog << resetiosflags(ios::fixed);

  closeFiles();

  _histf->cd();
  _histf->Write();
  _histf->Close();
  delete _histf;
}
bool MuonEfficiency::readJob(const string& jobFile, int& nFiles)
{
  AnaBase::readJob(jobFile, nFiles);

  static const int BUF_SIZE = 256;

  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }

  // note that you must use a pointer (reference!) to the cut map
  // in order to avoid scope related issues
  map<string, map<string, double>* > hmap;
  hmap.insert(pair<string, map<string, double>* >("evselCutList", &_evselCutMap));

  char buf[BUF_SIZE];
  vector<string> tokens;
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   

    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    AnaUtil::tokenize(line, tokens);
    int vsize = tokens.size();
    assert(vsize > 1);
    string key = tokens.at(0);

    if (key == "tagTrigPathList")
      AnaUtil::buildList(tokens, _triggerPathTagList);
    else if (key == "dumpEvent")
      _dumpEvent = atoi(tokens[1].c_str()) > 0 ? true : false;
    else if (key == "evselCutList")
      AnaUtil::storeCuts(tokens, hmap);

    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();

  return true;
}
void MuonEfficiency::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));

  // Trigger Path List
  if (_useTrigger) 
    AnaUtil::showList(_triggerPathTagList, ">>> INFO. Trigger Paths for Tag Leg:", os);

  AnaUtil::showCuts(hmap, os);
}
