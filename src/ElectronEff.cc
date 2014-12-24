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

#include "ElectronEff.h"
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
ElectronEff::ElectronEff()
  : AnaBase(),
    _dumpEvent(false)
{
  _triggerPathTagList.clear();  
}
// ----------
// Destructor
// ----------
ElectronEff::~ElectronEff() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool ElectronEff::beginJob() 
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
void ElectronEff::bookHistograms() 
{
  new TH1D("evcounter", "Selected event counter", 9, -0.5, 8.5);

  //// FILLING HISTOGRAMS FOR FINDING OUT EFFICIENCY for abs(eta) < 1.44
  new TH1F("mass_ptle10_pass_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle20_pass_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle30_pass_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle40_pass_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle50_pass_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle100_pass_centralEta", "Invariant mass of the two electron", 120, 20, 140);

  new TH1F("mass_ptle10_fail_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle20_fail_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle30_fail_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle40_fail_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle50_fail_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle100_fail_centralEta", "Invariant mass of the two electron", 120, 20, 140);

  new TH1F("mass_ptle10_passIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle20_passIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle30_passIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle40_passIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle50_passIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle100_passIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);

  new TH1F("mass_ptle10_failIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle20_failIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle30_failIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle40_failIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle50_failIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle100_failIso_centralEta", "Invariant mass of the two electron", 120, 20, 140);
  //// FILLING HISTOGRAMS FOR FINDING OUT EFFICIENCY for 1.44 < abs(eta) < 2.1
  //=========================================================================================
  new TH1F("mass_ptle10_pass_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle20_pass_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle30_pass_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle40_pass_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle50_pass_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle100_pass_endEta", "Invariant mass of the two electron", 120, 20, 140);

  new TH1F("mass_ptle10_fail_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle20_fail_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle30_fail_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle40_fail_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle50_fail_endEta", "Invariant mass of the two electron", 120, 20, 140);
  new TH1F("mass_ptle100_fail_endEta", "Invariant mass of the two electron", 120, 20, 140);

  new TH1F("mass_ptle10_passIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle20_passIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle30_passIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle40_passIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle50_passIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle100_passIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);

  new TH1F("mass_ptle10_failIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle20_failIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle30_failIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle40_failIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle50_failIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle100_failIso_endEta", "Invariant mass of the two electrons", 120, 20, 140);
  //========================================================================================
  new TH1F("mass_ptle10_pass", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle20_pass", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle30_pass", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle40_pass", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle50_pass", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle100_pass", "Invariant mass of the two electrons", 120, 20, 140);

  new TH1F("mass_ptle10_fail", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle20_fail", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle30_fail", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle40_fail", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle50_fail", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle100_fail", "Invariant mass of the two electrons", 120, 20, 140);

  new TH1F("mass_ptle10_passIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle20_passIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle30_passIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle40_passIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle50_passIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle100_passIso", "Invariant mass of the two electrons", 120, 20, 140);

  new TH1F("mass_ptle10_failIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle20_failIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle30_failIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle40_failIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle50_failIso", "Invariant mass of the two electrons", 120, 20, 140);
  new TH1F("mass_ptle100_failIso", "Invariant mass of the two electrons", 120, 20, 140);

  new TH1F("mass_outwindow", "Invariant mass of the two electrons outside of the fixed window", 120, 0, 120);
  new TH1F("mass_inwindow", "Invariant mass of the two electrons inside of the fixed window", 120, 0, 120);
  new TH1F("met_", "met befor puting total charge #neq 0 cut on tag electron", 100, 0, 100);
  new TH2F("TnP_correlation","pt correlation between tag and probe electron within acceptance cut", 100, 0, 100, 100, 0, 100);
  // new TProfile("probe_trkD0_pt_before","profile plot of trkD0 as a func of pt just after the acceptance cut", 100, 0, 100, 0, 0.1);
  new TProfile("probe_dB_pt_before","profile plot of dB as a func of pt just after the acceptance cut", 100, 0, 100, 0, 0.1);
  new TProfile("probe_vtxDistZ_pt_before","profile plot of vtxDistZ as a func of pt just after the acceptance cut", 100, 0, 100, 0, 0.8);
  new TProfile("pt_pfRelIso", "pfRelIso of probe electron after ID and acceptance cuts as a function of pt", 100, 0, 100, 0,0.8);
  new TProfile("pt_pfRelIso_afterAllCuts", "pfRelIso of probe electron after ID, acceptance and Iso cuts as a function of pt", 100, 0, 100, 0,0.8);
  new TProfile("pt_pfRelIso_noIsoCut", "pfRelIso of probe electron after ID and acceptance cuts as a function of pt", 100, 0, 100, 0,0.8);
  new TProfile("eta_pfRelIso_noIsoCut", "pfRelIso of probe electron after ID and acceptance cuts as a function of pt", 50, -2.5, 2.5, 0,0.8);
  new TProfile("eta_pfRelIso_afterAllCuts", "pfRelIso of probe electron after ID, acceptance and Iso cuts as a function of pt", 50, -2.5, 2.5, 0, 0.8);

  new TH1F("pfRelIsoDB04","distribution of RelIso after id cuts",1000, 0, 10);

  new TH1F("RelIsoDB04Barrel","distribution of RelIso after id cuts in Barrel",1000, 0, 10);
  new TH1F("pfRelIsoDB04EndCap","distribution of RelIso after id cuts in EndCap",1000, 0, 10);

  //new TH1F("pfRelIso","distribution of pfRelIso after id cuts",1000, 0, 10);
  new TH1F("probeEta_before","eta distribution of electron before id cuts",100,-3,3);
  new TH1F("probeEta_after","eta distribution of electron after id cuts",100,-3,3);

  new TH1F("tagpt","pt distribution of tag electron within acceptance cut) ", 100, 0, 100);

  new TH1F("probept_beforeEtaCut","pt distribution of electron before acceptance cut (ID cuts are applied)", 100, 0, 100);
  new TH1F("probept_before","pt distribution of electron before cuts (only #eta cut applied)", 100, 0, 100);
  new TH1F("probept_before_low","pt distribution of electron for #eta < 1.4 before cuts ", 100, 0, 100);
  new TH1F("probept_before_high","pt distribution of electron for 1.4 <= #eta < 2.1 before cuts", 100, 0, 100);
  new TH1F("probept_after","pt distribution of electron after cuts", 100, 0, 100);
  new TH1F("probept_after_low","pt distribution of electron for #eta < 1.4 after cuts", 100, 0, 100);
  new TH1F("probept_after_high","pt distribution of electron for 1.4 <= #eta < 2.1 after cuts", 100, 0, 100);

  new TH1F("probept_after_AllCuts","pt distribution of probe electron after all ID, acceptance and pfRelIso cuts", 100, 0, 100);
  new TH1F("probept_after_AllCuts_lowEta","pt distribution of probe electron after all ID, pfRelIso cuts for #eta < 1.4", 100, 0, 100);
  new TH1F("probept_after_AllCuts_highEta","pt distribution of probe electron after all Id and pfRelIso cuts for 1.4<= #eta < 2.1", 100, 0, 100);

  new TH1F("eEleId_central", "pt of probe electron after eEleId85cIso cut", 100, 0, 100);
  new TH1F("dz_central", "pt of probe electron after dz cut", 100, 0, 100);
  new TH1F("pixHits_central", "pt distribution of probe electron after trkHits cut", 100, 0, 100);
  new TH1F("trkHits_central", "pt distribution of probe electron after trkHits cut", 100, 0, 100);
  new TH1F("trkD0_central", "pt distribution of probe electron after trkD0 cut", 100, 0, 100);
  new TH1F("Iso_central", "pt distribution of probe electron afterrelIso cut", 100, 0, 100);

  new TH1F("eEleId_end", "pt of probe electron after eEleId85cIso cut", 100, 0, 100);
  new TH1F("dz_end", "pt of probe electron after dz cut", 100, 0, 100);
  new TH1F("pixHits_end", "pt distribution of probe electron after trkHits cut", 100, 0, 100);
  new TH1F("trkHits_end", "pt distribution of probe electron after trkHits cut", 100, 0, 100);
  new TH1F("trkD0_end", "pt distribution of probe electron after trkD0 cut", 100, 0, 100);
  new TH1F("Iso_end", "pt distribution of probe electron afterrelIso cut", 100, 0, 100);

  new TH1F("trkD0_distribution","istribution of trkD0 for tag electron", 40, -10, 10);
}

// -------------------
// The main event loop
// -------------------
void ElectronEff::clearLists() {
  vtxList.clear();
  muonList.clear();
  tauList.clear();
  bjetList.clear();
  trigObjList.clear();
  probeEleList.clear();
  tagEleList.clear();
}
void ElectronEff::eventLoop() 
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
    dumpTriggerPaths(_fLog, false);
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
    findMuonInfo(muonList, vz, op, _fLog);

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
void ElectronEff::selectEvent() 
{
  if (muonList.size()) return;
  AnaUtil::fillHist1D("evcounter", 3, _puevWt);
  if (tauList.size()) return;
  AnaUtil::fillHist1D("evcounter", 4, _puevWt);
  if (bjetList.size()) return;
  AnaUtil::fillHist1D("evcounter", 5, _puevWt);

  findEleIDInfo();
}
void ElectronEff::findEleIDInfo() {
  // Selection of tag muon
  if (n_electron < 2) return;
  AnaUtil::fillHist1D("evcounter", 6, _puevWt);
  if (_dumpEvent) {
    dumpTriggerPaths(_fLog);
    dumpEvent("11111", _fLog, true);
    dumpTriggerObjectInfo(trigObjList, _fLog);
  }

  size_t tagindex;
  vector<Electron> preTagEleList;
  for (int indx = 0; indx < n_electron; ++indx) {
    const Electron* ele = dynamic_cast<Electron*>(electronA->At(indx));
    if (ele) preTagEleList.push_back(*ele);
  }
  random_shuffle (preTagEleList.begin(), preTagEleList.end());

  double maxPtDiff = AnaUtil::cutValue(_evselCutMap, "maxPtDiff");
  double vz = vtxList.at(0).z; // assuming events are selected using good vtx condition

  int ntobj = trigObjList.size();
  int indx = 0;
  for (vector<Electron>::const_iterator it  = preTagEleList.begin(); 
                                        it != preTagEleList.end(); ++it,++indx) {
    const Electron& ele = (*it);
    if (abs(ele.eta) >= AnaUtil::cutValue(_electronCutMap, "eta")) continue;
    if (abs(ele.eta) >= 1.4442 && abs(ele.eta) <= 1.566) continue;
   
    double tgscEta = ele.scEta;
    double tgmva = ele.idMVA;
    bool mva1 = abs(tgscEta) <  1.0 && tgmva > 0.133;
    bool mva2 = abs(tgscEta) >=  1.0 && abs(tgscEta) < 1.5 && tgmva > 0.465;
    bool mva3 = abs(tgscEta) >=  1.5 && abs(tgscEta) < 2.5 && tgmva > 0.518;
    bool mvalowpt = ele.pt < 20 && (mva1 || mva2 || mva3);

    bool mva4 = abs(tgscEta) <  1.0 && tgmva > 0.942;
    bool mva5 = abs(tgscEta) >=  1.0 && abs(tgscEta) < 1.5 && tgmva > 0.947;
    bool mva6 = abs(tgscEta) >=  1.5 && abs(tgscEta) < 2.5 && tgmva > 0.878;
    bool mvahighpt = ele.pt > 20 && (mva4 || mva5 || mva6);

    bool ElectronId = mvalowpt || mvahighpt;

     
    int sbit = 0;
    if (!ElectronId)                                                              sbit |= (1<<0);
    //    if (ele.simpleEleId85cIso != AnaUtil::cutValue(_electronCutMap, "eleId"))     sbit |= (1 << 0);
    bool isGoodVtx;
    TVector3 vele = findLeptonVtx(ele.vtxIndex, isGoodVtx);
    double dz = vele.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(_electronCutMap, "dz"))         sbit |= (1 << 1);
    if (ele.pfRelIsoDB04v2 >= AnaUtil::cutValue(_electronCutMap, "relIso"))              sbit |= (1 << 2);
    if (ele.pixHits <= AnaUtil::cutValue(_electronCutMap,"pixHits"))               sbit |= (1 << 3);
    if (ele.trkHits <= AnaUtil::cutValue(_electronCutMap,"trkHits"))               sbit |= (1 << 4);
    if (abs(ele.trkD0) >= AnaUtil::cutValue(_electronCutMap,"trkD0"))              sbit |= (1 << 5);

    // +Now apply cuts
    if (sbit) continue;  

    // Trigger match
    TLorentzVector taglv;
    taglv.SetPtEtaPhiE(ele.pt, ele.eta, ele.phi, ele.energy);
    int tindx = -1;
    uint flag = 0;
    double drTag = matchTriggerObject(trigObjList, taglv, _triggerPathTagList, -1, maxPtDiff, tindx, flag);
    if(0)  cout << "=> indx = " << indx
                << " tindx = " << tindx
                << " drTag = " << drTag
                << " flag = " << flag
                << endl;
    if (tindx < 0 || tindx >= ntobj) continue;
    if (drTag >= AnaUtil::cutValue(_evselCutMap, "maxDr") || flag != 1) continue;

    tagindex = indx;
    tagEleList.push_back(ele);
    break;
  }
  if (tagEleList.size() < 1) return;
  AnaUtil::fillHist1D("evcounter", 7, _puevWt);
  const Electron& tagele = tagEleList.at(0);
    
  const MET* Met = dynamic_cast<MET*>(metA->At(0));
  if (!Met) return;
  AnaUtil::fillHist1D("met_", Met->met, _puevWt);
  if (Met->met > AnaUtil::cutValue(_evselCutMap, "maxMET")) return;
  AnaUtil::fillHist1D("evcounter", 8);
    
  TLorentzVector tag;
  tag.SetPtEtaPhiE(tagele.pt, tagele.eta, tagele.phi, tagele.energy);
  //  vector<Electron> probeAllList;
  size_t jndx = 0;
  for (vector<Electron>::const_iterator it  = preTagEleList.begin(); 
                                    it != preTagEleList.end(); ++it,++jndx) {
    if (jndx == tagindex) continue;
    const Electron& ele = (*it);
    TLorentzVector probe;
    probe.SetPtEtaPhiE(ele.pt, ele.eta, ele.phi, ele.energy); 
    TLorentzVector z = tag + probe;
    double mass = z.M();
    //probeAllList.push_back(ele);
    if (tagele.charge + ele.charge != 0) continue;
    AnaUtil::fillHist1D("mass_outwindow", mass, _puevWt);
    if (mass > AnaUtil::cutValue(_evselCutMap, "massLow") && mass < AnaUtil::cutValue(_evselCutMap, "massHigh")) {
      if (tagele.pt > AnaUtil::cutValue(_evselCutMap, "tagPt")) {
	AnaUtil::fillHist1D("tagpt", tagele.pt, _puevWt);
	AnaUtil::fillHist1D("trkD0_distribution", tagele.trkD0, _puevWt);
        probeEleList.push_back(ele);
	AnaUtil::fillHist1D("mass_inwindow", mass, _puevWt);
      }
    }
  }
  if (probeEleList.size() < 1) return;
  AnaUtil::fillHist1D("evcounter", 9);
  sort(probeEleList.begin(), probeEleList.end(), PtComparator<Electron>());

  double etaCut = AnaUtil::cutValue(_electronCutMap, "eta");
  if (abs(probeEleList.at(0).eta) < etaCut && abs(tagele.eta) < etaCut && tagele.pt > 15)
    AnaUtil::fillHist2D("TnP_correlation", tagele.pt, probeEleList.at(0).pt, _puevWt);

  computeEleEff();
}
void ElectronEff::computeEleEff() {
  TLorentzVector t, p;

  const Electron& tele = tagEleList.at(0);
  t.SetPtEtaPhiE(tele.pt, tele.eta, tele.phi, tele.energy);

  const Electron& pele = probeEleList.at(0);
  p.SetPtEtaPhiE(pele.pt, pele.eta, pele.phi, pele.energy);

  double prbeta = pele.eta;
  double prbpt  = pele.pt;

  double mass = (t+p).M();
  if (abs(prbeta) < 1.4442)
    AnaUtil::fillHist1D("probept_before_low", prbpt, _puevWt);
  else if (1.566 <= abs(prbeta) && abs(prbeta) < 2.5)
    AnaUtil::fillHist1D("probept_before_high", prbpt, _puevWt);

  const string HistAfterCentral[] = { // contains name of some histograms, used to calculate 
                               // the efficiency for each individual cuts
    "eEleId_central",
    "dz_central",
    "pixHits_central",
    "trkHits_central",
    "trkD0_central",
    "Iso_central"
  };

  const string HistAfterEnd[] = { // contains name of some histograms, used to calculate                                                        
    // the efficiency for each individual cuts                                                                           
    "eEleId_end",
    "dz_end",
    "pixHits_end",
    "trkHits_end",
    "trkD0_end",
    "Iso_end"
  };


  double vz = vtxList.at(0).z; // assuming events are selected using good vtx condition

  // Selection of probe electron 

  double prbscEta = pele.scEta;
  double prbmva = pele.idMVA;
  bool mva1 = abs(prbscEta) <  1.0 && prbmva > 0.133;
  bool mva2 = abs(prbscEta) >=  1.0 && abs(prbscEta) < 1.5 && prbmva > 0.465;
  bool mva3 = abs(prbscEta) >=  1.5 && abs(prbscEta) < 2.5 && prbmva > 0.518;
  bool mvalowpt = pele.pt < 20 && (mva1 || mva2 || mva3);

  bool mva4 = abs(prbscEta) <  1.0 && prbmva > 0.942;
  bool mva5 = abs(prbscEta) >=  1.0 && abs(prbscEta) < 1.5 && prbmva > 0.947;
  bool mva6 = abs(prbscEta) >=  1.5 && abs(prbscEta) < 2.5 && prbmva > 0.878;
  bool mvahighpt = pele.pt > 20 && (mva4 || mva5 || mva6);

  bool ElectronId = mvalowpt || mvahighpt;

  int sbit = 0;
  //  if (pele.simpleEleId85cIso != AnaUtil::cutValue(_electronCutMap, "eleId"))     sbit |= (1 << 0);
  if (!ElectronId)                                                                sbit |= (1 << 0);
  bool isGoodVtx;
  TVector3 vele = findLeptonVtx(pele.vtxIndex, isGoodVtx);
  double dz = vele.z() - vz;
  if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(_electronCutMap, "dz"))          sbit |= (1 << 1);
  if (pele.pixHits <= AnaUtil::cutValue(_electronCutMap,"pixHits"))               sbit |= (1 << 2);
  if (pele.trkHits <= AnaUtil::cutValue(_electronCutMap,"trkHits"))               sbit |= (1 << 3);
  if (abs(pele.trkD0) >= AnaUtil::cutValue(_electronCutMap,"trkD0"))              sbit |= (1 << 4);
    
    
  AnaUtil::fillHist1D("probeEta_before", prbeta, _puevWt);
    
  if (sbit == 0) {
    AnaUtil::fillHist1D("probeEta_after", prbeta, _puevWt);
    AnaUtil::fillHist1D("probept_beforeEtaCut", prbpt, _puevWt);
    if (abs(prbeta) < 1.4442)
      AnaUtil::fillHist1D("probept_after_low", prbpt, _puevWt);
    else if (1.566 <= abs(prbeta) && abs(prbeta) < 2.5)
      AnaUtil::fillHist1D("probept_after_high", prbpt, _puevWt);
  }

  if (abs(prbeta) >= AnaUtil::cutValue(_electronCutMap, "eta")) return;                //acceptance cut
  // Filling Histograms for finding out efficiency
  if (sbit == 0) {
    // Allowed pseudorapidity region
    if (8 < prbpt && prbpt <= 10)        AnaUtil::fillHist1D("mass_ptle10_pass", mass, _puevWt);
    else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_pass", mass, _puevWt);
    else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_pass", mass, _puevWt);
    else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_pass", mass, _puevWt);
    else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_pass", mass, _puevWt);
    else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_pass", mass, _puevWt);
    if (pele.pfRelIsoDB04v2 < AnaUtil::cutValue(_electronCutMap, "relIso")) {
      if (8 < prbpt && prbpt <= 10)        AnaUtil::fillHist1D("mass_ptle10_passIso", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_passIso", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_passIso", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_passIso", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_passIso", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_passIso", mass, _puevWt);
    }
    else {
      if (8 < prbpt && prbpt <= 10)        AnaUtil::fillHist1D("mass_ptle10_failIso", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_failIso", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_failIso", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_failIso", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_failIso", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_failIso", mass, _puevWt); 
    }
    // Central region
    if (abs(prbeta) < 1.4442) {
      if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_pass_centralEta", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_pass_centralEta", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_pass_centralEta", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_pass_centralEta", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_pass_centralEta", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_pass_centralEta", mass, _puevWt);
      // study isolation only when sbit = 0
      if (pele.pfRelIsoDB04v2 < AnaUtil::cutValue(_electronCutMap, "relIso")) {
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
    else if (1.566 <= abs(prbeta)) {
      if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_pass_endEta", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_pass_endEta", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_pass_endEta", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_pass_endEta", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_pass_endEta", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_pass_endEta", mass, _puevWt);
     
      if (pele.pfRelIsoDB04v2 < AnaUtil::cutValue(_electronCutMap, "relIso")) {
	if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_passIso_endEta", mass, _puevWt);
	else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_passIso_endEta", mass, _puevWt);
        else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_passIso_endEta", mass, _puevWt);
	else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_passIso_endEta", mass, _puevWt);
	else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_passIso_endEta", mass, _puevWt);
	else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_passIso_endEta", mass, _puevWt);
      }
      else {
	if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_failIso_endEta", mass, _puevWt);
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
    if (abs(prbeta) < 1.4442) {
      if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_fail_centralEta", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_fail_centralEta", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_fail_centralEta", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_fail_centralEta", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_fail_centralEta", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_fail_centralEta", mass, _puevWt);
    }
    // Forward region
    else if (1.566 <= abs(prbeta) && abs(prbeta) < 2.5) {
      if ( 8 < prbpt && prbpt <= 10)       AnaUtil::fillHist1D("mass_ptle10_fail_endEta", mass, _puevWt);
      else if (10 < prbpt && prbpt <= 20)  AnaUtil::fillHist1D("mass_ptle20_fail_endEta", mass, _puevWt);
      else if (20 < prbpt && prbpt <= 30)  AnaUtil::fillHist1D("mass_ptle30_fail_endEta", mass, _puevWt);
      else if (30 < prbpt && prbpt <= 40)  AnaUtil::fillHist1D("mass_ptle40_fail_endEta", mass, _puevWt);
      else if (40 < prbpt && prbpt <= 50)  AnaUtil::fillHist1D("mass_ptle50_fail_endEta", mass, _puevWt);
      else if (50 < prbpt && prbpt <= 100) AnaUtil::fillHist1D("mass_ptle100_fail_endEta", mass, _puevWt);
    }
  } 
  // =====================================================================================================
   ++nProbe[0];
  ++nSingleCut[0];
  AnaUtil::fillHist1D("probept_before", prbpt, _puevWt); // within total acceptance range 
 

  if (sbit == 0) {
    AnaUtil::fillHist1D("probept_after",prbpt,_puevWt);
    AnaUtil::fillProfile("pt_pfRelIso_noIsoCut", prbpt, pele.pfRelIsoDB04v2, _puevWt);
    AnaUtil::fillProfile("eta_pfRelIso_noIsoCut", prbeta, pele.pfRelIsoDB04v2, _puevWt);

    if (pele.pfRelIsoDB04v2 <= AnaUtil::cutValue(_electronCutMap, "relIso")) {
      AnaUtil::fillHist1D("probept_after_AllCuts", prbpt, _puevWt);
      AnaUtil::fillProfile("pt_pfRelIso_afterAllCuts", prbpt, pele.pfRelIsoDB04v2, _puevWt);
      AnaUtil::fillProfile("eta_pfRelIso_afterAllCuts", prbeta, pele.pfRelIsoDB04v2, _puevWt);

      if (abs(prbeta) < 1.4442)
  	AnaUtil::fillHist1D("probept_after_AllCuts_lowEta", prbpt, _puevWt);
      else if (1.566 <= abs(prbeta) && abs(prbeta) < 2.5)
  	AnaUtil::fillHist1D("probept_after_AllCuts_highEta", prbpt, _puevWt);

    }

  }
  AnaUtil::fillProfile("pt_pfRelIso", prbpt, pele.pfRelIso, _puevWt);
  AnaUtil::fillHist1D("pfRelIso", pele.pfRelIso, _puevWt);
  AnaUtil::fillHist1D("RelIso", pele.relIso, _puevWt);
  AnaUtil::fillHist1D("mvaIso", pele.isoMVA, _puevWt);
  if (abs(prbeta) < 1.4442){
    AnaUtil::fillHist1D("pfRelIsoBarrel", pele.pfRelIso, _puevWt);
    AnaUtil::fillHist1D("RelIsoBarrel", pele.relIso, _puevWt);
    AnaUtil::fillHist1D("mvaIsoBarrel", pele.isoMVA, _puevWt);
  }
  else if (1.566 <= abs(prbeta) && abs(prbeta) < 2.5){
    AnaUtil::fillHist1D("pfRelIsoEndCap", pele.pfRelIso, _puevWt);
    AnaUtil::fillHist1D("RelIsoEndCap", pele.relIso, _puevWt);
    AnaUtil::fillHist1D("mvaIsoEndCap", pele.isoMVA, _puevWt);
  }


  // Isolation Cut 
  if (pele.pfRelIso >= AnaUtil::cutValue(_electronCutMap, "relIso")) sbit |= (1 << 5);

  for (int i = 1; i < 7; ++i) {
    int j = pow(2,i) - 1;
    if ((sbit & j) == 0) 
     ++nProbe[i];
    if ((sbit & (1<< (i-1))) == 0) {
      ++nSingleCut[i];
      if (abs(pele.eta) < 1.4442)
        AnaUtil::fillHist1D(HistAfterCentral[i-1], prbpt, _puevWt);
      else if (abs(pele.eta) >= 1.566)
	AnaUtil::fillHist1D(HistAfterEnd[i-1], prbpt, _puevWt);
    }
  }

}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void ElectronEff::endJob() {  
  // Print Muon ID Efficiency Numbers
  const string tagnprobe[] =
  { 
    "eEleId85cIso",
    "dz",
    "pixHits",
    "trkHits",
    "trkD0"
    "Iso"                
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
bool ElectronEff::readJob(const string& jobFile, int& nFiles)
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
void ElectronEff::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));

  // Trigger Path List
  if (_useTrigger) 
    AnaUtil::showList(_triggerPathTagList, ">>> INFO. Trigger Paths for Tag Leg:", os);

  AnaUtil::showCuts(hmap, os);
}
