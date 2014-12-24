#include "configana.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <cmath>

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

#include "ELTau.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"
#include "MVASkim.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::abs;
using std::fabs;
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
ELTau::ELTau()
  : AnaBase(),
  _createMVATree(false),
  _readMVA(false),
  _mvaInputFile(""),
  _skimObj(0)

{}
// ----------
// Destructor
// ----------
ELTau::~ELTau() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool ELTau::beginJob() 
{
  AnaBase::beginJob();
  histf()->cd();
  bookHistograms();
  
  return true;
}
// ---------------
// Book histograms
// ---------------
void ELTau::bookHistograms() 
{
  new TH1D("evcounter", "Selected event counter", 43, -0.5, 42.5);
  new TH1D("elcounter", "electron Selection event counter", 6, -0.5, 5.5);
  new TH1D("mucounter", "muon Selection event counter", 10, -0.5, 9.5);
  new TH1D("taucounter", "tau Selection event counter", 6, -0.5, 5.5);

  new TH1D("counter_emt", "Selected event counter", 15, -0.5, 14.5);
  new TH1D("counter_eet", "Selected event counter", 14, -0.5, 13.5);

  new TH1D("Electron_Pt_emt", "#P_T of electron", 140, 0, 140);
  new TH1D("Muon_Pt_emt", "#P_T of Muon", 140, 0, 140);
  new TH1D("Tau_Pt_emt", "#P_T of Tau", 140, 0, 140);
  new TH1D("visMass_l2t_emtChannel", "visible mass of sub leading lepton and #tau", 30, 0, 300);
  new TH1D("MET_emtChannel", "met of event for emt channel", 100, 0, 100);

  new TH1D("Electron1_Pt_eet", "#P_T of electron", 140, 0, 140);
  new TH1D("Electron2_Pt_eet", "#P_T of Muon", 140, 0, 140);
  new TH1D("Tau_Pt_eet", "#P_T of Tau", 140, 0, 140);
  new TH1D("visMass_e2t_eetChannel", "visible mass of sub leading electron and #tau", 30, 0, 300);
  new TH1D("MET_eetChannel", "met of event for eet channel", 100, 0, 100);
}

// -------------------
// The main event loop
// -------------------
void ELTau::clearLists() {
  vtxList.clear();
  bjetList.clear();
  trigObjList.clear();
  genMuList.clear();
  genEleList.clear();
  genTauList.clear();
  genList.clear(); // for GenLevel Matching
  fEleList.clear();
  fMuList.clear();
  fTauList.clear();
}
void ELTau::eventLoop() 
{

  // Initialize analysis
  if (!beginJob()) return;
  int nPrint = max(10000, nEvents()/1000);

  Options op;
  op.verbose = false;
  op.usesbit = true;  // Crucial
  op.printselected = false;

  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  for (int ev = 0; ev < nEvents(); ++ev) {
    clearEvent();
    int lflag = chain()->LoadTree(ev);
    int nbytes = getEntry(lflag);    // returns total bytes read

    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));

    const Event& evt = eventColl()->at(0);

    histf()->cd();

    // PileUP weight
    puevWt_ = 1; // for data

#if 0
    if (isMC()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
    }
#endif

    AnaUtil::fillHist1D("evcounter", 0, puevWt_);
    AnaUtil::fillHist1D("counter_emt", 0, puevWt_);
    AnaUtil::fillHist1D("counter_eet", 0, puevWt_);


    // Show status of the run
    int run   = evt.run;
    int event = evt.event;
    int lumis = evt.lumis;

    if (currentFile != lastFile) 
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << currentFile 
	   << " <<< Run# " << run
	   << " Lumis# " << lumis
	   << " Event# " << setw(8) << event << " >>> " 
	   << " Events proc. " << setw(8) << ev 
	   << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << chain()->GetCurrentFile()->GetName()
	   << " <<< Run# " << run 
	   << " Lumis# " << lumis
	   << " Event# " << setw(8) << event << " >>> " 
	   << " Events proc. " << setw(8) << ev 
	   << endl;

    clearLists();
#if 0
    if (eventIdMap().size()) {
      std::ostringstream mkey;
      mkey << event << "-" << lumis << "-" << run;
      if (eventIdMap().find(mkey.str()) == eventIdMap().end()) continue;
    }
#endif

    if (isMC()) {

      findLLTGenInfo(genMuList, genEleList, genTauList);
      if(!emt && !eet) continue;
      if(emt)
	counter = "counter_emt";
      else if(eet)
	counter = "counter_eet" ;
      //findGenInfo(11, genEleList, genTauList);
      //bool genflag = (genEleList.size() == 1 && genTauList.size() == 2); // for test (gen level pdgid matching)
      //if (!genflag)  continue;
#if 0
      if (emt) {
	fLog() << "=> Event " << event
	       << " Lumis " << lumis
	       << " Run " << run
	       << endl;
	fLog() << "=> nEle: " << genEleList.size()
	       << " nMu: " << genMuList.size()
	       << " nTau: " << genTauList.size()
	       << endl;
	dumpGenInfo(fLog());
      }
#endif
    }

    AnaUtil::fillHist1D(counter.c_str(), 1, puevWt_);

    // Trigger selection
    //if (useTrigger() && !isTriggered()) continue;
    AnaUtil::fillHist1D(counter.c_str(), 2, puevWt_);
    //dumpTriggerPaths(fLog(), true);
    //dumpTriggerObjectInfo(trigObjList, fLog());

    op.verbose = (logOption() >> 1 & 0x1); 
    findVtxInfo(vtxList, op, fLog());
    // Event Selection Starts here .....
    // presence of > 1 good vertex
    if (vtxList.size() < 1) continue;
    AnaUtil::fillHist1D(counter.c_str(), 3, puevWt_);
    op.verbose = (logOption() >> 5 & 0x1);
    //findJetInfo(bjetList, op, fLog());
    //findTriggerObjectInfo(trigObjList);
    selectEvent();
  }  
  // Analysis is over
  endJob();
}
void ELTau::selectEvent() {

  selectElectron();
  selectMuon();
  selectTau();

  if(emt) calculateEMTEff();
  if(eet) calculateEETEff();
}

void ELTau::calculateEETEff() {
  if (fEleList.size() < 2) return;
  AnaUtil::fillHist1D("counter_eet", 4, puevWt_);
  
  if (fTauList.size() < 1) return;
  AnaUtil::fillHist1D("counter_eet", 5, puevWt_);

  sort(fEleList.begin(), fEleList.end(), PtComparator<Electron>());
  sort(fTauList.begin(), fTauList.end(), PtComparator<Tau>());

  const Electron& ele1 = fEleList.at(0);
  const Electron& ele2 = fEleList.at(1);
  const Tau& tau = fTauList.at(0);

  // b-tagged jet Veto
  //  if (bjetList.size() > 0) return;
  // mu veto
  if (AnaBase::vetoMuon(15, 0.2) > 0) return;
  AnaUtil::fillHist1D("counter_eet", 6, puevWt_);

  // e veto
  if (AnaBase::vetoElectron(15, 0.2) > 2) return;
  AnaUtil::fillHist1D("counter_eet", 7, puevWt_);

  double highestPt = std::max(ele1.pt, ele2.pt);
  if (highestPt < 20) return;
  AnaUtil::fillHist1D("counter_eet", 8, puevWt_);

  TLorentzVector E1,E2, T1;
  E1.SetPtEtaPhiE(ele1.pt, ele1.eta, ele1.phi, ele1.energy);
  E2.SetPtEtaPhiE(ele2.pt, ele2.eta, ele2.phi, ele2.energy);
  T1.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);

  if (ele1.charge + ele2.charge == 0) return;
  AnaUtil::fillHist1D("counter_eet", 9, puevWt_);
  if (ele1.charge + tau.charge != 0) return;
  AnaUtil::fillHist1D("counter_eet", 10, puevWt_);


  double dr_e1t = E1.DeltaR(T1);
  double dr_e2t = E2.DeltaR(T1);
  double dr_ee = E1.DeltaR(E2);

  bool dr = dr_e1t <= 0.5
         || dr_e2t <= 0.5
         || dr_ee <= 0.5;

  if(dr) return;
  AnaUtil::fillHist1D("counter_eet", 11, puevWt_);

  if (tau.againstElectronLooseMVA5 <= 0.5) return;
  AnaUtil::fillHist1D("counter_eet", 12, puevWt_);

  if(tau.againstMuonTight3 <= 0.5) return;
  AnaUtil::fillHist1D("counter_eet", 13, puevWt_);


  AnaUtil::fillHist1D("Electron1_Pt_eet", ele1.pt);
  AnaUtil::fillHist1D("Electron2_Pt_eet", ele2.pt);
  AnaUtil::fillHist1D("Tau_Pt_eet", tau.pt);
  AnaUtil::fillHist1D("visMass_e2t_eetChannel", (E2+T1).M());


}

void ELTau::calculateEMTEff() {
  if (fEleList.size() < 1) return;
  AnaUtil::fillHist1D("counter_emt", 4, puevWt_);
  if (fMuList.size() < 1) return;
  AnaUtil::fillHist1D("counter_emt", 5, puevWt_);
  if (fTauList.size() < 1) return;
  AnaUtil::fillHist1D("counter_emt", 6, puevWt_);

  sort(fEleList.begin(), fEleList.end(), PtComparator<Electron>());
  sort(fMuList.begin(), fMuList.end(), PtComparator<Muon>());
  sort(fTauList.begin(), fTauList.end(), PtComparator<Tau>());

  const Electron& ele = fEleList.at(0);
  const Muon& mu = fMuList.at(0);
  const Tau& tau = fTauList.at(0);

  // b-tagged jet Veto                                                                                                                 
  //  if (bjetList.size() > 0) return;
  // mu veto
  if (AnaBase::vetoMuon(15, 0.2) > 0) return;
  AnaUtil::fillHist1D("counter_emt", 7, puevWt_);

  // e veto
  if (AnaBase::vetoElectron(15, 0.2) > 1) return;
  AnaUtil::fillHist1D("counter_emt", 8, puevWt_);

  double highestPt = std::max(ele.pt, mu.pt);
  if (highestPt < 20) return;
  AnaUtil::fillHist1D("counter_emt", 9, puevWt_);
  
  TLorentzVector E1,M1, T1;
  E1.SetPtEtaPhiE(ele.pt, ele.eta, ele.phi, ele.energy);
  M1.SetPtEtaPhiE(mu.pt, mu.eta, mu.phi, mu.energy);
  T1.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
  
  if (ele.charge + mu.charge == 0) return;
  AnaUtil::fillHist1D("counter_emt", 10, puevWt_);
  if (ele.charge + tau.charge != 0) return;
  AnaUtil::fillHist1D("counter_emt", 11, puevWt_);
  

  double dr_et = E1.DeltaR(T1);
  double dr_mt = M1.DeltaR(T1);
  double dr_em = E1.DeltaR(M1);

  bool dr = dr_et <= 0.5
         || dr_mt <= 0.5
         || dr_em <= 0.5;

  if(dr) return;
  AnaUtil::fillHist1D("counter_emt", 12, puevWt_);

  double mass_ET = (E1+T1).M();
  double mass_MT = (M1+T1).M();

  if (fabs(mass_ET -91) < 20 ) {
    if(tau.againstElectronMediumMVA5 <= 0.5) return;
  }
  else if (tau.againstElectronLooseMVA5 <= 0.5) return;
  AnaUtil::fillHist1D("counter_emt", 13, puevWt_);  

  if (fabs(mass_MT -91) < 20) {
    if(tau.againstMuonTight3 <= 0.5) return;
  }
  
  else if (tau.againstMuonLoose3 <= 0.5) return;
  AnaUtil::fillHist1D("counter_emt", 14, puevWt_);


  AnaUtil::fillHist1D("Electron_Pt_emt", ele.pt);
  AnaUtil::fillHist1D("Muon_Pt_emt", mu.pt);
  AnaUtil::fillHist1D("Tau_Pt_emt", tau.pt);
  //  AnaUtil::fillHist1D("visMass_e2t_eetChannel", );

  // ensure no extra significant objects
  // DR angle between any two objects
  // plot various quantities

}

void ELTau::selectElectron() {
  int indx=0;
  int ecounters[] = {0,0,0,0,0,0};
  //vector<Electron> fEleList;
  for (auto it = electronColl()->begin(); it != electronColl()->end(); ++it,++indx) {
    const Electron&  electron = (*it);

    if (electron.pt <= AnaUtil::cutValue(electronCutMap(), "pt")) continue;
    ++ecounters[0];

    double eleta  = std::fabs(electron.eta);
    bool   etaCut = (eleta >= AnaUtil::cutValue(electronCutMap(), "etaLow") &&
		     eleta <= AnaUtil::cutValue(electronCutMap(), "etaUp")) ||
                    eleta  >= AnaUtil::cutValue(electronCutMap(), "eta");
    if (etaCut) continue;
    ++ecounters[1];

    if(!eleId(electron,25)) continue;
    ++ecounters[2];

    if (std::fabs(electron.trkD0) >= AnaUtil::cutValue(electronCutMap(),"trkD0")) continue;
    ++ecounters[3];

    if (std::fabs(electron.vtxDistZ) >= AnaUtil::cutValue(electronCutMap(), "dz")) continue;
    ++ecounters[4];

    if (electron.pfRelIso > 0.2) continue;
    ++ecounters[5];

    fEleList.push_back(electron);
    //electron_indx = indx;
    //break;
  }
  for (size_t i = 0; i < NEL(ecounters); ++i) {
    if (ecounters[i]) AnaUtil::fillHist1D("elcounter", i, puevWt_);
  }
}

void ELTau::selectMuon() {
  double vz = vtxList.at(0).z;
  int mcounters[] = {0,0,0,0,0,0,0,0,0,0};
  //vector<Muon> fMuList;
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);
    
    if (muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt")
	|| fabs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta")) continue;
    ++mcounters[0];
    
    if (!muon.isTrackerMuon || !muon.isGlobalMuonPromptTight) continue;
    ++mcounters[1];
    if (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) continue;
    ++mcounters[2];
    
    if (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches") ||
        muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) continue;
    ++mcounters[3];
    
    if (muon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) continue; //// ???
    ++mcounters[4];
    
    if (muon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) continue;
    ++mcounters[5];
    
    if (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) continue;
    ++mcounters[6];
    
    if (fabs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) continue;
    ++mcounters[7];
    
    if (abs(muon.vtxDistZ) >= AnaUtil::cutValue(muonCutMap(), "vtxDistZ")) continue;
    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;
    ++mcounters[8];
    
    if (fabs(muon.pfChargedIsoR04/muon.pt) >= AnaUtil::cutValue(muonCutMap(), "relIso")) continue; /// 0.15(0.1) at endcap(barrel)
    ++mcounters[9];
    
    fMuList.push_back(muon);
  }

  for (size_t i = 0; i < NEL(mcounters); ++i) {
    if (mcounters[i]) AnaUtil::fillHist1D("mucounter", i, puevWt_);
  }

}


void ELTau::selectTau() {
  int taucounters[] = {0,0,0,0,0,0};
  for (auto it = tauColl()->begin(); it != tauColl()->end(); ++it) {
    const Tau& tau = (*it);
    
    //    std::cout << "ptTau1 = " << AnaUtil::cutValue(_evselCutMap,"ptTau1")) << std::endl;
    if (tau.pt < AnaUtil::cutValue(_evselCutMap,"ptTau1")) continue;
    ++taucounters[0];
    if( std::fabs(tau.eta) >= AnaUtil::cutValue(_evselCutMap,"etaTau1")) continue;
    ++taucounters[1];
    
    if (tau.decayModeFinding != 1.0) continue;
    ++taucounters[2];
    
    if(tau.chargedIsoPtSum >= 2) continue; /// loose combined Isolation
    ++taucounters[3];
    
    if (tau.againstMuonLoose3 <= 0.5) continue;
    ++taucounters[4];
    
    if (tau.againstElectronLooseMVA5 <= 0.5) continue;
    ++taucounters[5];
    
    
    fTauList.push_back(tau);
  }

  for (size_t i = 0; i < NEL(taucounters); ++i) {
    if (taucounters[i]) AnaUtil::fillHist1D("taucounter", i, puevWt_);
  }

}


// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void ELTau::endJob() {  
  histf()->cd();

  TH1 *h = AnaUtil::getHist1D(counter.c_str());
  if (h) {
    fLog() << setprecision(3);
    fLog() << "Events Processed: " << int(h->GetBinContent(1)) << endl;
    for (int i = 2; i <= h->GetNbinsX(); ++i)
      fLog() << "Cut: " 
	     << setw(3) << i 
	     << setw(10) << int(h->GetBinContent(i))
	     << setw(10) << ((h->GetBinContent(i-1)>0) ? h->GetBinContent(i)/h->GetBinContent(i-1) : 0.0)
	     << setw(10) << ((h->GetBinContent(2)>0) ? h->GetBinContent(i)/h->GetBinContent(2) : 0.0)
	     << endl;
  }
  fLog() << resetiosflags(ios::fixed);
  histf()->Write();
  histf()->Close();
  delete  histf();
  //  if (_skimObj) _skimObj->close();
  closeFiles();
}
// -------------------------------------------------------------------------------
// Poor man's way of a datacard. Each line between the 'START' and 'END' tags
// is read in turn, split into words, where the first element is the 'key' and
// the rest the value(s). If more than one values are present they are 
// stored in a vector. No safety mechanism is in place. Any line with an unknown 
// key is skipped. Comments lines should start with either '#' or '//', preferably
// in the first column. Empty lines are skipped. The file containing the datacards 
// is passed as the only argument of the program, there is no default
// -------------------------------------------------------------------------------
bool ELTau::readJob(const string& jobFile, int& nFiles) {
  if (!AnaBase::readJob(jobFile, nFiles)) return false;
  
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
    
    // Split the line into words
    AnaUtil::tokenize(line, tokens);
    assert(tokens.size() > 1);
    string key = tokens[0];
    string value = tokens[1];
    if (key == "createMVATree")
      _createMVATree = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "mvaInputFile")
      _mvaInputFile = value;
    else if (key == "evselCutList")
      AnaUtil::storeCuts(tokens, hmap);
    else if (key == "readMVA")
      _readMVA = (atoi(value.c_str()) > 0) ? true : false;
    
    tokens.clear();
  }
  // Close the file
  fin.close();
  
  printJob();
  
  return true;
}

void ELTau::printJob(ostream& os) const {
  AnaBase::printJob(os);
  
  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));
  AnaUtil::showCuts(hmap, os);
}
