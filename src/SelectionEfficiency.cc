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
#include "SelectionEfficiency.h"
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
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
SelectionEfficiency::SelectionEfficiency()
  : AnaBase(),
    _createMVATree(false),
    _readMVA(false),
    _readMVAFK(false),
    _mvaInputFile(""),
    _MVAdisFile(""),
    _MVAFKdisFile(""),
    _skimObj(0)
{}
// ----------
// Destructor
// ----------
SelectionEfficiency::~SelectionEfficiency() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool SelectionEfficiency::beginJob() 
{ 
  AnaBase::beginJob();

   histf()->cd();
  bookHistograms();

  // Optionally write selected events to another tree                                                                                                                   
  if (_createMVATree) _skimObj = new MVASkim(_mvaInputFile);

  if (_readMVAFK) {
    reader1 = new TMVA::Reader("Silent");
    reader1->AddVariable("LeadPt", &leadPt);
    reader1->AddVariable("SubPt", &subPt);
    reader1->AddVariable("Met", &met);
    reader1->AddVariable("DeltaRDiTau", &deltaRDiTau);
    reader1->AddVariable("PtRatio", &ptRatio);

    //reader->BookMVA("MLP::MLPBNN", "TMVAClassification_MLPBNN.weights.xml");
    //reader->BookMVA("MLPBNN", _MVAdisFile);
    reader1->BookMVA("BDT8", _MVAFKdisFile);
  }
  if (_readMVA) {
    reader = new TMVA::Reader("Silent");
    //reader->AddVariable("muEta", &muEta);
    reader->AddVariable("muPt", &muPt);
    //reader->AddVariable("tau1Eta", &tau1Eta);
    reader->AddVariable("tau1Pt", &tau1Pt);
    //reader->AddVariable("tau2Eta", &tau2Eta);
    reader->AddVariable("tau2Pt", &tau2Pt);
    //reader->AddVariable("diTaudR", &diTaudR);
    //reader->AddVariable("dphiMuTau1", &dphiMuTau1);
    reader->AddVariable("dphiMuDiTau", &dphiMuDiTau);
    reader->AddVariable("met", &met);
    reader->AddVariable("diTauPt/(tau1Pt+tau2Pt)", &ptRatio);

    //reader->BookMVA("MLP::MLPBNN", "TMVAClassification_MLPBNN.weights.xml");
    //reader->BookMVA("MLPBNN", _MVAdisFile);
    reader->BookMVA("MLP", _MVAdisFile);
  }
  return true;
}
// ---------------
// Book histograms
// ---------------
void SelectionEfficiency::bookHistograms() 
{
  new TH1D("counter_gen", "Selected event counter", 49, -0.5, 48.5);
  new TH1F("muPt", "pt distribution of the muon after all cut", 140, 0, 140);
  new TH1F("tau1Pt", "pt distribution of the tau1 after all cut", 140, 0, 140);
  new TH1F("tau2Pt", "pt distribution of the tau2 after all cut", 140, 0, 140);

  new TH1F("mueta", "eta distribution of the muon after all cut", 100, -2.5, 2.5);
  new TH1F("tau1eta", "eta distribution of the tau1 after all cut", 100, -2.5, 2.5);
  new TH1F("tau2eta", "eta distribution of the tau2 after all cut", 100, -2.5, 2.5);

  new TH1F("mvamet", "met distribution of the event after all cut", 350, 0, 350);
  new TH1F("pfmet", "met distribution of the event after all cut", 350, 0, 350);
  new TH1F("mvamet_t2", "met distribution of the event after Tau2 selection", 350, 0, 350);
  new TH1F("pfmet_t2", "met distribution of the event after Tau2 selection", 350, 0, 350);
  new TH1F("MuTauMass", "vissible mass distribution of the muon and leading tau after all cut", 200, 0, 200);
  new TH1F("Mt_Mu_Met", "transverse mass distribution of the muon and met after all cut", 300, 0, 300);
  new TH1F("diTauPt", "Pt distribution of the diTau after all cut", 300, 0, 300);
  new TH1F("TauTauMass", "vissible mass distribution of the lead and sub-lead tau after mva5  cut", 320, 0, 320);
  new TH2D("elecIso_B", "electron pfRelIso as a function of elec pT in the barrel region", 150, 0, 150, 100, 0, 1);  
  new TH2D("elecIso_E", "electron pfRelIso as a function of elec pT in the endcap region", 150, 0, 150, 100, 0, 1);  
  
  new TH1F("mvaoutput_BL", "mvaoutput after BL", 100, -2, 2);
  new TH1F("mvaOutputFK_BL", "fake mvaoutput after BL", 100, -1, 1);
  new TH1F("TauTauMass_BL", "vissible mass distribution of the lead and sub-lead tau after all BL cut", 320, 0, 320);
  new TH1F("filter", "WH normal vs WH contamination", 5, 0.5, 5.5);

  new TH1F("TauTauMass_mva1", "vissible mass distribution of the lead and sub-lead tau after all BL cut+ mva1", 320, 0, 320);
  new TH1F("TauTauMass_mva2", "vissible mass distribution of the lead and sub-lead tau after all BL cut+ mva2", 320, 0, 320);
  new TH1F("TauTauMass_mva3", "vissible mass distribution of the lead and sub-lead tau after all BL cut+ mva3", 320, 0, 320);
  new TH1F("TauTauMass_mva4", "vissible mass distribution of the lead and sub-lead tau after all BL cut+ mva4", 320, 0, 320);
  new TH1F("TauTauMass_mva5", "vissible mass distribution of the lead and sub-lead tau after all BL cut+ mva5", 320, 0, 320);
  new TH1F("TauTauMass_mva6", "vissible mass distribution of the lead and sub-lead tau after all BL cut+ mva6", 320, 0, 320);
  new TH1F("TauTauMass_mva7", "vissible mass distribution of the lead and sub-lead tau after all BL cut+ mva7", 320, 0, 320);

  new TH1F("mvaoutput_mva1", "Irr  mvaoutput after BL & mva1", 100, -2, 2);
  new TH1F("mvaoutput_mva2", "Irr  mvaoutput after BL & mva2", 100, -2, 2);
  new TH1F("mvaoutput_mva3", "Irr  mvaoutput after BL & mva3", 100, -2, 2);
  new TH1F("mvaoutput_mva4", "Irr  mvaoutput after BL & mva4", 100, -2, 2);
  new TH1F("mvaoutput_mva5", "Irr  mvaoutput after BL & mva5", 100, -2, 2);
  new TH1F("mvaoutput_mva6", "Irr  mvaoutput after BL & mva6", 100, -2, 2);
  new TH1F("mvaoutput_mva7", "Irr  mvaoutput after BL & mva7", 100, -2, 2);

  new TProfile("testdeno", "profile of t1+t2 and pt ratio", 50, 0, 1, 0, 200);
  new TProfile("testnume", "profile of t1+t2 and pt ratio", 50, 0, 1, 0, 200);
  new TH1F("test1", "Pt of Tau1 in the desired pt ratio range", 200, 0, 200);
  new TH1F("test2", "Pt of Tau2 in the desired pt ratio range", 200, 0, 200);
  new TH1D("testnin", "nume in ratio range", 200, 0, 200);
  new TH1D("testdin", "deno in ratio range", 200, 0, 200);
  new TH1D("testnout", "nume in ratio", 200, 0, 200);
  new TH1D("testdout", "deno in ratio", 200, 0, 200);
  new TH1D("ratio", "ratio", 100, 0, 1);

  new TH1F("mupt_mu", "mupt after muon", 140, 0, 140);
  new TH1F("mueta_mu", "mueta after muon", 100, -2.5, 2.5);
  new TH1F("tau1pt_tau1", "tau1pt after tau1", 140, 0, 140);
  new TH1F("tau1eta_tau1", "tau1eta after tau1", 100, -2.5, 2.5); 
}
// -------------------
// The main event loop
// -------------------
void SelectionEfficiency::clearLists() {
  vtxList.clear();
  muoList.clear();
  eleList.clear();
  tauList.clear();
  bjetList.clear();

  genMuonList.clear();
  genTauList.clear();
}
void SelectionEfficiency::eventLoop() 
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

    /*
    if (isMC()) {
      int npu = 0;
      puevWt_ = wtPileUp(npu);
    }
    */
    AnaUtil::fillHist1D("counter_gen", 0, puevWt_);

    int run   = evt.run;
    int event = evt.event;
    int lumis = evt.lumis;

    // Show status of the run
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

    if (eventIdMap().size()) {
      std::ostringstream mkey;
      mkey << event << "-" << lumis << "-" << run;
      if (eventIdMap().find(mkey.str()) == eventIdMap().end()) continue;
      evLog() << "Event " << event
            << " Lumis " << lumis
            << " Run " << run << endl;
      //dumpGenInfo(evLog()); 
      //dumpEvent("1111", evLog(), false);
    }

    if (logOption() > 0) {
      fLog() << "Event " << event
            << " Lumis " << lumis
            << " Run " << run << endl;
      fLog() << "n_tau: "<< ntau()
	     << ", n_muon: "<< nmuon()
	     << ", n_jet: " << njet()
	     << ", n_vertex: " << nvertex()
	     << ", n_met: " << nmet()
             << ", n_electron: " << nelectron()
             << endl;
    }
    clearLists();
 
    //findGenInfo(genMuonList, genTauList);
    //if (genMuonList.size() != 1 || genTauList.size() != 2) continue;
    //if (logOption() >> 1 & 0x1) dumpGenInfo(fLog()); 
  
    //if (!WhHmhFilter()) continue;                                        //Put this off if not required::ATTENTION
    AnaUtil::fillHist1D("counter_gen", 1, puevWt_);

    // Trigger selection
    if (useTrigger() && !isTriggered()) continue;
    AnaUtil::fillHist1D("counter_gen", 2, puevWt_);

    op.verbose = (logOption() >> 2 & 0x1); 
    findVtxInfo(vtxList, op, fLog());
    int nvtx = vtxList.size();
    double vz = ( nvtx > 0 ) ? vtxList.at(0).z : -999;

    if (logOption() > 0)
    fLog() << "Event " << event
          << " Lumis " << evt.lumis
          << " Run " << run 
          << " n_vertex_good " << nvtx
          << endl;

    op.verbose = (logOption() >> 5 & 0x1);
    findJetInfo(bjetList, op, fLog());

    // Event Selection Starts here .....
    // presence of > 1 good vertex
    if (vtxList.size() < 1) continue;
    AnaUtil::fillHist1D("counter_gen", 3, puevWt_);

    selectEvent();

  }  
  // Analysis is over
  endJob();
}
void SelectionEfficiency::selectEvent() {
  study_eff();
}
void SelectionEfficiency::study_eff()
{
  double vz = vtxList.at(0).z;

  int mcounters[] = {0,0,0,0,0,0,0,0,0,0};
  vector<Muon> fMuoList;
  for (auto it = muonColl()->begin(); it != muonColl()->end(); ++it) {
    const Muon& muon = (*it);

    if (muon.pt <= AnaUtil::cutValue(muonCutMap(), "pt") 
      || abs(muon.eta) >= AnaUtil::cutValue(muonCutMap(), "eta")) continue;
    ++mcounters[0];

    if (!muon.isTrackerMuon || !muon.isGlobalMuonPromptTight) continue;
    ++mcounters[1];

    if (muon.nMatchedStations < AnaUtil::cutValue(muonCutMap(), "nMatchedStations")) continue;
    ++mcounters[2];

    if (muon.nMatches < AnaUtil::cutValue(muonCutMap(), "nMatches") ||
        muon.nChambers < AnaUtil::cutValue(muonCutMap(), "nChambers")) continue;
    ++mcounters[3];

    if (muon.trkHits <= AnaUtil::cutValue(muonCutMap(),"trkHits")) continue;
    ++mcounters[4];

    if (muon.pixHits < AnaUtil::cutValue(muonCutMap(),"pixHits")) continue;
    ++mcounters[5];

    if (muon.globalChi2 >= AnaUtil::cutValue(muonCutMap(),"globalChi2")) continue;
    ++mcounters[6];

    if (abs(muon.trkD0) >= AnaUtil::cutValue(muonCutMap(),"trkD0")) continue;
    ++mcounters[7];

    //if (abs(muon.vtxDistZ) >= AnaUtil::cutValue(muonCutMap(), "vtxDistZ")) continue;
    bool isGoodVtx;
    TVector3 vmu = findLeptonVtx(muon.vtxIndex, isGoodVtx);
    double dz = vmu.z() - vz;
    if (!isGoodVtx || abs(dz) >= AnaUtil::cutValue(muonCutMap(), "dz")) continue;
    ++mcounters[8];

    //if (abs(muon.pfRelIso) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso")) continue;
    if (abs(muon.relIso) >= AnaUtil::cutValue(muonCutMap(), "pfRelIso")) continue;
    ++mcounters[9];

    fMuoList.push_back(muon);
  }
  int ishift = 4;
  for (size_t i = 0; i < NEL(mcounters); ++i) {
    if (mcounters[i]) AnaUtil::fillHist1D("counter_gen", ishift + i, puevWt_);
  }
  
  // atleast 1 good muon
  if (fMuoList.size() < 1) return;
  AnaUtil::fillHist1D("counter_gen", 14, puevWt_);
  const Muon& muo = fMuoList.at(0);

  AnaUtil::fillHist1D("mupt_mu", muo.pt, puevWt_);
  AnaUtil::fillHist1D("mueta_mu", muo.eta, puevWt_);

  TLorentzVector M, T1;
  M.SetPtEtaPhiE(muo.pt, muo.eta, muo.phi, muo.energy);

  //---------------------------------------------
  //
  //              Tau1 SELECTION//OS to Muon
  //
  //---------------------------------------------
  int t1counters[] = {0,0,0,0,0,0,0,0};
  vector<Tau> fTau1List;
  //int tau1_indx;
  for (auto it = tauColl()->begin(); it != tauColl()->end(); ++it) {
    const Tau& tau = (*it);

    T1.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
    double dr = AnaUtil::deltaR(M, T1);
    if (dr < AnaUtil::cutValue(_tau1CutMap, "drMuTau")) continue;
    ++t1counters[0];

    if (tau.pt <= AnaUtil::cutValue(_tau2CutMap, "pt") || 
        abs(tau.eta) >= AnaUtil::cutValue(_tau1CutMap, "eta")) continue;      //Notice its Tau2 Pt //not wrong//minimum of the two tau's
    ++t1counters[1];

    if (tau.decayModeFinding <= 0.5) continue;
    ++t1counters[2];

    if (tau.byLooseCombinedIsolationDeltaBetaCorr3Hits <= 0.5) continue;                     //Inverted-Medium Iso for MVA test, otherwise put it TIGHT
    ++t1counters[3];
    
    if (tau.againstMuonTight <= 0.5) continue;                        
    ++t1counters[4];

    if (tau.againstElectronLoose <= 0.5) continue;                                            
    ++t1counters[5];

    if (abs(tau.zvertex - vz) >= AnaUtil::cutValue(_tau1CutMap, "dz")) continue;   
    ++t1counters[6];

    if ((muo.charge + tau.charge) != 0) continue;   
    ++t1counters[7];

    fTau1List.push_back(tau);
  }
  ishift = 15;
  for (size_t i = 0; i < NEL(t1counters); ++i) {
    if (t1counters[i]) AnaUtil::fillHist1D("counter_gen", ishift + i, puevWt_);
  }
  
  if (fTau1List.size() < 1) return; 
  AnaUtil::fillHist1D("counter_gen", 23, puevWt_);

  AnaUtil::fillHist1D("tau1pt_tau1", fTau1List.at(0).pt , puevWt_);
  AnaUtil::fillHist1D("tau1eta_tau1", fTau1List.at(0).eta , puevWt_);

  //---------------------------------------------
  //
  //              Tau2 SELECTION
  //
  //---------------------------------------------
  int t2counters[] = {0,0,0,0,0,0,0,0,0};
  vector<Tau> fTau2List;
  TLorentzVector T2;
  for (auto it = tauColl()->begin(); it != tauColl()->end(); ++it) {
    const Tau& tau = (*it);

    T2.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
    double dr = AnaUtil::deltaR(M, T2);
    //double dr1 = AnaUtil::deltaR(T1, T2);
    if (dr < AnaUtil::cutValue(_tau2CutMap, "drMuTau")) continue;
    ++t2counters[0];

    //if (dr1 < AnaUtil::cutValue(_tau2CutMap, "drTauTau")) continue;
    ++t2counters[1];

    if (tau.pt <= AnaUtil::cutValue(_tau2CutMap, "pt") || 
      abs(tau.eta) >= AnaUtil::cutValue(_tau2CutMap, "eta")) continue;
    ++t2counters[2];

    if (tau.decayModeFinding <= 0.5) continue;
    ++t2counters[3];
  
    if (tau.byMediumCombinedIsolationDeltaBetaCorr3Hits <= 0.5) continue;                        
    ++t2counters[4];

    if (tau.againstMuonTight <= 0.5) continue; 
    ++t2counters[5];

    if (tau.againstElectronLoose <= 0.5) continue;                         
    ++t2counters[6];

    if (abs(tau.zvertex - vz) >= AnaUtil::cutValue(_tau2CutMap, "dz")) continue;   
    ++t2counters[7];

    if ((muo.charge + tau.charge) == 0) continue;   
    ++t2counters[8];

    fTau2List.push_back(tau);
  }
  ishift = 24;
  for (size_t i = 0; i < NEL(t2counters); ++i) {
    if (t2counters[i]) AnaUtil::fillHist1D("counter_gen", ishift + i, puevWt_);
  }


  if (fTau2List.size() < 1) return; 
  AnaUtil::fillHist1D("counter_gen", 33, puevWt_);


  int tau1_indx, tau2_indx;
  bool taupair = false;
  for (unsigned int indx = 0; indx < fTau1List.size(); ++indx) {
    const Tau& tau1 = fTau1List[indx];     
    T1.SetPtEtaPhiE(tau1.pt, tau1.eta, tau1.phi, tau1.energy);
    for (unsigned int jndx = 0; jndx < fTau2List.size(); ++jndx) {
      const Tau& tau2 = fTau2List[jndx];     
      T2.SetPtEtaPhiE(tau2.pt, tau2.eta, tau2.phi, tau2.energy);
      if (AnaUtil::deltaR(T1, T2) >= 0.1 && (T1.Pt() > AnaUtil::cutValue(_tau1CutMap, "pt") || T2.Pt() > AnaUtil::cutValue(_tau1CutMap, "pt"))) {
        taupair = true;
        tau2_indx = jndx;
        break;
      } 
    }
    if (taupair) {
      tau1_indx = indx;
      break;
    }
  }

  if (!taupair) return;
  AnaUtil::fillHist1D("counter_gen", 34, puevWt_);

  // Here comes MET
  const MET& mvamet = metColl()->at(0);
  const MET& pfmet = metColl()->at(0);

  AnaUtil::fillHist1D("mvamet_t2", mvamet.met, puevWt_);
  AnaUtil::fillHist1D("pfmet_t2", pfmet.met, puevWt_);

  const Tau& taua = fTau1List.at(tau1_indx);
  const Tau& taub = fTau2List.at(tau2_indx);

  if (abs(taua.zvertex - taub.zvertex) >= 0.14) return;
  AnaUtil::fillHist1D("counter_gen", 35, puevWt_);

  bool isGoodVtx;
  TVector3 vtxmu = findLeptonVtx(muo.vtxIndex, isGoodVtx);
  double mutaudz = abs(taua.zvertex - vtxmu.z());
  if (!isGoodVtx || mutaudz >= 0.14) return;
  AnaUtil::fillHist1D("counter_gen", 36, puevWt_);

  //---------------------------------------------
  //
  //              Muon Veto
  //
  //---------------------------------------------
  //if (vetoMuon(taua.zvertex, 15, 0.14) != 1) return;
  AnaUtil::fillHist1D("counter_gen", 37, puevWt_);
  
  //---------------------------------------------
  //
  //              Electron Veto
  //
  //---------------------------------------------
  //if (vetoElectron(taua.zvertex, 10, 0.14) != 0) return;
  AnaUtil::fillHist1D("counter_gen", 38, puevWt_);

  //---------------------------------------------
  //
  //              b-tagged Jets Veto
  //
  //---------------------------------------------
  if (bjetList.size() > 0) return; 
  AnaUtil::fillHist1D("counter_gen", 39, puevWt_);


  TLorentzVector zMuTau, diTau;
  //T2.SetPtEtaPhiE(taub.pt, taub.eta, taub.phi, taub.energy);     // already defined earlier

  zMuTau = M + T1;
  diTau = T1 + T2;

  //////////////////////////////////
  //////////////////////////////////
  //if (mvamet.met < 20) return;                        // No MET Cut for MVA study
  //if (pfmet.met < 20) return;                        // No MET Cut for MVA study
  AnaUtil::fillHist1D("counter_gen", 40, puevWt_);

  double mass1 = sqrt(2*M.Pt()*mvamet.met*(1-cos(AnaUtil::deltaPhi(M.Phi(), mvamet.metphi))));
  //double mass1 = sqrt(2*M.Pt()*pfmet.met*(1-cos(AnaUtil::deltaPhi(M.Phi(), pfmet.metphi))));

  //if (mass1 <= 20) return;
  AnaUtil::fillHist1D("counter_gen", 41, puevWt_);


  //if (zMuTau.M() < 80 && diTau.Pt() < 50) return;     // to reduce Z->TauTau bkg
  AnaUtil::fillHist1D("counter_gen", 42, puevWt_);


/*
  //This part will be deleted
  double ratioo = diTau.Pt()/(T1.Pt() + T2.Pt());
  if (ratioo >= 0.6 && ratioo <= 0.7) { 
    AnaUtil::fillHist1D("test1", T1.Pt(), puevWt_);
    AnaUtil::fillHist1D("test2", T2.Pt(), puevWt_);
    AnaUtil::fillHist1D("testnin", diTau.Pt(), puevWt_);
    AnaUtil::fillHist1D("testdin", T1.Pt()+T2.Pt(), puevWt_);
  }
  AnaUtil::fillProfile("testnume", ratioo, diTau.Pt(), puevWt_);
  AnaUtil::fillProfile("testdeno", ratioo, T1.Pt()+T2.Pt(),  puevWt_);
  AnaUtil::fillHist1D("testnout", diTau.Pt(), puevWt_);
  AnaUtil::fillHist1D("testdout", T1.Pt()+T2.Pt(), puevWt_);
  AnaUtil::fillHist1D("ratio", ratioo, puevWt_);
*/

  double mvaOutputFK = -999;       
  if (_readMVAFK) {
    leadPt     = T1.Pt();
    subPt     = T2.Pt();
    met        = mvamet.met;
    deltaRDiTau   = AnaUtil::deltaR(T1, T2);
    ptRatio = diTau.Pt()/(T1.Pt() + T2.Pt());

    mvaOutputFK = reader1->EvaluateMVA("BDT8");
    //if (mvaOutputFK <= -0.106) return;     
  }
  AnaUtil::fillHist1D("counter_gen", 43, puevWt_);
  AnaUtil::fillHist1D("TauTauMass_BL", diTau.M(), puevWt_);
  AnaUtil::fillHist1D("mvaOutputFK_BL", mvaOutputFK, puevWt_);


  double mvaOutput = -999;      //Set at high negative value, as MvaOutput > -0.106 
  if (_readMVA) {
    //muEta      = M.Eta();
    muPt       = M.Pt();
    //tau1Eta    = T1.Eta();
    tau1Pt     = T1.Pt();
    //tau2Eta    = T2.Eta();
    tau2Pt     = T2.Pt();
    //diTaudR   = AnaUtil::deltaR(T1, T2);
    //dphiMuTau1 = AnaUtil::deltaPhi(M, T1);
    dphiMuDiTau = AnaUtil::deltaPhi(M, diTau);
    met        = mvamet.met;
    ptRatio = diTau.Pt()/(T1.Pt() + T2.Pt());

    //mvaOutput = reader->EvaluateMVA("MLPBNN");
    mvaOutput = reader->EvaluateMVA("MLP");

    //for the time being returning at this value
    //if (mvaOutput <= 0.75) return;     //This value has to be optimized
  }
  AnaUtil::fillHist1D("mvaoutput_BL", mvaOutput, puevWt_);

  if (mvaOutputFK > -0.40) {     
    AnaUtil::fillHist1D("counter_gen", 44, puevWt_);
    AnaUtil::fillHist1D("mvaoutput_mva1", mvaOutput, puevWt_);
    AnaUtil::fillHist1D("TauTauMass_mva1", diTau.M(), puevWt_);
  }
  if (mvaOutputFK > -0.35) {     
    AnaUtil::fillHist1D("counter_gen", 45, puevWt_);
    AnaUtil::fillHist1D("mvaoutput_mva2", mvaOutput, puevWt_);
    AnaUtil::fillHist1D("TauTauMass_mva2", diTau.M(), puevWt_);
  }
  if (mvaOutputFK > -0.30) {     
    AnaUtil::fillHist1D("counter_gen", 46, puevWt_);
    AnaUtil::fillHist1D("mvaoutput_mva3", mvaOutput, puevWt_);
    AnaUtil::fillHist1D("TauTauMass_mva3", diTau.M(), puevWt_);
  }
  if (mvaOutputFK > -0.25) {     
    AnaUtil::fillHist1D("counter_gen", 47, puevWt_);
    AnaUtil::fillHist1D("mvaoutput_mva4", mvaOutput, puevWt_);
    AnaUtil::fillHist1D("TauTauMass_mva4", diTau.M(), puevWt_);
  }
  if (mvaOutputFK > -0.20) {     
    AnaUtil::fillHist1D("counter_gen", 48, puevWt_);
    AnaUtil::fillHist1D("mvaoutput_mva5", mvaOutput, puevWt_);
    AnaUtil::fillHist1D("TauTauMass_mva5", diTau.M(), puevWt_);
  }
  if (mvaOutputFK > -0.15) {     
    AnaUtil::fillHist1D("counter_gen", 49, puevWt_);
    AnaUtil::fillHist1D("mvaoutput_mva6", mvaOutput, puevWt_);
    AnaUtil::fillHist1D("TauTauMass_mva6", diTau.M(), puevWt_);
  }
  if (mvaOutputFK > -0.10) {     
    AnaUtil::fillHist1D("counter_gen", 50, puevWt_);
    AnaUtil::fillHist1D("mvaoutput_mva7", mvaOutput, puevWt_);
    AnaUtil::fillHist1D("TauTauMass_mva7", diTau.M(), puevWt_);
  }

/*
  //Gen Level Filter Information::for SIGNAL MC ONLY  
    AnaUtil::fillHist1D("filter", 1, puevWt_);
  findGenInfo(genMuonList, genTauList);
  if (genMuonList.size() == 1 && genTauList.size() == 2) 
    AnaUtil::fillHist1D("filter", 2, puevWt_);
  if (WhHmhFilter())
    AnaUtil::fillHist1D("filter", 3, puevWt_);
  if (WtmHhhFilter())
    AnaUtil::fillHist1D("filter", 4, puevWt_);
  if (genMuonList.size() == 1 && genTauList.size() == 2 && WtmHhhFilter()) 
    AnaUtil::fillHist1D("filter", 5, puevWt_);
*/

  //After All cut
  AnaUtil::fillHist1D("tau2Pt", taub.pt, puevWt_);
  AnaUtil::fillHist1D("tau1Pt", taua.pt, puevWt_);
  AnaUtil::fillHist1D("muPt", muo.pt, puevWt_);
  AnaUtil::fillHist1D("tau2eta", taub.eta, puevWt_);
  AnaUtil::fillHist1D("tau1eta", taua.eta, puevWt_);
  AnaUtil::fillHist1D("mueta", muo.eta, puevWt_);
  AnaUtil::fillHist1D("mvamet", mvamet.met, puevWt_);
  AnaUtil::fillHist1D("pfmet", pfmet.met, puevWt_);
  AnaUtil::fillHist1D("MuTauMass", zMuTau.M(), puevWt_);
  AnaUtil::fillHist1D("Mt_Mu_Met", mass1, puevWt_);
  AnaUtil::fillHist1D("diTauPt", diTau.Pt(), puevWt_);
  AnaUtil::fillHist1D("TauTauMass", diTau.M(), puevWt_);

/*
  evLog() << "vetoMuon = " << vetoMuon(taua.zvertex, 15, 0.14) 
         << ", vetoEle = " << vetoElectron(taua.zvertex, 10, 0.14)
         << ", met = " << mvamet.met
         << ", diTauPt = " << diTau.Pt()
         << ", MtMuMet = " << mass1
         << ", MuTauM = " << zMuTau.M()
         << endl;

  dumpEvent("11111", evLog(), false);
*/
  if (_skimObj) {
    TreeVariables varList;
    varList.muEta      = M.Eta();
    varList.muPt       = M.Pt();
    varList.tau1Eta    = T1.Eta();
    varList.tau1Pt     = T1.Pt();
    varList.tau2Eta    = T2.Eta();
    varList.tau2Pt     = T2.Pt();
    varList.diTauEta   = diTau.Eta();
    varList.diTaudR   = AnaUtil::deltaR(T1, T2);
    varList.diTauPt    = diTau.Pt();
    varList.diTauMass    = diTau.M();
    varList.dphiMuTau1 = AnaUtil::deltaPhi(M, T1);
    varList.dphiTau1Tau2 = AnaUtil::deltaPhi(T1, T2);
    varList.met        = mvamet.met;
    TVector3 Mu_Alpha(M.Px(), M.Py(), M.Pz());
    TVector3 Z_Alpha(0, 0, 1.0);
    TVector3 d1(T1.Px(), T1.Py(), T2.Pz());
    TVector3 d2(T2.Px(), T2.Py(), T2.Pz());
    TVector3 ditauCrossProd = d1.Cross(d2);
    TVector3 muZ = Z_Alpha.Cross(Mu_Alpha);
    double angle = muZ.Angle(ditauCrossProd);
    double dphi = AnaUtil::deltaPhi(M.Phi(), mvamet.metphi); // acoplanarity
    varList.alpha       = angle; 
    varList.acop        = dphi;
    varList.dphiMuDiTau = AnaUtil::deltaPhi(M, diTau);
    varList.dzTau12Vtx = abs(taua.zvertex - taub.zvertex);
    _skimObj->fill(varList);
  }
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void SelectionEfficiency::endJob() 
{  
  histf()->cd();

  TH1 *h = AnaUtil::getHist1D("counter_gen");
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
  histf()->Write();
  histf()->Close();
  delete histf();

  fLog() << resetiosflags(ios::fixed);

  if (_skimObj) _skimObj->close();

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
bool SelectionEfficiency::readJob(const string& jobFile, int& nFiles)
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
  hmap.insert(pair<string, map<string, double>* >("tau1CutList", &_tau1CutMap));
  hmap.insert(pair<string, map<string, double>* >("tau2CutList", &_tau2CutMap));

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
    string key = tokens.at(0);
    string value = tokens.at(1);

    if (key == "tau1CutList" || key == "tau2CutList")
      AnaUtil::storeCuts(tokens, hmap);
    else if (key == "createMVATree")
      _createMVATree = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "mvaInputFile")
      _mvaInputFile = value;
    else if (key == "readMVAFK")
      _readMVAFK = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "readMVA")
      _readMVA = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "MVAdisFile")
      _MVAdisFile = value;
    else if (key == "MVAFKdisFile")
      _MVAFKdisFile = value;


    tokens.clear();
  }
  // Close the file
  fin.close();

  printJob();

  return true;
}
void SelectionEfficiency::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("tau1CutList", _tau1CutMap));
  hmap.insert(pair<string, map<string, double> >("tau2CutList", _tau2CutMap));
  AnaUtil::showCuts(hmap, os);
}
