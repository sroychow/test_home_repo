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

#include "ETauTau.h"
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
ETauTau::ETauTau()
  : AnaBase(),
  _createMVATree(false),
  _readMVA(false),
  _mvaInputFile(""),
  _skimObj(0)

{}
// ----------
// Destructor
// ----------
ETauTau::~ETauTau() 
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool ETauTau::beginJob() 
{
  AnaBase::beginJob();
  histf()->cd();
  bookHistograms();
  
  // Optionally write selected events to another tree
  if (_createMVATree) _skimObj = new MVASkim(_mvaInputFile);

  reader1 = new TMVA::Reader("Silent");
  reader1->AddVariable("LeadPt", &leadPt);
  reader1->AddVariable("SubPt", &subPt);
  reader1->AddVariable("Met", &met);
  reader1->AddVariable("DeltaRDiTau", &deltaRDiTau);
  reader1->AddVariable("PtRatio", &ptRatio);

  reader1->BookMVA("BDT8","/afs/cern.ch/work/k/kchatter/public/VHTauTau_miniAOD/CMSSW_7_0_7_patch1/src/15Sept2014/WH_BDT8.weights.xml");


  if (_readMVA) {
    reader = new TMVA::Reader("Silent");
    reader->AddVariable("lepPt", &lepPt);
    reader->AddVariable("tauOSPt", &tauOSPt);
    reader->AddVariable("tauSSPt", &tauSSPt);
    reader->AddVariable("met", &met);

    reader->BookMVA("MLP", "/afs/cern.ch/work/k/kchatter/public/VHTauTau_miniAOD/CMSSW_7_0_7_patch1/src/15Sept2014/TMVAClassification_MLP.weights.xml");
  }
  
  return true;
}
// ---------------
// BookB histograms
// ---------------
void ETauTau::bookHistograms() 
{
  new TH1D("evcounter", "Selected event counter", 43, -0.5, 42.5);
  //new TH1D("evcounter_clone", "Selected event counter for pfIsolation", 38, -0.5, 37.5);

  new TH1F("EGen_eta", " #eta of gen level electron if  ElectronColl() size is zero", 60, -3, 3);
  new TH1F("electronPt_afterEselct", " Pt of electron after Electron selection ", 24, 0, 240);
  new TH1F("electronEta_afterEselct", " #eta of electron after Electron selection ", 12, -3, 3);

  new TH1F("OSTauPt_afterTauSelect", " #p_T of OS #tau after Tau Selection", 240, 0, 240);
  new TH1F("SSTauPt_afterTauSelect", " #p_T of SS #tau after Tau Selection", 240, 0, 240);
  new TH1F("OSTauEta_afterTauSelect", " #eta of OS #tau after Tau Selection", 24, -3, 3);
  new TH1F("SSTauEta_afterTauSelect", " #eta of SS #tau after Tau Selection", 24, -3, 3);

  new TH1F("electronPt_final", " #p_T of electron after all cuts", 48, 0, 240);
  new TH1F("leadTauPt_final", " #p_T of leading #tau after all cuts", 48, 0, 240);
  new TH1F("subleadTauPt_final", "#p_T of subleading #tau after all cuts", 48, 0, 240);
  new TH1F("electronEta_final", " #eta of electron after after all cuts ", 24, -3, 3);
  new TH1F("leadTauEta_final", " #eta of leading #tau after all cuts", 24, -3, 3);
  new TH1F("subleadTauEta_final", " #eta of subleading #tau after all cuts ", 24, -3, 3);
  new TH1F("electronPhi_final", " #phi of electron after after all cuts ", 65, 0, 6.5);
  new TH1F("leadTauPhi_final", " #phi of leading #tau   after after all cuts ", 65, 0, 6.5);
  new TH1F("subleadTauPhi_final", " #phi of subleading #tau after all cuts ", 65, 0, 6.5);

  new TH1F("diTauMass_afterTau2select"," diTau mass after #tau selection ", 32, 0, 320);
  new TH1F("diTauMass_oppositeCharge"," diTau mass of opposite charge #tau s ", 32, 0, 320);
  new TH1F("diTauMass_afterZtoeeVeto"," diTau mass of #tau s after | #m_z - #m_{eTau}| cut ", 32, 0, 320);
  new TH1F("diTauMass_afterZtoTauTauVeto"," diTau mass of #tau s after Z->#tau#tau veto ", 32, 0, 320);
  new TH1F("diTauMass_afterbjetVeto"," diTau mass of #tau s after bjet veto ", 32, 0, 320);

  new TH1F("lt_afterJetVeto","lepton energy after bJet veto", 30, 0, 300);
  new TH1F("lt_final","lepton energy after all cuts", 30, 0, 300);
  new TH1F("mvaMet_final"," MVA MET after all cuts", 14, 0, 140);
  new TH1F("pfMet_final"," PF MET after all cuts", 14, 0, 140);

  new TH1F("MtEleMet_final","transverse mass plot for Tau1 and MET", 14, 0, 140);
  new TH1F("diTauMass_final","diTau mass after final selection", 32, 0, 320);
  new TH1F("diTauPt_final", " #p_T of diTau after all cuts", 48, 0, 240);
  new TH1F("diTauEta_final","#eta distribution of diTau  after final selection", 24, -3, 3);
  new TH1F("diTauPhi_final","#phi distribution  of diTau after final selection", 65, 0, 6.5);

  new TH1F("pfRelIsoDB","distribution of pfRelIso DeltaBeta Corrected", 100, 0, 0.5);
  new TH1F("pfIso","distribution of pfIso: electron.sumChargedHadronPt/electron.pt", 100, 0, 0.5);
  new TH1F("el_trkD0","distribution of trackD0 of electron", 1000, -5, 5);
  new TH1F("el_trkdz","distribution of trackdz of electron", 100, -5, 5);
  new TH1F("el_vtxDistZ","distribution of vtxDistZ  of electron", 100, -5, 5);

  /// for test only

  new TH1F("GenMatch_TOS_fail_DMF","pdg Id of #tau^{os} candidates, rejected due to DecayModeFinding cut", 52, -1.5, 50.5);
  new TH1F("GenMatch_TOS_fail_antiEle","pdg Id of #tau^{os} candidates, rejected due to AgainstElectronTight cut", 52, -1.5, 50.5);
  new TH1F("GenMatch_TOS_fail_charge","pdg Id of #tau^{os} candidates, rejected due to charge cut", 52, -1.5, 50.5);

  new TH1F("GenMatch_TSS_fail_DMF","pdg Id of #tau^{ss} candidates, rejected due to DecayModeFinding cut", 52, -1.5, 50.5);
  new TH1F("GenMatch_TSS_fail_antiEle","pdg Id of #tau^{ss} candidates, rejected due to AgainstElectronTight cut", 52, -1.5, 50.5);
  new TH1F("GenMatch_TSS_fail_charge","pdg Id of #tau^{ss} candidates, rejected due to charge cut", 52, -1.5, 50.5);



  new TH1F("GenMatch_TOS_pass_DMF","pdg Id of #tau^{os} candidates,pass the DecayModeFinding cut", 52, -1.5, 50.5);
  new TH1F("GenMatch_TOS_pass_antiEle","pdg Id of #tau^{os} candidates,pass the AgainstElectronTight cut", 52, -1.5, 50.5);
  new TH1F("GenMatch_TOS_pass_charge","pdg Id of #tau^{os} candidates,pass the charge cut", 52, -1.5, 50.5);

  new TH1F("GenMatch_TSS_pass_DMF","pdg Id of #tau^{ss} candidates,pass the DecayModeFinding cut", 52, -1.5, 50.5);
  new TH1F("GenMatch_TSS_pass_antiEle","pdg Id of #tau^{ss} candidates,pass the AgainstElectronTight cut", 52, -1.5, 50.5);
  new TH1F("GenMatch_TSS_pass_charge","pdg Id of #tau^{ss} candidates,pass the charge cut", 52, -1.5, 50.5);

  new TH1F("GenMatch_afterEselected","pdg Id of the object selected as electrons, after electron selection", 52, -1.5, 50.5);

  new TH1F("GenMatch_Electron_final","pdg Id of the object selected as electrons, after final event selection", 52, -1.5, 50.5);
  new TH1F("GenMatch_TOS_final","pdg Id of the object selected as #tau^{os}, after final event selection", 52, -1.5, 50.5);
  new TH1F("GenMatch_TSS_final","pdg Id of the object selected as #tau^{ss}, after final event selection", 52, -1.5, 50.5);

#if 0
  new TH1F("tauOSdmf_beforeCut","distribution DMF of OS Tau before dmf cut", 11, -0.05, 1.05);
  new TH1F("tauOSdmf_afterCut","distribution DMF of OS Tau after dmf cut", 11, -0.05, 1.05);
  new TH1F("tauSSdmf_beforeCut","distribution DMF of SS Tau before dmf cut", 11, -0.05, 1.05);
  new TH1F("tauSSdmf_afterCut","distribution DMF of SS Tau after dmf cut", 11, -0.05, 1.05);
#endif
  /////


  new TH1F("dr_EGen_Eselected","difference between two object", 1000, 0, 10);
  new TH1F("dr_TauOSGen_Eselected","difference between two object", 1000, 0, 10);
  new TH1F("dr_TauSSGen_Eselected","difference between two object", 1000, 0, 10);

  new TH1F("dr_E_EGen","difference between two object", 1000, 0, 10);
  new TH1F("dr_TSS_EGen","difference between two object", 1000, 0, 10);
  new TH1F("dr_TOS_EGen","difference between two object", 1000, 0, 10);
  new TH1F("dr_E_TauOSGen","difference between two object", 1000, 0, 10);
  new TH1F("dr_TSS_TauOSGen","difference between two object", 1000, 0, 10);
  new TH1F("dr_TOS_TauOSGen","difference between two object", 1000, 0, 10);
  new TH1F("dr_E_TauSSGen","difference between two object", 1000, 0, 10);
  new TH1F("dr_TSS_TauSSGen","difference between two object", 1000, 0, 10);
  new TH1F("dr_TOS_TauSSGen","difference between two object", 1000, 0, 10);


#if 0
  for(int i = 1; i <= 22; ++i){
    std::string histo_name = "mvaFk_" + itoa(i);
    std::string histo_title = "fake mva distribution after MVA " + itoa(i);

     new TH1F("mvaFk_1"," fake mva distribution after MVA 1", 100, -1, 1);
  }

  new TH1F("mvaFk_1"," fake mva distribution after MVA 1", 100, -1, 1);
  new TH1F("ditauMass_mva1"," ditau mass distributio after MVA 1", 320, 0, 320);

  new TH1F("mvaFk_2"," fake mva distribution after MVA 2", 100, -1, 1);
  new TH1F("ditauMass_mva2"," ditau mass distributio after MVA 2", 320, 0, 320);

  new TH1F("mvaFk_3"," fake mva distribution after MVA 3", 100, -1, 1);
  new TH1F("ditauMass_mva3"," ditau mass distributio after MVA 3", 320, 0, 320);

  new TH1F("mvaFk_4"," fake mva distribution after MVA 4", 100, -1, 1);
  new TH1F("ditauMass_mva4"," ditau mass distributio after MVA 4", 320, 0, 320);

  new TH1F("mvaFk_5"," fake mva distribution after MVA 5", 100, -1, 1);
  new TH1F("ditauMass_mva5"," ditau mass distributio after MVA 5", 320, 0, 320);

  new TH1F("mvaFk_6"," fake mva distribution after MVA 6", 100, -1, 1);
  new TH1F("ditauMass_mva6"," ditau mass distributio after MVA 6", 320, 0, 320);

  new TH1F("mvaFk_7"," fake mva distribution after MVA 7", 100, -1, 1);
  new TH1F("ditauMass_mva7"," ditau mass distributio after MVA 7", 320, 0, 320);

  new TH1F("mvaFk_8"," fake mva distribution after MVA 8", 100, -1, 1);
  new TH1F("ditauMass_mva8"," ditau mass distributio after MVA 7", 320, 0, 320);

  new TH1F("mvaFk_9"," fake mva distribution after MVA 9", 100, -1, 1);
  new TH1F("ditauMass_mva9"," ditau mass distributio after MVA 9", 320, 0, 320);

  new TH1F("mvaFk_10"," fake mva distribution after MVA 10", 100, -1, 1);
  new TH1F("ditauMass_mva10"," ditau mass distributio after MVA 10", 320, 0, 320);

  new TH1F("mvaFk_11"," fake mva distribution after MVA 11", 100, -1, 1);
  new TH1F("ditauMass_mva11"," ditau mass distributio after MVA 11", 320, 0, 320);

  new TH1F("mvaFk_12"," fake mva distribution after MVA 12", 100, -1, 1);
  new TH1F("ditauMass_mva12"," ditau mass distributio after MVA 12", 320, 0, 320);

  new TH1F("mvaFk_13"," fake mva distribution after MVA 13", 100, -1, 1);
  new TH1F("ditauMass_mva13"," ditau mass distributio after MVA 13", 320, 0, 320);

  new TH1F("mvaFk_14"," fake mva distribution after MVA 14", 100, -1, 1);
  new TH1F("ditauMass_mva14"," ditau mass distributio after MVA 14", 320, 0, 320);

  new TH1F("mvaFk_15"," fake mva distribution after MVA 15", 100, -1, 1);
  new TH1F("ditauMass_mva15"," ditau mass distributio after MVA 15", 320, 0, 320);

  new TH1F("mvaFk_16"," fake mva distribution after MVA 16", 100, -1, 1);
  new TH1F("ditauMass_mva16"," ditau mass distributio after MVA 16", 320, 0, 320);

  new TH1F("mvaFk_17"," fake mva distribution after MVA 17", 100, -1, 1);
  new TH1F("ditauMass_mva17"," ditau mass distributio after MVA 17", 320, 0, 320);

  new TH1F("mvaFk_18"," fake mva distribution after MVA 18", 100, -1, 1);
  new TH1F("ditauMass_mva18"," ditau mass distributio after MVA 18", 320, 0, 320);

  new TH1F("mvaFk_19"," fake mva distribution after MVA 19", 100, -1, 1);
  new TH1F("ditauMass_mva19"," ditau mass distributio after MVA 19", 320, 0, 320);

  new TH1F("mvaFk_20"," fake mva distribution after MVA 20", 100, -1, 1);
  new TH1F("ditauMass_mva20"," ditau mass distributio after MVA 20", 320, 0, 320);

  new TH1F("mvaFk_21"," fake mva distribution after MVA 21", 100, -1, 1);
  new TH1F("ditauMass_mva21"," ditau mass distributio after MVA 21", 320, 0, 320);

  new TH1F("mvaFk_22"," fake mva distribution after MVA 22", 100, -1, 1);
  new TH1F("ditauMass_mva22"," ditau mass distributio after MVA 22", 320, 0, 320);

  new TH1F("mvaFk_BL"," fake mva distribution after basline selection", 100, -1, 1);
  new TH1F("mvaVZ_BL"," VZ mva distribution after basline selection", 100, -1, 1);
  new TH1F("ditauMass_BL"," ditau mass distributio after BL ", 320, 0, 320);
  new TH1F("mva_evcounter"," event counter for each MVA CUT", 23, -0.5, 22.5);

#endif
}

// -------------------
// The main event loop
// -------------------
void ETauTau::clearLists() {
  vtxList.clear();
  bjetList.clear();
  trigObjList.clear();
  genEleList.clear(); // for Gen Level test
  genTauList.clear(); // for Gen Level test
  genList.clear(); // for GenLevel Matching
}
void ETauTau::eventLoop() 
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

      //////// should be used after calling findgenInfo 
      //      fLog() << "=> Event " << event
      //	     << " Lumis " << lumis
      //	     << " Run " << run
      //	     << endl;
      //dumpGenInfo(fLog());
      ///////

      findGenInfo(11, genEleList, genTauList);
      bool genflag = (genEleList.size() == 1 && genTauList.size() == 2); // for test (gen level pdgid matching)
      if (!genflag)  continue;
      genList.push_back(genEleList.at(0));
      genList.push_back(genTauList.at(0));
      genList.push_back(genTauList.at(1));

      float elecharge = genEleList.at(0).charge;
      EGen.SetPtEtaPhiE(genEleList.at(0).pt, genEleList.at(0).eta, genEleList.at(0).phi, genEleList.at(0).energy);
      if(genTauList.at(0).charge + elecharge == 0){
	TauOSGen.SetPtEtaPhiE(genTauList.at(0).pt, genTauList.at(0).eta, genTauList.at(0).phi, genTauList.at(0).energy);
	TauSSGen.SetPtEtaPhiE(genTauList.at(1).pt, genTauList.at(1).eta, genTauList.at(1).phi, genTauList.at(1).energy);
      }
      else {
	TauOSGen.SetPtEtaPhiE(genTauList.at(1).pt, genTauList.at(1).eta, genTauList.at(1).phi, genTauList.at(1).energy);
	TauSSGen.SetPtEtaPhiE(genTauList.at(0).pt, genTauList.at(0).eta, genTauList.at(0).phi, genTauList.at(0).energy);
      } /// for test (gen level pdgid matching)

       
#if 0
      fLog() << "=> Event " << event
	     << " Lumis " << lumis
	     << " Run " << run
	     << endl;
      fLog() << "=> nEle: " << genEleList.size()
	     << " nTau: " << genTauList.size()
	     << endl;
      dumpGenInfo(fLog());
#endif
    }


    AnaUtil::fillHist1D("evcounter", 1, puevWt_);


    // Trigger selection
    //if (useTrigger() && !isTriggered()) continue;
    AnaUtil::fillHist1D("evcounter", 2, puevWt_);
    //dumpTriggerPaths(fLog(), true);
    //dumpTriggerObjectInfo(trigObjList, fLog());

    op.verbose = (logOption() >> 1 & 0x1); 
    findVtxInfo(vtxList, op, fLog());
    // Event Selection Starts here .....
    // presence of > 1 good vertex
    if (vtxList.size() < 1) continue;
    AnaUtil::fillHist1D("evcounter", 3, puevWt_);

    op.verbose = (logOption() >> 5 & 0x1);
    //findJetInfo(bjetList, op, fLog());
    //findTriggerObjectInfo(trigObjList);
    selectEvent();
  }  
  // Analysis is over
  endJob();
}
void ETauTau::selectEvent() {
  study_eff();
}
void ETauTau::study_eff()
{
  double vz = vtxList.at(0).z;

  //////////////////////////////////////////////
  /////////      ELECTRON SELECTION ///////////
  /////////////////////////////////////////////
  //  TLorentzVector E1;
  int ecounters[] = {0,0,0,0,0,0,0,0};
  vector<Electron> fEleList;
  int electron_indx;
  int indx = 0;

  if(electronColl()->size()==0) {
    AnaUtil::fillHist1D("EGen_eta",EGen.Eta());
    return;
  }
  AnaUtil::fillHist1D("evcounter", 4, puevWt_);

  for (auto it = electronColl()->begin(); it != electronColl()->end(); ++it,++indx) {
    const Electron&  electron = (*it);

    if (electron.pt <= AnaUtil::cutValue(electronCutMap(), "pt")) continue;
    ++ecounters[0];

    double eleta = std::fabs(electron.eta);
    if ((eleta >= AnaUtil::cutValue(electronCutMap(), "etaLow") &&
         eleta <= AnaUtil::cutValue(electronCutMap(), "etaUp")) ||
         eleta >= AnaUtil::cutValue(electronCutMap(), "eta") ) continue;
    ++ecounters[1];

    //  if (electron.trkHits < AnaUtil::cutValue(electronCutMap(),"trkHits")) continue;
    ++ecounters[2];
    //if (electron.pixHits < AnaUtil::cutValue(electronCutMap(),"pixHits")) continue;
    ++ecounters[3];

    double electronIso = electron.sumChargedHadronPt/electron.pt;
    AnaUtil::fillHist1D("pfIso", electronIso, puevWt_);
    AnaUtil::fillHist1D("pfRelIsoDB", electron.pfRelIso, puevWt_);

#if 0
    bool quality_EB_loose = electron.hasGsfTrack
                         && (fabs(electron.eta) <= 1.4442
                         &&  electron.sigmaEtaEta < 0.01
                         &&  electron.deltaEtaTrkSC < 0.007
                         &&  electron.deltaPhiTrkSC < 0.8
			 &&  electron.hoe < 0.15);
    bool quality_EE_loose = electron.hasGsfTrack
                         && (fabs(electron.eta) >= 1.566
                         &&  electron.sigmaEtaEta < 0.03
                         &&  electron.deltaEtaTrkSC < 0.01
                         &&  electron.deltaPhiTrkSC < 0.7
			 &&  electron.hoe < 0.07);
    bool quality_loose = quality_EB_loose || quality_EE_loose;
    if(!quality_loose) continue;
#endif
    if(!eleId(electron,25)) continue;
    ++ecounters[4];
    //if (std::fabs(electron.dxyPV) >= AnaUtil::cutValue(electronCutMap(),"trkD0")) continue;
    //////////////////////////////////////////////////////////
    AnaUtil::fillHist1D("el_vtxDistZ", electron.vtxDistZ, puevWt_);
    AnaUtil::fillHist1D("el_trkdz", electron.trkDz, puevWt_);
    AnaUtil::fillHist1D("el_trkD0", electron.trkD0, puevWt_);
    //////////////////////////////////////////////////////////

    //if (std::fabs(dxy) >= AnaUtil::cutValue(electronCutMap(),"trkD0")) continue;

    if (std::fabs(electron.trkD0) >= AnaUtil::cutValue(electronCutMap(),"trkD0")) continue;
    ++ecounters[5];

#if 0
    bool isGoodVtx;
    TVector3 vele = findLeptonVtx(electron.vtxIndex, isGoodVtx);
    double dz = vele.z() - vz;
    if (!isGoodVtx || std::fabs(dz) >= AnaUtil::cutValue(electronCutMap(), "dz")) continue;
#endif
    double dz  = electron.vtxDistZ;
    if (std::fabs(dz) >= AnaUtil::cutValue(electronCutMap(), "dz")) continue;

    //if (std::fabs(electron.trkDz) >= AnaUtil::cutValue(electronCutMap(), "dz")) continue;
    ++ecounters[6];

    //if ( !AnaBase::electronMVA(electron) ) continue;
    // ++ecounters[6];

    //double electronIso = electron.sumChargedHadronPt/electron.pt;
    //AnaUtil::fillHist1D("pfIso", electronIso, puevWt_);
    //double electronIso = electron.sumChargedHadronPt;
    // bool electronIsoBarrel = ( eleta < AnaUtil::cutValue(electronCutMap(), "etaLow") &&
    //                               electronIso < AnaUtil::cutValue(electronCutMap(), "relIsoBarrel") );
  //bool electronIsoEndCap = ( eleta > AnaUtil::cutValue(electronCutMap(), "etaUp") && 
    //                         electronIso < AnaUtil::cutValue(electronCutMap(), "relIsoEndCap"));
    //if (!(electronIsoBarrel || electronIsoEndCap)) continue;
    if (electron.pfRelIso > 0.2) continue;
    ++ecounters[7];

    fEleList.push_back(electron);
    electron_indx = indx;
    break;
  }
  int ishift = 5;
  for (size_t i = 0; i < NEL(ecounters); ++i) {
    if (ecounters[i] > 0) AnaUtil::fillHist1D("evcounter", ishift + i, puevWt_);
  }

  // atleast 1 good electron
  if (fEleList.size() < 1) return;
  const Electron& ele = fEleList.at(0);
  AnaUtil::fillHist1D("evcounter", 13, puevWt_);
  AnaUtil::fillHist1D("electronPt_afterEselct", ele.pt, puevWt_);
  AnaUtil::fillHist1D("electronEta_afterEselct", ele.eta, puevWt_);

  TLorentzVector E1;
  E1.SetPtEtaPhiE(ele.pt, ele.eta, ele.phi, ele.energy);
  //  AnaUtil::fillHist1D("dr_EGen_Eselected", E1.DeltaR(EGen)); // for MC only
  //  AnaUtil::fillHist1D("dr_TauOSGen_Eselected", E1.DeltaR(TauOSGen)); // for MC only
  //  AnaUtil::fillHist1D("dr_TauSSGen_Eselected", E1.DeltaR(TauSSGen)); // for MC only

  AnaUtil::fillHist1D("GenMatch_afterEselected", GenLevelMatching(E1, genList));

  //---------------------------------------------
  //             OS Tau SELECTION
  //---------------------------------------------

  auto it = metColl()->begin();
  const MET& missET = (*it);
  double Met = missET.met;
  double MetPhi = missET.metphi;

  vector<Tau> OStauList;
  vector<Tau>  SStauList;

  int ostcounters[] = {0,0,0,0,0,0,0,0,0};
  for (auto it = tauColl()->begin(); it != tauColl()->end(); ++it) {
    const Tau& tau = (*it);
    //if (!tau) continue;
    TLorentzVector T;
    T.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
    int pdg_Id = GenLevelMatching(T, genList); // for gen level testing
    double dr = E1.DeltaR(T);

    if (dr <= AnaUtil::cutValue(_evselCutMap, "drEleTau")) continue;
    ++ostcounters[0];

    //if (taua.pt <= AnaUtil::cutValue(_evselCutMap, "ptTau1")) continue;
    if (tau.pt <= 20) continue;
    ++ostcounters[1];
     
    if( std::fabs(tau.eta) >= AnaUtil::cutValue(_evselCutMap,"etaTau1")) continue;
    ++ostcounters[2];

    if (tau.decayModeFinding != 1.0) {
      AnaUtil::fillHist1D("GenMatch_TOS_fail_DMF", pdg_Id);
      continue;
    }
    ++ostcounters[3];
    AnaUtil::fillHist1D("GenMatch_TOS_pass_DMF", pdg_Id);

    if(tau.chargedIsoPtSum >= 2) continue;
    ++ostcounters[4];

    //if (tau.againstMuonTight <= 0.5) continue;
    if (tau.againstMuonTight3 <= 0.5) continue;
    ++ostcounters[5];


    //    if (tau.againstElectronLoose <= 0.5) { 
    if (tau.againstElectronLooseMVA5 <= 0.5) { 
      AnaUtil::fillHist1D("GenMatch_TOS_fail_antiEle", pdg_Id);
      continue;
    }
    ++ostcounters[6];
    AnaUtil::fillHist1D("GenMatch_TOS_pass_antiEle", pdg_Id);

    //if (tau.vtxDz >= 0.2) continue;
    if (fabs(tau.zvertex - vz ) >= 0.2) continue;
    ++ostcounters[7];

    if((tau.charge + ele.charge) !=0 ) {
      AnaUtil::fillHist1D("GenMatch_TOS_fail_charge", pdg_Id);
      continue;  /// ACTUAL PLACE OF THIS CUT.. MOVE IT ABOVE FOR TESTING ONLYk
    }
    ++ostcounters[8];
    AnaUtil::fillHist1D("GenMatch_TOS_pass_charge", pdg_Id);

    OStauList.push_back(tau);
    AnaUtil::fillHist1D("OSTauPt_afterTauSelect", tau.pt);
    AnaUtil::fillHist1D("OSTauEta_afterTauSelect", tau.eta);
    break;
  }
  ishift = 14;
  for (size_t i = 0; i < NEL(ostcounters); ++i) {
    if (ostcounters[i]) AnaUtil::fillHist1D("evcounter", ishift + i, puevWt_);
  }

  if(OStauList.size()<1) return;
  AnaUtil::fillHist1D("evcounter", 23, puevWt_);
  const Tau& tauOS = OStauList.at(0);
  TLorentzVector TOS;
  TOS.SetPtEtaPhiE(tauOS.pt, tauOS.eta, tauOS.phi, tauOS.energy);


  //---------------------------------------------
  //             SS Tau SELECTION
  //---------------------------------------------
  int sstcounters[] = {0,0,0,0,0,0,0,0,0,0};
  for (auto it = tauColl()->begin(); it != tauColl()->end(); ++it) {
    const Tau& tau = (*it);
    TLorentzVector T;
    T.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
    int pdg_Id = GenLevelMatching(T, genList);
    double dr = E1.DeltaR(T);

    if (dr <= AnaUtil::cutValue(_evselCutMap, "drEleTau")) continue;
    ++sstcounters[0];

    double drtt = T.DeltaR(TOS);
    if (drtt <= AnaUtil::cutValue(_evselCutMap, "drTauTau")) continue;
    ++sstcounters[1];

    if (tau.pt <= 20) continue;
    ++sstcounters[2];

    if( std::fabs(tau.eta) >= AnaUtil::cutValue(_evselCutMap,"etaTau1")) continue;
    ++sstcounters[3];

    if (tau.decayModeFinding != 1.0) {
      AnaUtil::fillHist1D("GenMatch_TSS_fail_DMF", pdg_Id);
      continue;
    }
    ++sstcounters[4];
    AnaUtil::fillHist1D("GenMatch_TSS_pass_DMF", pdg_Id);
    

    if(tau.chargedIsoPtSum >= 1) continue;
    ++sstcounters[5];

    //if (tau.againstMuonTight <= 0.5) continue;
    if (tau.againstMuonTight3 <= 0.5) continue;
    ++sstcounters[6];


    //    if (tau.againstElectronLoose <= 0.5) {
    if (tau.againstElectronLooseMVA5 <= 0.5) {
          AnaUtil::fillHist1D("GenMatch_TSS_fail_antiEle", pdg_Id);
	  continue;
    }
    ++sstcounters[7];
    AnaUtil::fillHist1D("GenMatch_TSS_pass_antiEle", pdg_Id);

    if (fabs(tau.zvertex - vz) >= 0.2) continue;
    ++sstcounters[8];
    if((tau.charge + ele.charge) ==0 ) {
      AnaUtil::fillHist1D("GenMatch_TSS_fail_charge", pdg_Id);
      continue;
    }
    ++sstcounters[9];
    AnaUtil::fillHist1D("GenMatch_TSS_pass_charge", pdg_Id);


    SStauList.push_back(tau);
    AnaUtil::fillHist1D("SSTauPt_afterTauSelect", tau.pt);
    AnaUtil::fillHist1D("SSTauEta_afterTauSelect", tau.eta);
    break;
  }
  ishift = 24;
  for (size_t i = 0; i < NEL(sstcounters); ++i) {
    if (sstcounters[i]) AnaUtil::fillHist1D("evcounter", ishift + i, puevWt_);
  }
  
  /////////////////////////////////////////////////
  
  if (SStauList.size() < 1) return;
  AnaUtil::fillHist1D("evcounter", 34, puevWt_);
  
  
  const Tau& tauSS = SStauList.at(0);
  
  TLorentzVector TSS;
  TSS.SetPtEtaPhiE(tauSS.pt, tauSS.eta, tauSS.phi, tauSS.energy);
  TLorentzVector diTau = TOS + TSS;
  double diTauMass = diTau.M();
  
  ///// TOPOLOGICAL CUTS //////    
  AnaUtil::fillHist1D("diTauMass_oppositeCharge", diTauMass, puevWt_);
  
  double highestPt = std::max(tauOS.pt, tauSS.pt);
  if (highestPt <= 25) return;
  AnaUtil::fillHist1D("evcounter", 35, puevWt_);
  
  // b-tagged jet Veto
  //  if (bjetList.size() > 0) return;
  AnaUtil::fillHist1D("evcounter", 36, puevWt_);
  AnaUtil::fillHist1D("diTauMass_aftrbjetVeto", diTauMass, puevWt_);    


  AnaUtil::fillHist1D("GenMatch_Electron_final", GenLevelMatching(E1, genList));  /////////////////
  AnaUtil::fillHist1D("GenMatch_TOS_final", GenLevelMatching(TOS,genList));       ////////
  AnaUtil::fillHist1D("GenMatch_TSS_final", GenLevelMatching(TSS,genList));       ///////

  AnaUtil::fillHist1D("dr_E_EGen", E1.DeltaR(EGen), puevWt_);                     //////
  AnaUtil::fillHist1D("dr_TSS_EGen", TSS.DeltaR(EGen), puevWt_);
  AnaUtil::fillHist1D("dr_TOS_EGen", TOS.DeltaR(EGen), puevWt_);                  ////// for gen level testing
  AnaUtil::fillHist1D("dr_TSS_TauOSGen", TSS.DeltaR(TauOSGen), puevWt_);
  AnaUtil::fillHist1D("dr_TOS_TauOSGen", TOS.DeltaR(TauOSGen), puevWt_);         ///////
  AnaUtil::fillHist1D("dr_E_TauOSGen", E1.DeltaR(TauSSGen), puevWt_);
  AnaUtil::fillHist1D("dr_E_TauSSGen", E1.DeltaR(TauSSGen), puevWt_);            ////////
  AnaUtil::fillHist1D("dr_TSS_TauSSGen", TSS.DeltaR(TauSSGen), puevWt_);         ////////
  AnaUtil::fillHist1D("dr_TOS_TauSSGen", TOS.DeltaR(TauSSGen), puevWt_);         ///////////////////
  


  
#if 1
  if (_skimObj) {
    
    //double MtEleMet = std::sqrt(pow(E1.Pt() + mvaMet,2) - pow(E1.Pt()*cos(E1.Phi()) + mvaMet*cos(MetPhi),2) - 
    //                           pow(E1.Pt()*sin(E1.Phi()) + mvaMet*sin(MetPhi),2));
    //if (MtEleMet < 20) return;
    TreeVariables varList;
    varList.lepEta       = ele.eta;
    varList.lepPt        = ele.pt;
    varList.tauOSEta     = TOS.Eta();
    varList.tauOSPt      = TOS.Pt();
    varList.tauSSEta     = TSS.Eta();
    varList.tauSSPt      = TSS.Pt();
    varList.diTauEta     = diTau.Eta();
    varList.diTauPt      = diTau.Pt();
    varList.dphilepTauOS = AnaUtil::deltaPhi(E1, TOS);
    // varList.met       = mvaMet;
    varList.met          = Met;

    TVector3 lep(E1.Px(), E1.Py(), E1.Pz());
    TVector3 Z_Alpha(0, 0, 1.0);
    TVector3 d1(TOS.Px(), TOS.Py(), TOS.Pz());
    TVector3 d2(TSS.Px(), TSS.Py(), TSS.Pz());
    TVector3 ditauCrossProd = d1.Cross(d2);      
    TVector3 lepZ = Z_Alpha.Cross(lep);
    double   angle = lepZ.Angle(ditauCrossProd);

    varList.alpha        = angle;
    varList.alphatt      = d1.Angle(d2);
    varList.alphalditau  = lep.Angle(d1+d2);
    varList.alphalprod   = lep.Angle(ditauCrossProd);
    varList.acop         = AnaUtil::deltaPhi(E1.Phi(), MetPhi);
    varList.dphilepDiTau = AnaUtil::deltaPhi(E1, diTau);
    varList.DeltaRDiTau  = TOS.DeltaR(TSS);
    varList.PtRatio      = (diTau.Pt())/(TOS.Pt() + TSS.Pt());
    //varList.mass       = diTauMass;
    
    _skimObj->fill(varList);
    return;
  }
#endif
  
  // mu veto
  if (AnaBase::vetoMuon(15, 0.2) > 0) return;
  AnaUtil::fillHist1D("evcounter", 37, puevWt_);
  
  // e veto
  if (AnaBase::vetoElectron(15, 0.2) > 1) return;
  AnaUtil::fillHist1D("evcounter", 38, puevWt_);
  
#if 0
  // z -> ee veto
  double tauOSemfrac = tauOS.emFraction;
  bool ztoee = AnaBase::DYtoeeVeto(TOS, TSS, tauOSemfrac, ele, electron_indx);
  if (ztoee == true) return;
#endif
  AnaUtil::fillHist1D("evcounter", 39, puevWt_);
  
  AnaUtil::fillHist1D("diTauMass_afterZtoeeVeto", diTauMass, puevWt_);
  
  double massETauOS = ( E1 + TOS ).M();
  double massETauSS = ( E1 + TSS ).M();
  
  //if (Met < 40 && massETauOS < 90 && massETauSS < 90) return; // z to tau tau veto 
  AnaUtil::fillHist1D("diTauMass_afterZtoTauTauVeto", diTauMass, puevWt_);
  AnaUtil::fillHist1D("evcounter", 40, puevWt_);
  
  double MtEleMet = std::sqrt(pow(E1.Pt() + Met,2) - pow(E1.Pt()*cos(E1.Phi()) + Met*cos(MetPhi),2) - 
                                                     pow(E1.Pt()*sin(E1.Phi()) + Met*sin(MetPhi),2));   
  
  if (MtEleMet < 20) return;
  AnaUtil::fillHist1D("evcounter", 41, puevWt_);
  
  // dummy (was for mvaMet)
  AnaUtil::fillHist1D("evcounter", 42, puevWt_);
  
  //dumpEvent("01010", _evLog, false);
  //AnaUtil::fillHist1D("MtEleMet_final", MtEleMet, puevWt_);
  AnaUtil::fillHist1D("diTauMass_final", diTauMass,   puevWt_);
  AnaUtil::fillHist1D("diTauEta_final",  diTau.Eta(), puevWt_);
  AnaUtil::fillHist1D("diTauPt_final",   diTau.Pt(),  puevWt_);
  AnaUtil::fillHist1D("diTauPhi_final",  diTau.Phi(), puevWt_);
  //AnaUtil::fillHist1D("mvaMet_final",  mvaMet,      puevWt_);
  //AnaUtil::fillHist1D("pfMet_final",   pfMet,       puevWt_);
  // AnaUtil::fillHist1D("lt_final",     lt,          puevWt_);
  AnaUtil::fillHist1D("electronPt_final",    ele.pt,  puevWt_);
  AnaUtil::fillHist1D("electronEta_final",   ele.eta, puevWt_);
  AnaUtil::fillHist1D("electronPhi_final",   ele.phi, puevWt_);
  AnaUtil::fillHist1D("leadTauPt_final",     TOS.Pt(),  puevWt_);
  AnaUtil::fillHist1D("leadTauEta_final",    TOS.Eta(), puevWt_);
  AnaUtil::fillHist1D("leadTauPhi_final",    TOS.Phi(), puevWt_);
  AnaUtil::fillHist1D("subleadTauPt_final",  TSS.Pt(), puevWt_);
  AnaUtil::fillHist1D("subleadTauEta_final", TSS.Eta(), puevWt_);
  AnaUtil::fillHist1D("subleadTauPhi_final", TSS.Phi(), puevWt_);
  
  /*
  //Armin's MVA 
  double mvaOutputFK = -999; 
  leadPt      = TOS.Pt();
  subPt       = TSS.Pt();
  met         = mvaMet;
  deltaRDiTau = AnaUtil::deltaR(TOS, TSS);
  ptRatio     = diTau.Pt()/(TOS.Pt() + TSS.Pt());
  
  mvaOutputFK = reader1->EvaluateMVA("BDT8");
  
  double mvaOutput = -999;      //Set at high negative value, as MvaOutput > -0.106
  if (_readMVA){
  lepPt        = ele.pt;
  tauOSPt      = TOS.Pt();
  tauSSPt      = TSS.Pt();
  //dphilepDiTau = AnaUtil::deltaPhi(E1, diTau);
  //PtRatio      = (diTau.Pt())/(TOS.Pt() + TSS.Pt());
  //mass       = diTauMass;
  
  mvaOutput = reader->EvaluateMVA("MLP");
  }
  */
#if 0
  _evLog << setw(15) <<  "ditauMass: " << setw(8) << diTauMass 
	 << setw(15) << "diTauEta: "   << setw(8) << diTau.Eta() 
	 << setw(15) <<  "diTauPt: "   << setw(8) << diTau.Pt() 
	 << setw(15) << "diTauPhi: "   << setw(8) << diTau.Phi() << endl;
  
  double d3d = ele.vtxDist3D;
  double dz  = ele.vtxDistZ;
  double dxy = std::sqrt(pow(d3d,2) - pow(dz,2));
  
  _evLog << setw(8) << "electron Eta: " << setw(8) << ele.eta 
	 << setw(8) << "electron Pt: "  << setw(8) << ele.pt 
	 << setw(8) << "electron Phi: " << setw(8) << ele.phi 
	 << setw(8) << "electron MVA: " << setw(8) << ele.mvaPOGNonTrig
	 << setw(8) << "Iso" << setw(8) << ele.pfRelIso
	 << setw(8) << "dz" << setw(8)  << dz 
	 << setw(8) << "dxy" << setw(8) << dxy << endl;
  
  double eletau1dz = std::fabs(taua.vtxDz - ele.vtxDistZ);
  _evLog << setw(8) << "Tau1 Eta: " << setw(8) << taua.eta
	 << setw(8) << "Tau1 Pt: "  << setw(8) << taua.pt
	 << setw(8) << "Tau1 Phi: " << setw(8) << taua.phi
	 << setw(8) << "drtau1e: "  << setw(8) << E1.DeltaR(T1)
	 << setw(8) << "dztau1e: "  << setw(8) << eletau1dz << endl;
  
  double eletau2dz = std::fabs(taub.vtxDz - ele.vtxDistZ);
  _evLog << setw(8) << "Tau2 Eta: "   << setw(8) << taub.eta
	 << setw(8) << "Tau2 Pt: "    << setw(8) << taub.pt
	 << setw(8) << "Tau2 Phi: "   << setw(8) << taub.phi
	 << setw(8) << "drTau2E: "    << setw(8) << E1.DeltaR(T2)
	 << setw(8) << "drtau1Tau2: " << setw(8) << T1.DeltaR(T2)
	 << setw(8) << "dztau2e: "    << setw(8) << eletau2dz << endl;
#endif
  
  histf()->cd();
#if 0
  AnaUtil::fillHist1D("ditauMass_BL",  diTauMass,   puevWt_);
  AnaUtil::fillHist1D("mvaFk_BL",      mvaOutputFK, puevWt_);
  AnaUtil::fillHist1D("mvaVZ_BL",      mvaOutput,   puevWt_);
  AnaUtil::fillHist1D("mva_evcounter", 0.0,         puevWt_); // base line
  
  if(mvaOutput > 0.65){
    AnaUtil::fillHist1D("mvaFk_1", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva1", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 1.0, puevWt_);
  }
  if(mvaOutput > 0.70){
    AnaUtil::fillHist1D("mvaFk_2", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva2", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter",  2.0,       puevWt_);
  }
  if(mvaOutput > 0.75){
    AnaUtil::fillHist1D("mvaFk_3", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva3", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 3.0, puevWt_);
  }
  if(mvaOutput > 0.78){
    AnaUtil::fillHist1D("mvaFk_4", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva4",diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 4.0, puevWt_);
  }
  
  if(mvaOutput > 0.79){
    AnaUtil::fillHist1D("mvaFk_5", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva5", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 5.0, puevWt_);
  }
  if(mvaOutput > 0.80){
    AnaUtil::fillHist1D("mvaFk_6", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva6", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 6.0, puevWt_);
  }
  if(mvaOutput > 0.81){
    AnaUtil::fillHist1D("mvaFk_7", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva7", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 7.0, puevWt_);
  }
  if(mvaOutput > 0.82){
    AnaUtil::fillHist1D("mvaFk_8", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva8", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 8.0, puevWt_);
  }
  if(mvaOutput > 0.83){
    AnaUtil::fillHist1D("mvaFk_9", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva9", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 9.0, puevWt_);
  }
  if(mvaOutput > 0.84){
    AnaUtil::fillHist1D("mvaFk_10", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva10", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 10.0, puevWt_);
  }
  if(mvaOutput > 0.85){
    AnaUtil::fillHist1D("mvaFk_11", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva11", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 11.0, puevWt_);
  }
  if(mvaOutput > 0.86){
    AnaUtil::fillHist1D("mvaFk_12", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva12", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 12.0, puevWt_);
  }
  
  if(mvaOutput > 0.87){
    AnaUtil::fillHist1D("mvaFk_13", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva13", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 13.0, puevWt_);
  }
  if(mvaOutput > 0.88){
    AnaUtil::fillHist1D("mvaFk_14", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva14", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 14.0, puevWt_);
  }
  if(mvaOutput > 0.89){
    AnaUtil::fillHist1D("mvaFk_15", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva15", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 15.0, puevWt_);
  }
  if(mvaOutput > 0.90){
    AnaUtil::fillHist1D("mvaFk_16", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva16", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 16.0, puevWt_);
  }
  if(mvaOutput > 0.91){
    AnaUtil::fillHist1D("mvaFk_17", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva17", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 17.0, puevWt_);
  }
  if(mvaOutput > 0.92){
    AnaUtil::fillHist1D("mvaFk_18", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva18", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 18.0, puevWt_);
  }
  if(mvaOutput > 0.93){
    AnaUtil::fillHist1D("mvaFk_19", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva19", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 19.0, puevWt_);
  }
  if(mvaOutput > 0.94){
    AnaUtil::fillHist1D("mvaFk_20", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva20", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 20.0, puevWt_);
  }
  if(mvaOutput > 0.95){
    AnaUtil::fillHist1D("mvaFk_21", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva21", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 21.0, puevWt_);
  }
  if(mvaOutput > 0.96){
    AnaUtil::fillHist1D("mvaFk_22", mvaOutputFK, puevWt_);
    AnaUtil::fillHist1D("ditauMass_mva22", diTauMass, puevWt_);
    AnaUtil::fillHist1D("mva_evcounter", 22.0, puevWt_);
  }
#endif 
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void ETauTau::endJob() 
{  
  histf()->cd();
  
  TH1 *h = AnaUtil::getHist1D("evcounter");
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
bool ETauTau::readJob(const string& jobFile, int& nFiles)
{
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

void ETauTau::printJob(ostream& os) const
{
  AnaBase::printJob(os);
  
  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));
  AnaUtil::showCuts(hmap, os);
}
