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

#include "FakeRate.h"
#include "AnaUtil.h"
#include "PhysicsObjects.h"

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
using std::max;
using std::sqrt;
using std::sort;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
FakeRate::FakeRate()
  : AnaBase(),
    _readMVA(false)
{}
// ----------
// Destructor
// ----------
FakeRate::~FakeRate()
{}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool FakeRate::beginJob() 
{ 
  AnaBase::beginJob();

   histf()->cd();
  bookHistograms();



  reader1 = new TMVA::Reader("Silent");
  reader1->AddVariable("LeadPt", &leadPt);
  reader1->AddVariable("SubPt", &subPt);
  reader1->AddVariable("Met", &met);
  reader1->AddVariable("DeltaRDiTau", &deltaRDiTau);
  reader1->AddVariable("PtRatio", &ptRatio);
  reader1->BookMVA("BDT8","/cmshome/kchatter/CMSSW_5_3_7_patch4/src/VHTauTau/TreeAnalyzer/June2013/WH_BDT8.weights.xml");


  if (_readMVA){
    reader = new TMVA::Reader("Silent");
    reader->AddVariable("lepPt", &lepPt);
    reader->AddVariable("tauOSPt", &tauOSPt);
    reader->AddVariable("tauSSPt", &tauSSPt);
    //reader->AddVariable("dphilepDiTau", &dphilepDiTau);
    reader->AddVariable("met", &met);
    //reader->AddVariable("PtRatio", &PtRatio);
    //reader->AddVariable("mass", &mass);
    reader->BookMVA("MLP", "/cmshome/kchatter/CMSSW_5_3_7_patch4/src/VHTauTau/TreeAnalyzer/June2013/TMVAClassification_MLP.weights.xml");
  }

  return true;
}
// ---------------
// Book histograms
// ---------------
void FakeRate::bookHistograms() 
{

  new TH1F("evtcounter", " selct tau counter", 6, -0.5, 5.5);
  new TH1F("fakeElec_deno", "Denominator of the JetToMu Fake Function", 140, 0, 140);
  new TH1F("faketau1_deno_DY", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketau2_deno_DY", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketau1_deno_WJets", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketau2_deno_WJets", "Denominator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("fakeElec_nume", "Numerator of the JetToMu Fake Function", 140, 0, 140);
  new TH1F("faketau1_nume_DY", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketau2_nume_DY", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketau1_nume_WJets", "Tight WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("faketau2_nume_WJets", "Medium WP Numerator of the JetToTau Fake Function", 140, 0, 140);
  new TH1F("ZMass_e", "Z mass of the MuMu pair after selecting the Denominator of the JetToMu Fake function", 100, 0, 100);

  new TH1F("tau1fake_pt", "pt distribution of antiisolated tau1::same sign to the electron", 140, 0, 140);
  new TH1F("tau2fake_pt", "pt distribution of antiisolated tau2::same sign to the electron", 140, 0, 140);
  new TH1F("tau2fake_eta", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH1F("tautauMass_Signal", "tautau Mass in the signal region", 320, 0, 320);
  new TH1D("EventCount", "no of events after all cut, in signal region selection", 2, 0.5, 2.5);

  //new TH1D("mvaFK", "distribution of fake MVA ", 100, -1, 1);
  // new TH2D("mvaFKCR1", "distribution of fake MVA ",140, 0, 140,  100, -1, 1);
  //new TH2D("mvaFKCR2", "distribution of fake MVA ", 140, 0, 140, 100, -1, 1);

  new TH2D("mvaFK_FKCR_tau2FK", "distribution of fake MVA in jet to tau2 CR", 140, 0, 140, 100, -1, 1);
  new TH2D("mvaVZ_FKCR_tau2FK", "distribution of fake MVA in jet to tau2 CR ", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR ", 140, 0, 140, 320, 0, 320);
  new TH2D("tauOSpt_FKCR_tau2FK", "distribution of pt of OS #tau with Pt of #tauSS in jet to tau CR ", 140, 0, 140, 140, 0, 140);
  new TH2D("tauSSpt_FKCR_tau2FK", "distribution of pt of SS #tau with Pt of #tauSS in jet to tau CR ", 140, 0, 140, 140, 0, 140);
  new TH2D("tauOSeta_FKCR_tau2FK", "distribution of eta of OS #tau with Pt of #tauSS in jet to tau CR ", 140, 0, 140, 6, -3, 3);
  new TH2D("tauSSeta_FKCR_tau2FK", "distribution of eta of SS #tau with Pt of #tauSS in jet to tau CR ", 140, 0, 140, 6, -3, 3);

  new TH1F("tau2fake_pt1", "pt distribution of #tau1 in jet to tau2 CR after MVA1", 140, 0, 140);
  new TH1F("tau2fake_eta1", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR1", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_1", "distribution of fake MVA in jet to tau2 CR after MVA1", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK1", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA1 ", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt2", "pt distribution of #tau1 in jet to tau2 CR after MVA2", 140, 0, 140);
  new TH1F("tau2fake_eta2", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3, 3);
  new TH1F("tautauMass_fakeCR2", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_2", "distribution of fake MVA in jet to tau2 CR after MVA2", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK2", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA2 ", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt3", "pt distribution of #tau1 in jet to tau2 CR after MVA3", 140, 0, 140);
  new TH1F("tau2fake_eta3", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR3", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_3", "distribution of fake MVA in jet to tau2 CR after MVA3", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK3", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA3 ", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt4", "pt distribution of #tau1 in jet to tau2 CR after MVA4", 140, 0, 140);
  new TH1F("tau2fake_eta4", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR4", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_4", "distribution of fake MVA in jet to tau2 CR after MVA4", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK4", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA4", 140, 0, 140, 320, 0, 320);
  new TH2D("tauOSpt_FKCR_tau2FK4", "distribution of pt of OS #tau with Pt of #tauSS in jet to tau CR ", 140, 0, 140, 140, 0, 140);
  new TH2D("tauSSpt_FKCR_tau2FK4", "distribution of pt of SS #tau with Pt of #tauSS in jet to tau CR ", 140, 0, 140, 140, 0, 140);
  new TH2D("tauOSeta_FKCR_tau2FK4", "distribution of eta of OS #tau with Pt of #tauSS in jet to tau CR after mva4", 140, 0, 140, 6, -3, 3);
  new TH2D("tauSSeta_FKCR_tau2FK4", "distribution of eta of SS #tau with Pt of #tauSS in jet to tau CR after mva4 ", 140, 0, 140, 6, -3, 3);


  new TH1F("tau2fake_pt5", "pt distribution of #tau1 in jet to tau2 CR after MVA5", 140, 0, 140);
  new TH1F("tau2fake_eta5", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR5", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_5", "distribution of fake MVA in jet to tau2 CR after MVA5", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK5", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA5 ", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt6", "pt distribution of #tau1 in jet to tau2 CR after MVA6", 140, 0, 140);
  new TH1F("tau2fake_eta6", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR6", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_6", "distribution of fake MVA in jet to tau2 CR after MVA6", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK6", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA6 ", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt7", "pt distribution of #tau1 in jet to tau2 CR after MVA7", 140, 0, 140);
  new TH1F("tau2fake_eta7", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR7", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_7", "distribution of fake MVA in jet to tau2 CR after MVA7", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK7", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA7 ", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt8", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta8", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR8", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_8", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK8", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt9", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta9", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR9", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_9", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK9", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt10", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta10", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR10", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_10", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK10", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt11", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta11", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR11", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_11", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK11", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt12", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta12", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR12", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_12", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK12", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt13", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta13", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR13", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_13", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK13", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt14", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta14", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR14", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_14", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK14", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt15", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta15", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR15", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_15", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK15", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt16", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta16", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR16", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_16", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK16", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt17", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta17", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR17", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_17", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK17", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt18", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta18", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR18", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_18", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK18", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt19", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta19", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR19", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_19", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK19", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt20", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta20", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR20", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_20", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK20", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt21", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta21", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR21", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_21", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK21", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("tau2fake_pt22", "pt distribution of #tau1 in jet to tau2 CR after MVA", 140, 0, 140);
  new TH1F("tau2fake_eta22", "eta distribution of antiisolated tau2::same sign to the electron", 12, -3,3);
  new TH1F("tautauMass_fakeCR22", "tautau Mass in the SS-AntiIsolated control region", 320, 0, 320);
  new TH2D("mvaFK_FKCR_tau2FK_22", "distribution of fake MVA in jet to tau2 CR after MVA", 140, 0, 140, 100, -1, 1);
  new TH2D("tautauMass_FKCR_tau2FK22", "distribution of #tau, #tau mass with Pt of #tau2 in jet to tau1 CR after MVA", 140, 0, 140, 320, 0, 320);

  new TH1F("mvaFK_signal", "fake MVA in the signal region", 100, -1, 1);
  new TH1F("mvaVZ_signal", "VZ MVA in the signal region", 100, -1, 1);
  new TH1F("tauOSpt_signal", "pt of #tau^(os) in the signal region", 140, 0, 140);
  new TH1F("tauSSpt_signal", "pt of #tau^(ss) in the signal region", 140, 0, 140);
  new TH1F("tauOSeta_signal", "eta of #tau^(os) in the signal region", 6, -3, 3);
  new TH1F("tauSSeta_signal", "eta of #tau^(ss) in the signal region", 6, -3, 3);

  new TH1F("tautauMass_Signal_1", "#tau #tau mass in the signal region after MVA1", 320, 0, 320);
  new TH1F("mvaFK_signal_1", "fake MVA in the signal region after MVA1", 100, -1, 1);
  new TH1F("mvaVZ_signal_1", "VZ MVA in the signal region after MVA1", 100, -1, 1);

  new TH1F("tautauMass_Signal_2", "#tau #tau mass in the signal region after MVA2", 320, 0, 320);
  new TH1F("mvaFK_signal_2", "fake MVA in the signal region after MVA2", 100, -1, 1);
  new TH1F("mvaVZ_signal_2", "VZ MVA in the signal region after MVA2", 100, -1, 1);

  new TH1F("tautauMass_Signal_3", "#tau #tau mass in the signal region after MVA3", 320, 0, 320);
  new TH1F("mvaFK_signal_3", "fake MVA in the signal region after MVA3", 100, -1, 1);
  new TH1F("mvaVZ_signal_3", "VZ MVA in the signal region after MVA3", 100, -1, 1);

  new TH1F("tautauMass_Signal_4", "#tau #tau mass in the signal region after MVA4", 320, 0, 320);
  new TH1F("mvaFK_signal_4", "fake MVA in the signal region after MVA4", 100, -1, 1);
  new TH1F("mvaVZ_signal_4", "VZ MVA in the signal region after MVA4", 100, -1, 1);
  new TH1F("tauOSpt_signal_4", "pt of #tau^(os) in the signal region after mva4", 140, 0, 140);
  new TH1F("tauSSpt_signal_4", "pt of #tau^(ss) in the signal region after mva4", 140, 0, 140);
  new TH1F("tauOSeta_signal_4", "eta of #tau^(os) in the signal region after mva4", 6, -3, 3);
  new TH1F("tauSSeta_signal_4", "eta of #tau^(ss) in the signal region after mva4", 6, -3, 3);

  new TH1F("tautauMass_Signal_5", "#tau #tau mass in the signal region after MVA5", 320, 0, 320);
  new TH1F("mvaFK_signal_5", "fake MVA in the signal region after MVA5", 100, -1, 1);
  new TH1F("mvaVZ_signal_5", "VZ MVA in the signal region after MVA5", 100, -1, 1);

  new TH1F("tautauMass_Signal_6", "#tau #tau mass in the signal region after MVA6", 320, 0, 320);
  new TH1F("mvaFK_signal_6", "fake MVA in the signal region after MVA6", 100, -1, 1);
  new TH1F("mvaVZ_signal_6", "VZ MVA in the signal region after MVA6", 100, -1, 1);

  new TH1F("tautauMass_Signal_7", "#tau #tau mass in the signal region after MVA7", 320, 0, 320);
  new TH1F("mvaFK_signal_7", "fake MVA in the signal region after MVA7", 100, -1, 1);
  new TH1F("mvaVZ_signal_7", "VZ MVA in the signal region after MVA7", 100, -1, 1);

  new TH1F("tautauMass_Signal_8", "#tau #tau mass in the signal region after MVA8", 320, 0, 320);
  new TH1F("mvaFK_signal_8", "fake MVA in the signal region after MVA8", 100, -1, 1);
  new TH1F("mvaVZ_signal_8", "VZ MVA in the signal region after MVA8", 100, -1, 1);

  new TH1F("tautauMass_Signal_9", "#tau #tau mass in the signal region after MVA9", 320, 0, 320);
  new TH1F("mvaFK_signal_9", "fake MVA in the signal region after MVA9", 100, -1, 1);
  new TH1F("mvaVZ_signal_9", "VZ MVA in the signal region after MVA9", 100, -1, 1);

  new TH1F("tautauMass_Signal_10", "#tau #tau mass in the signal region after MVA10", 320, 0, 320);
  new TH1F("mvaFK_signal_10", "fake MVA in the signal region after MVA10", 100, -1, 1);
  new TH1F("mvaVZ_signal_10", "VZ MVA in the signal region after MVA10", 100, -1, 1);

  new TH1F("tautauMass_Signal_11", "#tau #tau mass in the signal region after MVA11", 320, 0, 320);
  new TH1F("mvaFK_signal_11", "fake MVA in the signal region after MVA11", 100, -1, 1);
  new TH1F("mvaVZ_signal_11", "VZ MVA in the signal region after MVA11", 100, -1, 1);

  new TH1F("tautauMass_Signal_12", "#tau #tau mass in the signal region after MVA12", 320, 0, 320);
  new TH1F("mvaFK_signal_12", "fake MVA in the signal region after MVA12", 100, -1, 1);
  new TH1F("mvaVZ_signal_12", "VZ MVA in the signal region after MVA12", 100, -1, 1);

  new TH1F("tautauMass_Signal_13", "#tau #tau mass in the signal region after MVA13", 320, 0, 320);
  new TH1F("mvaFK_signal_13", "fake MVA in the signal region after MVA13", 100, -1, 1);
  new TH1F("mvaVZ_signal_13", "VZ MVA in the signal region after MVA13", 100, -1, 1);

  new TH1F("tautauMass_Signal_14", "#tau #tau mass in the signal region after MVA14", 320, 0, 320);
  new TH1F("mvaFK_signal_14", "fake MVA in the signal region after MVA14", 100, -1, 1);
  new TH1F("mvaVZ_signal_14", "VZ MVA in the signal region after MVA14", 100, -1, 1);

  new TH1F("tautauMass_Signal_15", "#tau #tau mass in the signal region after MVA15", 320, 0, 320);
  new TH1F("mvaFK_signal_15", "fake MVA in the signal region after MVA15", 100, -1, 1);
  new TH1F("mvaVZ_signal_15", "VZ MVA in the signal region after MVA15", 100, -1, 1);

  new TH1F("tautauMass_Signal_16", "#tau #tau mass in the signal region after MVA16", 320, 0, 320);
  new TH1F("mvaFK_signal_16", "fake MVA in the signal region after MVA16", 100, -1, 1);
  new TH1F("mvaVZ_signal_16", "VZ MVA in the signal region after MVA16", 100, -1, 1);

  new TH1F("tautauMass_Signal_17", "#tau #tau mass in the signal region after MVA17", 320, 0, 320);
  new TH1F("mvaFK_signal_17", "fake MVA in the signal region after MVA17", 100, -1, 1);
  new TH1F("mvaVZ_signal_17", "VZ MVA in the signal region after MVA17", 100, -1, 1);

  new TH1F("tautauMass_Signal_18", "#tau #tau mass in the signal region after MVA18", 320, 0, 320);
  new TH1F("mvaFK_signal_18", "fake MVA in the signal region after MVA18", 100, -1, 1);
  new TH1F("mvaVZ_signal_18", "VZ MVA in the signal region after MVA18", 100, -1, 1);

  new TH1F("tautauMass_Signal_19", "#tau #tau mass in the signal region after MVA19", 320, 0, 320);
  new TH1F("mvaFK_signal_19", "fake MVA in the signal region after MVA19", 100, -1, 1);
  new TH1F("mvaVZ_signal_19", "VZ MVA in the signal region after MVA19", 100, -1, 1);

  new TH1F("tautauMass_Signal_20", "#tau #tau mass in the signal region after MVA20", 320, 0, 320);
  new TH1F("mvaFK_signal_20", "fake MVA in the signal region after MVA20", 100, -1, 1);
  new TH1F("mvaVZ_signal_20", "VZ MVA in the signal region after MVA20", 100, -1, 1);

  new TH1F("tautauMass_Signal_21", "#tau #tau mass in the signal region after MVA21", 320, 0, 320);
  new TH1F("mvaFK_signal_21", "fake MVA in the signal region after MVA21", 100, -1, 1);
  new TH1F("mvaVZ_signal_21", "VZ MVA in the signal region after MVA21", 100, -1, 1);

  new TH1F("tautauMass_Signal_22", "#tau #tau mass in the signal region after MVA22", 320, 0, 320);
  new TH1F("mvaFK_signal_22", "fake MVA in the signal region after MVA22", 100, -1, 1);
  new TH1F("mvaVZ_signal_22", "VZ MVA in the signal region after MVA22", 100, -1, 1);
}
// -------------------
// The main event loop
// -------------------
void FakeRate::clearLists() {
  vtxList.clear();
  muoList.clear();
  eleList.clear();
  // tauList.clear();
  bjetList.clear();
  selectTauList.clear();
  selectEList.clear();
  trigObjList.clear();
}
void FakeRate::eventLoop() 
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
  for (int ev = 0; ev < nEvents; ++ev) {
    clearEvent();
    int lflag = _chain->LoadTree(ev); 
    int nbytes = getEntry(lflag);    // returns total bytes read


    string currentFile(gSystem->BaseName(_chain->GetCurrentFile()->GetName())); 

    const Event* evt = dynamic_cast<Event*>(eventA->At(0));
    assert(evt);

    // PileUP weight
    _puevWt = 1; // for data

#if 0
    if (_isMC) {
      int npu = 0;
      _puevWt = wtPileUp(npu);
    }
#endif

    AnaUtil::fillHist1D("EventCount", 1.0, _puevWt);    
    
    int run   = evt->run;
    int event = evt->event;
    int lumis = evt->lumis;

    // Show status of the run
    if (currentFile != lastFile) 
    cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
         << " ==> " << currentFile 
         << " <<< Run# " << run
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;
    lastFile = currentFile;

    // Show the status 
    if (ev%nPrint == 0) 
    cout << "Tree# " << setw(4) << _chain->GetTreeNumber()  
         << " ==> " << _chain->GetCurrentFile()->GetName() 
         << " <<< Run# " << run 
         << " Lumis# " << lumis
         << " Event# " << setw(8) << event << " >>> " 
         << " Events proc. " << setw(8) << ev 
         << endl;
#if 0
    if (_eventIdMap.size()) {
      std::ostringstream mkey;
      mkey << event << "-" << lumis << "-" << run;
      if (_eventIdMap.find(mkey.str()) == _eventIdMap.end()) continue;
      _evLog << "Event " << event
            << " Lumis " << lumis
            << " Run " << run << endl;
      dumpGenInfo(_evLog); 
    }
#endif
    if (_logOption > 0) {
      _fLog << "Event " << event
            << " Lumis " << lumis
            << " Run " << run << endl;
      _fLog << "n_tau: "<< n_tau
            << ", n_muon: "<< n_muon
            << ", n_jet: " << n_jet
            << ", n_vertex: " << n_vertex
            << ", n_met: " << n_met
            << ", n_electron: " << n_electron 
            << endl;
    }
    clearLists();
 
    // Trigger selection
    if (_useTrigger && !isTriggered()) continue;

    op.verbose = (_logOption >> 2 & 0x1); 
    findVtxInfo(vtxList, op, _fLog);
    int nvtx = vtxList.size();
    double vz = ( nvtx > 0 ) ? vtxList.at(0).z : -999;

    if (_logOption > 0)
    _fLog << "Event " << event
          << " Lumis " << evt->lumis
          << " Run " << run 
          << " n_vertex_good " << nvtx
          << endl;

    op.verbose = (_logOption >> 5 & 0x1);
    findJetInfo(bjetList, op, _fLog);

    // Event Selection Starts here .....
    // presence of > 1 good vertex
    if (vtxList.size() < 1) continue;

    op.verbose = (_logOption >> 3 & 0x1);
    findMuonInfo(muoList, vz, op, _fLog);

    findTriggerObjectInfo(trigObjList);
    selectEvent();
  }  
  // Analysis is over
  endJob();
}
void FakeRate::selectEvent() {
  selectTau();
  selectElectron();
  JetToFakeDYJets();     //Here We Calculate Jet-> Mu, Tau1 and Tau2 fake rate in DYToMuMu + Jets Control region
  JetToTauFakeWJets();   // Here We Calculate Jet-> Tau1 and Tau2 fake function in W+Jets Control Region 
  ApplyFR();
}

void FakeRate::selectTau()                
{
  int tcounters[] = {0,0,0,0,0,0};
  for (int indx = 0; indx < n_tau; ++indx) {
    const Tau* tau = dynamic_cast<Tau*>(tauA->At(indx));
    if (!tau) continue;
    //    if ((abs(tau->eta) >= 2.3) ||
    //  (tau->decayModeFinding <= 0.5) ||
    //  (tau->againstMuonTight <= 0.5) ||
    //( tau->againstElectronLooseMVA3) ||
    //(tau->vtxDz <= 0.2)
        //(tau->againstElectronLoose <= 0.5)
    // ) continue; 
    ++tcounters[0];
    if (std::abs(tau->eta) >= 2.3) continue;
    ++tcounters[1]; 
    if (tau->decayModeFinding <= 0.5)  continue;
    ++tcounters[2];
    if (tau->againstMuonTight <= 0.5)  continue;
    ++tcounters[3];
    if ( tau->againstElectronLooseMVA3 <= 0.5)  continue;
    ++tcounters[4];
    if (std::abs(tau->vtxDz) >= 0.2) continue;
    ++tcounters[5];

    selectTauList.push_back(*tau);
  }
  if (selectTauList.size() > 1)
    sort(selectTauList.begin(), selectTauList.end(), PtComparator<Tau>());

  for (size_t i = 0; i < NEL(tcounters); ++i) {
    if (tcounters[i] > 0) AnaUtil::fillHist1D("evtcounter", i, _puevWt);
  }

}

void FakeRate::selectElectron()
{
  electron_indx = -1;
  for (int indx = 0; indx < n_electron; ++indx) {
    const Electron* electron = dynamic_cast<Electron*>(electronA->At(indx));
    if (!electron) continue;


    double eleta = std::abs(electron->eta);
    if (( eleta >= AnaUtil::cutValue(_electronCutMap, "etaLow")
       &&   eleta <= AnaUtil::cutValue(_electronCutMap, "etaUp"))
       || eleta >= AnaUtil::cutValue(_electronCutMap, "eta") ) continue;

    // double d3d = electron->vtxDist3D;
    double dz  = electron->vtxDistZ;
    //double dxy = std::sqrt(pow(d3d,2) - pow(dz,2));
    if (std::abs(electron->dxyPV) >= AnaUtil::cutValue(_electronCutMap,"trkD0")) continue;
    if (std::abs(dz) >= AnaUtil::cutValue(_electronCutMap, "dz")) continue;
    if ( !AnaBase::electronMVA(electron) ) continue;

    selectEList.push_back(*electron);
    double electronIso = electron->pfRelIsoDB04;
    bool electronIsoBarrel = ( electronIso < AnaUtil::cutValue(_electronCutMap, "relIsoBarrel") 
			       && eleta < AnaUtil::cutValue(_electronCutMap, "etaLow") ); //Iso < .15 for abs(eta) < 1.4442
    bool electronIsoEndCap = ( electronIso < AnaUtil::cutValue(_electronCutMap, "relIsoEndCap")
			       && eleta > AnaUtil::cutValue(_electronCutMap, "etaUp") ); //Iso < 0.10 for 1.556 < abs(eta) < 2.5
    if (!(electronIsoBarrel || electronIsoEndCap)) continue;
    if (electron->pt <= AnaUtil::cutValue(_electronCutMap, "pt")) continue;

    eleList.push_back(*electron);
    electron_indx = indx;
  }
}


void FakeRate::JetToFakeDYJets()
{
  if (muoList.size() > 0)  return;

  vector<Electron> fEleList, frEleList;
  for (unsigned int indx = 0; indx < selectEList.size(); ++indx){
    const Electron& ele = selectEList.at(indx);

    if (fEleList.size() >= 2)
      frEleList.push_back(ele);                      // This is to Catch the 3rd elctron // Jet->Electron Fake                                                                                 
    double electronIso = ele.pfRelIsoDB04;
    double eleta = std::abs(ele.eta);
    bool electronIsoBarrel = ( electronIso < AnaUtil::cutValue(_electronCutMap, "relIsoBarrel")
                            && eleta < AnaUtil::cutValue(_electronCutMap, "etaLow") );
    bool electronIsoEndCap = ( electronIso < AnaUtil::cutValue(_electronCutMap, "relIsoEndCap")
                            && eleta > AnaUtil::cutValue(_electronCutMap, "etaUp") );
    if (!(electronIsoBarrel || electronIsoEndCap)) continue;
    if (ele.pt <= AnaUtil::cutValue(_electronCutMap, "pt")) continue;

    // if (abs(ele.pfRelIsoDB04v2) >= AnaUtil::cutValue(_electronCutMap, "pfRelIso") || 
    fEleList.push_back(ele);
  }

  // atleast 2 good electron
  if (fEleList.size() < 2) return;
  const Electron& ele1 = fEleList.at(0);
  const Electron& ele2 = fEleList.at(1);

  TLorentzVector e1, e2,t1,sum;
  e1.SetPtEtaPhiE (ele1.pt, ele1.eta, ele1.phi, ele1.energy);
  e2.SetPtEtaPhiE (ele2.pt, ele2.eta, ele2.phi, ele2.energy);
  sum = e1 + e2;
  if (AnaUtil::deltaR(e1, e2) < 0.5) return;
  if ((ele1.charge + ele2.charge) != 0) return;
  if (sum.M() < 85 || sum.M() > 95) return;
  AnaUtil::fillHist1D("ZMass_e", sum.M(), _puevWt);

  // b-tagged Jets Veto   
  if (bjetList.size() > 0) return;
  if (ele1.pt <= 24 || ele2.pt <= 10) return;
  const MET* mt = dynamic_cast<MET*>(mvametA->At(0));
  assert(mt);
  // jet to e fake function
  if (frEleList.size() >= 1){
    for (unsigned int indx = 0; indx < frEleList.size(); ++indx){
      const Electron& fake = frEleList.at(indx);

      //double mass1 = sqrt(2*fake.pt*mt->met*(1-cos(AnaUtil::deltaPhi(fake.phi, mt->metphi))));                                                                                               
      //if (mass1 >= 20) continue;                                                          // To Remove WZ Contamination                                                                      
      AnaUtil::fillHist1D("fakeElec_deno", fake.pt, _puevWt);

      double fakeIso = fake.pfRelIsoDB04;
      double fakeeta = std::abs(fake.eta);
      bool fakeIsoBarrel = ( fakeIso < AnaUtil::cutValue(_electronCutMap, "relIsoBarrel")
				 && fakeeta < AnaUtil::cutValue(_electronCutMap, "etaLow") );
      bool fakeIsoEndCap = ( fakeIso < AnaUtil::cutValue(_electronCutMap, "relIsoEndCap")
				 && fakeeta > AnaUtil::cutValue(_electronCutMap, "etaUp") );

      if ((fakeIsoBarrel || fakeIsoEndCap))
	AnaUtil::fillHist1D("fakeElec_nume", fake.pt, _puevWt);
    }
  }

  // Tau1 & Tau2 Fake Function for DYJets Control Region
  if (fEleList.size() !=2) return;
  if (!selectTauList.size()) return;      // No Tau1 means no Tau2, so we can safely reject the event 
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx) {
    const Tau& fake = selectTauList.at(indx);
    t1.SetPtEtaPhiE(fake.pt, fake.eta, fake.phi, fake.energy);
    double dr1 = AnaUtil::deltaR(e1, t1);
    double dr2 = AnaUtil::deltaR(e2, t1);
    if ((dr1 < AnaUtil::cutValue(_evselCutMap, "drEleTau")) ||
        (dr2 < AnaUtil::cutValue(_evselCutMap, "drEleTau"))) continue;

    // as we have two electrons with opposite charges, we do not impose any charge sign
    // on the tau
    //double mass1 = sqrt(2*fake.pt*mt->met*(1-cos(AnaUtil::deltaPhi(fake.phi, mt->metphi)))); 
    //if (mass1 >= 20) continue;                                            // To Remove WZ Contamination  

    AnaUtil::fillHist1D("faketau1_deno_DY", fake.pt, _puevWt);
    //if (fake.byMediumCombinedIsolationDeltaBetaCorr >= 0.5)
    if (fake.byMediumCombinedIsolationDeltaBetaCorr3Hits >= 0.5)
      AnaUtil::fillHist1D("faketau1_nume_DY", fake.pt, _puevWt);
  }
}

void FakeRate::JetToTauFakeWJets()
{
  double vz = vtxList.at(0).z;

  // atleast 1 good electron
  if (eleList.size() != 1) return;
  const Electron& ele1 = eleList.at(0);
  TLorentzVector e1;
  e1.SetPtEtaPhiE (ele1.pt, ele1.eta, ele1.phi, ele1.energy);
  const MET* mt = dynamic_cast<MET*>(mvametA->At(0));
  assert(mt);
  double mass1 = sqrt(2*ele1.pt*mt->met*(1-cos(AnaUtil::deltaPhi(ele1.phi, mt->metphi))));
  if (mass1 <= 40) return;

  //  Muon Veto 
  if (muoList.size() > 0) return;

  // b-tagged Jets Veto  
  if (bjetList.size() > 0) return;

  TLorentzVector t1, t2;

  if (selectTauList.size() < 2) return;
  vector<Tau> taulist;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx){
    const Tau& fake = selectTauList.at(indx);
    t1.SetPtEtaPhiE(fake.pt, fake.eta, fake.phi, fake.energy);
    double dr1 = AnaUtil::deltaR(e1, t1);
    if (dr1 < AnaUtil::cutValue(_evselCutMap, "drEleTau")) continue;
    taulist.push_back(fake);
  }
  if (taulist.size() < 2) return;
  const Tau& ta1 = taulist.at(0);
  const Tau& ta2 = taulist.at(1);

  // all the 3 objects must have the same sign
  if (ta1.charge != ta2.charge || ta1.charge != ele1.charge) return;
  for (int indx = 0; indx < 2; ++indx) {
    const Tau& fake = taulist.at(indx);
    AnaUtil::fillHist1D("faketau1_deno_WJets", fake.pt, _puevWt);
    if (fake.byMediumCombinedIsolationDeltaBetaCorr3Hits >= 0.5)
      AnaUtil::fillHist1D("faketau1_nume_WJets", fake.pt, _puevWt);
  }
}

void FakeRate::ApplyFR()
{
  double vz = vtxList.at(0).z;

  // atleast 1 good electron
  if (eleList.size() < 1) return;
  const Electron& ele = eleList.at(0);
  TLorentzVector e1;
  e1.SetPtEtaPhiE (ele.pt, ele.eta, ele.phi, ele.energy);

  if (selectTauList.size() < 2) return;

  const MET* mt = dynamic_cast<MET*>(mvametA->At(0));
  assert(mt);
  double mvaMet = mt->met;
  double MetPhi = mt->metphi;

  ////// SEARCH FOR  TAUs /////////
  vector<Tau> tauSSList, tauOSList,antitauSSList;
  for (unsigned int indx = 0; indx < selectTauList.size(); ++indx){
    const Tau& tau = selectTauList.at(indx);
    TLorentzVector T;
    T.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
    double dr1 = AnaUtil::deltaR(e1, T);

    //double eletaudz = std::abs(tau.vtxDz - ele.vtxDistZ);
    bool decision1 = tau.pt > 20 &&
                    dr1 >= AnaUtil::cutValue(_evselCutMap, "drEleTau");
  
    bool isotauSS = (decision1 &&  
                     tau.byMediumCombinedIsolationDeltaBetaCorr3Hits >= 0.5 &&
                     (tau.charge + ele.charge) != 0 );

    bool isotauOS = (decision1 &&
                     tau.byMediumCombinedIsolationDeltaBetaCorr3Hits >= 0.5 &&
                     (tau.charge + ele.charge) == 0 &&
                     tau.againstElectronTightMVA3 >= 0.5);

    bool antiisotauSS = ( decision1 &&
                        tau.byMediumCombinedIsolationDeltaBetaCorr3Hits < 0.5 &&
                        (tau.charge + ele.charge) != 0 );

    // we apply Trigger object mactching on the OS tau canditate only to ensure that
    // the SS leg is unbiased, even at the Trigger level regarding isolation.
    if (isotauOS){
      TLorentzVector tauOS;
      tauOS.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
      double maxPtDiff = AnaUtil::cutValue(_evselCutMap, "maxPtDiff");
      int ntobj = trigObjList.size();
      int tindx = -1;
      uint flag = 0;
      double drTag = matchTriggerObject(trigObjList, tauOS, _trigPathList, -1, maxPtDiff, tindx, flag);
      if (tindx < 0 || tindx >= ntobj) continue;
      if (drTag >= AnaUtil::cutValue(_evselCutMap, "maxDr") || flag != 1) continue;
      tauOSList.push_back(tau);
    }
    else if (isotauSS)
       tauSSList.push_back(tau);
    else if (antiisotauSS) 
       antitauSSList.push_back(tau);
  }
  //  if(tau1List.size() <1 && tau2List.size() <1 ) return; // event should have atleast 1 Iso Tau
  bool Signal = (tauSSList.size() && tauOSList.size()) ? true : false;
  bool JetToTauFake = (!tauSSList.size() && antitauSSList.size() && tauOSList.size()) ? true : false;
  bool result = Signal || JetToTauFake;
  if (!result) return;

  const Tau& tauOS = tauOSList.at(0);
  TLorentzVector TSS, TOS;
  TOS.SetPtEtaPhiE(tauOS.pt, tauOS.eta, tauOS.phi, tauOS.energy);
  if (JetToTauFake)
    TSS.SetPtEtaPhiE(antitauSSList.at(0).pt, antitauSSList.at(0).eta, antitauSSList.at(0).phi, antitauSSList.at(0).energy);
  else if (Signal)
    TSS.SetPtEtaPhiE(tauSSList.at(0).pt, tauSSList.at(0).eta, tauSSList.at(0).phi, tauSSList.at(0).energy);
  //bool Signal = tau1List.size() && tauList.size();

  double highestPt = std::max(TSS.Pt(), TOS.Pt());
  if (highestPt <= 25) return;

  double drtt = TOS.DeltaR(TSS);
  if (drtt <= AnaUtil::cutValue(_evselCutMap, "drTauTau")) return;

  if (AnaBase::vetoMuon(15, 0.2) > 0) return;
  if (AnaBase::vetoElectron(15, 0.2) > 1) return;

  // z -> ee veto
  double tauOSemfrac = tauOS.emFraction;
  bool ztoee;
  ztoee = AnaBase::DYtoeeVeto( TOS, TSS, tauOSemfrac, ele, electron_indx);
  if (ztoee == true) return;

 // z to tau tau veto
  double massETauOS = ( e1 + TOS ).M();
  double massETauSS = ( e1 + TSS ).M();
  if (mvaMet < 40 && massETauOS <90 && massETauSS < 90) return;

  double MtEleMet = sqrt(pow(e1.Pt() + mvaMet,2) - pow(e1.Pt()*cos(e1.Phi()) + mvaMet*cos(MetPhi),2) -
                    pow(e1.Pt()*sin(e1.Phi()) + mvaMet*sin(MetPhi),2));

  if (MtEleMet < 20) return; // analysis overlap removing veto

  // if (mvaMet < 15) return;

  //No b-Jet Veto
  if (bjetList.size() > 0) return;


  TLorentzVector diTau;
  diTau = TSS + TOS;
  //////// READING FAKE MVA ////////
  double mvaOutputFK = -999; 
  leadPt     = TOS.Pt();
  subPt     = TSS.Pt();
  met        = mvaMet;
  deltaRDiTau   = AnaUtil::deltaR(TOS, TSS);
  ptRatio = diTau.Pt()/(TOS.Pt() + TSS.Pt());
  mvaOutputFK = reader1->EvaluateMVA("BDT8");

  ////////// READING IRREDUCIBLE MVA /////////
  double mvaOutput = -999;      //Set at high negative value, as MvaOutput > -0.106
  if (_readMVA){
    lepPt = ele.pt;
    tauOSPt = TOS.Pt();
    tauSSPt = TSS.Pt();
    //dphilepDiTau = AnaUtil::deltaPhi(e1, diTau);
    //PtRatio = (diTau.Pt())/(TOS.Pt() + TSS.Pt()); 
    //mass = diTau.M();
    mvaOutput = reader->EvaluateMVA("MLP");
  }


  if (JetToTauFake){
    AnaUtil::fillHist1D("tau2fake_pt", TSS.Pt(), _puevWt);
    AnaUtil::fillHist1D("tau2fake_eta", TSS.Eta(), _puevWt);
    AnaUtil::fillHist1D("tautauMass_fakeCR", diTau.M(), _puevWt);
    AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK",TSS.Pt(), mvaOutputFK, _puevWt);
    AnaUtil::fillHist2D("mvaVZ_FKCR_tau2FK",TSS.Pt(), mvaOutput, _puevWt);
    AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK",TSS.Pt(), diTau.M(), _puevWt);
    AnaUtil::fillHist2D("tauOSpt_FKCR_tau2FK",TSS.Pt(), TOS.Pt(), _puevWt);
    AnaUtil::fillHist2D("tauSSpt_FKCR_tau2FK",TSS.Pt(), TSS.Pt(), _puevWt);
    AnaUtil::fillHist2D("tauOSeta_FKCR_tau2FK",TSS.Pt(), TOS.Eta(), _puevWt);
    AnaUtil::fillHist2D("tauSSeta_FKCR_tau2FK",TSS.Pt(), TSS.Eta(), _puevWt);

    if (mvaOutput > 0.66){
      AnaUtil::fillHist1D("tau2fake_pt1", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta1", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR1", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_1",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK1",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.70){
      AnaUtil::fillHist1D("tau2fake_pt2", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta2", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR2", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_2",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK2",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.75){
      AnaUtil::fillHist1D("tau2fake_pt3", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta3", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR3", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_3",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK3",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.78){
      AnaUtil::fillHist1D("tau2fake_pt4", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta4", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR4", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_4",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK4",TSS.Pt(), diTau.M(), _puevWt);
      AnaUtil::fillHist2D("tauOSpt_FKCR_tau2FK4",TSS.Pt(), TOS.Pt(), _puevWt);
      AnaUtil::fillHist2D("tauSSpt_FKCR_tau2FK4",TSS.Pt(), TSS.Pt(), _puevWt);
      AnaUtil::fillHist2D("tauOSeta_FKCR_tau2FK4",TSS.Pt(), TOS.Eta(), _puevWt);
      AnaUtil::fillHist2D("tauSSeta_FKCR_tau2FK4",TSS.Pt(), TSS.Eta(), _puevWt);
    }
    if (mvaOutput > 0.79){
      AnaUtil::fillHist1D("tau2fake_pt5", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta5", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR5", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_5",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK5",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.80){
      AnaUtil::fillHist1D("tau2fake_pt6", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta6", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR6", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_6",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK6",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.81){
      AnaUtil::fillHist1D("tau2fake_pt7", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta7", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR7", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_7",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK7",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.82){
      AnaUtil::fillHist1D("tau2fake_pt8", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta8", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR8", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_8",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK8",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.83){
      AnaUtil::fillHist1D("tau2fake_pt9", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta9", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR9", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_9",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK9",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.84){
      AnaUtil::fillHist1D("tau2fake_pt10", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta10", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR10", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_10",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK10",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.85){
      AnaUtil::fillHist1D("tau2fake_pt11", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta11", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR11", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_11",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK11",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.86){
      AnaUtil::fillHist1D("tau2fake_pt12", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta12", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR12", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_12",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK12",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.87){
      AnaUtil::fillHist1D("tau2fake_pt13", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta13", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR13", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_13",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK13",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.88){
      AnaUtil::fillHist1D("tau2fake_pt14", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta14", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR14", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_14",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK14",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.89){
      AnaUtil::fillHist1D("tau2fake_pt15", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta15", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR15", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_15",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK15",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.90){
      AnaUtil::fillHist1D("tau2fake_pt16", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta16", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR16", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_16",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK16",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.91){
      AnaUtil::fillHist1D("tau2fake_pt17", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta17", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR17", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_17",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK17",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.92){
      AnaUtil::fillHist1D("tau2fake_pt18", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta18", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR18", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_18",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK18",TSS.Pt(), diTau.M(), _puevWt);
    }

    if (mvaOutput > 0.93){
      AnaUtil::fillHist1D("tau2fake_pt19", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta19", TSS.Eta(), _puevWt);  
      AnaUtil::fillHist1D("tautauMass_fakeCR19", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_19",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK19",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.94){
      AnaUtil::fillHist1D("tau2fake_pt20", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta20", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR20", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_20",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK20",TSS.Pt(), diTau.M(), _puevWt); 
    }
    if (mvaOutput > 0.95){
      AnaUtil::fillHist1D("tau2fake_pt21", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta21", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR21", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_21",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK21",TSS.Pt(), diTau.M(), _puevWt);
    }
    if (mvaOutput > 0.96){
      AnaUtil::fillHist1D("tau2fake_pt22", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tau2fake_eta22", TSS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tautauMass_fakeCR22", diTau.M(), _puevWt);
      AnaUtil::fillHist2D("mvaFK_FKCR_tau2FK_22",TSS.Pt(), mvaOutputFK, _puevWt);
      AnaUtil::fillHist2D("tautauMass_FKCR_tau2FK22",TSS.Pt(), diTau.M(), _puevWt);
   }

  }

  if (Signal){
    AnaUtil::fillHist1D("EventCount", 2.0, _puevWt);
    AnaUtil::fillHist1D("tautauMass_Signal", diTau.M(), _puevWt);
    AnaUtil::fillHist1D("mvaFK_signal", mvaOutputFK, _puevWt);
    AnaUtil::fillHist1D("mvaVZ_signal", mvaOutput, _puevWt);
    AnaUtil::fillHist1D("tauOSpt_signal", TOS.Pt(), _puevWt);
    AnaUtil::fillHist1D("tauSSpt_signal", TSS.Pt(), _puevWt);
    AnaUtil::fillHist1D("tauOSeta_signal", TOS.Eta(), _puevWt);
    AnaUtil::fillHist1D("tauSSeta_signal", TSS.Eta(), _puevWt);


    if (mvaOutput > 0.65){
      AnaUtil::fillHist1D("tautauMass_Signal_1", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_1", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_1", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.70){
      AnaUtil::fillHist1D("tautauMass_Signal_2", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_2", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_2", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.75){
      AnaUtil::fillHist1D("tautauMass_Signal_3", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_3", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_3", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.78){
      AnaUtil::fillHist1D("tautauMass_Signal_4", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_4", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_4", mvaOutput, _puevWt);
      AnaUtil::fillHist1D("tauOSpt_signal_4", TOS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tauSSpt_signal_4", TSS.Pt(), _puevWt);
      AnaUtil::fillHist1D("tauOSeta_signal_4", TOS.Eta(), _puevWt);
      AnaUtil::fillHist1D("tauSSeta_signal_4", TSS.Eta(), _puevWt);
    }
    if (mvaOutput > 0.79){
      AnaUtil::fillHist1D("tautauMass_Signal_5", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_5", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_5", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.80){
      AnaUtil::fillHist1D("tautauMass_Signal_6", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_6", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_6", mvaOutput, _puevWt);
    }

    if (mvaOutput > 0.81){
      AnaUtil::fillHist1D("tautauMass_Signal_7", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_7", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_7", mvaOutput, _puevWt);
    }

    if (mvaOutput > 0.82){
      AnaUtil::fillHist1D("tautauMass_Signal_8", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_8", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_8", mvaOutput, _puevWt);
    }

    if (mvaOutput > 0.83){
      AnaUtil::fillHist1D("tautauMass_Signal_9", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_9", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_9", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.84){
      AnaUtil::fillHist1D("tautauMass_Signal_10", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_10", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_10", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.85){
      AnaUtil::fillHist1D("tautauMass_Signal_11",diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_11", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_11", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.86){
      AnaUtil::fillHist1D("tautauMass_Signal_12", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_12", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_12", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.87){
      AnaUtil::fillHist1D("tautauMass_Signal_13", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_13", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_13", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.88){
      AnaUtil::fillHist1D("tautauMass_Signal_14", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_14", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_14", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.89){
      AnaUtil::fillHist1D("tautauMass_Signal_15", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_15", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_15", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.90){
      AnaUtil::fillHist1D("tautauMass_Signal_16", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_16", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_16", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.91){
      AnaUtil::fillHist1D("tautauMass_Signal_17", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_17", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_17", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.92){
      AnaUtil::fillHist1D("tautauMass_Signal_18", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_18", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_18", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.93){
      AnaUtil::fillHist1D("tautauMass_Signal_19", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_19", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_19", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.94){
      AnaUtil::fillHist1D("tautauMass_Signal_20", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_20", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_20", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.95){
      AnaUtil::fillHist1D("tautauMass_Signal_21", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_21", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_21", mvaOutput, _puevWt);
    }
    if (mvaOutput > 0.96){
      AnaUtil::fillHist1D("tautauMass_Signal_22", diTau.M(), _puevWt);
      AnaUtil::fillHist1D("mvaFK_signal_22", mvaOutputFK, _puevWt);
      AnaUtil::fillHist1D("mvaVZ_signal_22", mvaOutput, _puevWt);
    }

    //dumpEvent("01010", _evLog, false);
#if 0
    _evLog << setw(15) << "diTauEta: " << setw(8) << diTau.Eta()
           << setw(15) <<  "diTauPt: " << setw(8) << diTau.Pt()
           << setw(15) << "diTauPhi: " << setw(8) << diTau.Phi()
           << setw(15) <<  "ditauMass: " << setw(8) << diTau.M() << endl;

    double d3d = ele.vtxDist3D;
    double dz  = ele.vtxDistZ;
    double dxy = std::sqrt(pow(d3d,2) - pow(dz,2));

    _evLog << setw(8) << "electron Eta: " << setw(8) << ele.eta
           << setw(8) << "electron Pt: " << setw(8) << ele.pt
           << setw(8) << "electron Phi: " << setw(8) << ele.phi
           << setw(8) << "electron MVA: " << setw(8) << ele.mvaPOGNonTrig
           << setw(8) << "Iso" << setw(8) << ele.pfRelIso
           << setw(8) << "dz" << setw(8) << dz
           << setw(8) << "dxy" << setw(8) << dxy << endl;

    const Tau& taua = tau1List.at(0);
    const Tau& taub = tau2List.at(0);
    double eletau1dz = std::abs(taua.vtxDz - ele.vtxDistZ);
    _evLog << setw(8) << "Tau1 Eta: " << setw(8) << taua.eta
           << setw(8) << "Tau1 Pt: " << setw(8) << taua.pt
           << setw(8) << "Tau1 Phi: " << setw(8) << taua.phi
           << setw(8) << "drtau1e: " << setw(8) << e1.DeltaR(T1)
           << setw(8) << "dztau1e: " << setw(8) << eletau1dz << endl;

    double eletau2dz = std::abs(taub.vtxDz - ele.vtxDistZ);
    _evLog << setw(8) << "Tau2 Eta: " << setw(8) << taub.eta
	   << setw(8) << "Tau2 Pt: " << setw(8) << taub.pt
	   << setw(8) << "Tau2 Phi: " << setw(8) << taub.phi
	   << setw(8) << "drTau2E: " << setw(8) << e1.DeltaR(T)
	   << setw(8) << "drtau1Tau2: " << setw(8) << T1.DeltaR(T)
	   << setw(8) << "dztau2e: " << setw(8) << eletau2dz << endl;
#endif
  }
}
// ------------------------------------------------------------------
// Analysis is over, print summary, save histograms release resources
// ------------------------------------------------------------------
void FakeRate::endJob() 
{  
  histf()->cd();

  histf()->Write();
  histf()->Close();
  delete histf();

  _fLog << resetiosflags(ios::fixed);
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
bool FakeRate::readJob(const string& jobFile, int& nFiles)
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
  hmap.insert(pair<string, map<string, double>* >("electronCutList", &_electronCutMap));

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

    if (key == "evselCutList" || key == "electronCutList")
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
void FakeRate::printJob(ostream& os) const
{
  AnaBase::printJob(os);

  map<string, map<string, double> > hmap;
  hmap.insert(pair<string, map<string, double> >("evselCutList", _evselCutMap));
  // hmap.insert(pair<string, map<string, double> >("electronCutList", _electronCutMap));
  AnaUtil::showCuts(hmap, os);
}
