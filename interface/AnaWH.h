#ifndef __ANAWH__HH
#define __ANAWH__HH

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"

#include "PhysicsObjects.h"
#include "AnaBase.h"

using namespace vhtm;

class AnaWH : public AnaBase {
    
public:

  AnaWH();
  virtual ~AnaWH();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();

   // void setData(int val);  
  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout);
  
  void clearLists();
  //void findVtxInfo();
  void findGenInfo();
  //void findElectronInfo();
  //void findMuonInfo();
  //void findTauInfo();
  //void findJetInfo();
  void dumpGenInfo();
  
  void selectMuMuTauEvent(const std::vector<Muon>& muoList, const std::vector<Tau>& tauList, const std::vector<Jet>& bjetList);
  void selectMuTauTauEvent(const std::vector<Muon>& muoList, const std::vector<Tau>& tauList, const std::vector<Jet>& bjetList);
  void computeSignificance(const std::vector<Muon>& muoList, const std::vector<Tau>& tauList, const std::vector<Jet>& bjetList);
  //  void GenDetLevelMatch (const std::vector<Muon>& muoList, const std::vector<GenParticle>& genMuonList);
  void computeDeltaR(const std::vector<Muon>& muoList, const std::vector<Tau>& tauList );
  void selectMuMuTauTauEvent(const std::vector<Muon>& muoList, const std::vector<Tau>& tauList, const std::vector<Jet>& bjetList);
  void selectMuMuMuTauEvent(const std::vector<Muon>& muoList, const std::vector<Tau>& tauList, const std::vector<Jet>& bjetList);

  virtual bool selectEvent();
  virtual void bookHistograms();
  void fillTriggerHistograms(int run); 

  template <class T>
  bool fillHist1D(const std::string& hname,  T value, double w=1.0);
  template <class T1, class T2>
  bool fillHist2D(const std::string& hname,  T1 xvalue, T2 yvalue, double w=1.0);
  bool fillProfile(const std::string& hname, float xvalue, float yvalue, double w=1.0);
  static TH1* getHist1D(const std::string& hname);
  //  static TH2* getHist2D(const string& hname);
  static TProfile* getProfile(const std::string& hname);
 
public:
  int nEvtMuMuTau[11];
  int nEvtMuMuTauTau[11];
  int nEvtMuMuMuTau[10];
  int nEvtMuTauTau[9];
  int nSig[20];
  int count_sig_event;

  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::Muon> muoList;
  std::vector<vhtm::Electron> eleList;
  std::vector<vhtm::Tau> tauList;
  std::vector<vhtm::Jet> bjetList;
  std::vector<vhtm::GenParticle> genMuonList;
  std::vector<vhtm::GenParticle> genTauList;
  std::vector<vhtm::GenParticle> GEN_MUON_LIST;
  std::vector<vhtm::GenParticle> genMuoaList;
  std::vector<vhtm::GenParticle> genMuobList;
  std::vector<vhtm::GenParticle> genHList;
  std::vector<vhtm::GenParticle> genWList;
  std::vector<vhtm::GenParticle> genTau_H_List;
  std::vector<vhtm::GenParticle> genMuon_Z_List;

public:
  std::map<std::string, double> _vtxCutMap;
  std::map<std::string, double> _muonCutMap;
  std::map<std::string, double> _electronCutMap;
  std::map<std::string, double> _tauCutMap;
  std::map<std::string, double> _bjetCutMap;
};
#endif
