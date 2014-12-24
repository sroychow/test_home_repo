#ifndef __TAUEFFICIENCY__HH
#define __TAUEFFICIENCY__HH

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

class TauEfficiency : public AnaBase {
    
public:

  TauEfficiency(const std::string& filename="pippo.out");
  virtual ~TauEfficiency();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();

   // void setData(int val);  
  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout);
  
  void clearLists();
  //void findVtxInfo();
  //void findMuonIDInfo();
  
  //void computeMuonEff(const std::vector<Muon>& probemuoList);
  void computeTauEff(const std::vector<vhtm::Muon>& muoList,const std::vector<vhtm::Tau>& tauList);

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
  //int nProbe[14];
  //int nSingleCut[14];
  int nSteps[6];

  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::Muon> muoList;
  std::vector<vhtm::Tau> tauList;
  std::vector<vhtm::Jet> bjetList;
  std::vector<vhtm::Electron> eleList;
  //std::vector<Muon> probemuoList;
  //std::vector<Muon> tagmuoList;

public:
  std::map<std::string, double> _vtxCutMap;
  std::map<std::string, double> _muonCutMap;
  std::map<std::string, double> _electronCutMap;
  std::map<std::string, double> _tauCutMap;
  std::map<std::string, double> _bjetCutMap;
};
#endif
