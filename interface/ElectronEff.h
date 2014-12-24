#ifndef __ElectronEff__hh
#define __ElectronEff__hh

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


class ElectronEff : public AnaBase {
    
public:

  ElectronEff();
  virtual ~ElectronEff();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();

  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout) const;
  
  void clearLists();
  void findEleIDInfo();
  void computeEleEff();

  void selectEvent();
  virtual void bookHistograms();

public:
  int nProbe[6];
  int nSingleCut[6];

  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::Muon> muonList;
  std::vector<vhtm::Tau> tauList;
  std::vector<vhtm::Jet> bjetList;
  std::vector<vhtm::TriggerObject> trigObjList;
  std::vector<vhtm::Electron> probeEleList;
  std::vector<vhtm::Electron> tagEleList;

public:
  std::map<std::string, double> _evselCutMap;
  std::vector<std::string> _triggerPathTagList;
  bool _dumpEvent;
};
#endif
