#ifndef __SelectionEfficiency__HH
#define __SelectionEfficiency__HH

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

class MVASkim; // forward
namespace TMVA {
  class Reader;
}

class SelectionEfficiency : public AnaBase {
    
public:

  SelectionEfficiency();
  virtual ~SelectionEfficiency();
    
  void eventLoop();  // the main analysis 
  bool beginJob();
  void endJob();
  void selectEvent();

  bool readJob(const std::string& jobFile, int& nFiles);
  void printJob(std::ostream& os=std::cout) const;
  
  void clearLists();

  virtual void bookHistograms();
  void study_eff();

public:

  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::Muon> muoList;
  std::vector<vhtm::Electron> eleList;
  std::vector<vhtm::Tau> tauList;
  std::vector<vhtm::Jet> bjetList;

  std::vector<vhtm::GenParticle> genMuonList;
  std::vector<vhtm::GenParticle> genTauList;

public:
  std::map<std::string, double> _tau1CutMap;
  std::map<std::string, double> _tau2CutMap;
  bool _createMVATree;
  bool _readMVA;
  bool _readMVAFK;
  std::string _mvaInputFile;
  std::string _MVAdisFile;
  std::string _MVAFKdisFile;

  // MVA input variables while reading
  float muEta;
  float muPt;
  float tau1Eta;
  float tau1Pt;
  float tau2Eta;
  float tau2Pt;
  float diTaudR;
  float dphiMuTau1;
  float dphiMuDiTau;
  float met;
  float ptRatio;

  // MVA input variables while reading for fake mva
  float leadPt;
  float subPt;
  float deltaRDiTau;

  MVASkim* _skimObj;
  TMVA::Reader* reader;
  TMVA::Reader* reader1;
};
#endif
