#ifndef __ETauTau__HH
#define __ETauTau__HH

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

class MVASkim;
namespace TMVA {
  class Reader;
  }
class ETauTau : public AnaBase {
    
public:

  ETauTau();
  virtual ~ETauTau();
    
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
  TLorentzVector EGen, TauOSGen, TauSSGen;

  std::vector<vhtm::Vertex> vtxList;
  std::vector<vhtm::Jet> bjetList;
  std::vector<vhtm::TriggerObject> trigObjList;
  std::vector<vhtm::GenParticle> genEleList;
  std::vector<vhtm::GenParticle> genTauList;
  
  std::vector<vhtm::GenParticle> genList;
  std::vector<vhtm::GenParticle> genEleListtemp; //temp

public:
  std::map<std::string, double> _evselCutMap;

  //MVA input variables while reading
  float lepPt;
  float tauOSPt;
  float tauSSPt;
  float dphilepDiTau;
  float met;
  float PtRatio;
  float mass;

  bool _createMVATree;
  bool _readMVA;
  std::string _mvaInputFile;
  MVASkim* _skimObj;
  TMVA::Reader* reader;

  //MVA input variables while reading
  float leadPt;
  float subPt;
  //double met;
  float deltaRDiTau;
  float ptRatio;

  TMVA::Reader* reader1;
};
#endif
