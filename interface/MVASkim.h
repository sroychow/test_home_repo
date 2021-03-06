#ifndef __MVASkim__h
#define __MVASkim__h

#include <fstream>
#include <string>

class TTree;
class TFile;

typedef struct  
{
  float lepEta;
  float lepPt;
  float tauOSEta;
  float tauOSPt;
  float tauSSEta;
  float tauSSPt;
  float diTauEta;
  float diTauPt;
  float dphilepTauOS;
  float met;
  float alpha;
  float acop;
  float dphilepDiTau;
  float DeltaRDiTau;
  float PtRatio;

  float alphatt;
  float alphalditau;
  float alphalprod;

} TreeVariables;

class MVASkim {
    
public:

  MVASkim(const std::string& filename);
  virtual ~MVASkim();

  void fill(const TreeVariables& varList);
  void close();

  TFile* _mvaFile;
  TTree* _tree;

  TreeVariables _varList;
};
#endif
