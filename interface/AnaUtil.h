#ifndef __ANAUTIL__HH
#define __ANAUTIL__HH

#include <climits>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <string>
#include <vector>
#include <map>

#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"


using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::setw;

using std::string;
using std::vector;
using std::pair;
using std::map; 
using std::abs;
using std::sqrt;

namespace AnaUtil {

  //public:
  static inline 
  void tokenize(const string& str, vector<string>& tokens, const string& delimiters=" ") {
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);

    // Find first "non-delimiter".
    string::size_type pos = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)  {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));

      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delimiters, pos);

      // Find next "non-delimiter"
      pos = str.find_first_of(delimiters, lastPos);
    }
  }
  static inline 
  void bit_print(int value, int pos=32, ostream& os=cout) {
    static const int INT_BIT = 4*CHAR_BIT; 
    int i, mask = 1 << 31; 
     
    if (pos > INT_BIT) i = INT_BIT; 
    for (i = 1; i <= (INT_BIT - pos); ++i) { 
      value <<= 1; 
    } 
    os.put(' ');
    for (i = 1; i <= pos; ++i) { 
      os.put(((value & mask) == 0) ? '0' : '1'); 
      value <<= 1; 
      if ((INT_BIT - pos + i) % CHAR_BIT ==0 && i != INT_BIT) os.put(' '); 
    } 
    os << endl; 
  }
  static inline
  double deltaPhi(double phia, double phib) {
    double dphi = abs(phia - phib);
    if (dphi > TMath::Pi()) dphi = 2 * TMath::Pi() - dphi;
    return dphi;
  }
  static inline 
  double deltaPhi(const TLorentzVector &a, const TLorentzVector& b) {
    return deltaPhi(a.Phi(), b.Phi());
  }
  static inline 
  double deltaR(const TLorentzVector &a, const TLorentzVector& b) {
    double dphi = deltaPhi(a,b);
    double deta = a.Eta() - b.Eta();
    return sqrt(dphi * dphi + deta * deta);
  }
  static inline double cutValue(const map<string, double>& m, string cname) {
    if (m.find(cname) == m.end()) {
      cerr << ">>> key: " << cname << " not found in the map!" << endl;
      for (auto jt  = m.begin(); jt != m.end(); ++jt)  
        cerr << jt->first << ": " << setw(7) << jt->second << endl;
    }
    //assert(m.find(cname) != m.end());
    return m.find(cname)->second;
  }
  static inline void buildList(const vector<string>& tokens, vector<string>& list) {
    for (vector<string>::const_iterator it  = tokens.begin()+1; 
                                        it != tokens.end(); ++it) {
       list.push_back(*it);       
    }
  }
  static inline void buildMap(const vector<string>& tokens, map<string, int>& hmap) {
    string key = tokens.at(1) + "-" + tokens.at(2) + "-" + tokens.at(3);
    hmap.insert(pair<string, int>(key, 1));
  }
  static inline void storeCuts(const vector<string>& tokens, map<string, map<string, double>* >& hmap) {
    string key = tokens.at(0);
    map<string, map<string, double>* >::const_iterator pos = hmap.find(key);
    if (pos != hmap.end()) {
      map<string, double>* m = pos->second;        
      for (vector<string>::const_iterator it  = tokens.begin()+1; 
                                          it != tokens.end(); ++it) {
        // Split the line into words
        vector<string> cutstr;
        AnaUtil::tokenize(*it, cutstr, "=");
        if (cutstr.size() < 2) continue;
        m->insert( pair<string,double>(cutstr.at(0), atof(cutstr.at(1).c_str())));
      }
    }    
  }
  static inline void showCuts(const map<string, map<string, double> >& hmap, ostream& os=cout) {
    for (map<string, map<string, double> >::const_iterator it  = hmap.begin(); 
                                                           it != hmap.end(); ++it)  
    {
      os << ">>> " << it->first << endl; 
      map<string, double> m = it->second;
      os << std::setprecision(2);
      for (map<string,double>::const_iterator jt  = m.begin(); 
                                              jt != m.end(); ++jt)  
        os << setw(16) << jt->first << ": " 
           << setw(7) << jt->second << endl;
    }
  }
  // ------------------------------------------------------------------------
  // Convenience routine for filling 1D histograms. We rely on root to keep 
  // track of all the histograms that are booked all over so that we do not 
  // have to use any global variables to save the histogram pointers. Instead, 
  // we use the name of the histograms and gROOT to retrieve them from the 
  // Root object pool whenever necessary. This is the closest one can go to 
  // hbook and ID based histogramming
  // -------------------------------------------------------------------------
  static inline
  TH1* getHist1D(const char* hname) {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 
    if (!obj) {
      cerr << "**** getHist1D: Histogram for <" << hname << "> not found!" << endl;
      return 0;
    }
    TH1 *h = 0;
    if (obj->InheritsFrom("TH1D"))
      h = dynamic_cast<TH1D*>(obj);
    else if (obj->InheritsFrom("TH1C"))
      h = dynamic_cast<TH1C*>(obj);
    else if (obj->InheritsFrom("TH1K"))
      h = dynamic_cast<TH1K*>(obj);
    else if (obj->InheritsFrom("TH1S"))
      h = dynamic_cast<TH1S*>(obj);
    else if (obj->InheritsFrom("TH1I"))
      h = dynamic_cast<TH1I*>(obj);
    else
      h = dynamic_cast<TH1F*>(obj);

    if (!h) {
      cerr << "**** getHist1D: <" << hname << "> may not be a 1D Histogram" << endl;
      return 0;
    }
    return h;
  }
  static inline
  TH1* getHist1D(const string& hname) {
    return getHist1D(hname.c_str());
  }

  template <class T>
  static inline bool fillHist1D(const char* hname, T value, double w=1.0) 
  {
    TH1* h = getHist1D(hname);
    if (!h) return false;
    h->Fill(value, w);
    return true;
  }
  template <class T>
  static inline bool fillHist1D(const string& hname, T value, double w=1.0) {
    return fillHist1D(hname.c_str(), value, w);
  }
  // ---------------------------------------------
  // Convenience routine for filling 2D histograms
  // ---------------------------------------------
  static inline 
  TH2* getHist2D(const char* hname) 
  {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 
    if (!obj) {
      cerr << "**** getHist2D: Histogram for <" << hname << "> not found!" << endl;
      return 0;
    }

    TH2 *h = 0;
    if (obj->InheritsFrom("TH2D"))
      h = dynamic_cast<TH2D*>(obj);
    else if (obj->InheritsFrom("TH2C"))
      h = dynamic_cast<TH2C*>(obj);
    else if (obj->InheritsFrom("TH2S"))
      h = dynamic_cast<TH2S*>(obj);
    else if (obj->InheritsFrom("TH2I"))
      h = dynamic_cast<TH2I*>(obj);
    else
      h = dynamic_cast<TH2F*>(obj);

    if (!h) {
      cerr << "**** fillHist2D: <<" << hname << ">> may not be a 2D Histogram" << endl;
      return 0;
    }
    return h;
  }
  static inline
  TH2* getHist2D(const string& hname) {
    return getHist2D(hname.c_str());
  }
  template <class T1, class T2>
  static inline bool fillHist2D(const char* hname, T1 xvalue, T2 yvalue, double w=1.0) 
  {
    TH2* h = getHist2D(hname);
    if (!h) return false;
    h->Fill(xvalue, yvalue, w);
    return true;
  }
  template <class T1, class T2>
  static inline bool fillHist2D(const string& hname, T1 xvalue, T2 yvalue, double w=1.0) 
  {
    return fillHist2D(hname.c_str(), xvalue, yvalue, w);
  }
  // ---------------------------------------------
  // Convenience routine for filling 3D histograms
  // ---------------------------------------------

  static inline 
  TH3* getHist3D(const char* hname) 
  {
    TObject *obj = gDirectory->GetList()->FindObject(hname); 
    if (!obj) {
      cerr << "**** getHist3D: Histogram for <" << hname << "> not found!" << endl;
      return 0;
    }

    TH3* h = 0;
    if (obj->InheritsFrom("TH3D"))
      h = dynamic_cast<TH3D*>(obj);
    else if (obj->InheritsFrom("TH3C"))
      h = dynamic_cast<TH3C*>(obj);
    else if (obj->InheritsFrom("TH3S"))
      h = dynamic_cast<TH3S*>(obj);
    else if (obj->InheritsFrom("TH3I"))
      h = dynamic_cast<TH3I*>(obj);
    else
      h = dynamic_cast<TH3F*>(obj);

    if (!h) {
      cerr << "**** fillHist3D: <<" << hname << ">> may not be a 3D Histogram" << endl;
      return 0;
    }
    return h;
  }
  static inline
  TH3* getHist3D(const string& hname) {
    return getHist3D(hname.c_str());
  }
  template <class T1, class T2, class T3>
  static inline bool fillHist3D(const char* hname, T1 xvalue, T2 yvalue, T3 zvalue, double w=1.0) 
  {
    TH3* h = getHist3D(hname);
    if (!h) return false;
    h->Fill(xvalue, yvalue, zvalue, w);
    return true;
  }
  template <class T1, class T2, class T3>
  static inline bool fillHist3D(const string& hname, T1 xvalue, T2 yvalue, T3 zvalue, double w=1.0) 
  {
    return fillHist3D(hname.c_str(), xvalue, yvalue, zvalue, w);
  }

  // --------------------------------------------------
  // Convenience routine for filling profile histograms
  // --------------------------------------------------
  static inline
  TProfile* getProfile(const char* hname) 
  {
    TProfile *h = dynamic_cast<TProfile*>(gDirectory->GetList()->FindObject(hname));
    if (!h) {
      cerr << "**** getProfile: Profile Histogram <" << hname << "> not found" << endl;
      return 0;
    }
    return h;
  }
  static inline
  TProfile* getProfile(const string& hname) 
  {
    return getProfile(hname.c_str());
  }

  static inline 
  bool fillProfile(const char *hname, float xvalue, float yvalue, double w=1.0) 
  {
    TProfile *h = getProfile(hname);
    if (!h) return false;

    h->Fill(xvalue, yvalue, w);
    return true;
  }
  static inline 
  bool fillProfile(const string& hname, float xvalue, float yvalue, double w=1.0) 
  {
    return fillProfile(hname.c_str(), xvalue, yvalue, w);
  }
  /* print_list_elements()
   * - prints optional C-string optcstr followed by
   * - all elements of the collection coll
   */
  template <class T>
  static inline void showList(const T& coll, const char* optcstr="", ostream& os=cout) 
  {
    os << optcstr << ", Total # = " << coll.size() << ":" << endl;
    for (typename T::const_iterator pos = coll.begin(); pos != coll.end(); ++pos)
      os << *pos << endl;
  }

  template <class T1, class T2>
  static inline void showMap(const map<T1,T2>& m, const char* optcstr="", ostream& os=cout) 
  {
    os << optcstr << endl;
    for (typename map<T1,T2>::const_iterator it  = m.begin(); 
                                             it != m.end(); ++it)  
    {
      os << it->first << endl;
    }
  }
  /* copyList */
  template <class T>
  static inline void copyList (const T& sourceColl, T& destColl)  
  {
    destColl.clear();
    for (typename T::const_iterator pos = sourceColl.begin(); pos != sourceColl.end(); ++pos)   
      destColl.push_back(*pos); 
  }
};
#endif
