lltEff(){
  TFile *f = TFile::Open("test_eltau.root");
  TH1D *h1 = dynamic_cast<TH1D*>(f->Get("counter_emt"));
  const string cuts[] =
   {
     "GenFilter",
     "//Trigger",
     "vtx",
     "ne > 0",
     "nmu > 0",
     "ntau > 0",
     "mu veto",
     "e veto",
     "highestPt < 20",
     "same charge L",
     "opp charge tau",
     "dr",
     "anti el",
     "anti mu",
   };


  int nbins = h1->GetNbinsX();
  std::cout  << " nbins = " << nbins << std::endl;
  std::cout  << std::setw(30) << "Cuts" <<  std::setw(16) << "Survived Events" << std::setw(12) << "Rel Eff" << std::setw(12) << "Abs Eff" << std::endl;
  std::cout  <<  std::setw(30) << "Initial Eevents" << std::setw(16) << h1->GetBinContent(1) << std::endl;
  for (int i=2; i<= nbins; ++i){
    int binval = h1->GetBinContent(i);
    int binvalpre = h1->GetBinContent(i-1);
    if(binvalpre!=0)
      //std::cout  << std::setw(30) << cuts[i-2] << std::setw(16)
      //		 << resetiosflags(ios::fixed) << setprecision(3)
      //		 << binval << std::setw(12) 
      // 		 << (double)binval/binvalpre <<  std::setw(12) 
      // 		 << (double)binval/h1->GetBinContent(2) << std::endl;
      std::cout << cuts[i-2] <<  std::endl; // cuts
      //std::cout  << binval << std::endl; // survived events
      //std::cout  << resetiosflags(ios::fixed) << setprecision(3)  << (double)binval/binvalpre << std::endl; // rel eff
      //std::cout  << resetiosflags(ios::fixed) << setprecision(3)  << (double)binval/h1->GetBinContent(2) << std::endl; // abs eff
  }
}
