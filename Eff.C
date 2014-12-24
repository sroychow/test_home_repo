Eff(){
  TFile *f = TFile::Open("test.root");
  TH1D *h1 = dynamic_cast<TH1D*>(f->Get("evcounter"));
  const string cuts[] =
   {
     "GenFilter",
     "//Trigger",
     "vtx",
     "e collection size > 0",
     "e pt",
     "e eta",
     "trkHits",
     "pixHits",
     "loose e id",
     "trkD0",
     "dz",
     "Iso",
     "#e >1",
     "drEleTau",
     "tau pt",
     "tau eta",
     "dmf",
     "chargedIsoPtSum < 2",
     "againstMuonTight",
     "againstElectronLoose",
     "vtxDz",
     "OS Tau",
     "#osTau > 1",
     "drEleTau",
     "dr(tau,tau)",
     "tau pt",
     "tau eta",
     "dmf",
     "chargedIsoPtSum < 1",
     "againstMuonTight",
     "againstElectronLoose",
     "vtxDz",
     " SS Tau",
     "#ssTau > 1",
     "highestPt > 25",
     "//b-jet veto",
     "mu veto",
     "e veto",
     "//DYto ee veto",
     "//DY to taus veto",
     "mt(e, met)",
     "met"
   };


  int nbins = h1->GetNbinsX();
  std::cout  << " nbins = " << nbins << std::endl;
  std::cout  << std::setw(30) << "Cuts" <<  std::setw(16) << "Survived Events" << std::setw(12) << "Rel Eff" << std::setw(12) << "Abs Eff" << std::endl;
  std::cout  <<  std::setw(30) << "Initial Eevents" << std::setw(16) << h1->GetBinContent(1) << std::endl;
  for (int i=2; i<= nbins; ++i){
    int binval = h1->GetBinContent(i);
    int binvalpre = h1->GetBinContent(i-1);
    if(binvalpre!=0)
      std::cout  << std::setw(30) << cuts[i-2] << std::setw(16)
      		 << resetiosflags(ios::fixed) << setprecision(3)
      		 << binval << std::setw(12) 
      		 << (double)binval/binvalpre <<  std::setw(12) 
      		 << (double)binval/h1->GetBinContent(2) << std::endl;
      //std::cout << cuts[i-2] <<  std::endl; // cuts
      //std::cout  << binval << std::endl; // survived events
      //std::cout  << resetiosflags(ios::fixed) << setprecision(3)  << (double)binval/binvalpre << std::endl; // rel eff
      //std::cout  << resetiosflags(ios::fixed) << setprecision(3)  << (double)binval/h1->GetBinContent(2) << std::endl; // abs eff
  }
}
