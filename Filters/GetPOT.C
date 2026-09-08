double GetPOT(std::string file,bool make_rs_list=false){

  TFile* f_in = TFile::Open(file.c_str());

  // Pandora
  TTree* treein_pd = (TTree*)f_in->Get("nuselection/SubRun");

  Float_t         pot;
  Int_t           run;
  Int_t           subRun;
  treein_pd->SetBranchAddress("pot",&pot);
  treein_pd->SetBranchAddress("run",&run);
  treein_pd->SetBranchAddress("subRun",&subRun);

  std::ofstream rs_file;
  if(make_rs_list){
    std::string outname = file.substr(file.find_last_of("/\\")+1);
    outname = outname.substr(0,outname.find_last_of('.')) + "_runsubrun.txt";
    rs_file.open(outname);
  }

  
  double totalpot = 0.0;
  for(size_t ientry=0;ientry<treein_pd->GetEntries();ientry++){
    treein_pd->GetEntry(ientry);
    totalpot += pot;
    if(make_rs_list) rs_file << run << " " << subRun << "\n";
  }

  if(make_rs_list) rs_file.close();

  std::cout << "Pandora POT=" << totalpot << std::endl;

  // Wirecell
  TTree* treein_wc = (TTree*)f_in->Get("wcpselection/T_pot");
  Double_t pot_tor875good;

  treein_wc->SetBranchAddress("pot_tor875good",&pot_tor875good);

  totalpot = 0.0;
  for(size_t ientry=0;ientry<treein_wc->GetEntries();ientry++){
    treein_wc->GetEntry(ientry);
    totalpot += pot_tor875good;
  }

  std::cout << "Wirecell POT=" << totalpot << std::endl;

  f_in->Close();

  return totalpot;

}
