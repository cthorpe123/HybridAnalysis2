// Merge the three trees from surprise into one enabling
// event by event combining of information from each reconstruction

void SpeedyMergeNtuple(){

  bool is_data = false;  
  bool save_syst = false; 

  std::string dir_in = "/pnfs/uboone/persistent/users/uboonepro/surprise/detvar_test/";
  std::string filename = "DetVar_Run45_v10_04_07_15_BNB_nu_overlay_SCE_surprise_reco2_hist.root";
  std::string dir_out = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/detvar/";

  //std::string dir_in = "/exp/uboone/data/uboonepro/MCC9.10/liangliu/v10_04_07_09/";
  //std::string filename = "MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root";
  //std::string filename = "MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root";
  //std::string filename = "MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  //std::string filename = "MCC9.10_Run4b_v10_04_07_09_BNB_nue_overlay_surprise_reco2_hist.root";
  //std::string dir_out = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";

  //std::string dir_in = "/pnfs/uboone/scratch/users/bbogart/v10_04_07_16/";
  //std::string filename = "larpid_patch_smart_patch_test10_full_more.root";
  //std::string dir_out = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";

  // Open the two files and setup branches to read the RSE numbers
  TFile* f_in = new TFile((dir_in + filename).c_str());
  TTree* pd_t_in = static_cast<TTree*>(f_in->Get("nuselection/NeutrinoSelectionFilter"));
  TTree* wc_t_in = static_cast<TTree*>(f_in->Get("wcpselection/T_PFeval"));
  TTree* lt_t_in = static_cast<TTree*>(f_in->Get("lantern/EventTree"));

  int pd_rse[3];
  int wc_rse[3];
  int lt_rse[3];

  pd_t_in->SetBranchAddress("run",&pd_rse[0]);
  pd_t_in->SetBranchAddress("sub",&pd_rse[1]);
  pd_t_in->SetBranchAddress("evt",&pd_rse[2]);

  wc_t_in->SetBranchAddress("run",&wc_rse[0]);
  wc_t_in->SetBranchAddress("subrun",&wc_rse[1]);
  wc_t_in->SetBranchAddress("event",&wc_rse[2]);

  lt_t_in->SetBranchAddress("run",&lt_rse[0]);
  lt_t_in->SetBranchAddress("subrun",&lt_rse[1]);
  lt_t_in->SetBranchAddress("event",&lt_rse[2]);

  if(pd_t_in->GetEntries() != wc_t_in->GetEntries() || wc_t_in->GetEntries() != lt_t_in->GetEntries())
    throw std::invalid_argument("Can't use speedy merging - trees don't have equal numbers of entries");

  // Pandora branches 
  std::vector<int>* mc_pdg=0;
  std::vector<double>* mc_E=0;
  std::vector<double>* mc_px=0;
  std::vector<double>* mc_py=0;
  std::vector<double>* mc_pz=0;
  std::vector<double>* mc_endx=0;
  std::vector<double>* mc_endy=0;
  std::vector<double>* mc_endz=0;

  float true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z;
  int nu_pdg,ccnc,interaction;
  float nu_e; 
  std::vector<unsigned short>* weightsGenie=0; 
  std::vector<unsigned short>* weightsReint=0; 
  std::vector<unsigned short>* weightsFlux=0; 
  Float_t weightSpline;
  Float_t weightSplineTimesTune;
  Float_t weightTune;

  std::vector<float>* trk_len_v=0;
  std::vector<float>* trk_dir_x_v=0;
  std::vector<float>* trk_dir_y_v=0;
  std::vector<float>* trk_dir_z_v=0;
  std::vector<float>* trk_range_muon_mom_v=0;
  std::vector<float>* trk_mcs_muon_mom_v=0;
  std::vector<float>* trk_llr_pid_score_v=0;
  std::vector<float>* trk_start_x_v=0;
  std::vector<float>* trk_start_y_v=0;
  std::vector<float>* trk_start_z_v=0;
  std::vector<float>* trk_end_x_v=0;
  std::vector<float>* trk_end_y_v=0;
  std::vector<float>* trk_end_z_v=0;
  std::vector<int>* backtracked_pdg=0;
  std::vector<int>* pfng2semlabel=0; 

  std::vector<float>* shr_px_v=0;
  std::vector<float>* shr_py_v=0;
  std::vector<float>* shr_pz_v=0;
  std::vector<float>* shr_energy_y_v=0;

  float reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z;

  if(!is_data){
    pd_t_in->SetBranchAddress("mc_pdg",&mc_pdg);
    pd_t_in->SetBranchAddress("mc_E",&mc_E);
    pd_t_in->SetBranchAddress("mc_px",&mc_px);
    pd_t_in->SetBranchAddress("mc_py",&mc_py);
    pd_t_in->SetBranchAddress("mc_pz",&mc_pz);
    pd_t_in->SetBranchAddress("mc_endx",&mc_endx);
    pd_t_in->SetBranchAddress("mc_endy",&mc_endy);
    pd_t_in->SetBranchAddress("mc_endz",&mc_endz);
    pd_t_in->SetBranchAddress("true_nu_vtx_x",&true_nu_vtx_x);
    pd_t_in->SetBranchAddress("true_nu_vtx_y",&true_nu_vtx_y);
    pd_t_in->SetBranchAddress("true_nu_vtx_z",&true_nu_vtx_z);
    pd_t_in->SetBranchAddress("nu_pdg",&nu_pdg);
    pd_t_in->SetBranchAddress("interaction",&interaction);
    pd_t_in->SetBranchAddress("ccnc",&ccnc);
    pd_t_in->SetBranchAddress("nu_e",&nu_e); 
    if(save_syst){
      pd_t_in->SetBranchAddress("weightsGenie",&weightsGenie); 
      pd_t_in->SetBranchAddress("weightsReint",&weightsReint); 
      pd_t_in->SetBranchAddress("weightsFlux",&weightsFlux); 
    }
    pd_t_in->SetBranchAddress("weightSpline",&weightSpline); 
    pd_t_in->SetBranchAddress("weightTune",&weightTune); 
    pd_t_in->SetBranchAddress("weightSplineTimesTune",&weightSplineTimesTune); 

  }

  pd_t_in->SetBranchAddress("trk_len_v",&trk_len_v);
  pd_t_in->SetBranchAddress("trk_range_muon_mom_v",&trk_range_muon_mom_v);
  pd_t_in->SetBranchAddress("trk_mcs_muon_mom_v",&trk_mcs_muon_mom_v);
  pd_t_in->SetBranchAddress("trk_dir_x_v",&trk_dir_x_v);
  pd_t_in->SetBranchAddress("trk_dir_y_v",&trk_dir_y_v);
  pd_t_in->SetBranchAddress("trk_dir_z_v",&trk_dir_z_v);
  pd_t_in->SetBranchAddress("trk_start_x_v",&trk_start_x_v);
  pd_t_in->SetBranchAddress("trk_start_y_v",&trk_start_y_v);
  pd_t_in->SetBranchAddress("trk_start_z_v",&trk_start_z_v);
  pd_t_in->SetBranchAddress("trk_end_x_v",&trk_end_x_v);
  pd_t_in->SetBranchAddress("trk_end_y_v",&trk_end_y_v);
  pd_t_in->SetBranchAddress("trk_end_z_v",&trk_end_z_v);
  pd_t_in->SetBranchAddress("trk_llr_pid_score_v",&trk_llr_pid_score_v);
  pd_t_in->SetBranchAddress("reco_nu_vtx_x",&reco_nu_vtx_x);
  pd_t_in->SetBranchAddress("reco_nu_vtx_y",&reco_nu_vtx_y);
  pd_t_in->SetBranchAddress("reco_nu_vtx_z",&reco_nu_vtx_z);
  pd_t_in->SetBranchAddress("shr_px_v",&shr_px_v);
  pd_t_in->SetBranchAddress("shr_py_v",&shr_py_v);
  pd_t_in->SetBranchAddress("shr_pz_v",&shr_pz_v);
  pd_t_in->SetBranchAddress("shr_energy_y_v",&shr_energy_y_v);
  if(!is_data) pd_t_in->SetBranchAddress("backtracked_pdg",&backtracked_pdg);
  pd_t_in->SetBranchAddress("pfng2semlabel",&pfng2semlabel);

  // Wirecell branches
  Int_t mc_nu_pdg;
  Int_t mc_nu_ccnc;
  Int_t mc_nu_mode;
  Float_t mc_nu_pos[4];
  Float_t mc_nu_mom[4];

  Int_t truth_Ntrack;
  Int_t truth_pdg[2000];   //[truth_Ntrack]
  Int_t truth_mother[2000];   //[truth_Ntrack]
  Float_t truth_startMomentum[2000][4];   //[truth_Ntrack]
  Float_t truth_endMomentum[2000][4];   //[truth_Ntrack]
  Float_t truth_endXYZT[2000][4];   //[truth_Ntrack]
  Int_t truth_id[2000];   //[truth_Ntrack]

  Float_t reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ;
  Int_t reco_Ntrack;
  Int_t reco_pdg[1000];   //[reco_Ntrack]
  Int_t reco_mother[1000];   //[reco_Ntrack]
  Float_t reco_startMomentum[1000][4];   //[reco_Ntrack]
  Float_t reco_endMomentum[1000][4];   //[reco_Ntrack]
  Float_t reco_endXYZT[1000][4];
  Int_t reco_id[1000];   //[reco_Ntrack]
  Int_t reco_truthMatch_pdg[1000];   //[reco_Ntrack]
  Int_t reco_truthMatch_id[1000];   //[reco_Ntrack]
  Double_t mcs_mu_tracklen;
  Double_t mcs_emu_tracklen;
  Double_t mcs_emu_MCS;
  Double_t mcs_ambiguity_MCS;


  if(!is_data){
    wc_t_in->SetBranchAddress("mc_nu_pdg",&mc_nu_pdg); 
    wc_t_in->SetBranchAddress("mc_nu_ccnc",&mc_nu_ccnc); 
    wc_t_in->SetBranchAddress("mc_nu_mode",&mc_nu_mode); 
    wc_t_in->SetBranchAddress("mc_nu_pos",mc_nu_pos); 
    wc_t_in->SetBranchAddress("mc_nu_mom",mc_nu_mom); 
    wc_t_in->SetBranchAddress("truth_Ntrack",&truth_Ntrack); 
    wc_t_in->SetBranchAddress("truth_pdg",truth_pdg); 
    wc_t_in->SetBranchAddress("truth_mother",truth_mother); 
    wc_t_in->SetBranchAddress("truth_startMomentum",truth_startMomentum); 
    wc_t_in->SetBranchAddress("truth_endMomentum",truth_endMomentum); 
    wc_t_in->SetBranchAddress("truth_endXYZT",truth_endXYZT); 
    wc_t_in->SetBranchAddress("truth_id",&truth_id);
  }

  wc_t_in->SetBranchAddress("reco_nuvtxX",&reco_nuvtxX); 
  wc_t_in->SetBranchAddress("reco_nuvtxY",&reco_nuvtxY); 
  wc_t_in->SetBranchAddress("reco_nuvtxZ",&reco_nuvtxZ); 
  wc_t_in->SetBranchAddress("reco_Ntrack",&reco_Ntrack); 
  wc_t_in->SetBranchAddress("reco_pdg",reco_pdg); 
  wc_t_in->SetBranchAddress("reco_mother",reco_mother); 
  wc_t_in->SetBranchAddress("reco_startMomentum",reco_startMomentum); 
  wc_t_in->SetBranchAddress("reco_endMomentum",reco_endMomentum); 
  wc_t_in->SetBranchAddress("reco_endXYZT",reco_endXYZT); 
  wc_t_in->SetBranchAddress("reco_id",reco_id);
  if(!is_data){
    wc_t_in->SetBranchAddress("reco_truthMatch_pdg",reco_truthMatch_pdg);
    wc_t_in->SetBranchAddress("reco_truthMatch_id",reco_truthMatch_id);
  }
  wc_t_in->SetBranchAddress("mcs_mu_tracklen", &mcs_mu_tracklen);
  wc_t_in->SetBranchAddress("mcs_emu_tracklen", &mcs_emu_tracklen);
  wc_t_in->SetBranchAddress("mcs_emu_MCS", &mcs_emu_MCS);
  wc_t_in->SetBranchAddress("mcs_ambiguity_MCS", &mcs_ambiguity_MCS);

  // Lantern branches
  Float_t         trueNuE;
  Int_t           trueNuPDG;
  Int_t           trueNuCCNC;
  Float_t         trueVtxX;
  Float_t         trueVtxY;
  Float_t         trueVtxZ;
  Int_t           nTruePrimParts;
  Int_t           truePrimPartPDG[200];   //[nTruePrimParts]
  Float_t         truePrimPartPx[200];   //[nTruePrimParts]
  Float_t         truePrimPartPy[200];   //[nTruePrimParts]
  Float_t         truePrimPartPz[200];   //[nTruePrimParts]
  Float_t         truePrimPartE[200];   //[nTruePrimParts]

  Int_t           foundVertex;
  Float_t         vtxX;
  Float_t         vtxY;
  Float_t         vtxZ;
  Int_t           nTracks;
  Int_t           trackIsSecondary[100];   //[nTracks]
  Float_t         trackStartPosX[100];   //[nTracks]
  Float_t         trackStartPosY[100];   //[nTracks]
  Float_t         trackStartPosZ[100];   //[nTracks]
  Float_t         trackStartDirX[100];   //[nTracks]
  Float_t         trackStartDirY[100];   //[nTracks]
  Float_t         trackStartDirZ[100];   //[nTracks]
  Float_t         trackEndPosX[100];   //[nTracks]
  Float_t         trackEndPosY[100];   //[nTracks]
  Float_t         trackEndPosZ[100];   //[nTracks]
  Int_t           trackClassified[100];   //[nTracks]
  Int_t           trackNPlanesAbove[100];   //[nTracks]
  Int_t           trackPID[100];   //[nTracks]
  Float_t         trackRecoE[100];   //[nTracks]
  Int_t           trackTruePID[100];   //[nTracks]

  Int_t           nShowers;
  Int_t           showerIsSecondary[100];   //[nShowers]
  Float_t         showerStartPosX[100];   //[nShowers]
  Float_t         showerStartPosY[100];   //[nShowers]
  Float_t         showerStartPosZ[100];   //[nShowers]
  Float_t         showerStartDirX[100];   //[nShowers]
  Float_t         showerStartDirY[100];   //[nShowers]
  Float_t         showerStartDirZ[100];   //[nShowers]
  Int_t           showerClassified[100];   //[nShowers]
  Int_t           showerNPlanesAbove[100];   //[nShowers]
  Int_t           showerPID[100];   //[nShowers]
  Float_t         showerRecoE[100];   //[nShowers]
  Int_t           showerTruePID[100];   //[nShowers]

  if(!is_data){
    lt_t_in->SetBranchAddress("trueNuE", &trueNuE);
    lt_t_in->SetBranchAddress("trueNuPDG", &trueNuPDG);
    lt_t_in->SetBranchAddress("trueNuCCNC", &trueNuCCNC);
    lt_t_in->SetBranchAddress("nTruePrimParts", &nTruePrimParts);
    lt_t_in->SetBranchAddress("truePrimPartPDG", truePrimPartPDG);
    lt_t_in->SetBranchAddress("truePrimPartPx", truePrimPartPx);
    lt_t_in->SetBranchAddress("truePrimPartPy", truePrimPartPy);
    lt_t_in->SetBranchAddress("truePrimPartPz", truePrimPartPz);
    lt_t_in->SetBranchAddress("truePrimPartE", truePrimPartE);
    lt_t_in->SetBranchAddress("trueVtxX", &trueVtxX);
    lt_t_in->SetBranchAddress("trueVtxY", &trueVtxY);
    lt_t_in->SetBranchAddress("trueVtxZ", &trueVtxZ);
  }

  lt_t_in->SetBranchAddress("foundVertex", &foundVertex);
  lt_t_in->SetBranchAddress("vtxX", &vtxX);
  lt_t_in->SetBranchAddress("vtxY", &vtxY);
  lt_t_in->SetBranchAddress("vtxZ", &vtxZ);

  lt_t_in->SetBranchAddress("nTracks", &nTracks);
  lt_t_in->SetBranchAddress("trackIsSecondary", trackIsSecondary);
  lt_t_in->SetBranchAddress("trackStartPosX", trackStartPosX);
  lt_t_in->SetBranchAddress("trackStartPosY", trackStartPosY);
  lt_t_in->SetBranchAddress("trackStartPosZ", trackStartPosZ);
  lt_t_in->SetBranchAddress("trackStartDirX", trackStartDirX);
  lt_t_in->SetBranchAddress("trackStartDirY", trackStartDirY);
  lt_t_in->SetBranchAddress("trackStartDirZ", trackStartDirZ);
  lt_t_in->SetBranchAddress("trackEndPosX", trackEndPosX);
  lt_t_in->SetBranchAddress("trackEndPosY", trackEndPosY);
  lt_t_in->SetBranchAddress("trackEndPosZ", trackEndPosZ);
  lt_t_in->SetBranchAddress("trackClassified", trackClassified);
  lt_t_in->SetBranchAddress("trackNPlanesAbove", trackNPlanesAbove);
  lt_t_in->SetBranchAddress("trackPID", trackPID);
  lt_t_in->SetBranchAddress("trackRecoE", trackRecoE);
  if(!is_data) lt_t_in->SetBranchAddress("trackTruePID", trackTruePID);

  lt_t_in->SetBranchAddress("nShowers", &nShowers);
  lt_t_in->SetBranchAddress("showerIsSecondary", showerIsSecondary);
  lt_t_in->SetBranchAddress("showerStartPosX", showerStartPosX);
  lt_t_in->SetBranchAddress("showerStartPosY", showerStartPosY);
  lt_t_in->SetBranchAddress("showerStartPosZ", showerStartPosZ);
  lt_t_in->SetBranchAddress("showerStartDirX", showerStartDirX);
  lt_t_in->SetBranchAddress("showerStartDirY", showerStartDirY);
  lt_t_in->SetBranchAddress("showerStartDirZ", showerStartDirZ);
  lt_t_in->SetBranchAddress("showerClassified", showerClassified);
  lt_t_in->SetBranchAddress("showerNPlanesAbove", showerNPlanesAbove);
  lt_t_in->SetBranchAddress("showerPID", showerPID);
  lt_t_in->SetBranchAddress("showerRecoE", showerRecoE);
  if(!is_data) lt_t_in->SetBranchAddress("showerTruePID", showerTruePID);

  // Setup the tree storing the data
  gSystem->Exec(("mkdir -p "+dir_out).c_str());
  TFile* f_out = new TFile((dir_out+"/Merged_"+filename).c_str(),"RECREATE");
  TTree* t_out = new TTree("MergedNtuple","tree");

  t_out->Branch("run",&pd_rse[0]);
  t_out->Branch("sub",&pd_rse[1]);
  t_out->Branch("evt",&pd_rse[2]);

  // Pandora Branches
  if(!is_data){
    t_out->Branch("mc_pdg",&mc_pdg);
    t_out->Branch("mc_E",&mc_E);
    t_out->Branch("mc_px",&mc_px);
    t_out->Branch("mc_py",&mc_py);
    t_out->Branch("mc_pz",&mc_pz);
    t_out->Branch("mc_endx",&mc_endx);
    t_out->Branch("mc_endy",&mc_endy);
    t_out->Branch("mc_endz",&mc_endz);
    t_out->Branch("true_nu_vtx_x",&true_nu_vtx_x);
    t_out->Branch("true_nu_vtx_y",&true_nu_vtx_y);
    t_out->Branch("true_nu_vtx_z",&true_nu_vtx_z);
    t_out->Branch("nu_pdg",&nu_pdg);
    t_out->Branch("interaction",&interaction);
    t_out->Branch("ccnc",&ccnc);
    t_out->Branch("nu_e",&nu_e); 
    if(save_syst){
      t_out->Branch("weightsGenie",&weightsGenie); 
      t_out->Branch("weightsReint",&weightsReint); 
      t_out->Branch("weightsFlux",&weightsFlux); 
    }
    t_out->Branch("weightSpline",&weightSpline);
    t_out->Branch("weightTune",&weightTune);
    t_out->Branch("weightSplineTimesTune",&weightSplineTimesTune);
  }

  t_out->Branch("trk_len_v",&trk_len_v);
  t_out->Branch("trk_range_muon_mom_v",&trk_range_muon_mom_v);
  t_out->Branch("trk_mcs_muon_mom_v",&trk_mcs_muon_mom_v);
  t_out->Branch("trk_dir_x_v",&trk_dir_x_v);
  t_out->Branch("trk_dir_y_v",&trk_dir_y_v);
  t_out->Branch("trk_dir_z_v",&trk_dir_z_v);
  t_out->Branch("trk_start_x_v",&trk_start_x_v);
  t_out->Branch("trk_start_y_v",&trk_start_y_v);
  t_out->Branch("trk_start_z_v",&trk_start_z_v);
  t_out->Branch("trk_end_x_v",&trk_end_x_v);
  t_out->Branch("trk_end_y_v",&trk_end_y_v);
  t_out->Branch("trk_end_z_v",&trk_end_z_v);
  t_out->Branch("trk_llr_pid_score_v",&trk_llr_pid_score_v);
  t_out->Branch("reco_nu_vtx_x",&reco_nu_vtx_x);
  t_out->Branch("reco_nu_vtx_y",&reco_nu_vtx_y);
  t_out->Branch("reco_nu_vtx_z",&reco_nu_vtx_z);
  if(!is_data) t_out->Branch("backtracked_pdg",&backtracked_pdg);
  t_out->Branch("pfng2semlabel",&pfng2semlabel);

  t_out->Branch("shr_px_v",&shr_px_v);
  t_out->Branch("shr_py_v",&shr_py_v);
  t_out->Branch("shr_pz_v",&shr_pz_v);
  t_out->Branch("shr_energy_y_v",&shr_energy_y_v);


  // Wirecell Branches
  if(!is_data){
    t_out->Branch("mc_nu_pdg",&mc_nu_pdg); 
    t_out->Branch("mc_nu_ccnc",&mc_nu_ccnc); 
    t_out->Branch("mc_nu_mode",&mc_nu_mode); 
    t_out->Branch("mc_nu_pos",mc_nu_pos,"mc_nu_pos[4]/F"); 
    t_out->Branch("mc_nu_mom",mc_nu_mom,"mc_nu_mom[4]/F"); 
    t_out->Branch("truth_Ntrack",&truth_Ntrack); 
    t_out->Branch("truth_pdg",truth_pdg,"truth_pdg[2000]/I"); 
    t_out->Branch("truth_mother",truth_mother,"truth_mother[2000]/I"); 
    t_out->Branch("truth_startMomentum",truth_startMomentum,"truth_startMomentum[2000][4]/F"); 
    t_out->Branch("truth_endMomentum",truth_endMomentum,"truth_endMomentum[2000][4]/F"); 
    t_out->Branch("truth_endXYZT",truth_endXYZT,"truth_endXYZT[2000][4]/F"); 
    t_out->Branch("truth_id",&truth_id,"truth_id[2000]/I");
  }

  t_out->Branch("reco_nuvtxX",&reco_nuvtxX); 
  t_out->Branch("reco_nuvtxY",&reco_nuvtxY); 
  t_out->Branch("reco_nuvtxZ",&reco_nuvtxZ); 
  t_out->Branch("reco_Ntrack",&reco_Ntrack); 
  t_out->Branch("reco_pdg",reco_pdg,"reco_pdg[1000]/I"); 
  t_out->Branch("reco_mother",reco_mother,"reco_mother[1000]/I"); 
  t_out->Branch("reco_startMomentum",reco_startMomentum,"reco_startMomentum[1000][4]/F"); 
  t_out->Branch("reco_endMomentum",reco_endMomentum,"reco_endMomentum[1000][4]/F"); 
  t_out->Branch("reco_endXYZT",reco_endXYZT,"reco_endXYZT[1000][4]/F"); 
  t_out->Branch("reco_id",reco_id,"reco_id[1000]/I");
  if(!is_data){
    t_out->Branch("reco_truthMatch_pdg",reco_truthMatch_pdg,"reco_truthMatch_pdg[1000]/I");
    t_out->Branch("reco_truthMatch_id",reco_truthMatch_id,"reco_truthMatch_id[1000]/I");
  }
  t_out->Branch("mcs_mu_tracklen", &mcs_mu_tracklen);
  t_out->Branch("mcs_emu_tracklen", &mcs_emu_tracklen);
  t_out->Branch("mcs_emu_MCS", &mcs_emu_MCS);
  t_out->Branch("mcs_ambiguity_MCS", &mcs_ambiguity_MCS);

  // Lantern branches
  if(!is_data){
    t_out->Branch("trueNuE", &trueNuE);
    t_out->Branch("trueNuPDG", &trueNuPDG);
    t_out->Branch("trueNuCCNC", &trueNuCCNC);
    t_out->Branch("nTruePrimParts", &nTruePrimParts);
    t_out->Branch("truePrimPartPDG", truePrimPartPDG,"truePrimPartPDG[200]/I");
    t_out->Branch("truePrimPartPx", truePrimPartPx,"truePrimPartPx[200]/F");
    t_out->Branch("truePrimPartPy", truePrimPartPy,"truePrimPartPy[200]/F");
    t_out->Branch("truePrimPartPz", truePrimPartPz,"truePrimPartPz[200]/F");
    t_out->Branch("truePrimPartE", truePrimPartE,"truePrimPartE[200]/F");
    t_out->Branch("trueVtxX", &trueVtxX);
    t_out->Branch("trueVtxY", &trueVtxY);
    t_out->Branch("trueVtxZ", &trueVtxZ);
  }

  t_out->Branch("foundVertex", &foundVertex);
  t_out->Branch("vtxX", &vtxX);
  t_out->Branch("vtxY", &vtxY);
  t_out->Branch("vtxZ", &vtxZ);

  t_out->Branch("nTracks", &nTracks);
  t_out->Branch("trackIsSecondary", trackIsSecondary,"trackIsSecondary[100]/I");
  t_out->Branch("trackStartPosX", trackStartPosX,"trackStartPosX[100]/F");
  t_out->Branch("trackStartPosY", trackStartPosY,"trackStartPosY[100]/F");
  t_out->Branch("trackStartPosZ", trackStartPosZ,"trackStartPosZ[100]/F");
  t_out->Branch("trackStartDirX", trackStartDirX,"trackStartDirX[100]/F");
  t_out->Branch("trackStartDirY", trackStartDirY,"trackStartDirY[100]/F");
  t_out->Branch("trackStartDirZ", trackStartDirZ,"trackStartDirZ[100]/F");
  t_out->Branch("trackEndPosX", trackEndPosX,"trackEndPosX[100]/F");
  t_out->Branch("trackEndPosY", trackEndPosY,"trackEndPosY[100]/F");
  t_out->Branch("trackEndPosZ", trackEndPosZ,"trackEndPosZ[100]/F");
  t_out->Branch("trackClassified", trackClassified,"trackClassified[100]/I");
  t_out->Branch("trackNPlanesAbove", trackNPlanesAbove,"trackNPlanesAbove[100]/I");
  t_out->Branch("trackPID", trackPID,"trackPID[100]/I");
  t_out->Branch("trackRecoE", trackRecoE,"trackRecoE[100]/F");
  if(!is_data) t_out->Branch("trackTruePID", trackTruePID,"trackTruePID[100]/I");

  t_out->Branch("nShowers", &nShowers);
  t_out->Branch("showerIsSecondary", showerIsSecondary,"showerIsSecondary[100]/I");
  t_out->Branch("showerStartPosX", showerStartPosX,"showerStartPosX[100]/F");
  t_out->Branch("showerStartPosY", showerStartPosY,"showerStartPosY[100]/F");
  t_out->Branch("showerStartPosZ", showerStartPosZ,"showerStartPosZ[100]/F");
  t_out->Branch("showerStartDirX", showerStartDirX,"showerStartDirX[100]/F");
  t_out->Branch("showerStartDirY", showerStartDirY,"showerStartDirY[100]/F");
  t_out->Branch("showerStartDirZ", showerStartDirZ,"showerStartDirZ[100]/F");
  t_out->Branch("showerClassified", showerClassified,"showerClassified[100]/I");
  t_out->Branch("showerNPlanesAbove", showerNPlanesAbove,"showerNPlanesAbove[100]/I");
  t_out->Branch("showerPID", showerPID,"showerPID[100]/I");
  t_out->Branch("showerRecoE", showerRecoE,"showerRecoE[100]/F");
  if(!is_data) t_out->Branch("showerTruePID", showerTruePID,"showerTruePID[100]/I");

  for(Long64_t ientry=0;ientry<pd_t_in->GetEntries();ientry++){
    if(ientry % 20000 == 0) std::cout << ientry << "/" << pd_t_in->GetEntries() <<  std::endl;
    pd_t_in->GetEntry(ientry);
    wc_t_in->GetEntry(ientry);
    lt_t_in->GetEntry(ientry);

    if(pd_rse[0] != wc_rse[0] || wc_rse[0] != lt_rse[0]) 
      throw std::invalid_argument("Can't use speedy merging - trees don't have events in same order");
    if(pd_rse[1] != wc_rse[1] || wc_rse[1] != lt_rse[1]) 
      throw std::invalid_argument("Can't use speedy merging - trees don't have events in same order");
    if(pd_rse[2] != wc_rse[2] || wc_rse[2] != lt_rse[2]) 
      throw std::invalid_argument("Can't use speedy merging - trees don't have events in same order");

    t_out->Fill();

  }



  t_out->Write("MergedNtuple");
  f_out->Close();
} 

