#ifndef _BranchList_h_
#define _BranchList_h_

int run,subrun,event;

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
Float_t weightSpline = 1.0;
Float_t weightSplineTimesTune = 1.0;
Float_t weightTune = 1.0;


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

std::vector<float>* shr_px_v=0;
std::vector<float>* shr_py_v=0;
std::vector<float>* shr_pz_v=0;
std::vector<float>* shr_energy_y_v=0;

float reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z;

// Wirecell Branches
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

// Lantern Branches
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


void LoadTree(std::string filename,TFile*& f_in,TTree*& t_in,bool is_data,bool load_syst){

  f_in = TFile::Open(filename.c_str());
  t_in = static_cast<TTree*>(f_in->Get("MergedNtuple"));

  weightSpline = 1.0;
  weightSplineTimesTune = 1.0;
  weightTune = 1.0;

  t_in->SetBranchAddress("run",&run);
  t_in->SetBranchAddress("sub",&subrun);
  t_in->SetBranchAddress("evt",&event);

  // Pandora branches
  if(!is_data){
    t_in->SetBranchAddress("mc_pdg",&mc_pdg);
    t_in->SetBranchAddress("mc_E",&mc_E);
    t_in->SetBranchAddress("mc_px",&mc_px);
    t_in->SetBranchAddress("mc_py",&mc_py);
    t_in->SetBranchAddress("mc_pz",&mc_pz);
    t_in->SetBranchAddress("mc_endx",&mc_endx);
    t_in->SetBranchAddress("mc_endy",&mc_endy);
    t_in->SetBranchAddress("mc_endz",&mc_endz);
    t_in->SetBranchAddress("true_nu_vtx_x",&true_nu_vtx_x);
    t_in->SetBranchAddress("true_nu_vtx_y",&true_nu_vtx_y);
    t_in->SetBranchAddress("true_nu_vtx_z",&true_nu_vtx_z);
    t_in->SetBranchAddress("nu_pdg",&nu_pdg);
    t_in->SetBranchAddress("interaction",&interaction);
    t_in->SetBranchAddress("ccnc",&ccnc);
    t_in->SetBranchAddress("nu_e",&nu_e); 
    if(load_syst){
      t_in->SetBranchAddress("weightsGenie",&weightsGenie); 
      t_in->SetBranchAddress("weightsReint",&weightsReint); 
      t_in->SetBranchAddress("weightsFlux",&weightsFlux); 
    }
    t_in->SetBranchAddress("weightSpline",&weightSpline); 
    t_in->SetBranchAddress("weightTune",&weightTune); 
    t_in->SetBranchAddress("weightSplineTimesTune",&weightSplineTimesTune); 
  }

  t_in->SetBranchAddress("trk_len_v",&trk_len_v);
  t_in->SetBranchAddress("trk_range_muon_mom_v",&trk_range_muon_mom_v);
  t_in->SetBranchAddress("trk_mcs_muon_mom_v",&trk_mcs_muon_mom_v);
  t_in->SetBranchAddress("trk_dir_x_v",&trk_dir_x_v);
  t_in->SetBranchAddress("trk_dir_y_v",&trk_dir_y_v);
  t_in->SetBranchAddress("trk_dir_z_v",&trk_dir_z_v);
  t_in->SetBranchAddress("trk_start_x_v",&trk_start_x_v);
  t_in->SetBranchAddress("trk_start_y_v",&trk_start_y_v);
  t_in->SetBranchAddress("trk_start_z_v",&trk_start_z_v);
  t_in->SetBranchAddress("trk_end_x_v",&trk_end_x_v);
  t_in->SetBranchAddress("trk_end_y_v",&trk_end_y_v);
  t_in->SetBranchAddress("trk_end_z_v",&trk_end_z_v);
  t_in->SetBranchAddress("trk_llr_pid_score_v",&trk_llr_pid_score_v);
  t_in->SetBranchAddress("reco_nu_vtx_x",&reco_nu_vtx_x);
  t_in->SetBranchAddress("reco_nu_vtx_y",&reco_nu_vtx_y);
  t_in->SetBranchAddress("reco_nu_vtx_z",&reco_nu_vtx_z);
  t_in->SetBranchAddress("shr_px_v",&shr_px_v);
  t_in->SetBranchAddress("shr_py_v",&shr_py_v);
  t_in->SetBranchAddress("shr_pz_v",&shr_pz_v);
  t_in->SetBranchAddress("shr_energy_y_v",&shr_energy_y_v);
  if(!is_data) t_in->SetBranchAddress("backtracked_pdg",&backtracked_pdg);
  
  // Wirecell branches
  if(!is_data){
    t_in->SetBranchAddress("mc_nu_pdg",&mc_nu_pdg); 
    t_in->SetBranchAddress("mc_nu_ccnc",&mc_nu_ccnc); 
    t_in->SetBranchAddress("mc_nu_mode",&mc_nu_mode); 
    t_in->SetBranchAddress("mc_nu_pos",mc_nu_pos); 
    t_in->SetBranchAddress("mc_nu_mom",mc_nu_mom); 
    t_in->SetBranchAddress("truth_Ntrack",&truth_Ntrack); 
    t_in->SetBranchAddress("truth_pdg",truth_pdg); 
    t_in->SetBranchAddress("truth_mother",truth_mother); 
    t_in->SetBranchAddress("truth_startMomentum",truth_startMomentum); 
    t_in->SetBranchAddress("truth_endMomentum",truth_endMomentum); 
    t_in->SetBranchAddress("truth_endXYZT",truth_endXYZT); 
    t_in->SetBranchAddress("truth_id",&truth_id);
  }

  t_in->SetBranchAddress("reco_nuvtxX",&reco_nuvtxX); 
  t_in->SetBranchAddress("reco_nuvtxY",&reco_nuvtxY); 
  t_in->SetBranchAddress("reco_nuvtxZ",&reco_nuvtxZ); 
  t_in->SetBranchAddress("reco_Ntrack",&reco_Ntrack); 
  t_in->SetBranchAddress("reco_pdg",reco_pdg); 
  t_in->SetBranchAddress("reco_mother",reco_mother); 
  t_in->SetBranchAddress("reco_startMomentum",reco_startMomentum); 
  t_in->SetBranchAddress("reco_endMomentum",reco_endMomentum); 
  t_in->SetBranchAddress("reco_endXYZT",reco_endXYZT); 
  t_in->SetBranchAddress("reco_id",reco_id);
  if(!is_data){
    t_in->SetBranchAddress("reco_truthMatch_pdg",reco_truthMatch_pdg);
    t_in->SetBranchAddress("reco_truthMatch_id",reco_truthMatch_id);
  }
  t_in->SetBranchAddress("mcs_mu_tracklen", &mcs_mu_tracklen);
  t_in->SetBranchAddress("mcs_emu_tracklen", &mcs_emu_tracklen);
  t_in->SetBranchAddress("mcs_emu_MCS", &mcs_emu_MCS);
  t_in->SetBranchAddress("mcs_ambiguity_MCS", &mcs_ambiguity_MCS);


  // Lantern branches
  if(!is_data){
    t_in->SetBranchAddress("trueNuE", &trueNuE);
    t_in->SetBranchAddress("trueNuPDG", &trueNuPDG);
    t_in->SetBranchAddress("trueNuCCNC", &trueNuCCNC);
    t_in->SetBranchAddress("nTruePrimParts", &nTruePrimParts);
    t_in->SetBranchAddress("truePrimPartPDG", truePrimPartPDG);
    t_in->SetBranchAddress("truePrimPartPx", truePrimPartPx);
    t_in->SetBranchAddress("truePrimPartPy", truePrimPartPy);
    t_in->SetBranchAddress("truePrimPartPz", truePrimPartPz);
    t_in->SetBranchAddress("truePrimPartE", truePrimPartE);
    t_in->SetBranchAddress("trueVtxX", &trueVtxX);
    t_in->SetBranchAddress("trueVtxY", &trueVtxY);
    t_in->SetBranchAddress("trueVtxZ", &trueVtxZ);
  }

  t_in->SetBranchAddress("foundVertex", &foundVertex);
  t_in->SetBranchAddress("vtxX", &vtxX);
  t_in->SetBranchAddress("vtxY", &vtxY);
  t_in->SetBranchAddress("vtxZ", &vtxZ);

  t_in->SetBranchAddress("nTracks", &nTracks);
  t_in->SetBranchAddress("trackIsSecondary", trackIsSecondary);
  t_in->SetBranchAddress("trackStartPosX", trackStartPosX);
  t_in->SetBranchAddress("trackStartPosY", trackStartPosY);
  t_in->SetBranchAddress("trackStartPosZ", trackStartPosZ);
  t_in->SetBranchAddress("trackStartDirX", trackStartDirX);
  t_in->SetBranchAddress("trackStartDirY", trackStartDirY);
  t_in->SetBranchAddress("trackStartDirZ", trackStartDirZ);
  t_in->SetBranchAddress("trackEndPosX", trackEndPosX);
  t_in->SetBranchAddress("trackEndPosY", trackEndPosY);
  t_in->SetBranchAddress("trackEndPosZ", trackEndPosZ);
  t_in->SetBranchAddress("trackClassified", trackClassified);
  t_in->SetBranchAddress("trackNPlanesAbove", trackNPlanesAbove);
  t_in->SetBranchAddress("trackPID", trackPID);
  t_in->SetBranchAddress("trackRecoE", trackRecoE);
  if(!is_data) t_in->SetBranchAddress("trackTruePID", trackTruePID);

  t_in->SetBranchAddress("nShowers", &nShowers);
  t_in->SetBranchAddress("showerIsSecondary", showerIsSecondary);
  t_in->SetBranchAddress("showerStartPosX", showerStartPosX);
  t_in->SetBranchAddress("showerStartPosY", showerStartPosY);
  t_in->SetBranchAddress("showerStartPosZ", showerStartPosZ);
  t_in->SetBranchAddress("showerStartDirX", showerStartDirX);
  t_in->SetBranchAddress("showerStartDirY", showerStartDirY);
  t_in->SetBranchAddress("showerStartDirZ", showerStartDirZ);
  t_in->SetBranchAddress("showerClassified", showerClassified);
  t_in->SetBranchAddress("showerNPlanesAbove", showerNPlanesAbove);
  t_in->SetBranchAddress("showerPID", showerPID);
  t_in->SetBranchAddress("showerRecoE", showerRecoE);
  if(!is_data) t_in->SetBranchAddress("showerTruePID", showerTruePID);

}

// Branches used in filtered ntuple
Bool_t          is_data;
Bool_t          is_ext;
Bool_t          is_dirt;
Int_t           category;
Bool_t          is_signal_t;
Bool_t          in_tpc_t;
Bool_t          has_muon_t;
TVector3        *muon_mom_t=0;
Bool_t          muon_contained_t;
Double_t        W_t;
Int_t           nprot_t;
Int_t           npi_t;
Int_t           npi0_t;
TLorentzVector  *proton_p4_t=0;
TLorentzVector  *pion_p4_t=0;
TLorentzVector  *shower_p4_t=0;
std::vector<double>* est_nu_e_t=0;

Bool_t          sel_pd;
Bool_t          in_tpc_pd;
Bool_t          has_muon_pd;
TVector3        *muon_mom_pd=0;
TVector3        *muon_mom_mcs_pd=0;
Bool_t          muon_contained_pd;
Double_t        W_pd;
Int_t           nprot_pd;
Int_t           npi_pd;
Int_t           nsh_pd;
TLorentzVector  *proton_p4_pd=0;
TLorentzVector  *pion_p4_pd=0;
TLorentzVector  *shower_p4_pd=0;
std::vector<double>* est_nu_e_pd=0;

Bool_t          sel_wc;
Bool_t          in_tpc_wc;
Bool_t          has_muon_wc;
TVector3        *muon_mom_wc=0;
TVector3        *muon_mom_len_wc=0;
TVector3        *muon_mom_mcs_wc=0;
Bool_t          muon_contained_wc;
Double_t        W_wc;
Int_t           nprot_wc;
Int_t           npi_wc;
Int_t           nsh_wc;
TLorentzVector  *proton_p4_wc=0;
TLorentzVector  *pion_p4_wc=0;
TLorentzVector  *shower_p4_wc=0;
std::vector<double>* est_nu_e_wc=0;

Bool_t          sel_lt;
Bool_t          in_tpc_lt;
Bool_t          has_muon_lt;
TVector3        *muon_mom_lt=0;
Bool_t          muon_contained_lt;
Double_t        W_lt;
Int_t           nprot_lt;
Int_t           npi_lt;
Int_t           nsh_lt;
TLorentzVector  *proton_p4_lt=0;
TLorentzVector  *pion_p4_lt=0;
TLorentzVector  *shower_p4_lt=0;
std::vector<double>* est_nu_e_lt=0;

Bool_t          sel_h8;
Bool_t          in_tpc_h8;
Bool_t          has_muon_h8;
TVector3        *muon_mom_h8=0;
Bool_t          muon_contained_h8;
Double_t        W_h8;
Int_t           nprot_h8;
Int_t           npi_h8;
Int_t           nsh_h8;
TLorentzVector  *proton_p4_h8=0;
TLorentzVector  *pion_p4_h8=0;
TLorentzVector  *shower_p4_h8=0;
std::vector<double>* est_nu_e_h8=0;

void LoadTreeFiltered(std::string filename,TFile*& f_in,TTree*& t_in,bool is_data,bool load_syst){

  f_in = TFile::Open(filename.c_str());
  t_in = static_cast<TTree*>(f_in->Get("DISNtuple"));

  weightSpline = 1.0;
  weightSplineTimesTune = 1.0;
  weightTune = 1.0;

  t_in->SetBranchAddress("run",&run);
  t_in->SetBranchAddress("subrun",&subrun);
  t_in->SetBranchAddress("event",&event);
  t_in->SetBranchAddress("is_data",&is_data);
  t_in->SetBranchAddress("is_ext",&is_ext);
  t_in->SetBranchAddress("is_dirt",&is_dirt);
  t_in->SetBranchAddress("category",&category);

  if(!is_data){
    t_in->SetBranchAddress("mc_pdg",&mc_pdg);
    t_in->SetBranchAddress("mc_E",&mc_E);
    t_in->SetBranchAddress("mc_px",&mc_px);
    t_in->SetBranchAddress("mc_py",&mc_py);
    t_in->SetBranchAddress("mc_pz",&mc_pz);
    t_in->SetBranchAddress("mc_endx",&mc_endx);
    t_in->SetBranchAddress("mc_endy",&mc_endy);
    t_in->SetBranchAddress("mc_endz",&mc_endz);
    t_in->SetBranchAddress("true_nu_vtx_x",&true_nu_vtx_x);
    t_in->SetBranchAddress("true_nu_vtx_y",&true_nu_vtx_y);
    t_in->SetBranchAddress("true_nu_vtx_z",&true_nu_vtx_z);
    t_in->SetBranchAddress("nu_pdg",&nu_pdg);
    t_in->SetBranchAddress("interaction",&interaction);
    t_in->SetBranchAddress("ccnc",&ccnc);
    t_in->SetBranchAddress("nu_e",&nu_e); 
    t_in->SetBranchAddress("is_signal_t", &is_signal_t);
    t_in->SetBranchAddress("in_tpc_t", &in_tpc_t);
    t_in->SetBranchAddress("has_muon_t", &has_muon_t);
    t_in->SetBranchAddress("muon_mom_t", &muon_mom_t);
    t_in->SetBranchAddress("muon_contained_t", &muon_contained_t);
    t_in->SetBranchAddress("W_t", &W_t);
    t_in->SetBranchAddress("nprot_t", &nprot_t);
    t_in->SetBranchAddress("npi_t", &npi_t);
    t_in->SetBranchAddress("npi0_t", &npi0_t);
    t_in->SetBranchAddress("proton_p4_t", &proton_p4_t);
    t_in->SetBranchAddress("pion_p4_t", &pion_p4_t);
    t_in->SetBranchAddress("shower_p4_t", &shower_p4_t);
    t_in->SetBranchAddress("est_nu_e_t",&est_nu_e_t);
    if(load_syst){
      t_in->SetBranchAddress("weightsGenie",&weightsGenie); 
      t_in->SetBranchAddress("weightsReint",&weightsReint); 
      t_in->SetBranchAddress("weightsFlux",&weightsFlux); 
    }

    t_in->SetBranchAddress("weightSpline",&weightSpline); 
    t_in->SetBranchAddress("weightTune",&weightTune); 
    t_in->SetBranchAddress("weightSplineTimesTune",&weightSplineTimesTune); 

  }

  t_in->SetBranchAddress("sel_pd", &sel_pd);
  t_in->SetBranchAddress("in_tpc_pd", &in_tpc_pd);
  t_in->SetBranchAddress("has_muon_pd", &has_muon_pd);
  t_in->SetBranchAddress("muon_mom_pd", &muon_mom_pd);
  t_in->SetBranchAddress("muon_mom_mcs_pd", &muon_mom_mcs_pd);
  t_in->SetBranchAddress("muon_contained_pd", &muon_contained_pd);
  t_in->SetBranchAddress("W_pd", &W_pd);
  t_in->SetBranchAddress("nprot_pd", &nprot_pd);
  t_in->SetBranchAddress("npi_pd", &npi_pd);
  t_in->SetBranchAddress("nsh_pd", &nsh_pd);
  t_in->SetBranchAddress("proton_p4_pd", &proton_p4_pd);
  t_in->SetBranchAddress("pion_p4_pd", &pion_p4_pd);
  t_in->SetBranchAddress("shower_p4_pd", &shower_p4_pd);
  t_in->SetBranchAddress("est_nu_e_pd",&est_nu_e_pd);

  t_in->SetBranchAddress("sel_wc", &sel_wc);
  t_in->SetBranchAddress("in_tpc_wc", &in_tpc_wc);
  t_in->SetBranchAddress("has_muon_wc", &has_muon_wc);
  t_in->SetBranchAddress("muon_mom_wc", &muon_mom_wc);
  t_in->SetBranchAddress("muon_mom_mcs_wc", &muon_mom_mcs_wc);
  t_in->SetBranchAddress("muon_mom_len_wc", &muon_mom_len_wc);
  t_in->SetBranchAddress("muon_contained_wc", &muon_contained_wc);
  t_in->SetBranchAddress("W_wc", &W_wc);
  t_in->SetBranchAddress("nprot_wc", &nprot_wc);
  t_in->SetBranchAddress("npi_wc", &npi_wc);
  t_in->SetBranchAddress("nsh_wc", &nsh_wc);
  t_in->SetBranchAddress("proton_p4_wc", &proton_p4_wc);
  t_in->SetBranchAddress("pion_p4_wc", &pion_p4_wc);
  t_in->SetBranchAddress("shower_p4_wc", &shower_p4_wc);
  t_in->SetBranchAddress("est_nu_e_wc",&est_nu_e_wc);

  t_in->SetBranchAddress("sel_lt", &sel_lt);
  t_in->SetBranchAddress("in_tpc_lt", &in_tpc_lt);
  t_in->SetBranchAddress("has_muon_lt", &has_muon_lt);
  t_in->SetBranchAddress("muon_mom_lt", &muon_mom_lt);
  t_in->SetBranchAddress("muon_contained_lt", &muon_contained_lt);
  t_in->SetBranchAddress("W_lt", &W_lt);
  t_in->SetBranchAddress("nprot_lt", &nprot_lt);
  t_in->SetBranchAddress("npi_lt", &npi_lt);
  t_in->SetBranchAddress("nsh_lt", &nsh_lt);
  t_in->SetBranchAddress("proton_p4_lt", &proton_p4_lt);
  t_in->SetBranchAddress("pion_p4_lt", &pion_p4_lt);
  t_in->SetBranchAddress("shower_p4_lt", &shower_p4_lt);
  t_in->SetBranchAddress("est_nu_e_lt",&est_nu_e_lt);

  t_in->SetBranchAddress("sel_h8", &sel_h8);
  t_in->SetBranchAddress("in_tpc_h8", &in_tpc_h8);
  t_in->SetBranchAddress("has_muon_h8", &has_muon_h8);
  t_in->SetBranchAddress("muon_mom_h8", &muon_mom_h8);
  t_in->SetBranchAddress("muon_contained_h8", &muon_contained_h8);
  t_in->SetBranchAddress("W_h8", &W_h8);
  t_in->SetBranchAddress("nprot_h8", &nprot_h8);
  t_in->SetBranchAddress("npi_h8", &npi_h8);
  t_in->SetBranchAddress("nsh_h8", &nsh_h8);
  t_in->SetBranchAddress("proton_p4_h8", &proton_p4_h8);
  t_in->SetBranchAddress("pion_p4_h8", &pion_p4_h8);
  t_in->SetBranchAddress("shower_p4_h8", &shower_p4_h8);
  t_in->SetBranchAddress("est_nu_e_h8",&est_nu_e_h8);

}

const std::vector<std::string> categories = {"Signal","BG","Nue","OutFV","Dirt","EXT","Data"};
enum e_cat {kSignal,kBG,kNue,kOutFV,kDirt,kEXT,kData};

#endif
