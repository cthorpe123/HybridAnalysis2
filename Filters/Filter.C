#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "EnergyEstimatorFuncs.h"
#include "BranchList.h"
#include "Systematics.h"

void Filter(){

  is_data = false;
  is_ext = false;
  is_dirt = true;
  bool load_syst = false;

  //const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  //const std::string file = "Merged_larpid_patch_smart_patch_test10_full_more.root";
  //double POT_weight = 1.0; 
  //int detvar_univ = -1; // which detvar are you looking at

  // Main run4b files
  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";

  //const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_11_BNB_beam_on_surprise_reco2_hist.root";
  //double POT_weight = 1.0; 
  //int detvar_univ = -1; // which detvar are you looking at

  //const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  //double POT_weight = 1.332E+20/7.88166e+20; 
  //int detvar_univ = -1; // which detvar are you looking at

  //const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root";
  //double POT_weight = 31582916.0/88445969.0;
  //int detvar_univ = -1; // which detvar are you looking at

  const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root";
  double POT_weight = 1.332E+20/3.06E+20;
  int detvar_univ = -1; // which detvar are you looking at


  // Runs 4c -5, give same pot weight as run 4b for now so we can use to pad statistics
  //const std::string file = "Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4c.root";
  //double POT_weight = 1.332E+20/7.88166e+20; 
  //int detvar_univ = -1; // which detvar are you looking at
  
  //const std::string file = "Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4d.root";
  //double POT_weight = 1.332E+20/7.88166e+20; 
  //int detvar_univ = -1; // which detvar are you looking at

  //const std::string file = "Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_5.root";
  //double POT_weight = 1.332E+20/7.88166e+20; 
  //int detvar_univ = -1; // which detvar are you looking at

  // Run 4/5 detvars
/*
  //const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/detvar/";
  //const std::string file = "Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_cv_surprise_reco2_hist.root";
  //double POT_weight = 1.332E+20/5.42073e+20;
  //int detvar_univ = -1; // which detvar are you looking at
  
  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/detvar/";
  const std::string file = "Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_recomb2_surprise_reco2_hist.root";
  double POT_weight = 1.332E+20/1.14542e+20;
  int detvar_univ = syst::kRecomb2; // which detvar are you looking at
*/

  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(in_dir+file,f_in,t_in,(is_data || is_ext),load_syst);

  TFile* f_out = new TFile((in_dir+"Filtered_"+file).c_str(),"RECREATE");
  TTree* t_out = new TTree("DISNtuple","DISNtuple");

  if(is_data && is_ext || is_data && is_dirt || is_ext && is_dirt)
    throw std::invalid_argument("Sample can't be simultaneously flagged more than one of data/ext/dirt");

  // Truth branches
  t_out->Branch("run",&run);
  t_out->Branch("subrun",&subrun);
  t_out->Branch("event",&event);

  t_out->Branch("is_data",&is_data);
  t_out->Branch("is_ext",&is_ext);
  t_out->Branch("is_dirt",&is_dirt);
  t_out->Branch("category",&category);

  t_out->Branch("POT_weight",&POT_weight);
  t_out->Branch("detvar_univ",&detvar_univ);

  if(!is_data && !is_ext && !is_dirt){
    t_out->Branch("true_nu_vtx_x",&true_nu_vtx_x);
    t_out->Branch("true_nu_vtx_y",&true_nu_vtx_y);
    t_out->Branch("true_nu_vtx_z",&true_nu_vtx_z);
    t_out->Branch("nu_pdg",&nu_pdg);
    t_out->Branch("nu_e",&nu_e); 
    t_out->Branch("ccnc",&ccnc);
    t_out->Branch("interaction",&interaction);
    t_out->Branch("mc_pdg",&mc_pdg);
    t_out->Branch("mc_E",&mc_E);
    t_out->Branch("mc_px",&mc_px);
    t_out->Branch("mc_py",&mc_py);
    t_out->Branch("mc_pz",&mc_pz);
    t_out->Branch("mc_endx",&mc_endx);
    t_out->Branch("mc_endy",&mc_endy);
    t_out->Branch("mc_endz",&mc_endz);
  }

  // Analysis specific truth branches
  bool is_signal_t;
  bool in_tpc_t;
  TVector3 nu_vtx_t;
  bool has_muon_t;
  TVector3 muon_mom_t;
  bool muon_contained_t;
  double W_t;
  TLorentzVector proton_p4_t; 
  TLorentzVector pion_p4_t; 
  TLorentzVector pi0_p4_t; 
  TLorentzVector gamma_p4_t; 
  int nprot_t; 
  int npi_t; 
  int npi0_t; 
  int nsh_t; 
  std::vector<double> est_nu_e_t;

  if(!is_data && !is_ext && !is_dirt){
    t_out->Branch("is_signal_t",&is_signal_t);
    t_out->Branch("in_tpc_t",&in_tpc_t);
    t_out->Branch("nu_vtx_t",&nu_vtx_t);
    t_out->Branch("has_muon_t",&has_muon_t);
    t_out->Branch("muon_mom_t",&muon_mom_t);
    t_out->Branch("muon_contained_t",&muon_contained_t);
    t_out->Branch("W_t",&W_t);
    t_out->Branch("proton_p4_t",&proton_p4_t);
    t_out->Branch("pion_p4_t",&pion_p4_t);
    t_out->Branch("pi0_p4_t",&pi0_p4_t);
    t_out->Branch("gamma_p4_t",&gamma_p4_t);
    t_out->Branch("nprot_t",&nprot_t);
    t_out->Branch("npi_t",&npi_t);
    t_out->Branch("npi0_t",&npi0_t);
    t_out->Branch("nsh_t",&npi0_t);
    t_out->Branch("est_nu_e_t",&est_nu_e_t);
  }

  // PD Reco branches
  bool sel_pd;
  bool in_tpc_pd;
  TVector3 nu_vtx_pd;
  bool has_muon_pd;
  TVector3 muon_mom_pd;
  TVector3 muon_mom_mcs_pd;
  bool muon_contained_pd;
  double W_pd;
  int nprot_pd;
  int npi_pd; 
  int nsh_pd; 
  TLorentzVector proton_p4_pd; 
  TLorentzVector pion_p4_pd; 
  TLorentzVector gamma_p4_pd; 
  std::vector<double> est_nu_e_pd;

  t_out->Branch("sel_pd",&sel_pd);
  t_out->Branch("in_tpc_pd",&in_tpc_pd);
  t_out->Branch("nu_vtx_pd",&nu_vtx_pd);
  t_out->Branch("has_muon_pd",&has_muon_pd);
  t_out->Branch("muon_mom_pd",&muon_mom_pd);
  t_out->Branch("muon_mom_mcs_pd",&muon_mom_mcs_pd);
  t_out->Branch("muon_contained_pd",&muon_contained_pd);
  t_out->Branch("W_pd",&W_pd);
  t_out->Branch("nprot_pd",&nprot_pd);
  t_out->Branch("npi_pd",&npi_pd);
  t_out->Branch("nsh_pd",&nsh_pd);
  t_out->Branch("proton_p4_pd",&proton_p4_pd);
  t_out->Branch("pion_p4_pd",&pion_p4_pd);
  t_out->Branch("gamma_p4_pd",&gamma_p4_pd);
  t_out->Branch("est_nu_e_pd",&est_nu_e_pd);

  // WC Reco branches
  bool sel_wc;
  bool in_tpc_wc;
  TVector3 nu_vtx_wc;
  bool has_muon_wc;
  TVector3 muon_mom_wc;
  TVector3 muon_mom_len_wc;
  TVector3 muon_mom_mcs_wc;
  bool muon_contained_wc;
  double W_wc;
  int nprot_wc; 
  int npi_wc; 
  int nsh_wc; 
  TLorentzVector proton_p4_wc; 
  TLorentzVector pion_p4_wc; 
  TLorentzVector gamma_p4_wc; 
  std::vector<double> est_nu_e_wc;

  t_out->Branch("sel_wc",&sel_wc);
  t_out->Branch("in_tpc_wc",&in_tpc_wc);
  t_out->Branch("nu_vtx_wc",&nu_vtx_wc);
  t_out->Branch("has_muon_wc",&has_muon_wc);
  t_out->Branch("muon_mom_wc",&muon_mom_wc);
  t_out->Branch("muon_contained_wc",&muon_contained_wc);
  t_out->Branch("W_wc",&W_wc);
  t_out->Branch("nprot_wc",&nprot_wc);
  t_out->Branch("npi_wc",&npi_wc);
  t_out->Branch("nsh_wc",&nsh_wc);
  t_out->Branch("proton_p4_wc",&proton_p4_wc);
  t_out->Branch("pion_p4_wc",&pion_p4_wc);
  t_out->Branch("gamma_p4_wc",&gamma_p4_wc);
  t_out->Branch("est_nu_e_wc",&est_nu_e_wc);
  t_out->Branch("muon_mom_len_wc",&muon_mom_len_wc);
  t_out->Branch("muon_mom_mcs_wc",&muon_mom_mcs_wc);

  // LT Reco branches
  bool sel_lt;
  bool in_tpc_lt;
  TVector3 nu_vtx_lt;
  bool has_muon_lt;
  TVector3 muon_mom_lt;
  bool muon_contained_lt;
  double W_lt;
  int nprot_lt;
  int npi_lt; 
  int nsh_lt; 
  TLorentzVector proton_p4_lt; 
  TLorentzVector pion_p4_lt; 
  TLorentzVector gamma_p4_lt; 
  std::vector<double> est_nu_e_lt;

  t_out->Branch("sel_lt",&sel_lt);
  t_out->Branch("in_tpc_lt",&in_tpc_lt);
  t_out->Branch("nu_vtx_lt",&nu_vtx_lt);
  t_out->Branch("has_muon_lt",&has_muon_lt);
  t_out->Branch("muon_mom_lt",&muon_mom_lt);
  t_out->Branch("muon_contained_lt",&muon_contained_lt);
  t_out->Branch("W_lt",&W_lt);
  t_out->Branch("nprot_lt",&nprot_lt);
  t_out->Branch("npi_lt",&npi_lt);
  t_out->Branch("nsh_lt",&nsh_lt);
  t_out->Branch("proton_p4_lt",&proton_p4_lt);
  t_out->Branch("pion_p4_lt",&pion_p4_lt);
  t_out->Branch("gamma_p4_lt",&gamma_p4_lt);
  t_out->Branch("est_nu_e_lt",&est_nu_e_lt);

  // Hybrid methods 
  bool sel_h8;
  bool in_tpc_h8;
  TVector3 nu_vtx_h8;
  bool has_muon_h8;
  TVector3 muon_mom_h8;
  bool muon_contained_h8;
  double W_h8;
  int nprot_h8;
  int npi_h8; 
  int nsh_h8; 
  TLorentzVector proton_p4_h8; 
  TLorentzVector pion_p4_h8; 
  TLorentzVector gamma_p4_h8; 
  std::vector<double> est_nu_e_h8;

  t_out->Branch("sel_h8",&sel_h8);
  t_out->Branch("in_tpc_h8",&in_tpc_h8);
  t_out->Branch("nu_vtx_h8",&nu_vtx_h8);
  t_out->Branch("has_muon_h8",&has_muon_h8);
  t_out->Branch("muon_mom_h8",&muon_mom_h8);
  t_out->Branch("muon_contained_h8",&muon_contained_h8);
  t_out->Branch("W_h8",&W_h8);
  t_out->Branch("nprot_h8",&nprot_h8);
  t_out->Branch("npi_h8",&npi_h8);
  t_out->Branch("nsh_h8",&nsh_h8);
  t_out->Branch("proton_p4_h8",&proton_p4_h8);
  t_out->Branch("pion_p4_h8",&pion_p4_h8);
  t_out->Branch("gamma_p4_h8",&gamma_p4_h8);
  t_out->Branch("est_nu_e_h8",&est_nu_e_h8);

  // Systematics and tune 
  if(!is_data && !is_ext){
    t_out->Branch("weightSpline",&weightSpline); 
    t_out->Branch("weightTune",&weightTune); 
    t_out->Branch("weightSplineTimesTune",&weightSplineTimesTune); 
    if(load_syst){
      t_out->Branch("weightsGenie",&weightsGenie); 
      t_out->Branch("weightsReint",&weightsReint); 
      t_out->Branch("weightsFlux",&weightsFlux); 
      t_out->Branch("weightsUnisim",&weightsUnisim);
    }
  }

  for(int ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 1000) break;
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    category = -1;
    if(is_data) category = kData;
    else if(is_ext) category = kEXT;
    else if(is_dirt) category = kDirt;

    nprot_t = 0;
    npi_t = 0;
    npi0_t = 0;
    nsh_t = 0;
    has_muon_t = false;
    muon_mom_t = TVector3(0,0,0);
    proton_p4_t = TLorentzVector(0,0,0,0);
    pion_p4_t = TLorentzVector(0,0,0,0);
    pi0_p4_t = TLorentzVector(0,0,0,0);
    gamma_p4_t = TLorentzVector(0,0,0,0);
    est_nu_e_t = std::vector<double>(ee::kMAX,-1);

    if(!is_data && !is_ext && !is_dirt){

      // Select out signal events
      in_tpc_t = inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);
      nu_vtx_t = TVector3(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z);

      if(mc_pdg->size() && abs(mc_pdg->at(0)) == 13) muon_mom_t = TVector3(mc_px->at(0),mc_py->at(0),mc_pz->at(0));
      has_muon_t = mc_pdg->size() && abs(mc_pdg->at(0)) == 13 && muon_mom_t.Mag() > thresholds.at(13).first;
      muon_contained_t = has_muon_t && isContained(mc_endx->at(0),mc_endy->at(0),mc_endz->at(0));

      for(size_t i_p=0;i_p<mc_pdg->size();i_p++){
        TVector3 mom(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p));
        if(mc_pdg->at(i_p) == 2212 && mom.Mag() > thresholds.at(2212).first){
          nprot_t++;
          proton_p4_t += TLorentzVector(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p),mc_E->at(i_p));
        }
        if(abs(mc_pdg->at(i_p)) == 211 && mom.Mag() > thresholds.at(211).first){
          npi_t++;
          pion_p4_t += TLorentzVector(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p),mc_E->at(i_p));
        }
        if(mc_pdg->at(i_p) == 111){
          npi0_t++;
          pi0_p4_t += TLorentzVector(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p),mc_E->at(i_p));
        }
      } 

      std::vector<int> pi0_ids;
      for(int i=0;i<truth_Ntrack;i++){
        if(truth_pdg[i] == 111 && truth_mother[i] == 0)
          pi0_ids.push_back(truth_id[i]);
      }

      for(int i=0;i<truth_Ntrack;i++){
        if(truth_pdg[i] == 22 && std::find(pi0_ids.begin(),pi0_ids.end(),truth_mother[i]) != pi0_ids.end()){
          if(TVector3(truth_startMomentum[i][0],truth_startMomentum[i][1],truth_startMomentum[i][2]).Mag() > thresholds.at(22).first){
            gamma_p4_t += TLorentzVector(truth_startMomentum[i][0],truth_startMomentum[i][1],truth_startMomentum[i][2],truth_startMomentum[i][3]);               
            nsh_t++;
          }
        }
      }

      W_t = (proton_p4_t + pion_p4_t + gamma_p4_t).M();
      is_signal_t = abs(nu_pdg) == 14 && ccnc == 0 && in_tpc_t && muon_mom_t.Mag() > thresholds.at(13).first && nprot_t > 0; 

      // Energy estimators
      if(is_signal_t){
        TLorentzVector plepton_t(muon_mom_t.X(),muon_mom_t.Y(),muon_mom_t.Z(),sqrt(muon_mom_t.Mag()*muon_mom_t.Mag() + ml*ml));
        est_nu_e_t = ee::GetEnergyEst(plepton_t,W_t,proton_p4_t,nprot_t,pion_p4_t,npi_t,gamma_p4_t,npi0_t);
        category = kSignal;
      }
      else if(abs(nu_pdg) == 12 && in_tpc_t) category = kNue;
      else if(in_tpc_t) category = kBG;
      else category = kOutFV;
    }

    // Muon ID with each framework
    // PD Reco
    in_tpc_pd = inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z);
    nu_vtx_pd = TVector3(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z);
    int pd_muon = pd::SimpleMuonSelection(trk_llr_pid_score_v,trk_len_v);
    if(pd_muon != -1){
      muon_mom_pd = TVector3(trk_range_muon_mom_v->at(pd_muon)*trk_dir_x_v->at(pd_muon),trk_range_muon_mom_v->at(pd_muon)*trk_dir_y_v->at(pd_muon),trk_range_muon_mom_v->at(pd_muon)*trk_dir_z_v->at(pd_muon));
      muon_mom_mcs_pd = TVector3(trk_mcs_muon_mom_v->at(pd_muon)*trk_dir_x_v->at(pd_muon),trk_mcs_muon_mom_v->at(pd_muon)*trk_dir_y_v->at(pd_muon),trk_mcs_muon_mom_v->at(pd_muon)*trk_dir_z_v->at(pd_muon));
    }
    has_muon_pd = pd_muon != -1 && muon_mom_pd.Mag() > thresholds.at(13).first;
    muon_contained_pd = has_muon_pd && isContained(trk_end_x_v->at(pd_muon),trk_end_y_v->at(pd_muon),trk_end_z_v->at(pd_muon));

    std::vector<TLorentzVector> pd_proton_v = pd::RecoProton4MomV(trk_llr_pid_score_v,trk_len_v,trk_dir_x_v,trk_dir_y_v,trk_dir_z_v,pd_muon);
    nprot_pd = pd_proton_v.size();
    proton_p4_pd = SumTLorentzVector(pd_proton_v); 

    std::vector<TLorentzVector> pd_pion_v = pd::RecoPion4MomV(trk_llr_pid_score_v,trk_len_v,trk_dir_x_v,trk_dir_y_v,trk_dir_z_v,pd_muon);
    npi_pd = pd_pion_v.size();
    pion_p4_pd = SumTLorentzVector(pd_pion_v); 

    std::vector<TLorentzVector> pd_gamma_v = pd::RecoShower4MomV(shr_px_v,shr_py_v,shr_pz_v,shr_energy_y_v);
    nsh_pd = pd_gamma_v.size();
    gamma_p4_pd = SumTLorentzVector(pd_gamma_v);

    W_pd = (proton_p4_pd + pion_p4_pd + gamma_p4_pd).M();
    sel_pd = pd_muon != -1 && in_tpc_pd && nprot_pd > 0; 
    est_nu_e_pd = std::vector<double>(ee::kMAX,-1);
    TLorentzVector plepton_pd(muon_mom_pd.X(),muon_mom_pd.Y(),muon_mom_pd.Z(),sqrt(muon_mom_pd.Mag()*muon_mom_pd.Mag() + ml*ml));
    if(sel_pd) est_nu_e_pd = ee::GetEnergyEst(plepton_pd,W_pd,proton_p4_pd,nprot_pd,pion_p4_pd,npi_pd,gamma_p4_pd,nsh_pd);

    // WC Reco
    in_tpc_wc = inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ);
    nu_vtx_wc = TVector3(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ);
    int wc_muon = wc::SimpleMuonSelection(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT);
    if(wc_muon != -1){
      muon_mom_wc = TVector3(reco_startMomentum[wc_muon][0],reco_startMomentum[wc_muon][1],reco_startMomentum[wc_muon][2]);
      //TVector3 muon_dir = muon_mom_wc*(1.0/muon_mom_wc.Mag());
      muon_mom_mcs_wc = muon_mom_wc*(mcs_emu_MCS/muon_mom_wc.Mag());
      muon_mom_len_wc = muon_mom_wc*(mcs_emu_tracklen/muon_mom_wc.Mag()/1.04);
    } 
   has_muon_wc = wc_muon != -1 && muon_mom_wc.Mag() > thresholds.at(13).first;
    muon_contained_wc = has_muon_wc && isContained(reco_endXYZT[wc_muon][0],reco_endXYZT[wc_muon][1],reco_endXYZT[wc_muon][2]);
    

    std::vector<TLorentzVector> wc_proton_v = wc::RecoProton4MomV(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT,reco_id,wc_muon);
    nprot_wc = wc_proton_v.size();
    proton_p4_wc = SumTLorentzVector(wc_proton_v); 

    std::vector<TLorentzVector> wc_pion_v = wc::RecoPion4MomV(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT,reco_id,wc_muon);
    npi_wc = wc_pion_v.size();
    pion_p4_wc = SumTLorentzVector(wc_pion_v); 

    std::vector<TLorentzVector> wc_gamma_v = wc::RecoShower4MomV(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT,reco_id,wc_muon);
    nsh_wc = wc_gamma_v.size();
    gamma_p4_wc = SumTLorentzVector(wc_gamma_v);

    W_wc = (proton_p4_wc + pion_p4_wc + gamma_p4_wc).M();
    sel_wc = wc_muon != -1 && in_tpc_wc && nprot_wc > 0; 
    est_nu_e_wc = std::vector<double>(ee::kMAX,-1);
    TLorentzVector plepton_wc(muon_mom_wc.X(),muon_mom_wc.Y(),muon_mom_wc.Z(),sqrt(muon_mom_wc.Mag()*muon_mom_wc.Mag() + ml*ml));
    if(sel_wc) est_nu_e_wc = ee::GetEnergyEst(plepton_wc,W_wc,proton_p4_wc,nprot_wc,pion_p4_wc,npi_wc,gamma_p4_wc,nsh_wc);

    // LT Reco
    in_tpc_lt = inActiveTPC(vtxX,vtxY,vtxZ);
    nu_vtx_lt = TVector3(vtxX,vtxY,vtxZ);
    int lt_muon = lt::SimpleMuonSelection(nTracks,trackIsSecondary,trackPID);
    if(lt_muon != -1){
      double mom = sqrt(trackRecoE[lt_muon]*trackRecoE[lt_muon]/1e6 + 2*0.106*trackRecoE[lt_muon]/1e3);
      muon_mom_lt = TVector3(mom*trackStartDirX[lt_muon],mom*trackStartDirY[lt_muon],mom*trackStartDirZ[lt_muon]);
    }
    has_muon_lt = lt_muon != -1 && muon_mom_lt.Mag() > thresholds.at(13).first;
    muon_contained_lt = has_muon_lt && isContained(trackEndPosX[lt_muon],trackEndPosY[lt_muon],trackEndPosZ[lt_muon]);

    std::vector<TLorentzVector> lt_proton_v = lt::RecoProton4MomV(nTracks,trackIsSecondary,trackPID,trackRecoE,trackStartDirX,trackStartDirY,trackStartDirZ);
    nprot_lt = lt_proton_v.size();
    proton_p4_lt = SumTLorentzVector(lt_proton_v); 

    std::vector<TLorentzVector> lt_pion_v = lt::RecoPion4MomV(nTracks,trackIsSecondary,trackPID,trackRecoE,trackStartDirX,trackStartDirY,trackStartDirZ);
    npi_lt = lt_pion_v.size();
    pion_p4_lt = SumTLorentzVector(lt_pion_v); 

    std::vector<TLorentzVector> lt_gamma_v = lt::RecoShower4MomV(nShowers,showerIsSecondary,showerPID,showerRecoE,showerStartDirX,showerStartDirY,showerStartDirZ);
    nsh_lt = lt_gamma_v.size();
    gamma_p4_lt = SumTLorentzVector(lt_gamma_v);

    W_lt = (proton_p4_lt + pion_p4_lt + gamma_p4_lt).M();
    sel_lt = lt_muon != -1 && in_tpc_lt && nprot_lt > 0; 
    est_nu_e_lt = std::vector<double>(ee::kMAX,-1);
    TLorentzVector plepton_lt(muon_mom_lt.X(),muon_mom_lt.Y(),muon_mom_lt.Z(),sqrt(muon_mom_lt.Mag()*muon_mom_lt.Mag() + ml*ml));
    if(sel_lt) est_nu_e_lt = ee::GetEnergyEst(plepton_lt,W_lt,proton_p4_lt,nprot_lt,pion_p4_lt,npi_lt,gamma_p4_lt,nsh_lt);


    // H8 Reco - Lantern for the hadronic system, Pandora Range for contained muons, MCS for uncontained
    in_tpc_h8 = inActiveTPC(vtxX,vtxY,vtxZ);
    nu_vtx_h8 = TVector3(vtxX,vtxY,vtxZ);
    int h8_muon = pd_muon;
    if(h8_muon != -1)
      muon_mom_h8 = muon_contained_pd ? TVector3(muon_mom_pd.X(),muon_mom_pd.Y(),muon_mom_pd.Z()) : TVector3(muon_mom_mcs_pd.X(),muon_mom_mcs_pd.Y(),muon_mom_mcs_pd.Z());
    has_muon_h8 = pd_muon != -1 && muon_mom_pd.Mag() > thresholds.at(13).first;
    muon_contained_h8 = muon_contained_pd;

    nprot_h8 = lt_proton_v.size();
    proton_p4_h8 = SumTLorentzVector(lt_proton_v); 

    npi_h8 = lt_pion_v.size();
    pion_p4_h8 = SumTLorentzVector(lt_pion_v); 

    nsh_h8 = lt_gamma_v.size();
    gamma_p4_h8 = SumTLorentzVector(lt_gamma_v);

    W_h8 = (proton_p4_h8 + pion_p4_h8 + gamma_p4_h8).M();
    sel_h8 = has_muon_pd && sel_lt; 
    est_nu_e_h8 = std::vector<double>(ee::kMAX,-1);
    TLorentzVector plepton_h8(muon_mom_h8.X(),muon_mom_h8.Y(),muon_mom_h8.Z(),sqrt(muon_mom_h8.Mag()*muon_mom_h8.Mag() + ml*ml));
    if(sel_h8) est_nu_e_h8 = ee::GetEnergyEst(plepton_h8,W_h8,proton_p4_h8,nprot_h8,pion_p4_h8,npi_h8,gamma_p4_h8,nsh_h8);

    if(category == -1) std::cout << "Bad event" << std::endl;

    // Organise unisim systematics
    if(load_syst){
      weightsUnisim->clear();
      std::map<std::string,std::vector<double>>::iterator it;
      for(it = weights->begin();it != weights->end();it++){
        if(it->first != "flux_all" && it->first != "reint_all" && it->first != "All_UBGenie" && it->first != "xsr_scc_Fa3_SCC" && it->first != "xsr_scc_Fv3_SCC" && it->first != "ppfx_all"){
          weightsUnisim->insert(std::make_pair(it->first,it->second));
        } 
      }    
    }

    t_out->Fill();
  }

  t_out->Write();

  f_out->Close();
  f_in->Close();

}
