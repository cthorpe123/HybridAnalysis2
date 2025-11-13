#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"


void MuonEfficiency(){

  //const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  //const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";

  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  const std::string file = "Merged_larpid_patch_smart_patch_test10_full_more.root";

  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(in_dir+file,f_in,t_in,false,false);

  bool check_containment = false;

  TH1D* h_true_muon_mom = new TH1D("h_true_muon_mom",";True Muon Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_true_muon_costheta = new TH1D("h_true_muon_costheta",";True Muon Cos(#theta);Events",40,-1.0,1.0);

  TH1D* h_selected_true_muon_mom_pd = new TH1D("h_selected_true_muon_mom_pd",";True Muon Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_muon_costheta_pd = new TH1D("h_selected_true_muon_costheta_pd",";True Muon Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_muon_mom_reco_muon_mom_pd = new TH2D("h_selected_true_muon_mom_reco_muon_mom_pd",";True Muon Momentum (Gev);Reco Muon Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_muon_costheta_reco_muon_costheta_pd = new TH2D("h_selected_true_muon_costheta_reco_muon_costheta_pd",";True Muon Cos(#theta);Reco Muon Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_muon_mom_error_pd = new TH1D("h_muon_mom_error_pd",";(Reco - True)/True;Events",100,-1,1);


  TH1D* h_selected_true_muon_mom_wc = new TH1D("h_selected_true_muon_mom_wc",";True Muon Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_muon_costheta_wc = new TH1D("h_selected_true_muon_costheta_wc",";True Muon Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_muon_mom_reco_muon_mom_wc = new TH2D("h_selected_true_muon_mom_reco_muon_mom_wc",";True Muon Momentum (Gev);Reco Muon Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_muon_costheta_reco_muon_costheta_wc = new TH2D("h_selected_true_muon_costheta_reco_muon_costheta_wc",";True Muon Cos(#theta);Reco Muon Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_muon_mom_error_wc = new TH1D("h_muon_mom_error_wc",";(Reco - True)/True;Events",100,-1,1);


  TH1D* h_selected_true_muon_mom_lt = new TH1D("h_selected_true_muon_mom_lt",";True Muon Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_muon_costheta_lt = new TH1D("h_selected_true_muon_costheta_lt",";True Muon Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_muon_mom_reco_muon_mom_lt = new TH2D("h_selected_true_muon_mom_reco_muon_mom_lt",";True Muon Momentum (Gev);Reco Muon Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_muon_costheta_reco_muon_costheta_lt = new TH2D("h_selected_true_muon_costheta_reco_muon_costheta_lt",";True Muon Cos(#theta);Reco Muon Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_muon_mom_error_lt = new TH1D("h_muon_mom_error_lt",";(Reco - True)/True;Events",100,-1,1);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 50000) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);
    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 14) continue; 

    TVector3 plepton(mc_px->at(0),mc_py->at(0),mc_pz->at(0));
    double p = plepton.Mag();
    double theta = plepton.Theta();
    if(p < 0.1) continue;

    h_true_muon_mom->Fill(p);    
    h_true_muon_costheta->Fill(cos(theta));    

    // Pandora
    int pd_muon = pd::SimpleMuonSelection(trk_llr_pid_score_v,trk_len_v);
    if(inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z) && pd_muon != -1 && abs(backtracked_pdg->at(pd_muon)) == 13 && !(check_containment && !isContained(trk_end_x_v->at(pd_muon),trk_end_y_v->at(pd_muon),trk_end_z_v->at(pd_muon)))){

      double reco_p_range = trk_range_muon_mom_v->at(pd_muon); 
      double reco_p_mcs = trk_mcs_muon_mom_v->at(pd_muon); 
      double reco_costheta = TVector3(trk_dir_x_v->at(pd_muon),trk_dir_y_v->at(pd_muon),trk_dir_z_v->at(pd_muon)).CosTheta();

      h_selected_true_muon_mom_pd->Fill(p);
      h_selected_true_muon_costheta_pd->Fill(cos(theta));    

      h_selected_true_muon_mom_reco_muon_mom_pd->Fill(p,reco_p_range); 
      h_selected_true_muon_costheta_reco_muon_costheta_pd->Fill(cos(theta),reco_costheta);
      h_muon_mom_error_pd->Fill((reco_p_range-p)/p); 

    }

    // Wirecell
    int wc_muon = wc::SimpleMuonSelection(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT);
    if(inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ) && wc_muon != -1 && abs(reco_truthMatch_pdg[wc_muon]) == 13 && !(check_containment && !isContained(reco_endXYZT[wc_muon][0],reco_endXYZT[wc_muon][1],reco_endXYZT[wc_muon][2]))){

      TVector3 reco_plepton(reco_startMomentum[wc_muon][0],reco_startMomentum[wc_muon][1],reco_startMomentum[wc_muon][2]);

      h_selected_true_muon_mom_wc->Fill(plepton.Mag());
      h_selected_true_muon_costheta_wc->Fill(plepton.CosTheta());    

      h_selected_true_muon_mom_reco_muon_mom_wc->Fill(plepton.Mag(),reco_plepton.Mag()); 
      h_selected_true_muon_costheta_reco_muon_costheta_wc->Fill(plepton.CosTheta(),reco_plepton.CosTheta());
      h_muon_mom_error_wc->Fill((reco_plepton.Mag()-plepton.Mag())/plepton.Mag()); 

    }

    // Lantern
    int lt_muon = lt::SimpleMuonSelection(nTracks,trackIsSecondary,trackPID);
    if(foundVertex && inActiveTPC(vtxX,vtxY,vtxZ) && lt_muon != -1 && abs(trackTruePID[lt_muon]) == 13 && !(check_containment && !isContained(trackEndPosX[lt_muon],trackEndPosY[lt_muon],trackEndPosZ[lt_muon]))){

      double mom = sqrt(trackRecoE[lt_muon]*trackRecoE[lt_muon]/1e6 + 2*0.106*trackRecoE[lt_muon]/1e3);
      TVector3 reco_muon_mom(mom*trackStartDirX[lt_muon],mom*trackStartDirY[lt_muon],mom*trackStartDirZ[lt_muon]);

      h_selected_true_muon_mom_lt->Fill(plepton.Mag());
      h_selected_true_muon_costheta_lt->Fill(plepton.CosTheta());

      h_selected_true_muon_mom_reco_muon_mom_lt->Fill(plepton.Mag(),reco_muon_mom.Mag());
      h_selected_true_muon_costheta_reco_muon_costheta_lt->Fill(plepton.CosTheta(),reco_muon_mom.CosTheta());
      h_muon_mom_error_lt->Fill((reco_muon_mom.Mag()-plepton.Mag())/plepton.Mag()); 

    }

  }  


  gSystem->Exec("mkdir -p Plots/MuonEfficiency/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.8,0.8,0.95,0.95);

  h_selected_true_muon_mom_pd->Divide(h_true_muon_mom);
  h_selected_true_muon_costheta_pd->Divide(h_true_muon_costheta);
  h_selected_true_muon_mom_wc->Divide(h_true_muon_mom);
  h_selected_true_muon_costheta_wc->Divide(h_true_muon_costheta);
  h_selected_true_muon_mom_lt->Divide(h_true_muon_mom);
  h_selected_true_muon_costheta_lt->Divide(h_true_muon_costheta);

  THStack* hs_muon_mom_eff = new THStack("hs_muon_mom_eff",";True Muon Momentum (GeV);Efficiency");

  h_selected_true_muon_mom_pd->SetLineColor(1);
  h_selected_true_muon_mom_pd->SetLineWidth(2);
  hs_muon_mom_eff->Add(h_selected_true_muon_mom_pd);
  l->AddEntry(h_selected_true_muon_mom_pd,"PD","L");
    
  h_selected_true_muon_mom_wc->SetLineColor(2);
  h_selected_true_muon_mom_wc->SetLineWidth(2);
  hs_muon_mom_eff->Add(h_selected_true_muon_mom_wc);
  l->AddEntry(h_selected_true_muon_mom_wc,"WC","L");

  h_selected_true_muon_mom_lt->SetLineColor(3);
  h_selected_true_muon_mom_lt->SetLineWidth(2);
  hs_muon_mom_eff->Add(h_selected_true_muon_mom_lt);
  l->AddEntry(h_selected_true_muon_mom_lt,"LT","L");

  hs_muon_mom_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/MuonEfficiency/MomentumEfficiency.png");
  c->Clear();
  l->Clear();

  THStack* hs_muon_costheta_eff = new THStack("hs_muon_costheta_eff",";True Muon Cos(#theta);Efficiency");

  h_selected_true_muon_costheta_pd->SetLineColor(1);
  h_selected_true_muon_costheta_pd->SetLineWidth(2);
  hs_muon_costheta_eff->Add(h_selected_true_muon_costheta_pd);
  l->AddEntry(h_selected_true_muon_costheta_pd,"PD","L");
    
  h_selected_true_muon_costheta_wc->SetLineColor(2);
  h_selected_true_muon_costheta_wc->SetLineWidth(2);
  hs_muon_costheta_eff->Add(h_selected_true_muon_costheta_wc);
  l->AddEntry(h_selected_true_muon_costheta_wc,"WC","L");

  h_selected_true_muon_costheta_lt->SetLineColor(3);
  h_selected_true_muon_costheta_lt->SetLineWidth(2);
  hs_muon_costheta_eff->Add(h_selected_true_muon_costheta_lt);
  l->AddEntry(h_selected_true_muon_costheta_lt,"LT","L");

  hs_muon_costheta_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/MuonEfficiency/CosThetaEfficiency.png");
  c->Clear();
  l->Clear();

  h_selected_true_muon_mom_reco_muon_mom_pd->Draw("colz");
  h_selected_true_muon_mom_reco_muon_mom_pd->SetStats(0);
  c->Print("Plots/MuonEfficiency/MomentumReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_muon_mom_reco_muon_mom_pd);
  h_selected_true_muon_mom_reco_muon_mom_pd->Draw("colz");
  h_selected_true_muon_mom_reco_muon_mom_pd->SetStats(0);
  c->Print("Plots/MuonEfficiency/Normalised_MomentumReconstruction_PD.png");
  c->Clear();

  h_selected_true_muon_costheta_reco_muon_costheta_pd->Draw("colz");
  h_selected_true_muon_costheta_reco_muon_costheta_pd->SetStats(0);
  c->Print("Plots/MuonEfficiency/CosThetaReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_muon_costheta_reco_muon_costheta_pd);
  h_selected_true_muon_costheta_reco_muon_costheta_pd->Draw("colz");
  h_selected_true_muon_costheta_reco_muon_costheta_pd->SetStats(0);
  c->Print("Plots/MuonEfficiency/Normalised_CosThetaReconstruction_PD.png");
  c->Clear();

  h_selected_true_muon_mom_reco_muon_mom_wc->Draw("colz");
  h_selected_true_muon_mom_reco_muon_mom_wc->SetStats(0);
  c->Print("Plots/MuonEfficiency/MomentumReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_muon_mom_reco_muon_mom_wc);
  h_selected_true_muon_mom_reco_muon_mom_wc->Draw("colz");
  h_selected_true_muon_mom_reco_muon_mom_wc->SetStats(0);
  c->Print("Plots/MuonEfficiency/Normalised_MomentumReconstruction_WC.png");
  c->Clear();

  h_selected_true_muon_costheta_reco_muon_costheta_wc->Draw("colz");
  h_selected_true_muon_costheta_reco_muon_costheta_wc->SetStats(0);
  c->Print("Plots/MuonEfficiency/CosThetaReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_muon_costheta_reco_muon_costheta_wc);
  h_selected_true_muon_costheta_reco_muon_costheta_wc->Draw("colz");
  h_selected_true_muon_costheta_reco_muon_costheta_wc->SetStats(0);
  c->Print("Plots/MuonEfficiency/Normalised_CosThetaReconstruction_WC.png");
  c->Clear();


  h_selected_true_muon_mom_reco_muon_mom_lt->Draw("colz");
  h_selected_true_muon_mom_reco_muon_mom_lt->SetStats(0);
  c->Print("Plots/MuonEfficiency/MomentumReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_muon_mom_reco_muon_mom_lt);
  h_selected_true_muon_mom_reco_muon_mom_lt->Draw("colz");
  h_selected_true_muon_mom_reco_muon_mom_lt->SetStats(0);
  c->Print("Plots/MuonEfficiency/Normalised_MomentumReconstruction_LT.png");
  c->Clear();

  h_selected_true_muon_costheta_reco_muon_costheta_lt->Draw("colz");
  h_selected_true_muon_costheta_reco_muon_costheta_lt->SetStats(0);
  c->Print("Plots/MuonEfficiency/CosThetaReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_muon_costheta_reco_muon_costheta_lt);
  h_selected_true_muon_costheta_reco_muon_costheta_lt->Draw("colz");
  h_selected_true_muon_costheta_reco_muon_costheta_lt->SetStats(0);
  c->Print("Plots/MuonEfficiency/Normalised_CosThetaReconstruction_LT.png");
  c->Clear();

  THStack* hs_muon_mom_err = new THStack("hs_muon_mom_err",";Muon Momentum (Reco - True)/True;Events");

  h_muon_mom_error_pd->SetLineColor(1);
  h_muon_mom_error_pd->SetLineWidth(2);
  hs_muon_mom_err->Add(h_muon_mom_error_pd);
  l->AddEntry(h_muon_mom_error_pd,"PD","L");

  h_muon_mom_error_wc->SetLineColor(2);
  h_muon_mom_error_wc->SetLineWidth(2);
  hs_muon_mom_err->Add(h_muon_mom_error_wc);
  l->AddEntry(h_muon_mom_error_wc,"WC","L");

  h_muon_mom_error_lt->SetLineColor(3);
  h_muon_mom_error_lt->SetLineWidth(2);
  hs_muon_mom_err->Add(h_muon_mom_error_lt);
  l->AddEntry(h_muon_mom_error_lt,"LT","L");
    
  hs_muon_mom_err->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/MuonEfficiency/MomentumError.png");
  c->Clear();
  l->Clear();


}
