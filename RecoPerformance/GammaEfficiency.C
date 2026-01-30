#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

void GammaEfficiency(){

  //const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  //const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";

  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  const std::string file = "Merged_larpid_patch_smart_patch_test10_full_more.root";

  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(in_dir+file,f_in,t_in,false,false);

  bool check_containment = false;
  double FOM_t = 0.1; // Cut value in FOM calculation

  TH1D* h_true_leading_gamma_mom = new TH1D("h_true_leading_gamma_mom",";True Gamma Momentum (GeV);Events",40,0.0,0.6);
  TH1D* h_true_leading_gamma_costheta = new TH1D("h_true_leading_gamma_costheta",";True Gamma Cos(#theta);Events",40,-1.0,1.0);

  TH1D* h_true_subleading_gamma_mom = new TH1D("h_true_subleading_gamma_mom",";True Gamma Momentum (GeV);Events",40,0.0,0.6);
  TH1D* h_true_subleading_gamma_costheta = new TH1D("h_true_subleading_gamma_costheta",";True Gamma Cos(#theta);Events",40,-1.0,1.0);

  
  // Pandora histograms
  TH1D* h_selected_true_leading_gamma_mom_pd = new TH1D("h_selected_true_leading_gamma_mom_pd",";True Gamma Momentum (GeV);Events",40,0.0,0.6);
  TH1D* h_selected_true_leading_gamma_costheta_pd = new TH1D("h_selected_true_leading_gamma_costheta_pd",";True Gamma Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_pd = new TH2D("h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_pd",";True Gamma Momentum (Gev);Reco Gamma Momentum (Gev);",40,0.0,0.6,40,0.0,0.6);
  TH2D* h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_pd = new TH2D("h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_pd",";True Gamma Cos(#theta);Reco Gamma Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_leading_gamma_mom_error_pd = new TH1D("h_leading_gamma_mom_error_pd",";(Reco - True)/True;Events",100,-1,1);

  TH1D* h_selected_true_subleading_gamma_mom_pd = new TH1D("h_selected_true_subleading_gamma_mom_pd",";True Gamma Momentum (GeV);Events",40,0.0,0.6);
  TH1D* h_selected_true_subleading_gamma_costheta_pd = new TH1D("h_selected_true_subleading_gamma_costheta_pd",";True Gamma Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_pd = new TH2D("h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_pd",";True Gamma Momentum (Gev);Reco Gamma Momentum (Gev);",40,0.0,0.6,40,0.0,0.6);
  TH2D* h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_pd = new TH2D("h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_pd",";True Gamma Cos(#theta);Reco Gamma Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_subleading_gamma_mom_error_pd = new TH1D("h_subleading_gamma_mom_error_pd",";(Reco - True)/True;Events",100,-1,1);

  
  // Wirecell histograms
  TH1D* h_selected_true_leading_gamma_mom_wc = new TH1D("h_selected_true_leading_gamma_mom_wc",";True Gamma Momentum (GeV);Events",40,0.0,0.6);
  TH1D* h_selected_true_leading_gamma_costheta_wc = new TH1D("h_selected_true_leading_gamma_costheta_wc",";True Gamma Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_wc = new TH2D("h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_wc",";True Gamma Momentum (Gev);Reco Gamma Momentum (Gev);",40,0.0,0.6,40,0.0,0.6);
  TH2D* h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_wc = new TH2D("h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_wc",";True Gamma Cos(#theta);Reco Gamma Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_leading_gamma_mom_error_wc = new TH1D("h_leading_gamma_mom_error_wc",";(Reco - True)/True;Events",100,-1,1);

  TH1D* h_selected_true_subleading_gamma_mom_wc = new TH1D("h_selected_true_subleading_gamma_mom_wc",";True Gamma Momentum (GeV);Events",40,0.0,0.6);
  TH1D* h_selected_true_subleading_gamma_costheta_wc = new TH1D("h_selected_true_subleading_gamma_costheta_wc",";True Gamma Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_wc = new TH2D("h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_wc",";True Gamma Momentum (Gev);Reco Gamma Momentum (Gev);",40,0.0,0.6,40,0.0,0.6);
  TH2D* h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_wc = new TH2D("h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_wc",";True Gamma Cos(#theta);Reco Gamma Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_subleading_gamma_mom_error_wc = new TH1D("h_subleading_gamma_mom_error_wc",";(Reco - True)/True;Events",100,-1,1);

  
  // Lantern histograms
  TH1D* h_selected_true_leading_gamma_mom_lt = new TH1D("h_selected_true_leading_gamma_mom_lt",";True Gamma Momentum (GeV);Events",40,0.0,0.6);
  TH1D* h_selected_true_leading_gamma_costheta_lt = new TH1D("h_selected_true_leading_gamma_costheta_lt",";True Gamma Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_lt = new TH2D("h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_lt",";True Gamma Momentum (Gev);Reco Gamma Momentum (Gev);",40,0.0,0.6,40,0.0,0.6);
  TH2D* h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_lt = new TH2D("h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_lt",";True Gamma Cos(#theta);Reco Gamma Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_leading_gamma_mom_error_lt = new TH1D("h_leading_gamma_mom_error_lt",";(Reco - True)/True;Events",100,-1,1);

  TH1D* h_selected_true_subleading_gamma_mom_lt = new TH1D("h_selected_true_subleading_gamma_mom_lt",";True Gamma Momentum (GeV);Events",40,0.0,0.6);
  TH1D* h_selected_true_subleading_gamma_costheta_lt = new TH1D("h_selected_true_subleading_gamma_costheta_lt",";True Gamma Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_lt = new TH2D("h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_lt",";True Gamma Momentum (Gev);Reco Gamma Momentum (Gev);",40,0.0,0.6,40,0.0,0.6);
  TH2D* h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_lt = new TH2D("h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_lt",";True Gamma Cos(#theta);Reco Gamma Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_subleading_gamma_mom_error_lt = new TH1D("h_subleading_gamma_mom_error_lt",";(Reco - True)/True;Events",100,-1,1);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 100000) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);
    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 14) continue; 

    TVector3 plepton(mc_px->at(0),mc_py->at(0),mc_pz->at(0));
    if(plepton.Mag() < 0.1) continue;

    std::vector<int> pi0_ids;
    for(int i=0;i<truth_Ntrack;i++){
      if(truth_pdg[i] == 111 && truth_mother[i] == 0)
        pi0_ids.push_back(truth_id[i]);
    }

    int ngamma = 0;
    TVector3 p_leading_gamma(0,0,0);
    TVector3 p_subleading_gamma(0,0,0);
    for(int i=0;i<truth_Ntrack;i++){
      if(truth_pdg[i] == 22 && std::find(pi0_ids.begin(),pi0_ids.end(),truth_mother[i]) != pi0_ids.end()){
        TVector3 pgamma_tmp(truth_startMomentum[i][0],truth_startMomentum[i][1],truth_startMomentum[i][2]);
        ngamma++;
        if(pgamma_tmp.Mag() > p_leading_gamma.Mag()){
          p_subleading_gamma = p_leading_gamma;
          p_leading_gamma = pgamma_tmp;
        }
        else if(pgamma_tmp.Mag() > p_subleading_gamma.Mag()) p_subleading_gamma = pgamma_tmp;
      }
    }

    if(ngamma != 2) continue;

    double p_leading = p_leading_gamma.Mag();
    double costheta_leading = p_leading_gamma.CosTheta();
    h_true_leading_gamma_mom->Fill(p_leading);    
    h_true_leading_gamma_costheta->Fill(costheta_leading);    

    double p_subleading = p_subleading_gamma.Mag();
    double costheta_subleading = p_subleading_gamma.CosTheta();
    h_true_subleading_gamma_mom->Fill(p_subleading);    
    h_true_subleading_gamma_costheta->Fill(costheta_subleading);    

    // Pandora
    if(inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z)){

      std::vector<TLorentzVector> p4 = pd::RecoShower4MomV(shr_px_v,shr_py_v,shr_pz_v,shr_energy_y_v);
     
      if(p4.size() > 0){
        double p_leading_gamma_pd = p4.at(0).Vect().Mag();
        double costheta_leading_gamma_pd = p4.at(0).Vect().CosTheta();
        h_selected_true_leading_gamma_mom_pd->Fill(p_leading);
        h_selected_true_leading_gamma_costheta_pd->Fill(costheta_leading);    
        h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_pd->Fill(p_leading,p_leading_gamma_pd); 
        h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_pd->Fill(costheta_leading,costheta_leading_gamma_pd);
        h_leading_gamma_mom_error_pd->Fill((p_leading_gamma_pd-p_leading)/p_leading); 
      }
     
      if(p4.size() > 1){
        double p_subleading_gamma_pd = p4.at(1).Vect().Mag();
        double costheta_subleading_gamma_pd = p4.at(1).Vect().CosTheta();
        h_selected_true_subleading_gamma_mom_pd->Fill(p_subleading);
        h_selected_true_subleading_gamma_costheta_pd->Fill(costheta_subleading);    
        h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_pd->Fill(p_subleading,p_subleading_gamma_pd); 
        h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_pd->Fill(costheta_subleading,costheta_subleading_gamma_pd);
        h_subleading_gamma_mom_error_pd->Fill((p_subleading_gamma_pd-p_subleading)/p_subleading); 
      }

    }

    // Wirecell
    if(inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ)){
      std::vector<TLorentzVector> p4 = wc::RecoShower4MomV(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT,reco_id,-1);
     
      if(p4.size() > 0){
        double p_leading_gamma_wc = p4.at(0).Vect().Mag();
        double costheta_leading_gamma_wc = p4.at(0).Vect().CosTheta();
        h_selected_true_leading_gamma_mom_wc->Fill(p_leading);
        h_selected_true_leading_gamma_costheta_wc->Fill(costheta_leading);    
        h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_wc->Fill(p_leading,p_leading_gamma_wc); 
        h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_wc->Fill(costheta_leading,costheta_leading_gamma_wc);
        h_leading_gamma_mom_error_wc->Fill((p_leading_gamma_wc-p_leading)/p_leading); 
      }
     
      if(p4.size() > 1){
        double p_subleading_gamma_wc = p4.at(1).Vect().Mag();
        double costheta_subleading_gamma_wc = p4.at(1).Vect().CosTheta();
        h_selected_true_subleading_gamma_mom_wc->Fill(p_subleading);
        h_selected_true_subleading_gamma_costheta_wc->Fill(costheta_subleading);    
        h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_wc->Fill(p_subleading,p_subleading_gamma_wc); 
        h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_wc->Fill(costheta_subleading,costheta_subleading_gamma_wc);
        h_subleading_gamma_mom_error_wc->Fill((p_subleading_gamma_wc-p_subleading)/p_subleading); 
      }

    }

    // Lantern 
    if(foundVertex && inActiveTPC(vtxX,vtxY,vtxZ)){

      std::vector<TLorentzVector> p4 = lt::RecoShower4MomV(nShowers,showerIsSecondary,showerPID,showerRecoE,showerStartDirX,showerStartDirY,showerStartDirZ);
     
      if(p4.size() > 0){
        double p_leading_gamma_lt = p4.at(0).Vect().Mag();
        double costheta_leading_gamma_lt = p4.at(0).Vect().CosTheta();
        h_selected_true_leading_gamma_mom_lt->Fill(p_leading);
        h_selected_true_leading_gamma_costheta_lt->Fill(costheta_leading);    
        h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_lt->Fill(p_leading,p_leading_gamma_lt); 
        h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_lt->Fill(costheta_leading,costheta_leading_gamma_lt);
        h_leading_gamma_mom_error_lt->Fill((p_leading_gamma_lt-p_leading)/p_leading); 
      }
     
      if(p4.size() > 1){
        double p_subleading_gamma_lt = p4.at(1).Vect().Mag();
        double costheta_subleading_gamma_lt = p4.at(1).Vect().CosTheta();
        h_selected_true_subleading_gamma_mom_lt->Fill(p_subleading);
        h_selected_true_subleading_gamma_costheta_lt->Fill(costheta_subleading);    
        h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_lt->Fill(p_subleading,p_subleading_gamma_lt); 
        h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_lt->Fill(costheta_subleading,costheta_subleading_gamma_lt);
        h_subleading_gamma_mom_error_lt->Fill((p_subleading_gamma_lt-p_subleading)/p_subleading); 
      }

    }

  }  


  gSystem->Exec("mkdir -p Plots/GammaEfficiency/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.75,0.75,0.97,0.97);

  h_true_leading_gamma_mom->Draw("HIST");
  c->Print("Plots/GammaEfficiency/LeadingGammaMom.png");  
  c->Clear();

  h_true_leading_gamma_costheta->Draw("HIST");
  c->Print("Plots/GammaEfficiency/LeadingGammaCosTheta.png");  
  c->Clear();

  h_true_subleading_gamma_mom->Draw("HIST");
  c->Print("Plots/GammaEfficiency/SubleadingGammaMom.png");  
  c->Clear();

  h_true_subleading_gamma_costheta->Draw("HIST");
  c->Print("Plots/GammaEfficiency/SubleadingGammaCosTheta.png");  
  c->Clear();


  h_selected_true_leading_gamma_mom_pd->Divide(h_true_leading_gamma_mom);
  h_selected_true_leading_gamma_costheta_pd->Divide(h_true_leading_gamma_costheta);
  h_selected_true_subleading_gamma_mom_pd->Divide(h_true_subleading_gamma_mom);
  h_selected_true_subleading_gamma_costheta_pd->Divide(h_true_subleading_gamma_costheta);

  h_selected_true_leading_gamma_mom_wc->Divide(h_true_leading_gamma_mom);
  h_selected_true_leading_gamma_costheta_wc->Divide(h_true_leading_gamma_costheta);
  h_selected_true_subleading_gamma_mom_wc->Divide(h_true_subleading_gamma_mom);
  h_selected_true_subleading_gamma_costheta_wc->Divide(h_true_subleading_gamma_costheta);

  h_selected_true_leading_gamma_mom_lt->Divide(h_true_leading_gamma_mom);
  h_selected_true_leading_gamma_costheta_lt->Divide(h_true_leading_gamma_costheta);
  h_selected_true_subleading_gamma_mom_lt->Divide(h_true_subleading_gamma_mom);
  h_selected_true_subleading_gamma_costheta_lt->Divide(h_true_subleading_gamma_costheta);

  THStack* hs_leading_gamma_mom_eff = new THStack("hs_leading_gamma_mom_eff",";True Gamma Momentum (GeV);Efficiency");

  h_selected_true_leading_gamma_mom_pd->SetLineColor(1);
  h_selected_true_leading_gamma_mom_pd->SetLineWidth(2);
  hs_leading_gamma_mom_eff->Add(h_selected_true_leading_gamma_mom_pd);
  l->AddEntry(h_selected_true_leading_gamma_mom_pd,"PD","L");

  h_selected_true_leading_gamma_mom_wc->SetLineColor(2);
  h_selected_true_leading_gamma_mom_wc->SetLineWidth(2);
  hs_leading_gamma_mom_eff->Add(h_selected_true_leading_gamma_mom_wc);
  l->AddEntry(h_selected_true_leading_gamma_mom_wc,"WC","L");

  h_selected_true_leading_gamma_mom_lt->SetLineColor(3);
  h_selected_true_leading_gamma_mom_lt->SetLineWidth(2);
  hs_leading_gamma_mom_eff->Add(h_selected_true_leading_gamma_mom_lt);
  l->AddEntry(h_selected_true_leading_gamma_mom_lt,"LT","L");
    
  hs_leading_gamma_mom_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/GammaEfficiency/LeadingGammaMomentumEfficiency.png");
  c->Clear();
  l->Clear();

  THStack* hs_leading_gamma_costheta_eff = new THStack("hs_leading_gamma_costheta_eff",";True Gamma Cos(#theta);Efficiency");

  h_selected_true_leading_gamma_costheta_pd->SetLineColor(1);
  h_selected_true_leading_gamma_costheta_pd->SetLineWidth(2);
  hs_leading_gamma_costheta_eff->Add(h_selected_true_leading_gamma_costheta_pd);
  l->AddEntry(h_selected_true_leading_gamma_costheta_pd,"PD","L");

  h_selected_true_leading_gamma_costheta_wc->SetLineColor(2);
  h_selected_true_leading_gamma_costheta_wc->SetLineWidth(2);
  hs_leading_gamma_costheta_eff->Add(h_selected_true_leading_gamma_costheta_wc);
  l->AddEntry(h_selected_true_leading_gamma_costheta_wc,"WC","L");

  h_selected_true_leading_gamma_costheta_lt->SetLineColor(3);
  h_selected_true_leading_gamma_costheta_lt->SetLineWidth(2);
  hs_leading_gamma_costheta_eff->Add(h_selected_true_leading_gamma_costheta_lt);
  l->AddEntry(h_selected_true_leading_gamma_costheta_lt,"LT","L");

  hs_leading_gamma_costheta_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/GammaEfficiency/LeadingGammaCosThetaEfficiency.png");
  c->Clear();
  l->Clear();

  THStack* hs_subleading_gamma_mom_eff = new THStack("hs_subleading_gamma_mom_eff",";True Gamma Momentum (GeV);Efficiency");

  h_selected_true_subleading_gamma_mom_pd->SetLineColor(1);
  h_selected_true_subleading_gamma_mom_pd->SetLineWidth(2);
  hs_subleading_gamma_mom_eff->Add(h_selected_true_subleading_gamma_mom_pd);
  l->AddEntry(h_selected_true_subleading_gamma_mom_pd,"PD","L");

  h_selected_true_subleading_gamma_mom_wc->SetLineColor(2);
  h_selected_true_subleading_gamma_mom_wc->SetLineWidth(2);
  hs_subleading_gamma_mom_eff->Add(h_selected_true_subleading_gamma_mom_wc);
  l->AddEntry(h_selected_true_subleading_gamma_mom_wc,"WC","L");

  h_selected_true_subleading_gamma_mom_lt->SetLineColor(3);
  h_selected_true_subleading_gamma_mom_lt->SetLineWidth(2);
  hs_subleading_gamma_mom_eff->Add(h_selected_true_subleading_gamma_mom_lt);
  l->AddEntry(h_selected_true_subleading_gamma_mom_lt,"LT","L");
    
  hs_subleading_gamma_mom_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/GammaEfficiency/SubleadingGammaMomentumEfficiency.png");
  c->Clear();
  l->Clear();

  THStack* hs_subleading_gamma_costheta_eff = new THStack("hs_subleading_gamma_costheta_eff",";True Gamma Cos(#theta);Efficiency");

  h_selected_true_subleading_gamma_costheta_pd->SetLineColor(1);
  h_selected_true_subleading_gamma_costheta_pd->SetLineWidth(2);
  hs_subleading_gamma_costheta_eff->Add(h_selected_true_subleading_gamma_costheta_pd);
  l->AddEntry(h_selected_true_subleading_gamma_costheta_pd,"PD","L");

  h_selected_true_subleading_gamma_costheta_wc->SetLineColor(2);
  h_selected_true_subleading_gamma_costheta_wc->SetLineWidth(2);
  hs_subleading_gamma_costheta_eff->Add(h_selected_true_subleading_gamma_costheta_wc);
  l->AddEntry(h_selected_true_subleading_gamma_costheta_wc,"WC","L");

  h_selected_true_subleading_gamma_costheta_lt->SetLineColor(3);
  h_selected_true_subleading_gamma_costheta_lt->SetLineWidth(2);
  hs_subleading_gamma_costheta_eff->Add(h_selected_true_subleading_gamma_costheta_lt);
  l->AddEntry(h_selected_true_subleading_gamma_costheta_lt,"LT","L");

  hs_subleading_gamma_costheta_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/GammaEfficiency/SubleadingGammaCosThetaEfficiency.png");
  c->Clear();
  l->Clear();


  // Pandora 2D Plots
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_pd->Draw("colz");
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_pd->SetStats(0);
  c->Print("Plots/GammaEfficiency/LeadingGammaMomentumReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_pd);
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_pd->Draw("colz");
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_pd->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_LeadingGammaMomentumReconstruction_PD.png");
  c->Clear();

  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_pd->Draw("colz");
  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_pd->SetStats(0);
  c->Print("Plots/GammaEfficiency/LeadingGammaCosThetaReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_pd);
  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_pd->Draw("colz");
  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_pd->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_LeadingGammaCosThetaReconstruction_PD.png");
  c->Clear();


  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_pd->Draw("colz");
  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_pd->SetStats(0);
  c->Print("Plots/GammaEfficiency/SubleadingGammaMomentumReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_pd);
  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_pd->Draw("colz");
  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_pd->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_SubleadingGammaMomentumReconstruction_PD.png");
  c->Clear();

  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_pd->Draw("colz");
  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_pd->SetStats(0);
  c->Print("Plots/GammaEfficiency/SubleadingGammaCosThetaReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_pd);
  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_pd->Draw("colz");
  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_pd->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_SubleadingGammaCosThetaReconstruction_PD.png");
  c->Clear();


  // Wirecell 2D Plots
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_wc->Draw("colz");
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_wc->SetStats(0);
  c->Print("Plots/GammaEfficiency/LeadingGammaMomentumReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_wc);
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_wc->Draw("colz");
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_wc->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_LeadingGammaMomentumReconstruction_WC.png");
  c->Clear();

  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_wc->Draw("colz");
  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_wc->SetStats(0);
  c->Print("Plots/GammaEfficiency/LeadingGammaCosThetaReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_wc);
  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_wc->Draw("colz");
  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_wc->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_LeadingGammaCosThetaReconstruction_WC.png");
  c->Clear();


  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_wc->Draw("colz");
  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_wc->SetStats(0);
  c->Print("Plots/GammaEfficiency/SubleadingGammaMomentumReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_wc);
  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_wc->Draw("colz");
  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_wc->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_SubleadingGammaMomentumReconstruction_WC.png");
  c->Clear();

  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_wc->Draw("colz");
  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_wc->SetStats(0);
  c->Print("Plots/GammaEfficiency/SubleadingGammaCosThetaReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_wc);
  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_wc->Draw("colz");
  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_wc->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_SubleadingGammaCosThetaReconstruction_WC.png");
  c->Clear();


  // Lantern 2D Plots
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_lt->Draw("colz");
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_lt->SetStats(0);
  c->Print("Plots/GammaEfficiency/LeadingGammaMomentumReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_lt);
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_lt->Draw("colz");
  h_selected_true_leading_gamma_mom_reco_leading_gamma_mom_lt->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_LeadingGammaMomentumReconstruction_LT.png");
  c->Clear();

  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_lt->Draw("colz");
  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_lt->SetStats(0);
  c->Print("Plots/GammaEfficiency/LeadingGammaCosThetaReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_lt);
  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_lt->Draw("colz");
  h_selected_true_leading_gamma_costheta_reco_leading_gamma_costheta_lt->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_LeadingGammaCosThetaReconstruction_LT.png");
  c->Clear();


  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_lt->Draw("colz");
  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_lt->SetStats(0);
  c->Print("Plots/GammaEfficiency/SubleadingGammaMomentumReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_lt);
  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_lt->Draw("colz");
  h_selected_true_subleading_gamma_mom_reco_subleading_gamma_mom_lt->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_SubleadingGammaMomentumReconstruction_LT.png");
  c->Clear();

  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_lt->Draw("colz");
  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_lt->SetStats(0);
  c->Print("Plots/GammaEfficiency/SubleadingGammaCosThetaReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_lt);
  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_lt->Draw("colz");
  h_selected_true_subleading_gamma_costheta_reco_subleading_gamma_costheta_lt->SetStats(0);
  c->Print("Plots/GammaEfficiency/Normalised_SubleadingGammaCosThetaReconstruction_LT.png");
  c->Clear();

}
