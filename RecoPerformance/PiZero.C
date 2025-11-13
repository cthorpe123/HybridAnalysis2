#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"


void PiZero(){

 // const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
 // const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";

  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  const std::string file = "Merged_larpid_patch_smart_patch_test10_full_more.root";

  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(in_dir+file,f_in,t_in,false,false);

  TH1D* h_true_pizero_mom_pd = new TH1D("h_true_pizero_mom_pd",";True Pion Momentum (GeV);Events",30,0.0,1.0);
  TH1D* h_true_pizero_costheta_pd = new TH1D("h_true_pizero_costheta_pd",";True Pion Cos(#theta);Events",30,-1.0,1.0);
  TH1D* h_true_pizero_mom_wc = new TH1D("h_true_pizero_mom_wc",";True Pion Momentum (GeV);Events",30,0.0,1.0);
  TH1D* h_true_pizero_costheta_wc = new TH1D("h_true_pizero_costheta_wc",";True Pion Cos(#theta);Events",30,-1.0,1.0);
  TH1D* h_true_pizero_mom_lt = new TH1D("h_true_pizero_mom_lt",";True Pion Momentum (GeV);Events",30,0.0,1.0);
  TH1D* h_true_pizero_costheta_lt = new TH1D("h_true_pizero_costheta_lt",";True Pion Cos(#theta);Events",30,-1.0,1.0);

  TH1D* h_selected_true_pizero_mom_pd = new TH1D("h_selected_true_pizero_mom_pd",";True Muon Momentum (GeV);Events",30,0.0,1.0);
  TH1D* h_selected_true_pizero_costheta_pd = new TH1D("h_selected_true_pizero_costheta_pd",";True Muon Cos(#theta);Events",30,-1.0,1.0);
  TH2D* h_selected_true_pizero_mom_reco_pizero_mom_pd = new TH2D("h_selected_true_pizero_mom_reco_pizero_mom_pd",";True Muon Momentum (Gev);Reco Muon Momentum (Gev);",30,0.0,1.0,30,0.0,1.0);
  TH2D* h_selected_true_pizero_costheta_reco_pizero_costheta_pd = new TH2D("h_selected_true_pizero_costheta_reco_pizero_costheta_pd",";True Muon Cos(#theta);Reco Muon Cos(#theta);",30,-1.0,1.0,30,-1.0,1.0);
  TH1D* h_pizero_mom_error_pd = new TH1D("h_pizero_mom_error_pd",";(Reco - True)/True;Events",50,-1,1);
  TH1D* h_diphoton_mass_pd = new TH1D("h_diphoton_mass_pd",";W (GeV);Events",50,0.0,0.6);

  TH1D* h_selected_true_pizero_mom_wc = new TH1D("h_selected_true_pizero_mom_wc",";True Muon Momentum (GeV);Events",30,0.0,1.0);
  TH1D* h_selected_true_pizero_costheta_wc = new TH1D("h_selected_true_pizero_costheta_wc",";True Muon Cos(#theta);Events",30,-1.0,1.0);
  TH2D* h_selected_true_pizero_mom_reco_pizero_mom_wc = new TH2D("h_selected_true_pizero_mom_reco_pizero_mom_wc",";True Muon Momentum (Gev);Reco Muon Momentum (Gev);",30,0.0,1.0,30,0.0,1.0);
  TH2D* h_selected_true_pizero_costheta_reco_pizero_costheta_wc = new TH2D("h_selected_true_pizero_costheta_reco_pizero_costheta_wc",";True Muon Cos(#theta);Reco Muon Cos(#theta);",30,-1.0,1.0,30,-1.0,1.0);
  TH1D* h_pizero_mom_error_wc = new TH1D("h_pizero_mom_error_wc",";(Reco - True)/True;Events",50,-1,1);
  TH1D* h_diphoton_mass_wc = new TH1D("h_diphoton_mass_wc",";W (GeV);Events",50,0.0,0.6);


  TH1D* h_selected_true_pizero_mom_lt = new TH1D("h_selected_true_pizero_mom_lt",";True Muon Momentum (GeV);Events",30,0.0,1.0);
  TH1D* h_selected_true_pizero_costheta_lt = new TH1D("h_selected_true_pizero_costheta_lt",";True Muon Cos(#theta);Events",30,-1.0,1.0);
  TH2D* h_selected_true_pizero_mom_reco_pizero_mom_lt = new TH2D("h_selected_true_pizero_mom_reco_pizero_mom_lt",";True Muon Momentum (Gev);Reco Muon Momentum (Gev);",30,0.0,1.0,30,0.0,1.0);
  TH2D* h_selected_true_pizero_costheta_reco_pizero_costheta_lt = new TH2D("h_selected_true_pizero_costheta_reco_pizero_costheta_lt",";True Muon Cos(#theta);Reco Muon Cos(#theta);",30,-1.0,1.0,30,-1.0,1.0);
  TH1D* h_pizero_mom_error_lt = new TH1D("h_pizero_mom_error_lt",";(Reco - True)/True;Events",50,-1,1);
  TH1D* h_diphoton_mass_lt = new TH1D("h_diphoton_mass_lt",";W (GeV);Events",50,0.0,0.6);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 50000) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);

    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 0) continue; 

    // Calculating the true invariant mass
    double leading_pi0_p = 0.0;
    int n_pi0 = 0;
    TVector3 ppizero(0,0,0); 
    for(size_t i_p=0;i_p<mc_pdg->size();i_p++){
      if(mc_pdg->at(i_p) == 111){
        double  p = TVector3(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p)).Mag();
        if(p > leading_pi0_p){
          ppizero = TVector3(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p));
          leading_pi0_p = p;
          n_pi0++;
        }
      }
    }

    if(n_pi0 == 0) continue;

    if(inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z)){

      h_true_pizero_mom_pd->Fill(ppizero.Mag());    
      h_true_pizero_costheta_pd->Fill(ppizero.CosTheta());    

      // Pandora
      bool has_good_reco_pi0 = true;
      TLorentzVector p4_tot_reco(0,0,0,0);
      int showers = 0;
      for(size_t i_shr=0;i_shr<shr_px_v->size();i_shr++){
        double E = (shr_energy_y_v->at(i_shr))/1e3/0.83; 
        if(E < 0) continue;
        showers++;       
        p4_tot_reco += TLorentzVector(shr_px_v->at(i_shr)*E,shr_py_v->at(i_shr)*E,shr_pz_v->at(i_shr)*E,E);
      } 

      if(showers == 2){

        double reco_p = p4_tot_reco.Vect().Mag(); 
        double reco_costheta = p4_tot_reco.Vect().CosTheta();

        h_selected_true_pizero_mom_pd->Fill(ppizero.Mag());
        h_selected_true_pizero_costheta_pd->Fill(ppizero.CosTheta());    

        h_selected_true_pizero_mom_reco_pizero_mom_pd->Fill(ppizero.Mag(),reco_p); 
        h_selected_true_pizero_costheta_reco_pizero_costheta_pd->Fill(ppizero.CosTheta(),reco_costheta);
        h_pizero_mom_error_pd->Fill((reco_p-ppizero.Mag())/ppizero.Mag()); 
        h_diphoton_mass_pd->Fill(p4_tot_reco.M());

      }
    }


    // Wirecell
    if(inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ)){

      h_true_pizero_mom_wc->Fill(ppizero.Mag());    
      h_true_pizero_costheta_wc->Fill(ppizero.CosTheta());    

      int reco_pi0 = 0;
      bool has_reco_pi0 = false;
      double leading_reco_pi0_E = 0.0;
      int leading_reco_pi0_id = -1;
      std::map<int,int> part_by_id;
      for(size_t i_tr=0;i_tr<reco_Ntrack;i_tr++){
        part_by_id[reco_id[i_tr]] = reco_pdg[i_tr];
      }

      int showers = 0;
      TLorentzVector p4_tot_reco(0,0,0,0);
      for(size_t i_tr=0;i_tr<reco_Ntrack;i_tr++){
        if(reco_pdg[i_tr] == 22 && reco_mother[i_tr] != 0 && (part_by_id.at(reco_mother[i_tr]) != 2212 && part_by_id.at(reco_mother[i_tr]) != 211 && part_by_id.at(reco_mother[i_tr]) != 11 && part_by_id.at(reco_mother[i_tr]) != 13)){
          showers++;
          p4_tot_reco += TLorentzVector(reco_startMomentum[i_tr][0],reco_startMomentum[i_tr][1],reco_startMomentum[i_tr][2],reco_startMomentum[i_tr][3]);
        }
      }

      if(showers == 2){

        double reco_p = p4_tot_reco.Vect().Mag(); 
        double reco_costheta = p4_tot_reco.Vect().CosTheta();

        h_selected_true_pizero_mom_wc->Fill(ppizero.Mag());
        h_selected_true_pizero_costheta_wc->Fill(ppizero.CosTheta());    

        h_selected_true_pizero_mom_reco_pizero_mom_wc->Fill(ppizero.Mag(),reco_p); 
        h_selected_true_pizero_costheta_reco_pizero_costheta_wc->Fill(ppizero.CosTheta(),reco_costheta);
        h_pizero_mom_error_wc->Fill((reco_p-ppizero.Mag())/ppizero.Mag()); 
        h_diphoton_mass_wc->Fill(p4_tot_reco.M());

      }
    }





    // Lantern
    int lt_pizero = lt::SimpleMuonSelection(nTracks,trackIsSecondary,trackPID);
    if(foundVertex && inActiveTPC(vtxX,vtxY,vtxZ)){

      h_true_pizero_mom_lt->Fill(ppizero.Mag());    
      h_true_pizero_costheta_lt->Fill(ppizero.CosTheta());    

      int nphotons = 0;
      TLorentzVector p4_tot_reco(0,0,0,0); 
      for(int i=0;i<nShowers;i++){
        if(showerPID[i] == 22){
          p4_tot_reco += TLorentzVector(showerRecoE[i]*showerStartDirX[i]/1e3,showerRecoE[i]*showerStartDirY[i]/1e3,showerRecoE[i]*showerStartDirZ[i]/1e3,showerRecoE[i]/1e3);
          nphotons++;
        } 
      }    

      if(nphotons == 2){

        double reco_p = p4_tot_reco.Vect().Mag(); 
        double reco_costheta = p4_tot_reco.Vect().CosTheta();

        h_selected_true_pizero_mom_lt->Fill(ppizero.Mag());
        h_selected_true_pizero_costheta_lt->Fill(ppizero.CosTheta());    

        h_selected_true_pizero_mom_reco_pizero_mom_lt->Fill(ppizero.Mag(),reco_p); 
        h_selected_true_pizero_costheta_reco_pizero_costheta_lt->Fill(ppizero.CosTheta(),reco_costheta);
        h_pizero_mom_error_lt->Fill((reco_p-ppizero.Mag())/ppizero.Mag()); 
        h_diphoton_mass_lt->Fill(p4_tot_reco.M());

      }

    }

  }  


  gSystem->Exec("mkdir -p Plots/PiZero/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.8,0.8,0.95,0.95);

  h_selected_true_pizero_mom_pd->Divide(h_true_pizero_mom_pd);
  h_selected_true_pizero_costheta_pd->Divide(h_true_pizero_costheta_pd);
  h_selected_true_pizero_mom_wc->Divide(h_true_pizero_mom_wc);
  h_selected_true_pizero_costheta_wc->Divide(h_true_pizero_costheta_wc);
  h_selected_true_pizero_mom_lt->Divide(h_true_pizero_mom_lt);
  h_selected_true_pizero_costheta_lt->Divide(h_true_pizero_costheta_lt);

  THStack* hs_pizero_mom_eff = new THStack("hs_pizero_mom_eff",";True Muon Momentum (GeV);Efficiency");

  h_selected_true_pizero_mom_pd->SetLineColor(1);
  h_selected_true_pizero_mom_pd->SetLineWidth(2);
  hs_pizero_mom_eff->Add(h_selected_true_pizero_mom_pd);
  l->AddEntry(h_selected_true_pizero_mom_pd,"PD","L");

  h_selected_true_pizero_mom_wc->SetLineColor(2);
  h_selected_true_pizero_mom_wc->SetLineWidth(2);
  hs_pizero_mom_eff->Add(h_selected_true_pizero_mom_wc);
  l->AddEntry(h_selected_true_pizero_mom_wc,"WC","L");

  h_selected_true_pizero_mom_lt->SetLineColor(3);
  h_selected_true_pizero_mom_lt->SetLineWidth(2);
  hs_pizero_mom_eff->Add(h_selected_true_pizero_mom_lt);
  l->AddEntry(h_selected_true_pizero_mom_lt,"LT","L");

  hs_pizero_mom_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/PiZero/MomentumEfficiency.png");
  c->Clear();
  l->Clear();

  THStack* hs_pizero_costheta_eff = new THStack("hs_pizero_costheta_eff",";True Muon Cos(#theta);Efficiency");

  h_selected_true_pizero_costheta_pd->SetLineColor(1);
  h_selected_true_pizero_costheta_pd->SetLineWidth(2);
  hs_pizero_costheta_eff->Add(h_selected_true_pizero_costheta_pd);
  l->AddEntry(h_selected_true_pizero_costheta_pd,"PD","L");

  h_selected_true_pizero_costheta_wc->SetLineColor(2);
  h_selected_true_pizero_costheta_wc->SetLineWidth(2);
  hs_pizero_costheta_eff->Add(h_selected_true_pizero_costheta_wc);
  l->AddEntry(h_selected_true_pizero_costheta_wc,"WC","L");

  h_selected_true_pizero_costheta_lt->SetLineColor(3);
  h_selected_true_pizero_costheta_lt->SetLineWidth(2);
  hs_pizero_costheta_eff->Add(h_selected_true_pizero_costheta_lt);
  l->AddEntry(h_selected_true_pizero_costheta_lt,"LT","L");

  hs_pizero_costheta_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/PiZero/CosThetaEfficiency.png");
  c->Clear();
  l->Clear();

  h_selected_true_pizero_mom_reco_pizero_mom_pd->Draw("colz");
  h_selected_true_pizero_mom_reco_pizero_mom_pd->SetStats(0);
  c->Print("Plots/PiZero/MomentumReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_pizero_mom_reco_pizero_mom_pd);
  h_selected_true_pizero_mom_reco_pizero_mom_pd->Draw("colz");
  h_selected_true_pizero_mom_reco_pizero_mom_pd->SetStats(0);
  c->Print("Plots/PiZero/Normalised_MomentumReconstruction_PD.png");
  c->Clear();

  h_selected_true_pizero_costheta_reco_pizero_costheta_pd->Draw("colz");
  h_selected_true_pizero_costheta_reco_pizero_costheta_pd->SetStats(0);
  c->Print("Plots/PiZero/CosThetaReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_pizero_costheta_reco_pizero_costheta_pd);
  h_selected_true_pizero_costheta_reco_pizero_costheta_pd->Draw("colz");
  h_selected_true_pizero_costheta_reco_pizero_costheta_pd->SetStats(0);
  c->Print("Plots/PiZero/Normalised_CosThetaReconstruction_PD.png");
  c->Clear();

  h_selected_true_pizero_mom_reco_pizero_mom_wc->Draw("colz");
  h_selected_true_pizero_mom_reco_pizero_mom_wc->SetStats(0);
  c->Print("Plots/PiZero/MomentumReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_pizero_mom_reco_pizero_mom_wc);
  h_selected_true_pizero_mom_reco_pizero_mom_wc->Draw("colz");
  h_selected_true_pizero_mom_reco_pizero_mom_wc->SetStats(0);
  c->Print("Plots/PiZero/Normalised_MomentumReconstruction_WC.png");
  c->Clear();

  h_selected_true_pizero_costheta_reco_pizero_costheta_wc->Draw("colz");
  h_selected_true_pizero_costheta_reco_pizero_costheta_wc->SetStats(0);
  c->Print("Plots/PiZero/CosThetaReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_pizero_costheta_reco_pizero_costheta_wc);
  h_selected_true_pizero_costheta_reco_pizero_costheta_wc->Draw("colz");
  h_selected_true_pizero_costheta_reco_pizero_costheta_wc->SetStats(0);
  c->Print("Plots/PiZero/Normalised_CosThetaReconstruction_WC.png");
  c->Clear();


  h_selected_true_pizero_mom_reco_pizero_mom_lt->Draw("colz");
  h_selected_true_pizero_mom_reco_pizero_mom_lt->SetStats(0);
  c->Print("Plots/PiZero/MomentumReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_pizero_mom_reco_pizero_mom_lt);
  h_selected_true_pizero_mom_reco_pizero_mom_lt->Draw("colz");
  h_selected_true_pizero_mom_reco_pizero_mom_lt->SetStats(0);
  c->Print("Plots/PiZero/Normalised_MomentumReconstruction_LT.png");
  c->Clear();

  h_selected_true_pizero_costheta_reco_pizero_costheta_lt->Draw("colz");
  h_selected_true_pizero_costheta_reco_pizero_costheta_lt->SetStats(0);
  c->Print("Plots/PiZero/CosThetaReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_pizero_costheta_reco_pizero_costheta_lt);
  h_selected_true_pizero_costheta_reco_pizero_costheta_lt->Draw("colz");
  h_selected_true_pizero_costheta_reco_pizero_costheta_lt->SetStats(0);
  c->Print("Plots/PiZero/Normalised_CosThetaReconstruction_LT.png");
  c->Clear();

  THStack* hs_pizero_mom_err = new THStack("hs_pizero_mom_err",";#pi^{0} Momentum (Reco - True)/True;Events");

  h_pizero_mom_error_pd->SetLineColor(1);
  h_pizero_mom_error_pd->SetLineWidth(2);
  hs_pizero_mom_err->Add(h_pizero_mom_error_pd);
  l->AddEntry(h_pizero_mom_error_pd,"PD","L");

  h_pizero_mom_error_wc->SetLineColor(2);
  h_pizero_mom_error_wc->SetLineWidth(2);
  hs_pizero_mom_err->Add(h_pizero_mom_error_wc);
  l->AddEntry(h_pizero_mom_error_wc,"WC","L");

  h_pizero_mom_error_lt->SetLineColor(3);
  h_pizero_mom_error_lt->SetLineWidth(2);
  hs_pizero_mom_err->Add(h_pizero_mom_error_lt);
  l->AddEntry(h_pizero_mom_error_lt,"LT","L");

  hs_pizero_mom_err->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/PiZero/MomentumError.png");
  c->Clear();
  l->Clear();

  THStack* hs_diphoton_mass = new THStack("hs_diphoton_mass",";Diphoton Mass (GeV);Events");

  h_diphoton_mass_pd->SetLineColor(1);
  h_diphoton_mass_pd->SetLineWidth(2);
  hs_diphoton_mass->Add(h_diphoton_mass_pd);
  l->AddEntry(h_diphoton_mass_pd,"PD","L");

  h_diphoton_mass_wc->SetLineColor(2);
  h_diphoton_mass_wc->SetLineWidth(2);
  hs_diphoton_mass->Add(h_diphoton_mass_wc);
  l->AddEntry(h_diphoton_mass_wc,"WC","L");

  h_diphoton_mass_lt->SetLineColor(3);
  h_diphoton_mass_lt->SetLineWidth(2);
  hs_diphoton_mass->Add(h_diphoton_mass_lt);
  l->AddEntry(h_diphoton_mass_lt,"LT","L");

  hs_diphoton_mass->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/PiZero/DiphotonMass.png");
  c->Clear();
  l->Clear();


}
