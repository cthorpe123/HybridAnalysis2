#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

void DoubleProtons(){

  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";

  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(in_dir+file,f_in,t_in,false,false);

  TH2D* h_truth_2p = new TH2D("h_truth_2p",";Cos(#theta) Between Protons;2|p_{1}-p_{2}|/(p_{1}+p_{2});",100,0,3.14,100,0,1);
  TH2D* h_reco_2p = new TH2D("h_reco_2p",";Cos(#theta) Between Protons;2|p_{1}-p_{2}|/(p_{1}+p_{2});",100,0,3.14,100,0,1);
  TH2D* h_reco_1p = new TH2D("h_reco_1p",";Cos(#theta) Between Protons;2|p_{1}-p_{2}|/(p_{1}+p_{2});",100,0,3.14,100,0,1);

  TH1D* h_truth_2p_1D = new TH1D("h_truth_2p_1D",";2|p_{1}-p_{2}|/(p_{1}+p_{2});",100,0,2);
  TH1D* h_reco_2p_1D = new TH1D("h_reco_2p_1D",";2|p_{1}-p_{2}|/(p_{1}+p_{2});",100,0,2);
  TH1D* h_reco_1p_1D = new TH1D("h_reco_1p_1D",";2|p_{1}-p_{2}|/(p_{1}+p_{2});",100,0,2);
  
  TH1D* h_reco_1p_error_first = new TH1D("h_reco_1p_error_first","",100,-1,1);
  TH1D* h_reco_1p_error_second = new TH1D("h_reco_1p_error_second","",100,-1,1);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 100000) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);

    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 14) continue; 

    int nprot_t = 0;
    std::vector<TLorentzVector> p4_prot_t;
    for(size_t i_p=0;i_p<mc_pdg->size();i_p++){
      TVector3 mom(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p));
      if(abs(mc_pdg->at(i_p)) == 2212 && mom.Mag() > thresholds.at(2212).first){
        nprot_t++;
        p4_prot_t.push_back(TLorentzVector(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p),mc_E->at(i_p)));
      }
    } 

    // If leading proton has less than 0.3 GeV, event is not signal
    if(nprot_t == 0) continue;

    SortTLorentzVector(p4_prot_t);

    if(nprot_t == 2){
       double c = p4_prot_t.at(0).Vect().Angle(p4_prot_t.at(1).Vect());
       double frac = 2*(p4_prot_t.at(0).Vect().Mag() - p4_prot_t.at(1).Vect().Mag())/(p4_prot_t.at(0).Vect().Mag() + p4_prot_t.at(1).Vect().Mag());
       double d = (p4_prot_t.at(0).Vect() - p4_prot_t.at(1).Vect()).Mag()/(p4_prot_t.at(0).Vect() + p4_prot_t.at(1).Vect()).Mag();
       h_truth_2p->Fill(c,frac);
       h_truth_2p_1D->Fill(d);
    }

    // LT Reco

    std::vector<TLorentzVector> lt_proton_v = lt::RecoProton4MomV(nTracks,trackIsSecondary,trackPID,trackRecoE,trackStartDirX,trackStartDirY,trackStartDirZ);

    if(lt_proton_v.size() != 2) continue;
  
    double c = lt_proton_v.at(0).Vect().Angle(lt_proton_v.at(1).Vect());
    double frac = 2*(lt_proton_v.at(0).Vect().Mag() - lt_proton_v.at(1).Vect().Mag())/(lt_proton_v.at(0).Vect().Mag() + lt_proton_v.at(1).Vect().Mag());
    double d = (lt_proton_v.at(0).Vect() - lt_proton_v.at(1).Vect()).Mag()/(lt_proton_v.at(0).Vect() + lt_proton_v.at(1).Vect()).Mag();
    std::cout << d << std::endl;
  
    if(nprot_t == 2)  h_reco_2p->Fill(c,frac);
    if(nprot_t == 1)  h_reco_1p->Fill(c,frac);
  
    if(nprot_t == 2)  h_reco_2p_1D->Fill(d);
    if(nprot_t == 1)  h_reco_1p_1D->Fill(d);

    if(nprot_t == 1){
      double p_t = p4_prot_t.at(0).Vect().Mag();
      double p_r_first = lt_proton_v.at(0).Vect().Mag();
      double p_r_second = lt_proton_v.at(1).Vect().Mag();
      h_reco_1p_error_first->Fill((p_r_first - p_t)/p_t); 
      h_reco_1p_error_second->Fill((p_r_second - p_t)/p_t); 
    }

  }  

  gSystem->Exec("mkdir -p Plots/DoubleProtons/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.7,0.7,0.95,0.95);

  h_truth_2p->Draw("colz");
  c->Print("Plots/DoubleProtons/truth_2p.png");
  c->Clear();

  h_reco_2p->Draw("colz");
  c->Print("Plots/DoubleProtons/reco_2p.png");
  c->Clear();

  h_reco_1p->Draw("colz");
  c->Print("Plots/DoubleProtons/reco_1p.png");
  c->Clear();
  
  h_reco_1p_error_first->Draw("HIST");
  h_reco_1p_error_first->SetLineColor(1);
  h_reco_1p_error_first->SetLineWidth(2);
  h_reco_1p_error_first->SetStats(0);
  l->AddEntry(h_reco_1p_error_first,"Leading Reco Proton","L");
  
  h_reco_1p_error_second->Draw("HIST same");
  h_reco_1p_error_second->SetLineColor(2);
  h_reco_1p_error_second->SetLineWidth(2);
  h_reco_1p_error_second->SetStats(0);
  l->AddEntry(h_reco_1p_error_second,"Subleading Reco Proton","L");

  l->Draw();
 
  c->Print("Plots/DoubleProtons/reco_error.png"); 
  c->Clear(); 
  l->Clear();

  h_truth_2p_1D->Draw("HIST");
  h_truth_2p_1D->SetStats(0);
  h_truth_2p_1D->SetLineColor(1);
  h_truth_2p_1D->SetLineWidth(2);
  l->AddEntry(h_truth_2p_1D,"Truth 2p","L");

  h_reco_2p_1D->Draw("HIST same");
  h_reco_2p_1D->SetStats(0);
  h_reco_2p_1D->SetLineColor(2);
  h_reco_2p_1D->SetLineWidth(2);
  l->AddEntry(h_reco_2p_1D,"Reco 2p","L");

  h_reco_1p_1D->Draw("HIST same");
  h_reco_1p_1D->SetStats(0);
  h_reco_1p_1D->SetLineColor(3);
  h_reco_1p_1D->SetLineWidth(2);
  l->AddEntry(h_reco_1p_1D,"Reco 1p","L");

  c->Print("Plots/DoubleProtons/p_corr.png");

}
