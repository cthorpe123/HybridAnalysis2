#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

// Check if the muon reco'd by pandora is also reconstructed as something else by one of the
// other two frameworks

void CheckingForOverlap(){

  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  const std::string file = "Merged_larpid_patch_smart_patch_test10_full_more.root";

  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(in_dir+file,f_in,t_in,false,false);

  double wc_reco_protons = 0.0;
  double wc_bad_reco_protons = 0.0;
  double wc_reco_pions = 0.0;
  double wc_bad_reco_pions = 0.0;
  TH1D* h_bad_proton_angle_wc = new TH1D("h_bad_proton_angle_wc",";Muon - Proton Opening Angle;Events",50,0.0,180.0);
  TH1D* h_bad_pion_angle_wc = new TH1D("h_bad_pion_angle_wc",";Muon - Proton Opening Angle;Events",50,0.0,180.0);
  TH1D* h_all_proton_angle_wc = new TH1D("h_all_proton_angle_wc",";Muon - Proton Opening Angle;Events",50,0.0,180.0);
  TH1D* h_all_pion_angle_wc = new TH1D("h_all_pion_angle_wc",";Muon - Proton Opening Angle;Events",50,0.0,180.0);
  TH2D* h_bad_proton_momentum_wc = new TH2D("h_bad_proton_momentum_wc",";Muon Momentum (GeV);Proton Momentum (GeV);",50,0.0,2.0,50,0.0,2.0);
  TH2D* h_bad_pion_momentum_wc = new TH2D("h_bad_pion_momentum_wc",";Muon Momentum (GeV);Pion Momentum (GeV);",50,0.0,2.0,50,0.0,2.0);
  TH2D* h_all_proton_momentum_wc = new TH2D("h_all_proton_momentum_wc",";Muon Momentum (GeV);Proton Momentum (GeV);",50,0.0,2.0,50,0.0,2.0);
  TH2D* h_all_pion_momentum_wc = new TH2D("h_all_pion_momentum_wc",";Muon Momentum (GeV);Pion Momentum (GeV);",50,0.0,2.0,50,0.0,2.0);
  

  double lt_reco_protons = 0.0;
  double lt_bad_reco_protons = 0.0;
  double lt_reco_pions = 0.0;
  double lt_bad_reco_pions = 0.0;
  TH1D* h_bad_proton_angle_lt = new TH1D("h_bad_proton_angle_lt",";Muon - Proton Opening Angle;Events",50,0.0,180.0);
  TH1D* h_bad_pion_angle_lt = new TH1D("h_bad_pion_angle_lt",";Muon - Proton Opening Angle;Events",50,0.0,180.0);
  TH1D* h_all_proton_angle_lt = new TH1D("h_all_proton_angle_lt",";Muon - Proton Opening Angle;Events",50,0.0,180.0);
  TH1D* h_all_pion_angle_lt = new TH1D("h_all_pion_angle_lt",";Muon - Proton Opening Angle;Events",50,0.0,180.0);
  TH2D* h_bad_proton_momentum_lt = new TH2D("h_bad_proton_momentum_lt",";Muon Momentum (GeV);Proton Momentum (GeV);",50,0.0,2.0,50,0.0,2.0);
  TH2D* h_bad_pion_momentum_lt = new TH2D("h_bad_pion_momentum_lt",";Muon Momentum (GeV);Pion Momentum (GeV);",50,0.0,2.0,50,0.0,2.0);
  TH2D* h_all_proton_momentum_lt = new TH2D("h_all_proton_momentum_lt",";Muon Momentum (GeV);Proton Momentum (GeV);",50,0.0,2.0,50,0.0,2.0);
  TH2D* h_all_pion_momentum_lt = new TH2D("h_all_pion_momentum_lt",";Muon Momentum (GeV);Pion Momentum (GeV);",50,0.0,2.0,50,0.0,2.0);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 100) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);
    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 14) continue; 

    // Check there is a muon
    TVector3 plepton(mc_px->at(0),mc_py->at(0),mc_pz->at(0));
    double p = plepton.Mag();
    double theta = plepton.Theta();
    if(p < 0.1) continue;

    // Check there is a protons
    double leading_proton_mom = 0.0;
    double leading_proton_costheta = 0.0;
    TVector3 pproton(0,0,0);
    for(size_t i_p=0;i_p<mc_pdg->size();i_p++){
      double mom = TVector3(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p)).Mag();
      if(abs(mc_pdg->at(i_p)) == 2212 && mom > leading_proton_mom){
        leading_proton_mom = mom;
      }
    } 

    // If leading proton has less than 0.3 GeV, event is not signal
    if(leading_proton_mom < 0.3) continue;

    int pd_muon = pd::SimpleMuonSelection(trk_llr_pid_score_v,trk_len_v);
    bool good_muon = inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z) && pd_muon != -1 && abs(backtracked_pdg->at(pd_muon)) == 13;
    if(!good_muon) continue;

    TVector3 pmuon(trk_range_muon_mom_v->at(pd_muon)*trk_dir_x_v->at(pd_muon),trk_range_muon_mom_v->at(pd_muon)*trk_dir_y_v->at(pd_muon),trk_range_muon_mom_v->at(pd_muon)*trk_dir_z_v->at(pd_muon));

    // Get the list of hadrons from WC 
    for(int i_tr=0;i_tr<reco_Ntrack;i_tr++){
      if(reco_mother[i_tr] != 0 || reco_pdg[i_tr] == 13) continue;
      TVector3 dir(reco_startMomentum[i_tr][0],reco_startMomentum[i_tr][1],reco_startMomentum[i_tr][2]);
      if(reco_pdg[i_tr] == 2212 && abs(reco_truthMatch_pdg[i_tr]) == 2212){
        wc_reco_protons++;
        h_all_proton_angle_wc->Fill((180/3.14)*pmuon.Angle(dir));
        h_all_proton_momentum_wc->Fill(pmuon.Mag(),dir.Mag());    
      }
      if(reco_pdg[i_tr] == 2212 && abs(reco_truthMatch_pdg[i_tr]) == 13){
        wc_bad_reco_protons++; 
        h_bad_proton_angle_wc->Fill((180/3.14)*pmuon.Angle(dir));
        h_bad_proton_momentum_wc->Fill(pmuon.Mag(),dir.Mag());    
      }
      if(abs(reco_pdg[i_tr]) == 211 && abs(reco_truthMatch_pdg[i_tr]) == 211){
        wc_reco_pions++;
        h_all_pion_angle_wc->Fill((180/3.14)*pmuon.Angle(dir));
        h_all_pion_momentum_wc->Fill(pmuon.Mag(),dir.Mag());    
      }
      if(abs(reco_pdg[i_tr]) == 211 && abs(reco_truthMatch_pdg[i_tr]) == 13){
        wc_bad_reco_pions++; 
        h_bad_pion_angle_wc->Fill((180/3.14)*pmuon.Angle(dir));
        h_bad_pion_momentum_wc->Fill(pmuon.Mag(),dir.Mag());    
      }
    }

    // Lantern
    
    for(int i=0;i<nTracks;i++){
      if(trackIsSecondary[i]) continue;
      double mass =  trackPID[i] == 2212 ? Mp : mpi;
      double mom = sqrt(trackRecoE[i]*trackRecoE[i]/1e6 + 2*mass*trackRecoE[i]/1e3);
      TVector3 dir(mom*trackStartDirX[i],mom*trackStartDirY[i],mom*trackStartDirZ[i]);
      if(abs(trackPID[i]) == 2212 && abs(trackTruePID[i]) == 2212){
        lt_reco_protons++;
        h_all_proton_angle_lt->Fill((180/3.14)*pmuon.Angle(dir));
        h_all_proton_momentum_lt->Fill(pmuon.Mag(),dir.Mag());  
      } 
      if(abs(trackPID[i]) == 2212 && abs(trackTruePID[i]) == 13){
        lt_bad_reco_protons++;
        h_bad_proton_angle_lt->Fill((180/3.14)*pmuon.Angle(dir));
        h_bad_proton_momentum_lt->Fill(pmuon.Mag(),dir.Mag());  
      }
      if(abs(trackPID[i]) == 211 && abs(trackTruePID[i]) == 211){
        lt_reco_pions++;
        h_all_pion_angle_lt->Fill((180/3.14)*pmuon.Angle(dir));
        h_all_pion_momentum_lt->Fill(pmuon.Mag(),dir.Mag());  
      }
      if(abs(trackPID[i]) == 211 && abs(trackTruePID[i]) == 13){
        lt_bad_reco_pions++;
        h_bad_pion_angle_lt->Fill((180/3.14)*pmuon.Angle(dir));
        h_bad_pion_momentum_lt->Fill(pmuon.Mag(),dir.Mag());  
      }
    }

  }  

  std::cout << "wc_bad_reco_proton_frac = " << wc_bad_reco_protons/(wc_reco_protons+wc_bad_reco_protons) << std::endl;
  std::cout << "wc_bad_reco_pion_frac = " << wc_bad_reco_pions/(wc_reco_pions+wc_bad_reco_pions) << std::endl;
  std::cout << "lt_bad_reco_proton_frac = " << lt_bad_reco_protons/(lt_reco_protons+lt_bad_reco_protons) << std::endl;
  std::cout << "lt_bad_reco_pion_frac = " << lt_bad_reco_pions/(lt_reco_pions+lt_bad_reco_pions) << std::endl;

  gSystem->Exec("mkdir -p Plots/CheckingForOverlap/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.8,0.8,0.95,0.95);
   
  h_bad_proton_angle_wc->Draw("HIST");
  h_bad_proton_angle_wc->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_bad_proton_angle_wc.png");
  c->Clear();
   
  h_bad_pion_angle_wc->Draw("HIST");
  h_bad_pion_angle_wc->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_bad_pion_angle_wc.png");
  c->Clear();
   
  h_bad_proton_angle_lt->Draw("HIST");
  h_bad_proton_angle_lt->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_bad_proton_angle_lt.png");
  c->Clear();
   
  h_bad_pion_angle_lt->Draw("HIST");
  h_bad_pion_angle_lt->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_bad_pion_angle_lt.png");
  c->Clear();

  h_bad_proton_momentum_wc->Draw("colz");
  h_bad_proton_momentum_wc->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_bad_proton_momentum_wc.png");
  c->Clear(); 

  h_bad_pion_momentum_wc->Draw("colz");
  h_bad_pion_momentum_wc->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_bad_pion_momentum_wc.png");
  c->Clear(); 

  h_bad_proton_momentum_lt->Draw("colz");
  h_bad_proton_momentum_lt->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_bad_proton_momentum_lt.png");
  c->Clear(); 

  h_bad_pion_momentum_lt->Draw("colz");
  h_bad_pion_momentum_lt->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_bad_pion_momentum_lt.png");
  c->Clear(); 
   
  h_all_proton_angle_wc->Draw("HIST");
  h_all_proton_angle_wc->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_all_proton_angle_wc.png");
  c->Clear();
   
  h_all_pion_angle_wc->Draw("HIST");
  h_all_pion_angle_wc->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_all_pion_angle_wc.png");
  c->Clear();
   
  h_all_proton_angle_lt->Draw("HIST");
  h_all_proton_angle_lt->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_all_proton_angle_lt.png");
  c->Clear();
   
  h_all_pion_angle_lt->Draw("HIST");
  h_all_pion_angle_lt->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_all_pion_angle_lt.png");
  c->Clear();

  h_all_proton_momentum_wc->Draw("colz");
  h_all_proton_momentum_wc->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_all_proton_momentum_wc.png");
  c->Clear(); 

  h_all_pion_momentum_wc->Draw("colz");
  h_all_pion_momentum_wc->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_all_pion_momentum_wc.png");
  c->Clear(); 

  h_all_proton_momentum_lt->Draw("colz");
  h_all_proton_momentum_lt->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_all_proton_momentum_lt.png");
  c->Clear(); 

  h_all_pion_momentum_lt->Draw("colz");
  h_all_pion_momentum_lt->SetStats(0);
  c->Print("Plots/CheckingForOverlap/h_all_pion_momentum_lt.png");
  c->Clear(); 

}
