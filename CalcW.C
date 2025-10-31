#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

// Try using different combinations of cuts and frameworks to calculate W

void CalcW(){

  const double data_POT = 1.5e21;

  const std::string file = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(file,f_in,t_in,false);
  const double POT = 7.88166e+20;

  double signal = 0.0;

  TH1D* h_TrueW = new TH1D("h_TrueW",";True W (GeV);Events",50,0.9,5.0);

  // Try lots of combinations of reconstructions

  // All Pandora
  TH1D* h_RecoW_pd = new TH1D("h_RecoW_pd",";Reco W (GeV);Events",50,0.9,5.0);
  TH2D* h_TrueW_RecoW_pd = new TH2D("h_TrueW_RecoW_pd",";True W (GeV);Reco W (GeV);Events",50,0.9,5.0,50,0.9,5.0);
  TH1D* h_ErrorW_pd = new TH1D("h_ErrorW_pd",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ErrorW_HighW_pd = new TH1D("h_ErrorW_HighW_pd",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ShowerMass_pd = new TH1D("h_ShowerMass_pd",";Shower W (GeV);Events",50,-0.01,0.5);
  TH1D* h_PiZero_ErrorW_pd = new TH1D("h_PiZero_ErrorW_pd",";(Reco - True)/True;Events",51,-2,2);

  // All Wirecell
  TH1D* h_RecoW_wc = new TH1D("h_RecoW_wc",";Reco W (GeV);Events",50,0.9,5.0);
  TH2D* h_TrueW_RecoW_wc = new TH2D("h_TrueW_RecoW_wc",";True W (GeV);Reco W (GeV);Events",50,0.9,5.0,50,0.9,5.0);
  TH1D* h_ErrorW_wc = new TH1D("h_ErrorW_wc",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ErrorW_HighW_wc = new TH1D("h_ErrorW_HighW_wc",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ShowerMass_wc = new TH1D("h_ShowerMass_wc",";Shower W (GeV);Events",50,-0.01,0.5);
  TH1D* h_PiZero_ErrorW_wc = new TH1D("h_PiZero_ErrorW_wc",";(Reco - True)/True;Events",51,-2,2);

  // All Lantern 
  TH1D* h_RecoW_lt = new TH1D("h_RecoW_lt",";Reco W (GeV);Events",50,0.9,5.0);
  TH2D* h_TrueW_RecoW_lt = new TH2D("h_TrueW_RecoW_lt",";True W (GeV);Reco W (GeV);Events",50,0.9,5.0,50,0.9,5.0);
  TH1D* h_ErrorW_lt = new TH1D("h_ErrorW_lt",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ErrorW_HighW_lt = new TH1D("h_ErrorW_HighW_lt",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ShowerMass_lt = new TH1D("h_ShowerMass_lt",";Shower W (GeV);Events",50,-0.015,0.5);
  TH1D* h_PiZero_ErrorW_lt = new TH1D("h_PiZero_ErrorW_lt",";(Reco - True)/True;Events",51,-2,2);

  // Hybrid 1 - LT protons and pions, WC showers
  TH1D* h_RecoW_h1 = new TH1D("h_RecoW_h1",";Reco W (GeV);Events",50,0.9,5.0);
  TH2D* h_TrueW_RecoW_h1 = new TH2D("h_TrueW_RecoW_h1",";True W (GeV);Reco W (GeV);Events",50,0.9,5.0,50,0.9,5.0);
  TH1D* h_ErrorW_h1 = new TH1D("h_ErrorW_h1",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ErrorW_HighW_h1 = new TH1D("h_ErrorW_HighW_h1",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ShowerMass_h1 = new TH1D("h_ShowerMass_h1",";Shower W (GeV);Events",50,-0.01,0.5);
  TH1D* h_PiZero_ErrorW_h1 = new TH1D("h_PiZero_ErrorW_h1",";(Reco - True)/True;Events",51,-2,2);

  for(int ievent=0;ievent<t_in->GetEntries();ievent++){

    if(ievent > 10000) break;
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    // Select out signal events
    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 14) continue; 
   
    TVector3 plepton(mc_px->at(0),mc_py->at(0),mc_pz->at(0));
    double p = plepton.Mag();
    double theta = plepton.Theta();
    if(p < 0.1) continue;

    TLorentzVector true_proton(0,0,0,0); 
    TLorentzVector true_pion(0,0,0,0); 
    TLorentzVector true_shower(0,0,0,0); 
    int n_pi0 = 0;
    int n_proton = 0;
    int n_pion = 0;
    for(size_t i_p=0;i_p<mc_pdg->size();i_p++){
      TVector3 mom(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p));
      if(mc_pdg->at(i_p) == 2212 && mom.Mag() > 0.3){
        n_proton++;
        true_proton += TLorentzVector(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p),mc_E->at(i_p));
      }
      if(abs(mc_pdg->at(i_p)) == 211 && mom.Mag() > 0.1){
        n_pion++;
        true_pion += TLorentzVector(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p),mc_E->at(i_p));
      }
      if(mc_pdg->at(i_p) == 111){
        true_shower += TLorentzVector(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p),mc_E->at(i_p));
        n_pi0++;
      }
    } 
    if(n_proton == 0) continue;

    // Count number of signal events
    signal++;

    // Calculate the true hadronic invariant mass
    double TrueW = pd::CalcTrueW(mc_pdg,mc_px,mc_py,mc_pz,mc_E);
    h_TrueW->Fill(TrueW); 

    bool pd_vertex = inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z);
    bool wc_vertex = inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ);
    bool lt_vertex = inActiveTPC(vtxX,vtxY,vtxZ);

    // Muon ID with each framework
    int pd_muon = pd::SimpleMuonSelection(trk_llr_pid_score_v,trk_len_v);
    int wc_muon = wc::SimpleMuonSelection(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT);
    int lt_muon = lt::SimpleMuonSelection(nTracks,trackIsSecondary,trackPID);

    TLorentzVector pd_proton = pd::RecoProton4Mom(trk_llr_pid_score_v,trk_len_v,trk_dir_x_v,trk_dir_y_v,trk_dir_z_v,pd_muon);
    TLorentzVector pd_pion = pd::RecoPion4Mom(trk_llr_pid_score_v,trk_len_v,trk_dir_x_v,trk_dir_y_v,trk_dir_z_v,pd_muon);
    TLorentzVector pd_shower = pd::RecoShower4Mom(shr_px_v,shr_py_v,shr_pz_v,shr_energy_y_v);

    TLorentzVector wc_proton = wc::RecoProton4Mom(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT,reco_id,wc_muon);
    TLorentzVector wc_pion = wc::RecoPion4Mom(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT,reco_id,wc_muon);
    TLorentzVector wc_shower = wc::RecoShower4Mom(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT,reco_id,wc_muon);

    TLorentzVector lt_proton = lt::RecoProton4Mom(nTracks,trackIsSecondary,trackPID,trackRecoE,trackStartDirX,trackStartDirY,trackStartDirZ);
    TLorentzVector lt_pion = lt::RecoPion4Mom(nTracks,trackIsSecondary,trackPID,trackRecoE,trackStartDirX,trackStartDirY,trackStartDirZ);
    TLorentzVector lt_shower = lt::RecoShower4Mom(nShowers,showerIsSecondary,showerPID,showerRecoE,showerStartDirX,showerStartDirY,showerStartDirZ);

    // Calculate the hadronic invariant mass with each framework first
    double pd_RecoW = (pd_proton+pd_pion+pd_shower).M();
    double wc_RecoW = (wc_proton+wc_pion+wc_shower).M();
    double lt_RecoW = (lt_proton+lt_pion+lt_shower).M();
    double h1_RecoW = (lt_proton+lt_pion+wc_shower).M();

    if(pd_vertex && pd_muon != -1 && pd_RecoW > Mp-0.01){
      h_RecoW_pd->Fill(pd_RecoW);
      h_TrueW_RecoW_pd->Fill(TrueW,pd_RecoW);                
      h_ErrorW_pd->Fill((pd_RecoW - TrueW)/TrueW);
      if(TrueW > 1) h_ErrorW_HighW_pd->Fill((pd_RecoW - TrueW)/TrueW);
      h_ShowerMass_pd->Fill(pd_shower.M());
      if(n_pi0 == 1) h_PiZero_ErrorW_pd->Fill((pd_RecoW - TrueW)/TrueW);
    }

    if(wc_vertex && wc_muon != -1 && wc_RecoW > Mp-0.01){
      h_RecoW_wc->Fill(wc_RecoW);
      h_TrueW_RecoW_wc->Fill(TrueW,wc_RecoW);                
      h_ErrorW_wc->Fill((wc_RecoW - TrueW)/TrueW);
      if(TrueW > 1) h_ErrorW_HighW_wc->Fill((wc_RecoW - TrueW)/TrueW);
      h_ShowerMass_wc->Fill(wc_shower.M());
      if(n_pi0 == 1) h_PiZero_ErrorW_wc->Fill((wc_RecoW - TrueW)/TrueW);
    }

    if(lt_vertex && lt_muon != -1 && lt_RecoW > Mp-0.01){
      h_RecoW_lt->Fill(lt_RecoW);
      h_TrueW_RecoW_lt->Fill(TrueW,lt_RecoW);                
      h_ErrorW_lt->Fill((lt_RecoW - TrueW)/TrueW);
      if(TrueW > 1) h_ErrorW_HighW_lt->Fill((lt_RecoW - TrueW)/TrueW);
      h_ShowerMass_lt->Fill(lt_shower.M());
      if(n_pi0 == 1) h_PiZero_ErrorW_lt->Fill((lt_RecoW - TrueW)/TrueW);
    }

    // Hybrid methods
    if(wc_vertex && lt_vertex && h1_RecoW > Mp-0.01){
      h_RecoW_h1->Fill(h1_RecoW);
      h_TrueW_RecoW_h1->Fill(TrueW,h1_RecoW);                
      h_ErrorW_h1->Fill((h1_RecoW - TrueW)/TrueW);
      if(TrueW > 1) h_ErrorW_HighW_h1->Fill((h1_RecoW - TrueW)/TrueW);
      h_ShowerMass_h1->Fill(wc_shower.M());
      if(n_pi0 == 1) h_PiZero_ErrorW_h1->Fill((h1_RecoW - TrueW)/TrueW);
    } 

  }

  std::cout << "Signal = " << signal*data_POT/POT << std::endl; 

  gSystem->Exec("mkdir -p Plots/CalcW/"); 
  TLegend* l = new TLegend(0.8,0.8,0.95,0.95);
  TCanvas* c = new TCanvas("c","c");

  h_TrueW->Scale(data_POT/POT);
  h_TrueW->SetLineWidth(2);
  h_TrueW->SetLineColor(1);
  h_TrueW->SetStats(0);
  h_TrueW->Draw("HIST");   
  c->Print("Plots/CalcW/TrueE.png");
  c->Clear();
 
  // Draw the W dist 
  THStack* hs_RecoW = new THStack("hs_RecoW",";Reco W (GeV);");

  h_RecoW_pd->Scale(data_POT/POT);
  h_RecoW_pd->SetLineWidth(2);
  h_RecoW_pd->SetLineColor(2);
  hs_RecoW->Add(h_RecoW_pd); 
  l->AddEntry(h_RecoW_pd,"pd","L");

  h_RecoW_wc->Scale(data_POT/POT);
  h_RecoW_wc->SetLineWidth(2);
  h_RecoW_wc->SetLineColor(3);
  hs_RecoW->Add(h_RecoW_wc); 
  l->AddEntry(h_RecoW_wc,"wc","L");

  h_RecoW_lt->Scale(data_POT/POT);
  h_RecoW_lt->SetLineWidth(2);
  h_RecoW_lt->SetLineColor(4);
  hs_RecoW->Add(h_RecoW_lt); 
  l->AddEntry(h_RecoW_lt,"lt","L");

  h_RecoW_h1->Scale(data_POT/POT);
  h_RecoW_h1->SetLineWidth(2);
  h_RecoW_h1->SetLineColor(5);
  hs_RecoW->Add(h_RecoW_h1); 
  l->AddEntry(h_RecoW_h1,"h1","L");

  hs_RecoW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW/RecoW.png");
  c->Clear();
  l->Clear();

  // Draw the W Error 
  THStack* hs_Error = new THStack("hs_Error",";(Reco - True)/True;");

  h_ErrorW_pd->Scale(data_POT/POT);
  h_ErrorW_pd->SetLineWidth(2);
  h_ErrorW_pd->SetLineColor(2);
  hs_Error->Add(h_ErrorW_pd); 
  l->AddEntry(h_ErrorW_pd,"pd","L");

  h_ErrorW_wc->Scale(data_POT/POT);
  h_ErrorW_wc->SetLineWidth(2);
  h_ErrorW_wc->SetLineColor(3);
  hs_Error->Add(h_ErrorW_wc); 
  l->AddEntry(h_ErrorW_wc,"wc","L");

  h_ErrorW_lt->Scale(data_POT/POT);
  h_ErrorW_lt->SetLineWidth(2);
  h_ErrorW_lt->SetLineColor(4);
  hs_Error->Add(h_ErrorW_lt); 
  l->AddEntry(h_ErrorW_lt,"lt","L");

  h_ErrorW_h1->Scale(data_POT/POT);
  h_ErrorW_h1->SetLineWidth(2);
  h_ErrorW_h1->SetLineColor(5);
  hs_Error->Add(h_ErrorW_h1); 
  l->AddEntry(h_ErrorW_h1,"h1","L");

  hs_Error->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW/ErrorW.png");
  c->Clear();
  l->Clear();

  // Draw the W Error 
  THStack* hs_Error_HighW = new THStack("hs_Error_HighW",";(Reco - True)/True;");

  h_ErrorW_HighW_pd->Scale(data_POT/POT);
  h_ErrorW_HighW_pd->SetLineWidth(2);
  h_ErrorW_HighW_pd->SetLineColor(2);
  hs_Error_HighW->Add(h_ErrorW_HighW_pd); 
  l->AddEntry(h_ErrorW_HighW_pd,"pd","L");

  h_ErrorW_HighW_wc->Scale(data_POT/POT);
  h_ErrorW_HighW_wc->SetLineWidth(2);
  h_ErrorW_HighW_wc->SetLineColor(3);
  hs_Error_HighW->Add(h_ErrorW_HighW_wc); 
  l->AddEntry(h_ErrorW_HighW_wc,"wc","L");

  h_ErrorW_HighW_lt->Scale(data_POT/POT);
  h_ErrorW_HighW_lt->SetLineWidth(2);
  h_ErrorW_HighW_lt->SetLineColor(4);
  hs_Error_HighW->Add(h_ErrorW_HighW_lt); 
  l->AddEntry(h_ErrorW_HighW_lt,"lt","L");

  h_ErrorW_HighW_h1->Scale(data_POT/POT);
  h_ErrorW_HighW_h1->SetLineWidth(2);
  h_ErrorW_HighW_h1->SetLineColor(5);
  hs_Error_HighW->Add(h_ErrorW_HighW_h1); 
  l->AddEntry(h_ErrorW_HighW_h1,"h1","L");

  hs_Error_HighW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW/ErrorW_HighW.png");
  c->Clear();
  l->Clear();

  Normalise(h_TrueW_RecoW_pd);
  h_TrueW_RecoW_pd->Draw("colz");
  h_TrueW_RecoW_pd->SetStats(0);
  c->Print("Plots/CalcW/TrueW_RecoW_pd.png");
  c->Clear();

  Normalise(h_TrueW_RecoW_wc);
  h_TrueW_RecoW_wc->Draw("colz");
  h_TrueW_RecoW_wc->SetStats(0);
  c->Print("Plots/CalcW/TrueW_RecoW_wc.png");
  c->Clear();

  Normalise(h_TrueW_RecoW_lt);
  h_TrueW_RecoW_lt->Draw("colz");
  h_TrueW_RecoW_lt->SetStats(0);
  c->Print("Plots/CalcW/TrueW_RecoW_lt.png");
  c->Clear();

  Normalise(h_TrueW_RecoW_h1);
  h_TrueW_RecoW_h1->Draw("colz");
  h_TrueW_RecoW_h1->SetStats(0);
  c->Print("Plots/CalcW/TrueW_RecoW_h1.png");
  c->Clear();


  // Invariant mass of showers
  THStack* hs_ShowerMass = new THStack("hs_ShowerMass",";Reco W (GeV);");

  h_ShowerMass_pd->Scale(data_POT/POT);
  h_ShowerMass_pd->SetLineWidth(2);
  h_ShowerMass_pd->SetLineColor(2);
  hs_ShowerMass->Add(h_ShowerMass_pd); 
  l->AddEntry(h_ShowerMass_pd,"pd","L");

  h_ShowerMass_wc->Scale(data_POT/POT);
  h_ShowerMass_wc->SetLineWidth(2);
  h_ShowerMass_wc->SetLineColor(3);
  hs_ShowerMass->Add(h_ShowerMass_wc); 
  l->AddEntry(h_ShowerMass_wc,"wc","L");

  h_ShowerMass_lt->Scale(data_POT/POT);
  h_ShowerMass_lt->SetLineWidth(2);
  h_ShowerMass_lt->SetLineColor(4);
  hs_ShowerMass->Add(h_ShowerMass_lt); 
  l->AddEntry(h_ShowerMass_lt,"lt","L");

  h_ShowerMass_h1->Scale(data_POT/POT);
  h_ShowerMass_h1->SetLineWidth(2);
  h_ShowerMass_h1->SetLineColor(5);
  hs_ShowerMass->Add(h_ShowerMass_h1); 
  l->AddEntry(h_ShowerMass_h1,"h1","L");

  hs_ShowerMass->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW/ShowerMass.png");
  c->Clear();
  l->Clear();

  // Draw the W dist 
  THStack* hs_PiZero_ErrorW = new THStack("hs_PiZero_ErrorW",";(Reco - True)/True;");

  h_PiZero_ErrorW_pd->Scale(data_POT/POT);
  h_PiZero_ErrorW_pd->SetLineWidth(2);
  h_PiZero_ErrorW_pd->SetLineColor(2);
  hs_PiZero_ErrorW->Add(h_PiZero_ErrorW_pd); 
  l->AddEntry(h_PiZero_ErrorW_pd,"pd","L");

  h_PiZero_ErrorW_wc->Scale(data_POT/POT);
  h_PiZero_ErrorW_wc->SetLineWidth(2);
  h_PiZero_ErrorW_wc->SetLineColor(3);
  hs_PiZero_ErrorW->Add(h_PiZero_ErrorW_wc); 
  l->AddEntry(h_PiZero_ErrorW_wc,"wc","L");

  h_PiZero_ErrorW_lt->Scale(data_POT/POT);
  h_PiZero_ErrorW_lt->SetLineWidth(2);
  h_PiZero_ErrorW_lt->SetLineColor(4);
  hs_PiZero_ErrorW->Add(h_PiZero_ErrorW_lt); 
  l->AddEntry(h_PiZero_ErrorW_lt,"lt","L");

  h_PiZero_ErrorW_h1->Scale(data_POT/POT);
  h_PiZero_ErrorW_h1->SetLineWidth(2);
  h_PiZero_ErrorW_h1->SetLineColor(5);
  hs_PiZero_ErrorW->Add(h_PiZero_ErrorW_h1); 
  l->AddEntry(h_PiZero_ErrorW_h1,"h1","L");

  hs_PiZero_ErrorW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW/PiZero_ErrorW.png");
  c->Clear();
  l->Clear();


}
