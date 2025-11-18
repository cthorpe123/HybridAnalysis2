#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

// Try reconstructing the momentum carried by pi0s in signal events 

void CheckPi0(){

  const double data_POT = 1.5e21;

  const std::string file = "Filters/test.root";
  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTreeFiltered(file,f_in,t_in,false);
  const double POT = 7.88166e+20;

  // All Pandora
  TH1D* h_ShowerMom_Error_pd = new TH1D("h_ErrorW_pd",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ShowerMom_pd = new TH1D("h_ShowerMass_pd",";Shower Mom (GeV);Events",50,-0.01,0.5);
  TH1D* h_Showers_pd = new TH1D("h_Showers_pd",";N Showers;Events",6,-0.5,5.5);

  // All Wirecell
  TH1D* h_ShowerMom_Error_wc = new TH1D("h_ErrorW_wc",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ShowerMom_wc = new TH1D("h_ShowerMass_wc",";Shower Mom (GeV);Events",50,-0.01,0.5);
  TH1D* h_Showers_wc = new TH1D("h_Showers_wc",";N Showers;Events",6,-0.5,5.5);

  // All Lantern 
  TH1D* h_ShowerMom_Error_lt = new TH1D("h_ErrorW_lt",";(Reco - True)/True;Events",51,-2,2);
  TH1D* h_ShowerMom_lt = new TH1D("h_ShowerMass_lt",";Shower Mom (GeV);Events",50,-0.01,0.5);
  TH1D* h_Showers_lt = new TH1D("h_Showers_lt",";N Showers;Events",6,-0.5,5.5);

  for(int ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 10000) break;
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    if(!is_signal_t) continue;

    if(npi0_t != 1) continue;

    if(in_tpc_pd && has_muon_pd){
      h_ShowerMom_Error_pd->Fill((shower_p4_pd->Vect().Mag() - shower_p4_t->Vect().Mag())/shower_p4_t->Vect().Mag());
      h_ShowerMom_pd->Fill(shower_p4_pd->Vect().Mag());
      h_Showers_pd->Fill(nsh_pd);
    }

    if(in_tpc_wc && has_muon_wc){
      h_ShowerMom_Error_wc->Fill((shower_p4_wc->Vect().Mag() - shower_p4_t->Vect().Mag())/shower_p4_t->Vect().Mag());
      h_ShowerMom_wc->Fill(shower_p4_wc->Vect().Mag());
      h_Showers_wc->Fill(nsh_wc);
    }

    if(in_tpc_lt && has_muon_lt){
      h_ShowerMom_Error_lt->Fill((shower_p4_lt->Vect().Mag() - shower_p4_t->Vect().Mag())/shower_p4_t->Vect().Mag());
      h_ShowerMom_lt->Fill(shower_p4_lt->Vect().Mag());
      h_Showers_lt->Fill(nsh_lt);
    }

  }

  gSystem->Exec("mkdir -p Plots/CheckPi0/"); 
  TLegend* l = new TLegend(0.8,0.8,0.95,0.95);
  TCanvas* c = new TCanvas("c","c");
 
  // Draw the W dist 
  THStack* hs_ShowerMom = new THStack("hs_ShowerMom",";Shower Mom (GeV);");

  h_ShowerMom_pd->Scale(data_POT/POT);
  h_ShowerMom_pd->SetLineWidth(2);
  h_ShowerMom_pd->SetLineColor(2);
  hs_ShowerMom->Add(h_ShowerMom_pd); 
  l->AddEntry(h_ShowerMom_pd,"pd","L");

  h_ShowerMom_wc->Scale(data_POT/POT);
  h_ShowerMom_wc->SetLineWidth(2);
  h_ShowerMom_wc->SetLineColor(3);
  hs_ShowerMom->Add(h_ShowerMom_wc); 
  l->AddEntry(h_ShowerMom_wc,"wc","L");

  h_ShowerMom_lt->Scale(data_POT/POT);
  h_ShowerMom_lt->SetLineWidth(2);
  h_ShowerMom_lt->SetLineColor(4);
  hs_ShowerMom->Add(h_ShowerMom_lt); 
  l->AddEntry(h_ShowerMom_lt,"lt","L");

  hs_ShowerMom->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CheckPi0/ShowerMom.png");
  c->Clear();
  l->Clear();

  THStack* hs_ShowerMom_Error = new THStack("hs_ShowerMom_Error",";(Reco - True)/True;Events");

  h_ShowerMom_Error_pd->Scale(data_POT/POT);
  h_ShowerMom_Error_pd->SetLineWidth(2);
  h_ShowerMom_Error_pd->SetLineColor(2);
  hs_ShowerMom_Error->Add(h_ShowerMom_Error_pd); 
  l->AddEntry(h_ShowerMom_Error_pd,"pd","L");

  h_ShowerMom_Error_wc->Scale(data_POT/POT);
  h_ShowerMom_Error_wc->SetLineWidth(2);
  h_ShowerMom_Error_wc->SetLineColor(3);
  hs_ShowerMom_Error->Add(h_ShowerMom_Error_wc); 
  l->AddEntry(h_ShowerMom_Error_wc,"wc","L");

  h_ShowerMom_Error_lt->Scale(data_POT/POT);
  h_ShowerMom_Error_lt->SetLineWidth(2);
  h_ShowerMom_Error_lt->SetLineColor(4);
  hs_ShowerMom_Error->Add(h_ShowerMom_Error_lt); 
  l->AddEntry(h_ShowerMom_Error_lt,"lt","L");

  hs_ShowerMom_Error->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CheckPi0/ShowerMom_Error.png");
  c->Clear();
  l->Clear();

  // Draw the W dist 
  THStack* hs_Showers = new THStack("hs_Showers",";Showers;Events");

  h_Showers_pd->Scale(data_POT/POT);
  h_Showers_pd->SetLineWidth(2);
  h_Showers_pd->SetLineColor(2);
  hs_Showers->Add(h_Showers_pd); 
  l->AddEntry(h_Showers_pd,"pd","L");

  h_Showers_wc->Scale(data_POT/POT);
  h_Showers_wc->SetLineWidth(2);
  h_Showers_wc->SetLineColor(3);
  hs_Showers->Add(h_Showers_wc); 
  l->AddEntry(h_Showers_wc,"wc","L");

  h_Showers_lt->Scale(data_POT/POT);
  h_Showers_lt->SetLineWidth(2);
  h_Showers_lt->SetLineColor(4);
  hs_Showers->Add(h_Showers_lt); 
  l->AddEntry(h_Showers_lt,"lt","L");

  hs_Showers->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CheckPi0/Showers.png");
  c->Clear();
  l->Clear();


}
