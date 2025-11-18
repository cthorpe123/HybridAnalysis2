#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "EnergyEstimatorFuncs.h"

// Hybrid 7 - muon kinematics from pandora and hadrons from LT
// Hybrid 8 - muon kinematics from pandora with MCS for uncontained muons, and hadrons from LT
// Hybrid 9 - muon kinematics from wirecell, with MCS for uncontained muons, with hadrons from LT

const std::vector<std::string> methods_str = {"PD","WC","LT","H7","H8","H9"};
enum methds_e {kpd,kwc,klt,h7,h8,h9,kMethMAX};

// Try using different combinations of cuts and frameworks to calculate W

void CalcE(){

  const double data_POT = 1.5e21;

  const std::string file = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTreeFiltered(file,f_in,t_in,false,false);
  const double POT = 7.88166e+20;

  TH1D* h_TrueE = new TH1D("h_TrueE",";True E (GeV);Events",50,0.0,3.0);
  std::vector<TH1D*> h_EstE_v;
  for(std::string est : ee::estimators_str){
    h_EstE_v.push_back(new TH1D(("h_EstE_"+est).c_str(),";True Est E_{#nu} (GeV);Events",50,0.0,3.0)) ;
  } 

  std::vector<std::vector<TH1D*>> h_RecoE_v_v;
  std::vector<std::vector<TH1D*>> h_ErrorEstE_v_v;
  std::vector<std::vector<TH2D*>> h_EstE_RecoE_v_v;
  std::vector<std::vector<TH2D*>> h_TrueE_RecoE_v_v;
  std::vector<std::vector<TH1D*>> h_ErrorTrueE_v_v;

  std::vector<std::vector<TH1D*>> h_Selected_RecoE_v_v;
  std::vector<std::vector<TH1D*>> h_SelectedCorr_RecoE_v_v;
  std::vector<std::vector<TH1D*>> h_SelectedBG_RecoE_v_v;

  for(std::string est : ee::estimators_str){
    h_RecoE_v_v.push_back(std::vector<TH1D*>());
    h_ErrorEstE_v_v.push_back(std::vector<TH1D*>());
    h_EstE_RecoE_v_v.push_back(std::vector<TH2D*>());
    h_TrueE_RecoE_v_v.push_back(std::vector<TH2D*>());
    h_ErrorTrueE_v_v.push_back(std::vector<TH1D*>());

    h_Selected_RecoE_v_v.push_back(std::vector<TH1D*>());
    h_SelectedCorr_RecoE_v_v.push_back(std::vector<TH1D*>());
    h_SelectedBG_RecoE_v_v.push_back(std::vector<TH1D*>());

    for(std::string meth : methods_str){
      h_RecoE_v_v.back().push_back(new TH1D(("h_RecoE_"+est+"_"+meth).c_str(),";Reco E_{#nu};Events",50,0.0,3.0));   
      h_EstE_RecoE_v_v.back().push_back(new TH2D(("h_EstE_RecoE_"+est+"_"+meth).c_str(),";Est E_{#nu} (GeV);Reco E_{#nu};Events",50,0.0,3.0,50,0.0,3.0));   
      h_ErrorEstE_v_v.back().push_back(new TH1D(("h_ErrorEstE_"+est+"_"+meth).c_str(),";(Reco E_{#nu} - Est E_{#nu})/Est E_{#nu};Events",50,-1,1));  
      h_TrueE_RecoE_v_v.back().push_back(new TH2D(("h_TrueE_RecoE_"+est+"_"+meth).c_str(),";True E_{#nu} (GeV);Rec E_{#nu};Events",50,0.0,3.0,50,0.0,3.0));   
      h_ErrorTrueE_v_v.back().push_back(new TH1D(("h_ErrorTrueE_"+est+"_"+meth).c_str(),";(Reco E_{#nu} - True E_{#nu})/True E_{#nu};Events",50,-1,1));  
      h_Selected_RecoE_v_v.back().push_back(new TH1D(("h_Selected_RecoE_"+est+"_"+meth).c_str(),";Reco E_{#nu};Events",50,0.0,3.0));
      h_SelectedCorr_RecoE_v_v.back().push_back(new TH1D(("h_SelectedCorr_RecoE_"+est+"_"+meth).c_str(),";Reco E_{#nu};Events",50,0.0,3.0));
      h_SelectedBG_RecoE_v_v.back().push_back(new TH1D(("h_SelectedBG_RecoE_"+est+"_"+meth).c_str(),";Reco E_{#nu};Events",50,0.0,3.0));
    }
  }

  for(int ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 50000) break;
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    //if(!is_signal_t) continue;
    //if(!muon_contained_t) continue;

    h_TrueE->Fill(nu_e); 

    for(int i_e=0;i_e<ee::kMAX;i_e++) h_EstE_v.at(i_e)->Fill(est_nu_e_t->at(i_e));

    bool sel_h7 = has_muon_pd && sel_lt; 
    TLorentzVector plepton_pd(muon_mom_pd->X(),muon_mom_pd->Y(),muon_mom_pd->Z(),sqrt(muon_mom_pd->Mag()*muon_mom_pd->Mag() + ml*ml));
    std::vector<double> est_nu_e_h7 = ee::GetEnergyEst(plepton_pd,W_lt,*proton_p4_lt,nprot_lt,*pion_p4_lt,npi_lt,*shower_p4_lt,nsh_lt);

    bool sel_h9 = has_muon_wc && sel_lt; 
    TVector3 muon_mom_h9 = muon_contained_wc ? *muon_mom_len_wc : *muon_mom_mcs_wc;
    TLorentzVector plepton_h9 = TLorentzVector(muon_mom_h9.X(),muon_mom_h9.Y(),muon_mom_h9.Z(),sqrt(ml*ml+muon_mom_h9.Mag()*muon_mom_h9.Mag()));
    std::vector<double> est_nu_e_h9 = ee::GetEnergyEst(plepton_h9,W_lt,*proton_p4_lt,nprot_lt,*pion_p4_lt,npi_lt,*shower_p4_lt,nsh_lt);

    std::vector<bool> sel_v = {sel_pd,sel_wc,sel_lt,sel_h7,sel_h8,sel_h9};    
    std::vector<bool> cont_muon_v = {muon_contained_pd,muon_contained_wc,muon_contained_lt,muon_contained_pd,muon_contained_pd,muon_contained_wc};
    std::vector<double> W_v = {W_pd,W_wc,W_lt,W_lt,W_lt,W_lt};    
    std::vector<std::vector<double>*> e_v = {est_nu_e_pd,est_nu_e_wc,est_nu_e_lt,&est_nu_e_h7,est_nu_e_h8,&est_nu_e_h9}; 

    for(int i=0;i<kMethMAX;i++){
      if(!sel_v.at(i)) continue;
      //if(!cont_muon_v.at(i)) continue;
      for(int i_e=0;i_e<ee::kMAX;i_e++){
        if(is_signal_t){
          if(est_nu_e_t->at(i_e) > 0){
            h_RecoE_v_v.at(i_e).at(i)->Fill(e_v.at(i)->at(i_e));
            h_EstE_RecoE_v_v.at(i_e).at(i)->Fill(est_nu_e_t->at(i_e),e_v.at(i)->at(i_e));
            h_ErrorEstE_v_v.at(i_e).at(i)->Fill((e_v.at(i)->at(i_e)-est_nu_e_t->at(i_e))/est_nu_e_t->at(i_e));
            h_TrueE_RecoE_v_v.at(i_e).at(i)->Fill(nu_e,e_v.at(i)->at(i_e));
            h_ErrorTrueE_v_v.at(i_e).at(i)->Fill((e_v.at(i)->at(i_e)-nu_e)/nu_e);
          }
        }

        if(is_signal_t) h_Selected_RecoE_v_v.at(i_e).at(i)->Fill(e_v.at(i)->at(i_e));
        else h_SelectedBG_RecoE_v_v.at(i_e).at(i)->Fill(e_v.at(i)->at(i_e));

        if(is_signal_t && abs(e_v.at(i)->at(i_e)-est_nu_e_t->at(i_e))/est_nu_e_t->at(i_e) < 0.2) h_SelectedCorr_RecoE_v_v.at(i_e).at(i)->Fill(e_v.at(i)->at(i_e));
      }

    }

  }

  gSystem->Exec("mkdir -p Plots/CalcE/"); 
  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");
  l->SetNColumns(2);

  h_TrueE->Scale(data_POT/POT);
  h_TrueE->SetLineWidth(2);
  h_TrueE->SetLineColor(1);
  h_TrueE->SetStats(0);
  h_TrueE->Draw("HIST");   
  c->Print("Plots/CalcE/TrueE.png");
  c->Clear();

  THStack* hs_EstE = new THStack("hs_EstE",";True Est. E_{#nu};Events"); 
  for(int i_e=0;i_e<ee::kMAX;i_e++){
    h_EstE_v.at(i_e)->Scale(data_POT/POT);
    h_EstE_v.at(i_e)->SetLineWidth(2);
    h_EstE_v.at(i_e)->SetLineColor(ee::colors.at(i_e));
    hs_EstE->Add(h_EstE_v.at(i_e)); 
    l->AddEntry(h_EstE_v.at(i_e),ee::estimators_leg.at(i_e).c_str(),"L");
  }

  hs_EstE->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcE/EstE.png");
  c->Clear();
  l->Clear();

  // Draw the estimated energy spectra using each estimator and reconstruction combination
  for(int i_e=0;i_e<ee::kMAX;i_e++){

    THStack* hs_RecoE = new THStack(("hs_RecoE_"+ee::estimators_str.at(i_e)).c_str(),";Reco E_{#nu};Events"); 

    for(int i=0;i<kMethMAX;i++){
      h_RecoE_v_v.at(i_e).at(i)->Scale(data_POT/POT);
      h_RecoE_v_v.at(i_e).at(i)->SetLineWidth(2);
      h_RecoE_v_v.at(i_e).at(i)->SetLineColor(i+1);
      hs_RecoE->Add(h_RecoE_v_v.at(i_e).at(i)); 
      l->AddEntry(h_RecoE_v_v.at(i_e).at(i),methods_str.at(i).c_str(),"L");
    }

    hs_RecoE->Draw("nostack HIST");
    l->Draw();
    c->Print(("Plots/CalcE/RecoE_"+ee::estimators_str.at(i_e)+".png").c_str());
    c->Clear();
    l->Clear();

    delete hs_RecoE;

  }

  // Draw the estimated energy spectra using each estimator and reconstruction combination
  for(int i_e=0;i_e<ee::kMAX;i_e++){
    for(int i=0;i<kMethMAX;i++){
      h_EstE_RecoE_v_v.at(i_e).at(i)->Draw("colz");
      h_EstE_RecoE_v_v.at(i_e).at(i)->SetStats(0);
      c->Print(("Plots/CalcE/EstE_RecoE_"+ee::estimators_str.at(i_e)+"_"+methods_str.at(i)+".png").c_str());
      c->Clear();
      Normalise(h_EstE_RecoE_v_v.at(i_e).at(i));
      h_EstE_RecoE_v_v.at(i_e).at(i)->Draw("colz");
      h_EstE_RecoE_v_v.at(i_e).at(i)->SetStats(0);
      c->Print(("Plots/CalcE/Normalised_EstE_RecoE_"+ee::estimators_str.at(i_e)+"_"+methods_str.at(i)+".png").c_str());
      c->Clear();
    }
  }

  // Draw the estimated energy spectra using each estimator and reconstruction combination
  for(int i_e=0;i_e<ee::kMAX;i_e++){

    THStack* hs_ErrorEstE = new THStack(("hs_ErrorEstE_"+ee::estimators_str.at(i_e)).c_str(),";(Reco E_{#nu} - Est. E_{#nu})/Est. E_{#nu};Events"); 

    for(int i=0;i<kMethMAX;i++){
      h_ErrorEstE_v_v.at(i_e).at(i)->Scale(data_POT/POT);
      h_ErrorEstE_v_v.at(i_e).at(i)->SetLineWidth(2);
      h_ErrorEstE_v_v.at(i_e).at(i)->SetLineColor(i+1);
      hs_ErrorEstE->Add(h_ErrorEstE_v_v.at(i_e).at(i)); 
      l->AddEntry(h_ErrorEstE_v_v.at(i_e).at(i),methods_str.at(i).c_str(),"L");
    }

    hs_ErrorEstE->Draw("nostack HIST");
    l->Draw();
    c->Print(("Plots/CalcE/ErrorEstE_"+ee::estimators_str.at(i_e)+".png").c_str());
    c->Clear();
    l->Clear();

    delete hs_ErrorEstE;

  }

  for(int i_e=0;i_e<ee::kMAX;i_e++){
    for(int i=0;i<kMethMAX;i++){
      h_TrueE_RecoE_v_v.at(i_e).at(i)->Draw("colz");
      h_TrueE_RecoE_v_v.at(i_e).at(i)->SetStats(0);
      c->Print(("Plots/CalcE/TrueE_RecoE_"+ee::estimators_str.at(i_e)+"_"+methods_str.at(i)+".png").c_str());
      c->Clear();
      Normalise(h_TrueE_RecoE_v_v.at(i_e).at(i));
      h_TrueE_RecoE_v_v.at(i_e).at(i)->Draw("colz");
      h_TrueE_RecoE_v_v.at(i_e).at(i)->SetStats(0);
      c->Print(("Plots/CalcE/Normalised_TrueE_RecoE_"+ee::estimators_str.at(i_e)+"_"+methods_str.at(i)+".png").c_str());
      c->Clear();
    }
  }

  for(int i_e=0;i_e<ee::kMAX;i_e++){

    THStack* hs_ErrorTrueE = new THStack(("hs_ErrorTrueE_"+ee::estimators_str.at(i_e)).c_str(),";(Reco E_{#nu} - True. E_{#nu})/True. E_{#nu};Events"); 

    for(int i=0;i<kMethMAX;i++){
      h_ErrorTrueE_v_v.at(i_e).at(i)->Scale(data_POT/POT);
      h_ErrorTrueE_v_v.at(i_e).at(i)->SetLineWidth(2);
      h_ErrorTrueE_v_v.at(i_e).at(i)->SetLineColor(i+1);
      hs_ErrorTrueE->Add(h_ErrorTrueE_v_v.at(i_e).at(i)); 
      l->AddEntry(h_ErrorTrueE_v_v.at(i_e).at(i),methods_str.at(i).c_str(),"L");
    }

    hs_ErrorTrueE->Draw("nostack HIST");
    l->Draw();
    c->Print(("Plots/CalcE/ErrorTrueE_"+ee::estimators_str.at(i_e)+".png").c_str());
    c->Clear();
    l->Clear();

    delete hs_ErrorTrueE;

  }

  for(int i_e=0;i_e<ee::kMAX;i_e++){

    // Calculate signal/sqrt(signal + BG) in reco space
    THStack* hs_SSB_RecoE = new THStack("hs_SSB_RecoE",";Reco E (GeV);S/#sqrt{S+B}");

    for(int i=0;i<kMethMAX;i++){
      h_Selected_RecoE_v_v.at(i_e).at(i)->Scale(data_POT/POT);
      h_SelectedCorr_RecoE_v_v.at(i_e).at(i)->Scale(data_POT/POT);
      h_SelectedBG_RecoE_v_v.at(i_e).at(i)->Scale(data_POT/POT);     
      h_Selected_RecoE_v_v.at(i_e).at(i)->Add(h_SelectedBG_RecoE_v_v.at(i_e).at(i));
      for(int b=1;b<h_Selected_RecoE_v_v.at(i_e).at(i)->GetNbinsX()+1;b++){
        if(h_Selected_RecoE_v_v.at(i_e).at(i)->GetBinContent(b) > 0){
          h_SelectedCorr_RecoE_v_v.at(i_e).at(i)->SetBinContent(b,h_SelectedCorr_RecoE_v_v.at(i_e).at(i)->GetBinContent(b)/sqrt(h_Selected_RecoE_v_v.at(i_e).at(i)->GetBinContent(b))); 
        }
        else h_SelectedCorr_RecoE_v_v.at(i_e).at(i)->SetBinContent(b,0);
      }
      h_SelectedCorr_RecoE_v_v.at(i_e).at(i)->SetLineColor(i+1);
      h_SelectedCorr_RecoE_v_v.at(i_e).at(i)->SetLineWidth(2);
      hs_SSB_RecoE->Add(h_SelectedCorr_RecoE_v_v.at(i_e).at(i));
      l->AddEntry(h_SelectedCorr_RecoE_v_v.at(i_e).at(i),methods_str.at(i).c_str(),"L");
    }

    hs_SSB_RecoE->Draw("nostack HIST");
    l->Draw();
    c->Print(("Plots/CalcE/SSB_"+ee::estimators_str.at(i_e)+".png").c_str());
    c->Clear();
    l->Clear();

    delete hs_SSB_RecoE;

  }

}
