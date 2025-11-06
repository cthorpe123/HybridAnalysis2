#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

// Hybrid 1 - LT protons and pions, WC showers
// Hybrid 2 - LT protons and pions, WC showers, only count events even numbers of showers
// Hybrid 3 - Agreement cut + WC for showers 
// Hybrid 4 - Agreement cut between PD and LT, use LT for kinematics 
// Hybrid 5 - Agreement cut between WC and LT, use WC for kinematics 
// Hybrid 6 - Agreement cut between WC and LT, use LT for kinematics 

const std::vector<std::string> methods_str = {"PD","WC","LT","H1","H2","H3","H4","H5","H6"};
enum methds_e {kpd,kwc,klt,kh1,kh2,kh3,kh4,kh5,kh6,kMethMAX};

// Try using different combinations of cuts and frameworks to calculate W

void CalcW2(){

  const double FOM_Cut = 0.15;

  const double data_POT = 1.5e21;

  const std::string file = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTreeFiltered(file,f_in,t_in,false,false);
  const double POT = 7.88166e+20;

  TH1D* h_TrueW = new TH1D("h_TrueW",";True W (GeV);Events",40,0.9,5.0);

  // Try lots of combinations of reconstructions

  std::vector<TH1D*> h_RecoW_v;
  std::vector<TH2D*> h_TrueW_RecoW_v;
  std::vector<TH2D*> h_TrueW_RecoW_AltBin_v;
  std::vector<TH1D*> h_ErrorW_v;
  std::vector<TH1D*> h_ErrorW_HighW_v;
  std::vector<TH1D*> h_ShowerMass_v;
  std::vector<TH1D*> h_PiZero_ErrorW_v;
  std::vector<TH1D*> h_Selected_v;
  std::vector<TH1D*> h_SelectedCorr_v;

  std::vector<TH1D*> h_Selected_RecoW_v;
  std::vector<TH1D*> h_SelectedCorr_RecoW_v;
  std::vector<TH1D*> h_SelectedBG_RecoW_v;

  for(std::string meth : methods_str){
    h_RecoW_v.push_back(new TH1D(("h_RecoW_"+meth).c_str(),";Reco W (GeV);Events",40,0.9,5.0));
    h_TrueW_RecoW_v.push_back(new TH2D(("h_TrueW_RecoW_"+meth).c_str(),";True W (GeV);Reco W (GeV);Events",40,0.9,5.0,40,0.9,5.0));
    h_TrueW_RecoW_AltBin_v.push_back(new TH2D(("h_TrueW_RecoW_AltBin_"+meth).c_str(),";True W (GeV);Reco W (GeV);Events",40,0.9,5.0,400,0.0,10.0));
    h_ErrorW_v.push_back(new TH1D(("h_ErrorW_"+meth).c_str(),";(Reco - True)/True;Events",51,-2,2));
    h_ErrorW_HighW_v.push_back(new TH1D(("h_ErrorW_HighW_"+meth).c_str(),";(Reco - True)/True;Events",51,-2,2));
    h_ShowerMass_v.push_back(new TH1D(("h_ShowerMass_"+meth).c_str(),";Shower W (GeV);Events",40,-0.01,0.5));
    h_PiZero_ErrorW_v.push_back(new TH1D(("h_PiZero_ErrorW_"+meth).c_str(),";(Reco - True)/True;Events",51,-2,2));

    h_Selected_v.push_back(new TH1D(("h_Selected_"+meth).c_str(),";True W (GeV);Events",40,0.9,5.0));
    h_SelectedCorr_v.push_back(new TH1D(("h_SelectedCorr_"+meth).c_str(),";True W (GeV);Events",40,0.9,5.0));

    h_Selected_RecoW_v.push_back(new TH1D(("h_Selected_RecoW_"+meth).c_str(),";Reco W (GeV);Events",40,0.9,5.0));
    h_SelectedCorr_RecoW_v.push_back(new TH1D(("h_SelectedCorr_RecoW_"+meth).c_str(),";Reco W (GeV);Events",40,0.9,5.0));
    h_SelectedBG_RecoW_v.push_back(new TH1D(("h_SelectedBG_RecoW_"+meth).c_str(),";Reco W (GeV);Events",40,0.9,5.0));

  }

  for(int ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 20000) break;
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    h_TrueW->Fill(W_t); 

    bool sel_h1 = sel_lt && in_tpc_wc; 
    double W_h1 = (*proton_p4_lt+*pion_p4_lt+*shower_p4_wc).M();

    bool sel_h2 = sel_lt && in_tpc_wc && nsh_wc % 2 == 0; 
    double W_h2 = (*proton_p4_lt+*pion_p4_lt+*shower_p4_wc).M();

    bool sel_h3 = sel_pd && sel_wc && nprot_pd == nprot_wc && npi_pd == npi_wc; 
    double W_h3 = (*proton_p4_wc+*pion_p4_wc+*shower_p4_wc).M();

    bool sel_h4 = sel_pd && sel_lt && nprot_pd == nprot_lt && npi_pd == npi_lt; 
    double W_h4 = (*proton_p4_lt+*pion_p4_lt+*shower_p4_lt).M();

    bool sel_h5 = sel_wc && sel_lt && nprot_wc == nprot_lt && npi_wc == npi_lt; 
    double W_h5 = (*proton_p4_wc+*pion_p4_wc+*shower_p4_wc).M();

    bool sel_h6 = sel_wc && sel_lt && nprot_wc == nprot_lt && npi_wc == npi_lt; 
    double W_h6 = (*proton_p4_lt+*pion_p4_lt+*shower_p4_lt).M();

    std::vector<bool> sel_v = {sel_pd,sel_wc,sel_lt,sel_h1,sel_h2,sel_h3,sel_h4,sel_h5,sel_h6};    
    std::vector<double> W_v = {W_pd,W_wc,W_lt,W_h1,W_h2,W_h3,W_h4,W_h5,W_h6};    
    std::vector<TLorentzVector*> shower_p4_v = {shower_p4_pd,shower_p4_wc,shower_p4_lt,shower_p4_wc,shower_p4_wc,shower_p4_wc,shower_p4_lt,shower_p4_wc,shower_p4_lt};

    for(int i=0;i<kMethMAX;i++){
      if(!sel_v.at(i)) continue;
      if(is_signal_t){
        h_RecoW_v.at(i)->Fill(W_v.at(i));
        h_TrueW_RecoW_v.at(i)->Fill(W_t,W_v.at(i));
        h_TrueW_RecoW_AltBin_v.at(i)->Fill(W_t,W_v.at(i));
        h_ErrorW_v.at(i)->Fill((W_v.at(i) - W_t)/W_t);
        if(W_t > 1.0) h_ErrorW_HighW_v.at(i)->Fill((W_v.at(i) - W_t)/W_t);
        h_ShowerMass_v.at(i)->Fill(shower_p4_v.at(i)->M());
        if(npi0_t == 1) h_PiZero_ErrorW_v.at(i)->Fill((W_v.at(i) - W_t)/W_t);   
        h_Selected_v.at(i)->Fill(W_t);
        if(abs(W_v.at(i) - W_t)/W_t < FOM_Cut) h_SelectedCorr_v.at(i)->Fill(W_t);
        h_Selected_RecoW_v.at(i)->Fill(W_v.at(i));
        if(abs(W_v.at(i) - W_t)/W_t < FOM_Cut) h_SelectedCorr_RecoW_v.at(i)->Fill(W_v.at(i));
      }
      else {
        h_SelectedBG_RecoW_v.at(i)->Fill(W_v.at(i));
      }
    }

  }

  gSystem->Exec("mkdir -p Plots/CalcW2/"); 
  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");
  l->SetNColumns(3);

  h_TrueW->Scale(data_POT/POT);
  h_TrueW->SetLineWidth(2);
  h_TrueW->SetLineColor(1);
  h_TrueW->SetStats(0);
  h_TrueW->Draw("HIST");   
  c->Print("Plots/CalcW2/TrueE.png");
  c->Clear();

  // Draw the W dist 
  THStack* hs_RecoW = new THStack("hs_RecoW",";Reco W (GeV);");

  for(int i=0;i<kMethMAX;i++){
    h_RecoW_v.at(i)->Scale(data_POT/POT);
    h_RecoW_v.at(i)->SetLineWidth(2);
    h_RecoW_v.at(i)->SetLineColor(i+1);
    hs_RecoW->Add(h_RecoW_v.at(i)); 
    l->AddEntry(h_RecoW_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_RecoW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/RecoW.png");
  c->Clear();
  l->Clear();

  // Draw the W Error 
  THStack* hs_ErrorW = new THStack("hs_ErrorW",";(Reco - True)/True;Events");

  for(int i=0;i<kMethMAX;i++){
    if(i == kh1 || i == kh3 || i == kh5 || i == kh6) continue;
    h_ErrorW_v.at(i)->Scale(data_POT/POT);
    h_ErrorW_v.at(i)->SetLineWidth(2);
    h_ErrorW_v.at(i)->SetLineColor(i+1);
    hs_ErrorW->Add(h_ErrorW_v.at(i)); 
    l->AddEntry(h_ErrorW_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_ErrorW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/ErrorW.png");
  c->Clear();
  l->Clear();

  // Draw the W Error 
  THStack* hs_Normalised_ErrorW = new THStack("hs_Normalised_ErrorW",";(Reco - True)/True;");

  for(int i=0;i<kMethMAX;i++){
    if(i == kh1 || i == kh3 || i == kh5 || i == kh6) continue;
    h_ErrorW_v.at(i)->Scale(1.0/h_ErrorW_v.at(i)->Integral());
    h_ErrorW_v.at(i)->SetLineWidth(2);
    h_ErrorW_v.at(i)->SetLineColor(i+1);
    hs_Normalised_ErrorW->Add(h_ErrorW_v.at(i)); 
    l->AddEntry(h_ErrorW_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_Normalised_ErrorW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/Normalised_ErrorW.png");
  c->Clear();
  l->Clear();

  // Draw the W Error in high W only 
  THStack* hs_ErrorW_HighW = new THStack("hs_ErrorW_HighW",";(Reco - True)/True;");

  for(int i=0;i<kMethMAX;i++){
    if(i == kh1 || i == kh3 || i == kh5 || i == kh6) continue;
    h_ErrorW_HighW_v.at(i)->Scale(data_POT/POT);
    h_ErrorW_HighW_v.at(i)->SetLineWidth(2);
    h_ErrorW_HighW_v.at(i)->SetLineColor(i+1);
    hs_ErrorW_HighW->Add(h_ErrorW_HighW_v.at(i)); 
    l->AddEntry(h_ErrorW_HighW_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_ErrorW_HighW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/ErrorW_HighW.png");
  c->Clear();
  l->Clear();

  THStack* hs_Normalised_ErrorW_HighW = new THStack("hs_Normalised_ErrorW_HighW",";(Reco - True)/True;");

  for(int i=0;i<kMethMAX;i++){
    if(i == kh1 || i == kh3 || i == kh5 || i == kh6) continue;
    h_ErrorW_HighW_v.at(i)->Scale(1.0/h_ErrorW_HighW_v.at(i)->Integral());
    h_ErrorW_HighW_v.at(i)->SetLineWidth(2);
    h_ErrorW_HighW_v.at(i)->SetLineColor(i+1);
    hs_Normalised_ErrorW_HighW->Add(h_ErrorW_HighW_v.at(i)); 
    l->AddEntry(h_ErrorW_HighW_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_Normalised_ErrorW_HighW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/Normalised_ErrorW_HighW.png");
  c->Clear();
  l->Clear();

  // Draw the W Error 
  THStack* hs_ShowerMass = new THStack("hs_ShowerMass",";(Reco - True)/True;");

  for(int i=0;i<kMethMAX;i++){
    h_ShowerMass_v.at(i)->Scale(data_POT/POT);
    h_ShowerMass_v.at(i)->SetLineWidth(2);
    h_ShowerMass_v.at(i)->SetLineColor(i+1);
    hs_ShowerMass->Add(h_ShowerMass_v.at(i)); 
    l->AddEntry(h_ShowerMass_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_ShowerMass->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/ShowerMass.png");
  c->Clear();
  l->Clear();

  // Draw the W Error 
  THStack* hs_PiZero_ErrorW = new THStack("hs_PiZero_ErrorW",";(Reco - True)/True;");

  for(int i=0;i<kMethMAX;i++){
    h_PiZero_ErrorW_v.at(i)->Scale(data_POT/POT);
    h_PiZero_ErrorW_v.at(i)->SetLineWidth(2);
    h_PiZero_ErrorW_v.at(i)->SetLineColor(i+1);
    hs_PiZero_ErrorW->Add(h_PiZero_ErrorW_v.at(i)); 
    l->AddEntry(h_PiZero_ErrorW_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_PiZero_ErrorW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/PiZero_ErrorW.png");
  c->Clear();
  l->Clear();

  for(int i=0;i<kMethMAX;i++){
    TH2D* h = h_TrueW_RecoW_v.at(i);
    Normalise(h);
    h->Draw("colz");
    h->SetStats(0);
    c->Print(("Plots/CalcW2/TrueW_RecoW_"+methods_str.at(i)+".png").c_str());
    c->Clear();
  }

  // Calculate bias and variance AFO W
  THStack* hs_Bias = new THStack("hs_Bias",";True W (GeV);Frac. Bias");  
  THStack* hs_Variance = new THStack("hs_Variance",";True W (GeV);Frac. Variance");  
  std::vector<TH1D*> h_Bias_v(kMethMAX);
  std::vector<TH1D*> h_Variance_v(kMethMAX);

  for(int i=0;i<kMethMAX;i++){
    TH2D* h = h_TrueW_RecoW_AltBin_v.at(i);
    int nbins = h->GetNbinsX();
    double low = h->GetXaxis()->GetBinLowEdge(1);
    double high = h->GetXaxis()->GetBinLowEdge(nbins+1);
    h_Bias_v.at(i) = new TH1D(("h_Bias_"+methods_str.at(i)).c_str(),";True W (GeV);Bias",nbins,low,high);
    h_Variance_v.at(i) = new TH1D(("h_Variance_"+methods_str.at(i)).c_str(),";True W (GeV);Variance",nbins,low,high);
    GetBiasVariance(h,h_Bias_v.at(i),h_Variance_v.at(i));  

    h_Bias_v.at(i)->SetLineWidth(2);
    h_Bias_v.at(i)->SetLineColor(i+1);
    hs_Bias->Add(h_Bias_v.at(i));
    l->AddEntry(h_Bias_v.at(i),methods_str.at(i).c_str(),"L");

    h_Variance_v.at(i)->SetLineWidth(2);
    h_Variance_v.at(i)->SetLineColor(i+1);
    hs_Variance->Add(h_Variance_v.at(i));

  }

  hs_Bias->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/Bias.png");
  c->Clear();

  hs_Variance->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/Variance.png");
  c->Clear();

  l->Clear();

  // Efficiency AFO True W
  THStack* hs_Eff = new THStack("hs_Eff",";True W (GeV);");
  std::vector<TH1D*> h_Eff_v(kMethMAX);
  for(int i=0;i<kMethMAX;i++){
    if(i == kh1 || i == kh3 || i == kh5 || i == kh6) continue;
    h_Eff_v.at(i) = static_cast<TH1D*>(h_Selected_v.at(i)->Clone(("h_Eff_"+methods_str.at(i)).c_str()));
    h_Eff_v.at(i)->Scale(data_POT/POT);
    h_Eff_v.at(i)->Divide(h_TrueW);
    h_Eff_v.at(i)->SetLineWidth(2);
    h_Eff_v.at(i)->SetLineColor(i+1);
    hs_Eff->Add(h_Eff_v.at(i)); 
    l->AddEntry(h_Eff_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_Eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/Eff.png");
  c->Clear();
  l->Clear();

  // Efficiency AFO True W, only accepting events within 10% of 
  // correct value
  THStack* hs_EffCorr = new THStack("hs_EffCorr",";True W (GeV);");
  std::vector<TH1D*> h_EffCorr_v(kMethMAX);

  for(int i=0;i<kMethMAX;i++){
    h_EffCorr_v.at(i) = static_cast<TH1D*>(h_SelectedCorr_v.at(i)->Clone(("h_EffCorr_"+methods_str.at(i)).c_str()));
    h_EffCorr_v.at(i)->Scale(data_POT/POT);
    h_EffCorr_v.at(i)->Divide(h_TrueW);
    h_EffCorr_v.at(i)->SetLineWidth(2);
    h_EffCorr_v.at(i)->SetLineColor(i+1);
    hs_EffCorr->Add(h_EffCorr_v.at(i)); 
    l->AddEntry(h_EffCorr_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_EffCorr->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/EffCorr.png");
  c->Clear();
  l->Clear();

  // Calculate the purity (fraction of selected events with W within
  // 10% of the correct value)
  THStack* hs_Pur = new THStack("hs_Pur",";True W (GeV);");
  std::vector<TH1D*> h_Pur_v(kMethMAX);
  std::vector<TH1D*> h_Denom_v(kMethMAX);

  for(int i=0;i<kMethMAX;i++){
    h_Denom_v.at(i) = static_cast<TH1D*>(h_Selected_v.at(i)->Clone(("h_Denom_"+methods_str.at(i)).c_str()));
    h_Pur_v.at(i) = static_cast<TH1D*>(h_SelectedCorr_v.at(i)->Clone(("h_Pur_"+methods_str.at(i)).c_str()));
    h_Pur_v.at(i)->Divide(h_Denom_v.at(i));
    h_Pur_v.at(i)->SetLineWidth(2);
    h_Pur_v.at(i)->SetLineColor(i+1);
    hs_Pur->Add(h_Pur_v.at(i)); 
    l->AddEntry(h_Pur_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_Pur->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/Pur.png");
  c->Clear();
  l->Clear();

  // Efficiency x Purity
  THStack* hs_ExP = new THStack("hs_ExP",";True W (GeV);");

  for(int i=0;i<kMethMAX;i++){
    h_EffCorr_v.at(i)->Multiply(h_Pur_v.at(i));     
    hs_ExP->Add(h_EffCorr_v.at(i));
    l->AddEntry(h_EffCorr_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_ExP->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/ExP.png");
  c->Clear();
  l->Clear();
   
  // Calculate signal/sqrt(signal + BG) in reco space
  //THStack* hs_SSB_RecoW = new THStack("hs_SSB_RecoW",";Reco W (GeV);S/#sqrt{S+B}");
  THStack* hs_SSB_RecoW = new THStack("hs_SSB_RecoW",";Reco W (GeV);FOM");
  
  for(int i=0;i<kMethMAX;i++){
    if(i == kh1 || i == kh3 || i == kh5 || i == kh6) continue;
    h_Selected_RecoW_v.at(i)->Scale(data_POT/POT);
    h_SelectedCorr_RecoW_v.at(i)->Scale(data_POT/POT);
    h_SelectedBG_RecoW_v.at(i)->Scale(data_POT/POT);     
    h_Selected_RecoW_v.at(i)->Add(h_SelectedBG_RecoW_v.at(i));
    for(int b=1;b<h_Selected_RecoW_v.at(i)->GetNbinsX()+1;b++){
      if(h_Selected_RecoW_v.at(i)->GetBinContent(b) > 0){

        h_SelectedCorr_RecoW_v.at(i)->SetBinContent(b,h_SelectedCorr_RecoW_v.at(i)->GetBinContent(b)/sqrt(h_Selected_RecoW_v.at(i)->GetBinContent(b))); 
      }
      else h_SelectedCorr_RecoW_v.at(i)->SetBinContent(b,0);
    }
    h_SelectedCorr_RecoW_v.at(i)->SetLineColor(i+1);
    h_SelectedCorr_RecoW_v.at(i)->SetLineWidth(2);
    hs_SSB_RecoW->Add(h_SelectedCorr_RecoW_v.at(i));
    l->AddEntry(h_SelectedCorr_RecoW_v.at(i),methods_str.at(i).c_str(),"L");
  }

  hs_SSB_RecoW->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/CalcW2/SSB_RecoW.png");
  c->Clear();
  l->Clear();



}
