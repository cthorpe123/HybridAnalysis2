#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "EnergyEstimatorFuncs.h"
#include "Systematics.h"

using namespace syst;

void MakePrediction(){

  bool add_detvars = false;
  bool blinded = true;

  std::vector<std::string> label_v = {"MuonMom","MuonCosTheta","ProtonE","PionE","PiZeroE","NProt","NPi","NShr","W"};
  std::vector<bool> draw_underflow_v = {false,false,false,true,true,false,false,false,true};
  std::vector<bool> draw_overflow_v = {false,false,false,false,false,false,false,false,false};

  for(int i_e=0;i_e<ee::kMAX;i_e++){
    label_v.push_back(ee::estimators_str.at(i_e));
    draw_underflow_v.push_back(false);
    draw_overflow_v.push_back(false);
  }

  for(size_t i_f=0;i_f<label_v.size();i_f++){
    std::string label = label_v.at(i_f);
    bool draw_underflow = draw_underflow_v.at(i_f);
    bool draw_overflow = draw_overflow_v.at(i_f);

    TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
    TFile* f_in_detvar = add_detvars ? TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str()) : nullptr;

    // Build the CV prediction
    TH1D* h_CV_Tot = static_cast<TH1D*>(f_in_hist->Get("h_CV_Tot"));
    std::string axis_title = ";" + string(h_CV_Tot->GetXaxis()->GetTitle()) + ";" + h_CV_Tot->GetYaxis()->GetTitle();

    std::vector<TH1D*> h_CV;
    for(size_t i_c=0;i_c<categories.size();i_c++){
      if(i_c == kData) continue;
      h_CV.push_back(static_cast<TH1D*>(f_in_hist->Get(("h_CV_"+categories.at(i_c)).c_str())));
    }

    // Build total covariance 
    TH2D* h_Cov_Sum = static_cast<TH2D*>(f_in_hist->Get(("Cov_"+sys_str.at(0)).c_str())->Clone("h_Cov_Sum"));
    h_Cov_Sum->Reset();
    for(int i_s=0;i_s<kSystMAX;i_s++)
      h_Cov_Sum->Add(static_cast<TH2D*>(f_in_hist->Get(("Cov_"+sys_str.at(i_s)).c_str())));

    if(add_detvars){
      TH2D* h_Cov_Sum_Detvar = static_cast<TH2D*>(f_in_detvar->Get(("Cov_"+detvar_str.at(0)).c_str())->Clone("h_Cov_Sum_Detvar"));
      h_Cov_Sum_Detvar->Reset();
      for(int i_s=0;i_s<kDetvarMAX;i_s++) h_Cov_Sum_Detvar->Add(static_cast<TH2D*>(f_in_detvar->Get(("Cov_"+detvar_str.at(i_s)).c_str())));
      h_Cov_Sum->Add(h_Cov_Sum_Detvar);
    }

    for(int i=1;i<h_CV_Tot->GetNbinsX()+1;i++) h_CV_Tot->SetBinError(i,sqrt(h_Cov_Sum->GetBinContent(i,i)));

    // If not blindded, load the data
    TH1D* h_Data = !blinded ? (TH1D*)f_in_hist->Get("h_CV_Data") : nullptr; 

    std::string plot_dir = "Analysis/"+label+"/Plots/MakePrediction/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
    TCanvas* c = new TCanvas("c","c");
    l->SetNColumns(2);

    THStack* hs = new THStack("hs",axis_title.c_str());

    for(size_t i_c=0;i_c<categories.size();i_c++){
      if(i_c == kData) continue;
      h_CV.at(i_c)->SetFillColor(colors[i_c]);
      l->AddEntry(h_CV.at(i_c),categories.at(i_c).c_str(),"F");
      hs->Add(h_CV.at(i_c));
    }

    h_CV_Tot->SetFillStyle(3253);
    h_CV_Tot->SetFillColor(1);

    if(!blinded){
      h_Data->SetMarkerStyle(20);
      h_Data->SetMarkerSize(0.8);
      h_Data->SetLineColor(1);
    }

    hs->Draw("HIST");
    h_CV_Tot->Draw("same e2");
    if(!blinded) h_Data->Draw("same e1");
    l->Draw();
    c->Print((plot_dir+"Prediction.png").c_str());
    c->Clear();

    delete c;

    if(draw_underflow || draw_overflow){

      TLegend* l2 = new TLegend(0.75,0.75,0.98,0.98);
      l2->SetNColumns(2);
      THStack* hs_middle = new THStack("hs_middle",axis_title.c_str());
      for(size_t i_c=0;i_c<categories.size();i_c++){
        if(i_c == kData) continue;
        h_CV.at(i_c)->SetFillColor(colors[i_c]);
        l2->AddEntry(h_CV.at(i_c),categories.at(i_c).c_str(),"F");
        hs_middle->Add(h_CV.at(i_c));
      }
      h_CV_Tot->SetFillStyle(3253);
      h_CV_Tot->SetFillColor(1);

      TH1D *h_CV_O,*h_CV_U;
      MakeOU(label,h_CV_Tot,h_CV_U,h_CV_O,"","",h_Cov_Sum);
      h_CV_U->SetFillStyle(3253);
      h_CV_U->SetFillColor(1);
      h_CV_O->SetFillStyle(3253);
      h_CV_O->SetFillColor(1);

      std::vector<TH1D*> h_O(categories.size());
      std::vector<TH1D*> h_U(categories.size());
      THStack* hs_U = new THStack("hs_U","");
      THStack* hs_O = new THStack("hs_O","");
      for(size_t i_c=0;i_c<categories.size();i_c++){
        if(i_c == kData) continue;
        MakeOU(label+"_"+categories.at(i_c),h_CV.at(i_c),h_U.at(i_c),h_O.at(i_c));
        h_U.at(i_c)->SetFillColor(colors[i_c]);
        h_O.at(i_c)->SetFillColor(colors[i_c]);
        hs_U->Add(h_U.at(i_c));
        hs_O->Add(h_O.at(i_c));
      }

      TH1D *h_Data_U,*h_Data_O;
      if(!blinded) MakeOU(label+"_data",h_Data,h_Data_U,h_Data_O);

      if(!blinded){
        h_Data->SetMarkerStyle(20);
        h_Data->SetMarkerSize(0.8);
        h_Data->SetLineColor(1);
      }

      TCanvas* c2 = new TCanvas("c2","c2",1000,600);
      double  split_low = draw_underflow ? 0.18 : 0.0;
      double  split_high = draw_overflow ? 0.82 : 1.0;
      TPad* p_U = draw_underflow ? new TPad("p_U","p_U",0.0,0.0,split_low,1.0) : nullptr;
      TPad* p_middle = new TPad("p_middle","p_middle",split_low,0.0,split_high,1.0);
      TPad* p_O = draw_overflow ? new TPad("p_O","p_O",split_high,0.0,1.0,1.0) : nullptr;
      
      c2->cd();
      p_middle->Draw();
      p_middle->SetRightMargin(0.03);
      if(draw_underflow){
        p_U->Draw();
        p_U->SetLeftMargin(0.37);
        p_U->SetRightMargin(0.04);
        p_U->cd();
        hs_U->Draw("HIST");       
        h_CV_U->Draw("same e2");
        hs_U->SetMaximum(GetMax(h_CV_U)*1.1);
        if(h_CV_U->GetBinContent(1) < 1) hs_U->GetYaxis()->SetTitle("PMF");
        else hs_U->GetYaxis()->SetTitle("Events");
        if(!blinded) h_Data_U->Draw("same e1");
        hs_U->GetXaxis()->SetLabelSize(0.16);
        hs_U->GetYaxis()->SetTitleSize(0.1);
        hs_U->GetYaxis()->SetLabelSize(0.11);
        hs_U->GetYaxis()->SetLabelOffset(0.04);
        hs_U->GetYaxis()->SetTitleOffset(1.8);
        c2->cd();
      }
      if(draw_overflow){
        p_O->Draw(); 
        p_O->SetRightMargin(0.37);
        p_O->SetLeftMargin(0.04);
        p_O->cd();
        hs_O->Draw("HIST Y+");       
        h_CV_O->Draw("same e2");
        hs_O->SetMaximum(GetMax(h_CV_O)*1.1);
        if(h_CV_O->GetBinContent(1) < 1) hs_O->GetYaxis()->SetTitle("PMF");
        else hs_O->GetYaxis()->SetTitle("Events");
        hs_O->GetYaxis()->SetTitleSize(0.1);
        hs_O->GetXaxis()->SetLabelSize(0.16);
        hs_O->GetYaxis()->SetTitleOffset(1.8);
        hs_O->GetYaxis()->SetLabelSize(0.11);
        hs_O->GetYaxis()->SetLabelOffset(0.04);
        if(!blinded) h_Data_O->Draw("same e1");
        c2->cd();
       }    

      p_middle->cd(); 
      hs_middle->Draw("HIST");
      h_CV_Tot->Draw("same e2");
      hs_middle->SetMaximum(GetMax(h_CV_Tot)*1.1);
      if(!blinded) h_Data->Draw("same e1");

      c2->cd();
      l2->Draw();
      c2->Print((plot_dir+"test.png").c_str());

      delete c2;
    }

    f_in_hist->Close();
    if(add_detvars) f_in_detvar->Close();

  }

}
