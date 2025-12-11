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

void SysTest(){

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

    std::string plot_dir = "Analysis/"+label+"/Plots/SysTest/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());
    TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
    TCanvas* c = new TCanvas("c","c");

    // Print total covariance and total fractional covariance
    TH2D* h_Cov = static_cast<TH2D*>(f_in_hist->Get("Cov"));
    h_Cov->Draw("colz");
    h_Cov->SetStats(0);
    c->Print((plot_dir+"Cov.png").c_str());
    c->Clear(); 

    TH2D* h_FCov = static_cast<TH2D*>(f_in_hist->Get("FCov"));
    h_FCov->Draw("colz");
    h_FCov->SetStats(0);
    c->Print((plot_dir+"FCov.png").c_str());
    c->Clear(); 

    TH2D* h_Corr = CalcCorrelationMatrix("",h_Cov);
    h_Corr->Draw("colz");
    h_Corr->SetStats(0);
    c->Print((plot_dir+"Corr.png").c_str());
    c->Clear(); 

    // Print Genie, Flux and Reint covariance and fractional covariance
    for(int i_s=0;i_s<kSystMAX;i_s++){
      TH2D* h = static_cast<TH2D*>(f_in_hist->Get(("Cov_"+sys_str.at(i_s)).c_str()));
      h->Draw("colz");
      h->SetStats(0);
      c->Print((plot_dir+"Cov_"+sys_str.at(i_s)+".png").c_str());
      c->Clear();
      TH2D* h2 = CalcCorrelationMatrix(sys_str.at(i_s),h);
      h2->Draw("colz");
      c->Print((plot_dir+"Corr_"+sys_str.at(i_s)+".png").c_str());
      c->Clear(); 
      delete h;
      delete h2;
    }

    TH2D* h_Cov_MCStat = static_cast<TH2D*>(f_in_hist->Get("Cov_MCStat"));
    h_Cov_MCStat->Draw("colz");
    h_Cov_MCStat->SetStats(0);
    c->Print((plot_dir+"Cov_MCStat.png").c_str());
    c->Clear(); 

    TH2D* h_FCov_MCStat = static_cast<TH2D*>(f_in_hist->Get("FCov_MCStat"));
    h_FCov_MCStat->Draw("colz");
    h_FCov_MCStat->SetStats(0);
    c->Print((plot_dir+"FCov_MCStat.png").c_str());
    c->Clear(); 

    TH2D* h_Cov_EstDataStat = static_cast<TH2D*>(f_in_hist->Get("Cov_EstDataStat"));
    h_Cov_EstDataStat->Draw("colz");
    h_Cov_EstDataStat->SetStats(0);
    c->Print((plot_dir+"Cov_EstDataStat.png").c_str());
    c->Clear(); 

    TH2D* h_FCov_EstDataStat = static_cast<TH2D*>(f_in_hist->Get("FCov_EstDataStat"));
    h_FCov_EstDataStat->Draw("colz");
    h_FCov_EstDataStat->SetStats(0);
    c->Print((plot_dir+"FCov_EstDataStat.png").c_str());
    c->Clear(); 

    for(int i_s=0;i_s<kSystMAX;i_s++){
      TH2D* h = static_cast<TH2D*>(f_in_hist->Get(("FCov_"+sys_str.at(i_s)).c_str()));
      h->Draw("colz");
      h->SetStats(0);
      c->Print((plot_dir+"FCov_"+sys_str.at(i_s)+".png").c_str());
      delete h;
    }


    // Print the detector covarianc and fractional covariance  
    if(add_detvars){
      for(int i_s=0;i_s<kDetvarMAX;i_s++){
        TH2D* h = static_cast<TH2D*>(f_in_detvar->Get(("Cov_"+detvar_str.at(i_s)).c_str()));
        h->Draw("colz");
        h->SetStats(0);
        c->Print((plot_dir+"Cov_"+detvar_str.at(i_s)+".png").c_str());
        TH2D* h2 = CalcCorrelationMatrix(detvar_str.at(i_s),h);
        h2->Draw("colz");
        h2->SetStats(0);
        c->Print((plot_dir+"Corr_"+detvar_str.at(i_s)+".png").c_str());
        delete h;
        delete h2;
      }

      for(int i_s=0;i_s<kDetvarMAX;i_s++){
        TH2D* h = static_cast<TH2D*>(f_in_detvar->Get(("FCov_"+detvar_str.at(i_s)).c_str()));
        h->Draw("colz");
        h->SetStats(0);
        c->Print((plot_dir+"FCov_"+detvar_str.at(i_s)+".png").c_str());
        delete h;
      }
    }

    std::string axis_title = h_Cov->GetXaxis()->GetTitle();

    TH1D* h_FE_Tot = (TH1D*)f_in_hist->Get("h_CV_Tot")->Clone("h_FE_Tot");
    TH2D* h_Tot = (TH2D*)f_in_hist->Get("FCov")->Clone("FCov");

    THStack* hs_FE = new THStack("hs_FE",h_FE_Tot->GetTitle());

    // Calculate total frac error from all detvars combined
    TH2D* h_Tot_Detvar = nullptr;
    TH1D* h_FE_Tot_Detvar = nullptr;
    if(add_detvars){
      h_Tot_Detvar = (TH2D*)f_in_detvar->Get("FCov")->Clone("FCov_Detvar");
      h_FE_Tot_Detvar = (TH1D*)f_in_detvar->Get("h_Detvar_CV_Tot")->Clone("h_Detvar_FE_Tot");
      for(int i=0;i<h_FE_Tot_Detvar->GetNbinsX()+2;i++) h_FE_Tot_Detvar->SetBinContent(i,sqrt(h_Tot_Detvar->GetBinContent(i,i)));
      h_FE_Tot_Detvar->SetLineColor(detvar_color);
      h_FE_Tot_Detvar->SetLineWidth(2);
      l->AddEntry(h_FE_Tot_Detvar,"Detector","L");
      hs_FE->Add(h_FE_Tot_Detvar);
    }

    std::vector<TH1D*> h_FE;
    std::vector<std::string> legs;
    for(int i_s=0;i_s<kSystMAX;i_s++){
      h_FE.push_back((TH1D*)f_in_hist->Get("h_CV_Tot")->Clone(("h_FE_"+sys_str.at(i_s)).c_str()));
      TH2D* h = static_cast<TH2D*>(f_in_hist->Get(("FCov_"+sys_str.at(i_s)).c_str()));
      for(int i=0;i<h_FE.back()->GetNbinsX()+2;i++) h_FE.back()->SetBinContent(i,sqrt(h->GetBinContent(i,i)));
      h_FE.back()->SetLineColor(sys_color.at(i_s));
      h_FE.back()->SetLineWidth(2);
      delete h;
      hs_FE->Add(h_FE.back());
      l->AddEntry(h_FE.back(),sys_str.at(i_s).c_str(),"L");
      legs.push_back(sys_str.at(i_s));
    }

    h_FE.push_back((TH1D*)f_in_hist->Get("h_CV_Tot")->Clone("h_FE_MCStat"));
    for(int i=0;i<h_FE.back()->GetNbinsX()+2;i++) h_FE.back()->SetBinContent(i,sqrt(h_FCov_MCStat->GetBinContent(i,i)));
    h_FE.back()->SetLineColor(stat_color.at(kMCStat));
    h_FE.back()->SetLineWidth(2);
    hs_FE->Add(h_FE.back());
    l->AddEntry(h_FE.back(),"MCStat","L");
    legs.push_back("MCStat");

    h_FE.push_back((TH1D*)f_in_hist->Get("h_CV_Tot")->Clone("h_FE_EstDataStat"));
    for(int i=0;i<h_FE.back()->GetNbinsX()+2;i++) h_FE.back()->SetBinContent(i,sqrt(h_FCov_EstDataStat->GetBinContent(i,i)));
    h_FE.back()->SetLineColor(special_color.at(kEstDataStat));
    h_FE.back()->SetLineWidth(2);
    hs_FE->Add(h_FE.back());
    l->AddEntry( h_FE.back(),"EstDataStat","L");
    legs.push_back("EstDataStat");

    // Calculate total frac error from all uncertainties
    for(int i=0;i<h_FE_Tot->GetNbinsX()+2;i++){
      double total_cov = h_Tot->GetBinContent(i,i) + h_FCov_EstDataStat->GetBinContent(i,i);
      if(h_Tot_Detvar != nullptr) total_cov += h_Tot_Detvar->GetBinContent(i,i);
      h_FE_Tot->SetBinContent(i,sqrt(total_cov));
    }
    h_FE_Tot->SetLineColor(1);
    h_FE_Tot->SetLineWidth(3);
    l->AddEntry(h_FE_Tot,"Total","L");
    hs_FE->Add(h_FE_Tot);

    hs_FE->Draw("nostack HIST");
    hs_FE->GetXaxis()->SetTitle(axis_title.c_str());
    l->Draw();
    c->Print((plot_dir+"FE.png").c_str());
    c->Clear();
    l->Clear(); 

    if(add_detvars){
      THStack* hs_FE_Detvar = new THStack("hs_FE_Detvar",h_FE_Tot_Detvar->GetTitle());

      h_FE_Tot_Detvar->SetLineColor(detvar_color);
      hs_FE_Detvar->Add(h_FE_Tot_Detvar);
      l->AddEntry(h_FE_Tot_Detvar,"Total","L"); 

      std::vector<TH1D*> h_FE_Detvar;
      for(int i_s=0;i_s<kDetvarMAX;i_s++){
        h_FE_Detvar.push_back((TH1D*)f_in_detvar->Get("h_Detvar_CV_Tot")->Clone(("h_FE_"+detvar_str.at(i_s)).c_str()));
        TH2D* h = static_cast<TH2D*>(f_in_detvar->Get(("FCov_"+detvar_str.at(i_s)).c_str()));
        for(int i=0;i<h_FE_Detvar.back()->GetNbinsX()+2;i++) h_FE_Detvar.back()->SetBinContent(i,sqrt(h->GetBinContent(i,i)));
        delete h; 
        h_FE_Detvar.back()->SetLineColor(i_s+2);
        h_FE_Detvar.back()->SetLineWidth(2);
        l->AddEntry(h_FE_Detvar.back(),detvar_str.at(i_s).c_str(),"L");
        hs_FE_Detvar->Add(h_FE_Detvar.back());    
      }

      hs_FE_Detvar->Draw("nostack HIST");
      hs_FE_Detvar->GetXaxis()->SetTitle(axis_title.c_str());
      hs_FE_Detvar->GetYaxis()->SetTitle("Frac. Uncertainty");
      l->Draw();
      c->Print((plot_dir+"FE_Detvar.png").c_str());
      c->Clear();
      l->Clear(); 
    }

    delete c;

    if(draw_underflow || draw_overflow){

      THStack* hs_middle = new THStack("hs_middle",(";"+axis_title+";FE").c_str());
      THStack* hs_U = new THStack("hs_U",";;FE");
      THStack* hs_O = new THStack("hs_O",";;FE");
      TLegend* l2 = new TLegend(0.75,0.75,0.98,0.98);

      TH1D *h_FE_Tot_U,*h_FE_Tot_O;
      MakeOU(label,h_FE_Tot,h_FE_Tot_U,h_FE_Tot_O); 
      h_FE_Tot_U->SetLineWidth(2);
      h_FE_Tot_U->SetLineColor(1);
      h_FE_Tot_O->SetLineWidth(2);
      h_FE_Tot_O->SetLineColor(1);
      hs_middle->Add(h_FE_Tot);
      hs_U->Add(h_FE_Tot_U);
      hs_O->Add(h_FE_Tot_O);
      l2->AddEntry(h_FE_Tot,"Total","L");

      std::vector<TH1D*> h_O(h_FE.size());
      std::vector<TH1D*> h_U(h_FE.size());
      for(size_t i_s=0;i_s<h_FE.size();i_s++){
        MakeOU(h_FE.at(i_s)->GetName(),h_FE.at(i_s),h_U.at(i_s),h_O.at(i_s));
        h_U.at(i_s)->SetLineColor(h_FE.at(i_s)->GetLineColor());
        h_O.at(i_s)->SetLineColor(h_FE.at(i_s)->GetLineColor());
        h_U.at(i_s)->SetLineWidth(2);
        h_O.at(i_s)->SetLineWidth(2);
        hs_U->Add(h_U.at(i_s));
        hs_O->Add(h_O.at(i_s));
        hs_middle->Add(h_FE.at(i_s));
        l2->AddEntry(h_FE.at(i_s),legs.at(i_s).c_str(),"L");
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
        hs_U->Draw("nostack HIST");       
        //hs_U->SetMaximum(GetMax(h_FE_Tot_U)*1.1);
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
        hs_O->Draw("nostack HIST Y+");       
        //hs_O->SetMaximum(GetMax(h_FE_Tot_O)*1.1);
        hs_O->GetYaxis()->SetTitleSize(0.1);
        hs_O->GetXaxis()->SetLabelSize(0.16);
        hs_O->GetYaxis()->SetTitleOffset(1.8);
        hs_O->GetYaxis()->SetLabelSize(0.11);
        hs_O->GetYaxis()->SetLabelOffset(0.04);
        c2->cd();
      }    

      p_middle->cd(); 
      hs_middle->Draw("nostack HIST");
      //hs_middle->SetMaximum(GetMax(h_FE_Tot)*1.1);
      l2->Draw();

      c2->cd();
      c2->Print((plot_dir+"test.png").c_str());
      delete c2;

    }

    f_in_hist->Close();
    if(add_detvars) f_in_detvar->Close();




  }

}
