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

void SensitivityTest(){

  bool blinded = true;

  std::vector<TH1D*> h_chi2_v(ee::kMAX);
  for(int i_e=0;i_e<ee::kMAX;i_e++){

    std::string label = "RecoE_"+ee::estimators_str.at(i_e);

    TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
    TFile* f_in_detvar = TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str());

    // Build the CV prediction
    TH1D* h_CV_Tot = static_cast<TH1D*>(f_in_hist->Get("h_CV_Tot"));

    std::vector<TH1D*> h_CV;
    for(size_t i_c=0;i_c<categories.size();i_c++){
      if(i_c == kData) continue;
      h_CV.push_back(static_cast<TH1D*>(f_in_hist->Get(("h_CV_"+categories.at(i_c)).c_str())));
    }

    // Build total covariance without flux  
    TH2D* h_Cov = static_cast<TH2D*>(f_in_hist->Get("Cov"));
    h_Cov->Reset();

    h_Cov->Add((TH2D*)f_in_hist->Get("Cov_MCStat"));


    for(int i_s=0;i_s<kSystMAX;i_s++){
      if(i_s == kFlux) continue;
      h_Cov->Add(static_cast<TH2D*>(f_in_hist->Get(("Cov_"+sys_str.at(i_s)).c_str())));
    }  

    TH2D* h_Cov_Detvar = static_cast<TH2D*>(f_in_detvar->Get("Cov"));
    //h_Cov->Add(h_Cov_Detvar);

    for(int i=1;i<h_CV_Tot->GetNbinsX()+1;i++) h_CV_Tot->SetBinError(i,sqrt(h_Cov->GetBinContent(i,i)));

    std::string plot_dir = "Analysis/"+label+"/Plots/SensitivityTest/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
    TCanvas* c = new TCanvas("c","c");
    l->SetNColumns(3);

    THStack* hs = new THStack("hs",h_CV_Tot->GetTitle());

    for(size_t i_c=0;i_c<categories.size();i_c++){
      if(i_c == kData) continue;
      h_CV.at(i_c)->SetFillColor(i_c+2);
      l->AddEntry(h_CV.at(i_c),categories.at(i_c).c_str(),"F");
      hs->Add(h_CV.at(i_c));
    }

    h_CV_Tot->SetFillStyle(3253);
    h_CV_Tot->SetFillColor(1);

    hs->Draw("HIST");
    h_CV_Tot->Draw("same e2");
    l->Draw();
    c->Print((plot_dir+"Prediction.png").c_str());
    c->Clear();

    TMatrixDSym m_Cov = MakeCovMat(h_Cov); 
    TMatrixDSym m_CovInv = m_Cov;
    m_CovInv.Invert();

    // Calculate the chi2 WRT to the CV for every flux universe, make histogram of values
    h_chi2_v.at(i_e) = new TH1D(("h_chi2_"+ee::estimators_str.at(i_e)).c_str(),";#chi^{2}/ndof;Toys",25,0.0,1.0);
    for(int i_u=0;i_u<nuniv_Flux;i_u++){
      TH1D* h_alt = (TH1D*)f_in_hist->Get(("h_Vars_Tot_Flux_"+std::to_string(i_u)).c_str());
      double chi2 = 0.0;
      for(int i=1;i<h_CV_Tot->GetNbinsX()+1;i++)
        for(int j=1;j<h_CV_Tot->GetNbinsX()+1;j++)
          chi2 += (h_CV_Tot->GetBinContent(i) - h_alt->GetBinContent(i))*m_CovInv[i-1][j-1]*(h_CV_Tot->GetBinContent(j) - h_alt->GetBinContent(j));

      h_chi2_v.at(i_e)->Fill(chi2/h_CV_Tot->GetNbinsX());

    } 

    h_chi2_v.at(i_e)->Draw("HIST");
    c->Print((plot_dir+"Chi2.png").c_str());
    c->Clear();

   c->Close();

  }


  std::string plot_dir = "Analysis/Plots/SensitivityTest/";
  gSystem->Exec(("mkdir -p "+plot_dir).c_str());
  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");
  l->SetNColumns(2);

  THStack* hs_chi2 = new THStack("hs_chi2",";#chi^{2}/ndof;Toys");

  for(int i_e=0;i_e<ee::kMAX;i_e++){
    h_chi2_v.at(i_e)->SetLineWidth(2);
    h_chi2_v.at(i_e)->SetLineColor(ee::colors.at(i_e));
    l->AddEntry(h_chi2_v.at(i_e),ee::estimators_str.at(i_e).c_str(),"L");
    hs_chi2->Add(h_chi2_v.at(i_e)); 
  }

  hs_chi2->Draw("nostack HIST");
  l->Draw();
  c->Print((plot_dir+"Chi2.png").c_str());
  c->Clear();
  l->Clear();

}
