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

void ResponseCheck(){

  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");

  std::vector<std::string> label_v = {"W"};
  std::vector<std::string> special_univ_v = {"ExtraGamma","ExtraProt","ExtraPi"};

  for(size_t i_f=0;i_f<label_v.size();i_f++){

    std::string label = label_v.at(i_f);
    std::string plot_dir = "Analysis/"+label+"/Plots/ResponseCheck/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    TFile* f_in_res = TFile::Open(("Analysis/"+label+"/rootfiles/Response.root").c_str());

    TH1D* h_CV_Joint_Unrol = (TH1D*)f_in_res->Get("h_CV_Joint_Unrolled"); 
    TH2D* h_Cov = (TH2D*)f_in_res->Get("Cov_Tot"); 

    for(int i=1;i<h_Cov->GetNbinsX()+1;i++)
      if(!(h_CV_Joint_Unrol->GetBinContent(i) > 0)) 
        std::cout << "Bad response" << std::endl;

    for(int i=1;i<h_CV_Joint_Unrol->GetNbinsX()+1;i++)
      h_CV_Joint_Unrol->SetBinError(i,sqrt(h_Cov->GetBinContent(i,i)));

    h_CV_Joint_Unrol->SetLineColor(1);
    h_CV_Joint_Unrol->SetLineWidth(2);

    TMatrixDSym m_Cov = syst::MakeCovMat(h_Cov); 
    m_Cov.Invert();

    for(std::string su : special_univ_v){

       TH1D* h_cv_joint_unrol = (TH1D*)h_CV_Joint_Unrol->Clone("h_cv_joint_unrol");
       TH1D* h_su_joint_unrol = (TH1D*)f_in_res->Get(("h_Special_Joint_"+su+"_Unrolled").c_str());
        
       double sq_md = 0.0;
       for(int i=1;i<h_CV_Joint_Unrol->GetNbinsX()+1;i++)
         for(int j=1;j<h_CV_Joint_Unrol->GetNbinsX()+1;j++)
           sq_md += (h_su_joint_unrol->GetBinContent(i) - h_cv_joint_unrol->GetBinContent(i))*m_Cov[i-1][j-1]*(h_su_joint_unrol->GetBinContent(j) - h_cv_joint_unrol->GetBinContent(j));
       std::cout << label << " " << su << " sq_md/ndof = " << sq_md << "/" << h_cv_joint_unrol->GetNbinsX() << " = " << sq_md/h_cv_joint_unrol->GetNbinsX() << std::endl;
      
       THStack* hs = new THStack("hs","hs");

       hs->Add(h_cv_joint_unrol,"e1");
       l->AddEntry(h_cv_joint_unrol,"CV","M");

       h_su_joint_unrol->SetLineColor(2);
       h_su_joint_unrol->SetLineWidth(2);
       hs->Add(h_su_joint_unrol,"HIST");
       l->AddEntry(h_su_joint_unrol,su.c_str(),"L");

       hs->Draw("nostack");
       l->Draw();
       c->Print((plot_dir+"UnrolledResponses_"+su+".png").c_str()); 
       c->Clear();

       h_su_joint_unrol->Divide(h_cv_joint_unrol);
       h_cv_joint_unrol->Divide(h_cv_joint_unrol); 

       hs->Draw("nostack");
       l->Draw();
       c->Print((plot_dir+"UnrolledResponsesRatio_"+su+".png").c_str()); 
       c->Clear();

       delete hs;
       delete h_cv_joint_unrol;

        l->Clear();
    }

  }

} 

