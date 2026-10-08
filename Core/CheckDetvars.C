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
#include "DetvarHistograms.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"

using namespace syst;

void CheckDetvars(){

  // Label and set the branches defining the selection and systematics
  bool draw_truth = true;
  bool draw_hist = true; // Grab the CV from the non-detvar file
  bool draw_o=true,draw_u=true;

  std::vector<std::string> vars = {"MuonMom"/*,"MuonCosTheta","NProt","NPi","NSh","ProtonKE","PionE","PiZeroE","W"*/};

  for(std::string label : vars){

    std::string plot_dir = AnalysisDir()+"/"+label+"/Plots/CheckDetvars/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    TFile* f_in = TFile::Open((AnalysisDir()+"/"+label+"/rootfiles/Detvars.root").c_str());
    TFile* f_in_hist = draw_hist ? TFile::Open((AnalysisDir()+"/"+label+"/rootfiles/Histograms.root").c_str()) : nullptr;
 
    hist::MultiChannelHistogramManager mchm(label,true);
    mchm.LoadTemplates();

    std::vector<TH1D*> h_v;
    std::vector<int> fill_colors;
    std::vector<std::string> legs;

    h_v.push_back((TH1D*)f_in->Get("Reco/CV/h_Tot"));
    mchm.Restore(h_v.back());
    fill_colors.push_back(1);
    legs.push_back("CV");

    for(size_t i_s=0;i_s<detvar_str.size();i_s++){
      std::string var = detvar_str.at(i_s);
      h_v.push_back((TH1D*)f_in->Get(("Reco/Vars/"+var+"/h_Tot").c_str()));
      mchm.Restore(h_v.back());
      fill_colors.push_back(i_s+2);
      legs.push_back(var);
    }

    if(draw_hist){
      h_v.push_back((TH1D*)f_in_hist->Get("Reco/CV/h_Tot")->Clone("h_CV_Reco_Tot_Hist"));
      mchm.Restore(h_v.back());
      h_v.back()->SetLineStyle(2);
      fill_colors.push_back(1);
      legs.push_back("CV HIST");
    }

    pfs::DrawUnstacked(h_v,fill_colors,legs,draw_o,draw_u,false,true,plot_dir+"Reco.png"); 

    h_v.clear();

    // Draw the shape of the selected background
    h_v.push_back((TH1D*)f_in->Get("Reco/CV/h_AllBG"));
    mchm.Restore(h_v.back());
    h_v.back()->Scale(1.0/IntegralWithOU(h_v.back()));

    for(size_t i_s=0;i_s<detvar_str.size();i_s++){
      std::string var = detvar_str.at(i_s);
      h_v.push_back((TH1D*)f_in->Get(("Reco/Vars/"+var+"/h_AllBG").c_str()));
      mchm.Restore(h_v.back());
      h_v.back()->Scale(1.0/IntegralWithOU(h_v.back()));
    }

    if(draw_hist){
      h_v.push_back((TH1D*)f_in_hist->Get("Reco/CV/h_AllBG")->Clone("h_CV_Reco_AllBG_Hist"));
      mchm.Restore(h_v.back());
      h_v.back()->SetLineStyle(2);
      h_v.back()->Scale(1.0/IntegralWithOU(h_v.back()));
    }

    pfs::DrawUnstacked(h_v,fill_colors,legs,draw_o,draw_u,false,true,plot_dir+"AllBGShape.png");

    if(draw_truth){
      h_v.clear();
      h_v.push_back((TH1D*)f_in->Get("Truth/CV/h_Signal"));
      mchm.Restore(h_v.back(),"All",true);
      for(size_t i_s=0;i_s<detvar_str.size();i_s++){
        std::string var = detvar_str.at(i_s);
        h_v.push_back((TH1D*)f_in->Get(("Truth/Vars/"+var+"/h_Signal").c_str()));
        mchm.Restore(h_v.back(),"All",true);
      } 

      if(draw_hist){
        h_v.push_back((TH1D*)f_in_hist->Get("Truth/CV/h_Signal")->Clone("h_CV_Truth_Signal_Hist"));
        h_v.back()->SetLineStyle(2);
        mchm.Restore(h_v.back(),"All",true);
      }

      pfs::DrawUnstacked(h_v,fill_colors,legs,draw_o,draw_u,false,false,plot_dir+"Truth.png"); 

      // Draw the efficiency afo truth variable
      h_v.clear();

      TH2D* h_CV_Response_Signal = (TH2D*)f_in->Get("Response/CV/h_Signal");
      int bins = h_CV_Response_Signal->GetNbinsY()+1;
      h_v.push_back(h_CV_Response_Signal->ProjectionX("h_CV_Eff",0,bins)); 
      mchm.Restore(h_v.back(),"All",true);

      for(size_t i_s=0;i_s<detvar_str.size();i_s++){
        std::string var = detvar_str.at(i_s);
        TH2D* h_Var_Response_Signal = (TH2D*)f_in->Get(("Response/Vars/"+var+"/h_Signal").c_str());
        h_v.push_back(h_Var_Response_Signal->ProjectionX(("h_"+var+"_Eff").c_str(),0,bins));
        mchm.Restore(h_v.back(),"All",true);
      } 

      if(draw_hist){
        TH2D* h_CV_Hist_Response_Signal = (TH2D*)f_in_hist->Get("Response/CV/h_Signal");
        h_v.push_back(h_CV_Hist_Response_Signal->ProjectionX("h_CV_Hist_Eff",0,bins));
        mchm.Restore(h_v.back(),"All",true);
        h_v.back()->SetLineStyle(2);
      }

      pfs::DrawUnstacked(h_v,fill_colors,legs,draw_o,draw_u,false,false,plot_dir+"Efficiency.png"); 

    }


  }

  }
