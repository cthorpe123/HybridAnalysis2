#include "Funcs.h"
#include "BranchList.h"
#include "EnergyEstimatorFuncs.h"
#include "Systematics.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"
#include "WeightFuncs.h"

using namespace syst;

// Try forward folding the CV truth through the response
// calulated in the CV and special universes, calculate chi2
// between the CV and each special prediction

void FFTest_Truth(){

  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");

  bool add_detvars = true;
  const bool include_data_stat = true;
  const bool draw_underflow = false;
  const bool draw_overflow = false;
  const bool dbbw = true;
  const bool draw_chi2_curve = true;
  const bool diag_only = false;
  const bool draw_cov = false;

  std::vector<std::string> vars = {"MuonMom"};
  //std::vector<std::string> vars = var_names;
  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  weight::SetWeightFuncs();
  std::vector<std::string> special_univs;
  //special_univs.push_back("MuonMomShape_RH_Bias");
  for(const auto &item : weight::r_m)
    special_univs.push_back(item.first);

  for(size_t i_f=0;i_f<vars.size();i_f++){

    std::string label = vars.at(i_f);
    std::cout << label << std::endl;

    std::string plot_dir = "Analysis/"+label+"/Plots/FFTest/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    TFile* f_in = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
    TFile* f_in_detvar = add_detvars ? TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str()) : nullptr;

    // Load binning templates
    hist::MultiChannelHistogramManager mchm(label,true);
    mchm.SetTrueChannelList(channels_t);
    mchm.SetRecoChannelList(channels_r);
    mchm.LoadTemplates();
 
    TH1D* h_CV_Truth = (TH1D*)f_in->Get("Truth/CV/h_Signal");

    for(std::string s : special_univs){

      gSystem->Exec(("mkdir -p "+plot_dir+"/"+s).c_str());
      std::vector<std::pair<double,int>> spec_chi2;

      for(int i=0;i<weight::spline_pts;i++){

        std::string spec = s + "_" + std::to_string(i); 
        std::cout << spec << std::endl;

        TH1D* h_Spec_Truth = (TH1D*)f_in->Get(("Truth/Special/"+spec+"/h_Signal").c_str());

        // Draw the truth distribution in the special universe, useful for showing how the spline reweighting is modifying the distribution
        std::vector<TH1D*> h_Truth_v;
        std::vector<int> colors_ch;
        std::vector<std::string> legs_ch;
        for(size_t i_ch=0;i_ch<channels_t.size();i_ch++){ 
          std::string ch = channels_t.at(i_ch);
          h_Truth_v.push_back((TH1D*)h_Spec_Truth->Clone(("h_Spec_Truth_"+ch).c_str()));
          mchm.Restore(h_Truth_v.back(),ch,true);
          if(dbbw) DivideByBinWidth(h_Truth_v.back());
          h_Truth_v.back()->SetLineStyle(2);
          colors_ch.push_back(i_ch+1);
          legs_ch.push_back(channels_t.at(i_ch)+" "+spec);

          h_Truth_v.push_back((TH1D*)h_CV_Truth->Clone(("h_CV_Truth_"+ch).c_str()));
          mchm.Restore(h_Truth_v.back(),ch,true);
          if(dbbw) DivideByBinWidth(h_Truth_v.back());
          colors_ch.push_back(i_ch+1);
          legs_ch.push_back(channels_t.at(i_ch)+" CV");

        }
        pfs::DrawUnstacked2(h_Truth_v,colors_ch,legs_ch,plot_dir+"/"+s+"/"+spec+"_Truth.png",false);
        for(TH1D* hh : h_Truth_v) delete hh;
        legs_ch.clear();
        colors_ch.clear();

      }

    }
    
    f_in->Close();

  }

}

