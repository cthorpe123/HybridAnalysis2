#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "BranchList.h"
#include "Systematics.h"
#include "EnergyEstimatorFuncs.h"
#include "BinningFuncs.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"

using namespace syst;

// Strip a trailing "_<digits>" universe index off a histogram name, so
// e.g. "Pred_0" and "Truth_12" group as "Pred" and "Truth" respectively.
// A name with no such suffix (the unisim case) is returned unchanged.

std::string StripUniverseSuffix(const std::string& name){
  size_t pos = name.rfind('_');
  if(pos == std::string::npos) return name;
  std::string suffix = name.substr(pos+1);
  if(suffix.empty() || !std::all_of(suffix.begin(),suffix.end(),::isdigit)) return name;
  return name.substr(0,pos);
}

// Read the FFGenerators.root file written by MakeFoldingIngredients.C and,
// for every sub-subfolder of Vars/ (e.g. Vars/Genie/BG, Vars/Flux/NuWro...),
// draw the first 5 histograms it contains overlaid on one canvas. Folders
// that mix histogram types (e.g. Vars/Flux/<gen> holds both "Pred_i" and
// "Truth_i") are split by type, each getting its own canvas. Also draws,
// for each systematic, the response matrix in every universe (read from
// Histograms.root, since MakeFoldingIngredients.C doesn't copy those in).

void PlotFoldingIngredients(){

  std::vector<std::string> vars = {"MuonMom","MuonCosTheta","LeadProtonKE","ProtonKE"};
  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro","GiBUU"};

  const size_t n_max = 5;
  const std::vector<int> colors = {kBlack,kRed,kBlue,kGreen+2,kMagenta+1};
  const bool draw_o = false;
  const bool draw_u = false;

  for(const std::string& var : vars){
    std::cout << var << std::endl;

    std::string plot_dir = "Analysis/"+var+"/Plots/FoldingIngredients/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    TFile* f_in = TFile::Open(("Analysis/"+var+"/rootfiles/FFGenerators.root").c_str());
    if(!f_in || f_in->IsZombie()){
      std::cout << "Could not open FFGenerators.root for " << var << ", run MakeFoldingIngredients.C first" << std::endl;
      continue;
    }

    // Response matrices themselves aren't copied into FFGenerators.root by
    // MakeFoldingIngredients.C, so pull them from Histograms.root instead
    TFile* f_hist = TFile::Open(("Analysis/"+var+"/rootfiles/Histograms.root").c_str());
    if(!f_hist || f_hist->IsZombie())
      std::cout << "Could not open Histograms.root for " << var << ", won't be able to plot response matrices" << std::endl;

    // Response matrices are stored in the raw (template) binning, so need
    // to be run through RestoreRes to get the physical binning back, same
    // as MakeFoldingIngredients.C does for the regular histograms
    hist::MultiChannelHistogramManager mchm(var,true);
    mchm.LoadTemplates();

    // Draw the CV response matrix, also read from Histograms.root
    if(f_hist && !f_hist->IsZombie()){
      TH2D* h_res_cv = (TH2D*)f_hist->Get("Response/CV/h_Signal");
      if(h_res_cv){
        TH2D* h_res_cv_r = (TH2D*)h_res_cv->Clone("h_res_cv");
        mchm.RestoreRes(h_res_cv_r);
        std::string name = plot_dir+"CV_Response.png";
        pfs::Draw2DHist(h_res_cv_r,name);
        std::cout << "Wrote " << name << std::endl;
        delete h_res_cv_r;
      }
    }

    TDirectory* d_vars = f_in->GetDirectory("Vars");
    if(!d_vars){
      std::cout << "No Vars directory found in " << f_in->GetName() << std::endl;
      f_in->Close();
      continue;
    }

    TIter next_sys(d_vars->GetListOfKeys());
    TKey* k_sys;
    while((k_sys = (TKey*)next_sys())){

      std::string sys = k_sys->GetName();
      TDirectory* d_sys = d_vars->GetDirectory(sys.c_str());
      if(!d_sys) continue;

      std::string plot_dir_sys = plot_dir+sys+"/";
      gSystem->Exec(("mkdir -p " + plot_dir_sys).c_str());

      // Draw the response matrix for this systematic. It lives in
      // Histograms.root rather than FFGenerators.root, and doesn't depend
      // on subdir (BG/Data/generator), so only needs drawing once per sys:
      // one matrix per universe for a multisim, or a single one for a unisim.
      if(f_hist && !f_hist->IsZombie()){
        auto it_ms = std::find(sys_str.begin(),sys_str.end(),sys);
        if(it_ms != sys_str.end()){
          int i_s = std::distance(sys_str.begin(),it_ms);
          for(int i_u=0;i_u<std::min<int>(sys_nuniv.at(i_s),n_max);i_u++){
            TH2D* h_res = (TH2D*)f_hist->Get(("Response/Vars/"+sys+"/h_Signal_"+std::to_string(i_u)).c_str());
            if(!h_res) continue;
            TH2D* h_res_r = (TH2D*)h_res->Clone("h_res");
            mchm.RestoreRes(h_res_r);
            std::string name = plot_dir_sys+"Response_"+std::to_string(i_u)+".png";
            pfs::Draw2DHist(h_res_r,name);
            std::cout << "Wrote " << name << std::endl;
            delete h_res_r;
          }
        } else {
          TH2D* h_res = (TH2D*)f_hist->Get(("Response/Vars/"+sys+"/h_Signal").c_str());
          if(h_res){
            TH2D* h_res_r = (TH2D*)h_res->Clone("h_res");
            mchm.RestoreRes(h_res_r);
            std::string name = plot_dir_sys+"Response.png";
            pfs::Draw2DHist(h_res_r,name);
            std::cout << "Wrote " << name << std::endl;
            delete h_res_r;
          }
        }
      }

      TIter next_sub(d_sys->GetListOfKeys());
      TKey* k_sub;
      while((k_sub = (TKey*)next_sub())){

        std::string subdir = k_sub->GetName();
        TDirectory* d_sub = d_sys->GetDirectory(subdir.c_str());
        if(!d_sub) continue;

        // Group the histograms by their name with the universe index
        // stripped off, keeping the first n_max of each group, in the
        // order they were written.
        std::vector<std::string> group_order;
        std::map<std::string,std::vector<TH1D*>> h_m;
        std::map<std::string,std::vector<std::string>> legs_m;
        TIter next_h(d_sub->GetListOfKeys());
        TKey* k_h;
        while((k_h = (TKey*)next_h())){
          std::string hname = k_h->GetName();
          std::string prefix = StripUniverseSuffix(hname);
          std::vector<TH1D*>& h_v = h_m[prefix];
          if(h_v.size() >= n_max) continue;
          TH1D* h = dynamic_cast<TH1D*>(k_h->ReadObj());
          if(!h) continue;
          if(h_v.empty()) group_order.push_back(prefix);
          h_v.push_back(h);
          legs_m[prefix].push_back(hname);
        }

        for(const std::string& prefix : group_order){
          std::vector<TH1D*>& h_v = h_m[prefix];
          std::vector<std::string>& legs = legs_m[prefix];

          std::vector<int> cols;
          for(size_t i=0;i<h_v.size();i++) cols.push_back(colors.at(i % colors.size()));

          std::string name = group_order.size() > 1 ? plot_dir_sys+subdir+"_"+prefix+".png"
                                                     : plot_dir_sys+subdir+".png";
          pfs::DrawUnstacked(h_v,cols,legs,draw_o,draw_u,false,true,name);
          std::cout << "Wrote " << name << std::endl;

          for(TH1D* h : h_v) delete h;
        }

        // Draw the covariance matrix that goes with this Vars/<sys>/<subdir>
        // folder, stored under Cov/<sys>/<subdir> with a name that depends
        // on what subdir holds (a background, data, or a generator prediction)
        TDirectory* d_cov_sub = f_in->GetDirectory(("Cov/"+sys+"/"+subdir).c_str());
        if(d_cov_sub){
          std::string cov_key = subdir=="BG"      ? "Cov_BG"
                               : subdir=="Data"    ? "Cov_Data"
                               : subdir=="BGSData" ? "Cov_BGSData"
                                                   : "Cov_Pred";
          TH2D* h_cov = dynamic_cast<TH2D*>(d_cov_sub->Get(cov_key.c_str()));
          if(h_cov){
            std::string name = plot_dir_sys+subdir+"_Cov.png";
            pfs::Draw2DHist(h_cov,name);
            std::cout << "Wrote " << name << std::endl;
          }
        }
      }
    }

    // Also plot, for the CV, the difference between the data (and the
    // background-subtracted data) and each generator's folded prediction
    TDirectory* d_cv = f_in->GetDirectory("CV");
    if(d_cv){
      auto MakeDiffPlot = [&](const std::string& ref_name){
        TH1D* h_ref = dynamic_cast<TH1D*>(d_cv->Get(ref_name.c_str()));
        if(!h_ref) return;

        std::vector<TH1D*> h_diff_v;
        std::vector<std::string> legs;
        for(const std::string& gen : generators){
          TH1D* h_gen = dynamic_cast<TH1D*>(d_cv->Get(gen.c_str()));
          if(!h_gen) continue;
          TH1D* h_diff = (TH1D*)h_ref->Clone((ref_name+"_minus_"+gen).c_str());
          h_diff->Add(h_gen,-1);
          h_diff_v.push_back(h_diff);
          legs.push_back(gen);
        }
        if(h_diff_v.empty()) return;

        std::vector<int> cols;
        for(size_t i=0;i<h_diff_v.size();i++) cols.push_back(colors.at(i % colors.size()));

        std::string name = plot_dir+"CV_"+ref_name+"MinusGenerators.png";
        pfs::DrawUnstacked(h_diff_v,cols,legs,draw_o,draw_u,false,true,name);
        std::cout << "Wrote " << name << std::endl;

        for(TH1D* h : h_diff_v) delete h;
      };

      MakeDiffPlot("Data");
      MakeDiffPlot("BGSData");
    }

    f_in->Close();
    if(f_hist) f_hist->Close();
  }

}
