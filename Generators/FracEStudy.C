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
using namespace binning;

void FracEStudy(){

  bool load_asimov = true;
  //std::vector<std::string> vars = {"Norm","Enu","MuonMom","MuonCosTheta"};
  std::vector<std::string> vars = var_names;
  vars.push_back("Enu");
  vars.push_back("Norm");
  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro","GiBUU"};
  
  bool blinded = true;
  bool restore = true;
  bool dbbw = true;
  bool add_detvars = false;
  bool draw_o = false;
  bool draw_u = false;

  for(const std::string& var : vars){
    
    std::cout << var << std::endl;

    std::string plot_dir = "Analysis/"+var+"/Plots/FracEStudy/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    hist::MultiChannelHistogramManager mchm(var);
    mchm.LoadTemplates();    

    std::vector<std::string> legs;
    std::vector<int> cols;
    std::vector<std::string> chi2s;

    TFile* f_hist = TFile::Open(("Analysis/"+var+"/rootfiles/Histograms.root").c_str());
    TFile* f_gen = TFile::Open(("Analysis/"+var+"/rootfiles/GeneratorXSec.root").c_str());
    const double POT = ((TH1D*)f_hist->Get("Meta/POT"))->GetBinContent(1);

    TH2D* h_res_cv = (TH2D*)f_hist->Get("Response/CV/h_Signal");
    TH2D* h_cov_template = (TH2D*)f_hist->Get("Reco/Cov/Total/Cov_Tot");
    h_cov_template->Reset();

    std::vector<TH1D*> h_gen_ff_cv_v;

    for(size_t i_g=0;i_g<generators.size();i_g++){
      std::string gen = generators.at(i_g);
      h_gen_ff_cv_v.push_back(Multiply((TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str()),h_res_cv,"h_xsec_ff_"+gen)); 
      legs.push_back(gen);
      cols.push_back(i_g+2);
    }

    std::map<std::string,std::vector<TH2D*>> h_fcov_m;
    for(int i_s=0;i_s<kSystMAX;i_s++) h_fcov_m[sys_str.at(i_s)] = std::vector<TH2D*>();
    for(int i_s=0;i_s<kUnisimMAX;i_s++) h_fcov_m[unisims_str.at(i_s)] = std::vector<TH2D*>();
    h_fcov_m["Unisims"] = std::vector<TH2D*>();

    // Multiply the generator predictions by the response in each genie, g4 and detvar universe
    for(size_t i_g=0;i_g<generators.size();i_g++){

      std::string gen = generators.at(i_g);
      std::cout << gen << std::endl;

      // load the generator prediction in truth space
      const TH1D* h_gen_truth = (TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str());
      
      // Multisims
      for(int i_s=0;i_s<kSystMAX;i_s++){
        std::vector<TH1D*> h;
        for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
          std::string name = "h_"+gen+"_ff_"+sys_str.at(i_s)+"_"+std::to_string(i_u);
          h.push_back(Multiply(h_gen_truth,(TH2D*)f_hist->Get(("Response/Vars/"+sys_str.at(i_s)+"/h_Signal_"+std::to_string(i_u)).c_str()),name.c_str()));
        }
        TH2D *c,*fc;
        CalcCovMultisim(sys_str.at(i_s),h,c,fc);
        h_fcov_m.at(sys_str.at(i_s)).push_back((TH2D*)fc->Clone(("h_cov_"+sys_str.at(i_s)+"_"+gen).c_str()));
        h_fcov_m.at(sys_str.at(i_s)).back()->SetDirectory(0);
        for(TH1D*& hh : h) mchm.Restore(hh);
        pfs::DrawUnstacked(h,std::vector<int>(sys_nuniv.at(i_s),2),std::vector<std::string>(sys_nuniv.at(i_s),""),draw_o,draw_u,false,dbbw,plot_dir+"Univ_"+sys_str.at(i_s)+"_"+gen+".png");
        for(TH1D* hh : h) delete hh;
        h.clear();
        mchm.Restore(c);
        mchm.Restore(fc);
        pfs::Draw2DHist(c,plot_dir+"Cov_"+sys_str.at(i_s)+"_"+gen+".png");
        pfs::Draw2DHist(fc,plot_dir+"FCov_"+sys_str.at(i_s)+"_"+gen+".png");
        delete c;
        delete fc;
      }

      
      // Unisims - group together into one systematic
      TH1D* h_cv = (TH1D*)h_gen_ff_cv_v.at(i_g)->Clone("h_cv");
      
      TH2D* h_cov_unisim = (TH2D*)h_cov_template->Clone("h_cov_unisim");
      TH2D* h_fcov_unisim = (TH2D*)h_cov_template->Clone("h_fcov_unisim");
      h_cov_unisim->Reset();
      h_fcov_unisim->Reset();

      for(int i_s=0;i_s<kUnisimMAX;i_s++){

        std::string name = "h_"+gen+"_ff_"+unisims_str.at(i_s);
      
        // FF generator pred in this universe
        TH1D* h = Multiply(h_gen_truth,(TH2D*)f_hist->Get(("Response/Vars/"+unisims_str.at(i_s)+"/h_Signal").c_str()),name.c_str());

        TH2D *c,*fc;
        CalcCovUnisim(unisims_str.at(i_s),h_cv,h,c,fc);
        h_cov_unisim->Add(c);
        h_fcov_unisim->Add(fc);
        h_fcov_m.at(unisims_str.at(i_s)).push_back((TH2D*)fc->Clone(("h_fcov_"+unisims_str.at(i_s)+"_"+gen).c_str()));
        h_fcov_m.at(unisims_str.at(i_s)).back()->SetDirectory(0);
        delete h;
        delete c;
        delete fc;
      }

      delete h_cv; 
      h_fcov_m.at("Unisims").push_back((TH2D*)h_fcov_unisim->Clone(("h_fcov_Unisims_"+gen).c_str()));
      h_fcov_m.at("Unisims").back()->SetDirectory(0);
      delete h_cov_unisim;
      delete h_fcov_unisim;
      
    }

    for(auto& item : h_fcov_m){
      std::vector<TH1D*> h_fe_v;
      for(TH2D* h : item.second){
        h_fe_v.push_back((TH1D*)h_gen_ff_cv_v.back()->Clone("fe"));
        for(int i_b=0;i_b<h_gen_ff_cv_v.back()->GetNbinsX()+2;i_b++)
          h_fe_v.back()->SetBinContent(i_b,sqrt(h->GetBinContent(i_b,i_b)));
        mchm.Restore(h_fe_v.back());
        h_fe_v.back()->GetYaxis()->SetTitle("Frac. Unc.");
      }
      pfs::DrawUnstacked(h_fe_v,cols,legs,draw_o,draw_u,false,false,plot_dir+"Breakdown_"+item.first+".png");
    }

    f_hist->Close();
    f_gen->Close();

  }

}
