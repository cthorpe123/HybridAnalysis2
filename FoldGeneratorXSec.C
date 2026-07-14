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

// Fold generator cross sections through the detector response.
// Do the "correct way" - assume generator has access to full stack of 
// response matrices in every systematic universe

void FoldGeneratorXSec(){

  bool load_asimov = true;
  std::vector<std::string> vars = {"Enu","MuonMom","MuonCosTheta","Norm"};
  std::vector<std::string> generators = {"v3.0.6","NuWro"};
  bool blinded = true;
  bool add_detvars = false;

  for(const std::string& var : vars){
    gSystem->Exec(("mkdir -p Analysis/"+var+"/Plots/FoldGeneratorXSec").c_str());

    std::vector<TH1D*> h_v;
    std::vector<std::string> legs;
    std::vector<int> cols;

    TFile* f_hist = TFile::Open(("Analysis/"+var+"/rootfiles/Histograms.root").c_str());
    const double POT = ((TH1D*)f_hist->Get("Meta/POT"))->GetBinContent(1);
    
    // Load the reconstructed data with stat errors
    TH1D* h_reco_data = blinded ? (TH1D*)f_hist->Get("Reco/CV/h_Tot") : (TH1D*)f_hist->Get("Reco/CV/h_Data");

    // Subtract the predicted background from the data
    TH1D* h_reco_bg = (TH1D*)f_hist->Get("Reco/CV/h_AllBG");
    TH1D* h_bg_sub_data = (TH1D*)h_reco_data->Clone("h_bg_sub_data");
    h_bg_sub_data->Add(h_reco_bg,-1);

    // Data just gets stat error for now, plus MC stat error from BG subtracted
    TH2D* h_cov_data_stat = blinded ? (TH2D*)f_hist->Get("Reco/Cov/EstDataStat/Cov_Tot") : (TH2D*)f_hist->Get("Reco/Cov/MCStat/Cov_Data");
    h_cov_data_stat->Add((TH1D*)f_hist->Get("Reco/Cov/MCStat/Cov_AllBG"));

    // Add errors to the background subtracted data
    for(int i=0;i<h_bg_sub_data->GetNbinsX()+2;i++) 
      h_bg_sub_data->SetBinError(i,sqrt(h_cov_data_stat->GetBinContent(i,i)));

    // Convert to cross section
    CrossSectionH(h_bg_sub_data,POT);

    // Now take the generator prediction and fold through the CV detector resposne
    TH2D* h_res_cv = (TH2D*)f_hist->Get("Response/CV/h_Signal");
    TFile* f_gen = TFile::Open(("Analysis/"+var+"/rootfiles/GeneratorXSec.root").c_str());
    std::vector<TH1D*> h_gen_ff_cv_v;
    for(size_t i=0;i<generators.size();i++){
      std::string gen = generators.at(i);
      h_gen_ff_cv_v.push_back(Multiply((TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str()),h_res_cv,"h_xsec_ff_"+gen)); 
      legs.push_back(gen);
      cols.push_back(i+2);
    }

    // Plot the FF gen predictions along with BG subtracted data
    std::vector<TH1D*> h_ff_tmp_v = h_gen_ff_cv_v;
    h_ff_tmp_v.push_back(h_bg_sub_data);
    cols.push_back(1);
    legs.push_back("BG Subtracted Data");
    pfs::DrawUnstacked(h_ff_tmp_v,cols,legs,true,true,false,"Analysis/"+var+"/Plots/FoldGeneratorXSec/CV.png");

    // Multiply the generator predictions by the response in each genie, g4 and detvar universe
    std::vector<TH2D*> h_gen_cov_v;
    for(size_t i=0;i<generators.size();i++){

      std::string gen = generators.at(i);
      const TH1D* h_gen_truth = (TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str());
      h_gen_cov_v.push_back((TH2D*)h_cov_data_stat->Clone(("h_gen_cov_"+gen).c_str()));
      h_gen_cov_v.back()->Reset();

      // Multisims
      for(int i_s=0;i_s<kSystMAX;i_s++){
        if(i_s == kFlux) continue;
        std::vector<TH1D*> h;
        for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
          std::string name = "h_"+gen+"_ff_"+sys_str.at(i_s)+"_"+std::to_string(i_u);
          h.push_back(Multiply(h_gen_truth,(TH2D*)f_hist->Get(("Response/Vars/"+sys_str.at(i_s)+"/h_Signal_"+std::to_string(i_u)).c_str()),name.c_str()));
          TH1D* h_bg = (TH1D*)f_hist->Get(("Reco/Vars/"+sys_str.at(i_s)+"/h_AllBG_"+std::to_string(i_u)).c_str());
          CrossSectionH(h_bg,POT);
          ForceAddTH1D(h.back(),h_bg);
        }
        TH2D *c,*fc; 
        CalcCovMultisim(sys_str.at(i_s),h,c,fc); 
        h_gen_cov_v.back()->Add(c);
        for(TH1D* hh : h) delete hh;
        h.clear();
        delete c;
        delete fc;
      }
       
      // FF generator in CV universe, plus CV background, need this for unisims
      TH1D* h_cv = (TH1D*)h_gen_ff_cv_v.at(i)->Clone("h_cv");
      TH1D* h_cv_bg = (TH1D*)f_hist->Get("Reco/CV/h_AllBG");
      CrossSectionH(h_cv_bg,POT);
      h_cv->Add(h_cv_bg);

      // Unisims
      for(int i_s=0;i_s<kUnisimMAX;i_s++){
        std::cout << "Unisim: " << unisims_str.at(i_s) << std::endl;
        std::string name = "h_"+gen+"_ff_"+unisims_str.at(i_s);
      
        // FF generator pred in this universe
        TH1D* h = Multiply(h_gen_truth,(TH2D*)f_hist->Get(("Response/Vars/"+unisims_str.at(i_s)+"/h_Signal").c_str()),name.c_str());
       
        // Selected BG in this universe
        TH1D* h_bg = (TH1D*)f_hist->Get(("Reco/Vars/"+unisims_str.at(i_s)+"/h_AllBG").c_str());
        CrossSectionH(h_bg,POT);
        ForceAddTH1D(h,h_bg);

        TH2D *c,*fc; 
        CalcCovUnisim(unisims_str.at(i_s),h_cv,h,c,fc); 
        h_gen_cov_v.back()->Add(c);
        delete h;
        delete c;
        delete fc;
      }

      // Flux - use the 2D distribution and multiply by the flux ratios, then multiply by the response
      TFile* f_flux_ratios = TFile::Open("../Flux/FluxRatios.root");
      const TH2D* h_gen_truth_2d = (TH2D*)f_gen->Get(("h_xsec_2D_"+var+"_"+gen).c_str());
      std::vector<TH1D*> h_gen_ff_flux_v;
      for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++){
        h_gen_ff_flux_v.push_back(Multiply(Multiply((TH1D*)f_flux_ratios->Get(Form("ShapeRatios/NuMu/h_NuMu_FluxRatio_%i",i_u)),h_gen_truth_2d,Form("h_gen_truth_f%i",i_u)),h_res_cv,Form("h_gen_f%i",i_u)));
      }

      TH2D *c,*fc; 
      CalcCovMultisim(sys_str.at(kFlux),h_gen_ff_flux_v,c,fc); 
      pfs::Draw2DHist(fc,"Analysis/"+var+"/Plots/FoldGeneratorXSec/Flux_FCov_"+gen+".png");
      pfs::Draw2DHist(c,"Analysis/"+var+"/Plots/FoldGeneratorXSec/Flux_Cov_"+gen+".png");
      pfs::Draw2DHist(h_gen_cov_v.at(i),"Analysis/"+var+"/Plots/FoldGeneratorXSec/Cov_"+gen+".png");
      h_gen_cov_v.back()->Add(c);     
      for(TH1D* h : h_gen_ff_flux_v) delete h;
      h_gen_ff_flux_v.clear(); 

      for(int i_b=0;i_b<h_gen_ff_cv_v.at(i)->GetNbinsX()+2;i_b++)
        h_gen_ff_cv_v.at(i)->SetBinError(i_b,sqrt(h_gen_cov_v.at(i)->GetBinContent(i_b,i_b)));

      // Calculate the chi2 of each generator with the data
      std::pair<double,int> chi2 = Chi2(h_gen_ff_cv_v.at(i),h_bg_sub_data,h_gen_cov_v.at(i),true,true);
      std::cout << var << " " << gen << ", chi2/ndof = " << chi2.first << "/" << chi2.second << " = " << chi2.first/chi2.second << std::endl;

    }

    pfs::DrawUnstacked(h_ff_tmp_v,cols,legs,true,true,true,"Analysis/"+var+"/Plots/FoldGeneratorXSec/Test.png");

    // Clean up
    delete h_bg_sub_data;
    for(TH1D* h : h_gen_ff_cv_v) delete h;

    f_hist->Close();
    f_gen->Close();

  }

}
