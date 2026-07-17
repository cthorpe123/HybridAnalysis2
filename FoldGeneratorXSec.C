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
// Do all uncertainty calculations in "cross section space"

void FoldGeneratorXSec(){

  bool load_asimov = true;
  std::vector<std::string> vars = {"Norm","Enu","MuonMom","MuonCosTheta"};
  //std::vector<std::string> vars = {"Enu"};

  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro"};
  //std::vector<std::string> generators = {"Untunedv3.0.6"};
  bool blinded = true;
  bool add_detvars = false;

  for(const std::string& var : vars){
    std::cout << var << std::endl;
    gSystem->Exec(("mkdir -p Analysis/"+var+"/Plots/FoldGeneratorXSec").c_str());

    hist::MultiChannelHistogramManager mchm(var);
    mchm.LoadTemplates();    

    std::vector<TH1D*> h_v;
    std::vector<std::string> legs;
    std::vector<int> cols;
    std::vector<std::string> chi2s;

    TFile* f_hist = TFile::Open(("Analysis/"+var+"/rootfiles/Histograms.root").c_str());
    const double POT = ((TH1D*)f_hist->Get("Meta/POT"))->GetBinContent(1);

    // Store each individual covariance matrix contributing to the total
    // uncertainty budget, for comparison later
    TH2D* h_cov_data = (TH2D*)f_hist->Get("Reco/Cov/Total/Cov_Tot")->Clone("h_cov_data");
    h_cov_data->Reset();
    std::map<std::string,TH2D*> h_cov_breakdown_v;

    // Load the reconstructed data with stat errors
    TH1D* h_reco_data = blinded ? (TH1D*)f_hist->Get("Reco/CV/h_Tot")->Clone("h_reco_data") : (TH1D*)f_hist->Get("Reco/CV/h_Data")->Clone("h_reco_data");

    // Subtract the predicted background from the data
    h_reco_data->Add((TH1D*)f_hist->Get("Reco/CV/h_AllBG"),-1);

    // Convert BG subtracted data to cross section space
    CrossSectionH(h_reco_data,POT);

    // Data just gets stat error for now, plus MC stat error from BG subtracted
    TH2D* h_cov_data_stat = blinded ? (TH2D*)f_hist->Get("Reco/Cov/EstDataStat/Cov_Tot") : (TH2D*)f_hist->Get("Reco/Cov/MCStat/Cov_Data");
    CrossSectionCov(h_cov_data_stat,POT);
    h_cov_breakdown_v["DataStat"] = h_cov_data_stat;   
    h_cov_data->Add(h_cov_data_stat); 

    TH2D* h_cov_bg_mc_stat = (TH2D*)f_hist->Get("Reco/Cov/MCStat/Cov_AllBG");
    CrossSectionCov(h_cov_bg_mc_stat,POT);
    h_cov_breakdown_v["BGMCStat"] = h_cov_bg_mc_stat;
    h_cov_data->Add(h_cov_bg_mc_stat); 

    // Add errors to the background subtracted data
    for(int i=0;i<h_reco_data->GetNbinsX()+2;i++) 
      h_reco_data->SetBinError(i,sqrt(h_cov_data->GetBinContent(i,i)));

    // Now take the generator prediction and fold through the CV detector resposne,
    // and create the total covariance matrix for each
    TH2D* h_res_cv = (TH2D*)f_hist->Get("Response/CV/h_Signal");
    TFile* f_gen = TFile::Open(("Analysis/"+var+"/rootfiles/GeneratorXSec.root").c_str());
    std::vector<TH1D*> h_gen_ff_cv_v;
    std::vector<TH2D*> h_gen_cov_v;
    for(size_t i=0;i<generators.size();i++){
      std::string gen = generators.at(i);
      h_gen_ff_cv_v.push_back(Multiply((TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str()),h_res_cv,"h_xsec_ff_"+gen)); 
      h_gen_cov_v.push_back((TH2D*)h_cov_data->Clone(("h_gen_cov_"+gen).c_str()));
      h_gen_cov_v.back()->Reset();
      legs.push_back(gen);
      cols.push_back(i+2);
    }

    // Multiply the generator predictions by the response in each genie, g4 and detvar universe
    for(size_t i_g=0;i_g<generators.size();i_g++){

      std::string gen = generators.at(i_g);
      std::cout << gen << std::endl;

      // load the generator prediction in truth space
      const TH1D* h_gen_truth = (TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str());

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
        h_gen_cov_v.at(i_g)->Add(c);
        h_cov_breakdown_v[sys_str.at(i_s)] = (TH2D*)c->Clone(("h_cov_"+sys_str.at(i_s)).c_str());
        h_cov_breakdown_v.at(sys_str.at(i_s))->SetDirectory(0);
        for(TH1D* hh : h) delete hh;
        h.clear();
        delete c;
        delete fc;
      }
      
      // For the unisims, we need to compare the unisim universe prediction, 
      // including the background, to the forward folded generator prediction,
      // therefore we need the FF generator prediction + BG 

      // CV background in xsec space
      TH1D* h_cv_bg = (TH1D*)f_hist->Get("Reco/CV/h_AllBG");
      CrossSectionH(h_cv_bg,POT);

      TH1D* h_cv = (TH1D*)h_gen_ff_cv_v.at(i_g)->Clone("h_cv");
      ForceAddTH1D(h_cv,h_cv_bg);

      // Unisims
      for(int i_s=0;i_s<kUnisimMAX;i_s++){

        std::string name = "h_"+gen+"_ff_"+unisims_str.at(i_s);
      
        // FF generator pred in this universe
        TH1D* h = Multiply(h_gen_truth,(TH2D*)f_hist->Get(("Response/Vars/"+unisims_str.at(i_s)+"/h_Signal").c_str()),name.c_str());
       
        // Selected BG in this universe
        TH1D* h_bg = (TH1D*)f_hist->Get(("Reco/Vars/"+unisims_str.at(i_s)+"/h_AllBG").c_str());
        CrossSectionH(h_bg,POT);

        // FF generator pred in this universe + background
        ForceAddTH1D(h,h_bg);

        TH2D *c,*fc;
        CalcCovUnisim(unisims_str.at(i_s),h_cv,h,c,fc);
        h_gen_cov_v.at(i_g)->Add(c);
        h_cov_breakdown_v[unisims_str.at(i_s)] = (TH2D*)c->Clone(("h_cov_"+unisims_str.at(i_s)).c_str());
        h_cov_breakdown_v.at(unisims_str.at(i_s))->SetDirectory(0);
        delete h;
        delete c;
        delete fc;
      }

      delete h_cv; 

      
      // Flux - use the 2D distribution and multiply by the flux ratios, then multiply by the response
      TFile* f_flux_ratios = TFile::Open("../Flux/FluxRatios.root");
      const TH2D* h_gen_truth_2d = (TH2D*)f_gen->Get(("h_xsec_2D_"+var+"_"+gen).c_str());
      const TH1D* h_integrals = (TH1D*)f_flux_ratios->Get("IntegralRatios/h_NuMu_IntegralRatio");

      // Calculate FF generator prediction in each flux universe + bg in each flux universe
      // Generator truth is built by taking Generator CV in 2d (flux and var) and multiplying
      // it by hist holding flux syst/flux cv  
      std::vector<TH1D*> h_gen_ff_flux_v;
      for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++){
        
        double numu = h_integrals->GetBinContent(i_u+1); 
        TH1D* h_flux_ratio = (TH1D*)f_flux_ratios->Get(Form("ShapeRatios/NuMu/h_NuMu_FluxRatio_%i",i_u)); // Ratio of alt flux to CV in true nu_e
        TH1D* h_gen_truth_flux = Multiply(h_flux_ratio,h_gen_truth_2d,Form("h_gen_truth_f%i",i_u)); // Gen truth in flux universe
        h_gen_ff_flux_v.push_back(Multiply(h_gen_truth_flux,h_res_cv,Form("h_gen_f%i",i_u))); // Gen in reco space in flux universe

        // Calculate the cross section in the data using the flux in this universe
        TH1D* h_data = blinded ? (TH1D*)f_hist->Get("Reco/CV/h_Tot")->Clone("h_data") : (TH1D*)f_hist->Get("Reco/CV/h_Data")->Clone("h_data");
        TH1D* h_bg = (TH1D*)f_hist->Get(("Reco/Vars/"+sys_str.at(kFlux)+"/h_AllBG_"+std::to_string(i_u)).c_str());
        h_data->Add(h_bg,-1.0);

        CrossSectionH(h_data,POT*numu);           
        h_data->Scale(-1.0);

        ForceAddTH1D(h_gen_ff_flux_v.back(),h_data);

        delete h_gen_truth_flux;
        delete h_data;

      }
      
      TH2D *c,*fc;
      CalcCovMultisim(sys_str.at(kFlux),h_gen_ff_flux_v,c,fc);
      h_gen_cov_v.at(i_g)->Add(c);
      pfs::Draw2DHist(c,"Analysis/"+var+"/Plots/FoldGeneratorXSec/Cov_Flux_"+gen+".png");
      h_cov_breakdown_v[sys_str.at(kFlux)] = (TH2D*)c->Clone(("h_cov_"+sys_str.at(kFlux)).c_str());
      h_cov_breakdown_v.at(sys_str.at(kFlux))->SetDirectory(0);

      for(TH1D* h : h_gen_ff_flux_v) delete h;
      h_gen_ff_flux_v.clear();
      delete c; 
      delete fc;
      f_flux_ratios->Close();

      for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++)
        h_gen_ff_cv_v.at(i_g)->SetBinError(i_b,sqrt(h_gen_cov_v.at(i_g)->GetBinContent(i_b,i_b)));

      // For the chi2, include the data statistics
      TH2D* h_cov = (TH2D*)h_gen_cov_v.at(i_g)->Clone("h_cov"); 
      h_cov->Add(h_cov_data);

      // Calculate the chi2 of each generator with the data
      std::pair<double,int> chi2 = Chi2(h_gen_ff_cv_v.at(i_g),h_reco_data,h_cov,false,false);
      std::cout << var << " " << gen << ", chi2/ndof = " << chi2.first << "/" << chi2.second << " = " << chi2.first/chi2.second << std::endl;
      chi2s.push_back(to_string_with_precision(chi2.first/chi2.second,2));
      legs.at(i_g) += ", #chi^{2}/n = " + chi2s.at(i_g);

      // Draw breakdown of systematics budget
      std::vector<TH1D*> h_sys;
      std::vector<std::string> sys_legs;
      std::vector<int> sys_cols;
      int i_s=0;
      for(auto const& item : h_cov_breakdown_v){
        TH1D* h = (TH1D*)h_reco_data->Clone(("h_sys_"+item.first).c_str()); 
        for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++) h->SetBinContent(i_b,sqrt(item.second->GetBinContent(i_b,i_b))/h_gen_ff_cv_v.at(i_g)->GetBinContent(i_b));
        mchm.Restore(h);
        h_sys.push_back(h);
        sys_legs.push_back(item.first);
        sys_cols.push_back(i_s+2);
        i_s++;
      } 
      
      TH1D* h = (TH1D*)h_reco_data->Clone("h_sys_tot"); 
      for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++) h->SetBinContent(i_b,sqrt(h_cov->GetBinContent(i_b,i_b))/h_gen_ff_cv_v.at(i_g)->GetBinContent(i_b));
      mchm.Restore(h);
      h_sys.push_back(h); 
      sys_legs.push_back("Total");
      sys_cols.push_back(1);

      /*
      for(size_t i_s=0;i_s<h_sys.size();i_s++){
        std::cout << sys_legs.at(i_s) << std::endl;
        for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++) std::cout << h_sys.at(i_s)->GetBinContent(i_b) << std::endl;
      }
      */

      pfs::DrawUnstacked(h_sys,sys_cols,sys_legs,true,true,false,"Analysis/"+var+"/Plots/FoldGeneratorXSec/SysBreakdown_"+gen+".png");

      delete h_cov;
      for(auto const& item : h_cov_breakdown_v){
        if(item.first == "DataStat" || item.first == "BGMCStat") continue;
        delete item.second;
        h_cov_breakdown_v.erase(item.first);
      }

      std::cout << "Finished " << gen << std::endl;

    }
    

    std::string name = "Analysis/"+var+"/Plots/FoldGeneratorXSec/Test";
   
    h_gen_ff_cv_v.push_back(h_reco_data);
    cols.push_back(1);
    legs.push_back(!blinded ? "uboone Data" : "Asimov Data");
 
    for(TH1D*& h : h_gen_ff_cv_v) mchm.Restore(h);
      
    pfs::DrawUnstacked(h_gen_ff_cv_v,cols,legs,true,true,true,name+".png");

    // Clean up
    std::cout << "Cleaning up" << std::endl;
    //delete h_reco_data;
    //for(TH1D* h : h_gen_ff_cv_v) delete h;

    f_hist->Close();
    f_gen->Close();

    std::cout << "Finished cleaning up" << std::endl;
  }

}
