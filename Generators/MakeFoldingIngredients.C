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

// Convert everything to "cross section space", restore binnings
// and generate covariance matrices needed for studies on different
// data release recipes

void MakeFoldingIngredients(){

  bool load_asimov = true;
  std::vector<std::string> vars = {"Norm","Enu","MuonMom","MuonCosTheta"};
  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro","GiBUU"};

  bool blinded = true;
  bool add_detvars = false;

  for(const std::string& var : vars){
    std::cout << var << std::endl;

    std::string plot_dir = "Analysis/"+var+"/Plots/FoldGeneratorXSec/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    hist::MultiChannelHistogramManager mchm(var);
    mchm.LoadTemplates();    

    std::vector<TH1D*> h_v;
    std::vector<std::string> legs;
    std::vector<int> cols;
    std::vector<std::string> chi2s;


    // Open the file containing the histograms
    TFile* f_hist = TFile::Open(("Analysis/"+var+"/rootfiles/Histograms.root").c_str());
    const double POT = ((TH1D*)f_hist->Get("Meta/POT"))->GetBinContent(1);

    // Open the file containing the generator predictions
    TFile* f_gen = TFile::Open(("Analysis/"+var+"/rootfiles/GeneratorXSec.root").c_str());

    // Open a file to write the outputs
    TFile* f_out = new TFile(("Analysis/"+var+"/rootfiles/FFGenerators.root").c_str(),"RECREATE");

    // Open the file containing the flux variations
    TFile* f_flux_ratios = TFile::Open("../Flux/FluxRatios.root");
    //const TH2D* h_gen_truth_2d = (TH2D*)f_gen->Get(("h_xsec_2D_"+var+"_"+gen).c_str());
    const TH1D* h_numu_integrals = (TH1D*)f_flux_ratios->Get("IntegralRatios/h_NuMu_IntegralRatio");
    const TH1D* h_numubar_integrals = (TH1D*)f_flux_ratios->Get("IntegralRatios/h_NuMuBar_IntegralRatio");

    // Fold the generators through the CV response and write
    f_out->mkdir("CV");
    f_out->cd("CV");

    TH2D* h_res_cv = (TH2D*)f_hist->Get("Response/CV/h_Signal");
    std::vector<TH1D*> h_gen_ff_cv_v;
    for(std::string gen : generators){
      h_gen_ff_cv_v.push_back(Multiply((TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str()),h_res_cv,"h_xsec_ff_"+gen)); 
      mchm.Restore(h_gen_ff_cv_v.back());
      h_gen_ff_cv_v.back()->Write(gen.c_str()); 
    }

    // Also write the data and background
    TH1D* h_bg_cv = (TH1D*)f_hist->Get("Reco/CV/h_AllBG");
    CrossSectionH(h_bg_cv,POT);
    mchm.Restore(h_bg_cv);
    h_bg_cv->Write("BG");

    TH1D* h_reco_data = blinded ? (TH1D*)f_hist->Get("Reco/CV/h_Tot")->Clone("h_reco_data") 
                                : (TH1D*)f_hist->Get("Reco/CV/h_Data")->Clone("h_reco_data");

    CrossSectionH(h_reco_data,POT);
    mchm.Restore(h_reco_data);
    h_reco_data->Write("Data");

    TH1D* h_bgs_data = (TH1D*)h_reco_data->Clone("BGSData");
    h_bgs_data->Add(h_bg_cv,-1);
    h_bgs_data->Write("BGSData");

    f_out->cd();
  
    // Calculate and store the data stat and BG MC stat errors 
    f_out->mkdir("Cov");
    f_out->mkdir("Vars");

    f_out->mkdir("Cov/DataStat");
    f_out->cd("Cov/DataStat"); 
    
    TH2D* h_cov_data_stat = blinded ? (TH2D*)f_hist->Get("Reco/Cov/EstDataStat/Cov_Tot") 
                                    : (TH2D*)f_hist->Get("Reco/Cov/MCStat/Cov_Data");
    CrossSectionCov(h_cov_data_stat,POT);
    mchm.Restore(h_cov_data_stat); 
    h_cov_data_stat->Write("h_Cov");

    f_out->cd();
    f_out->mkdir("Cov/BGMCStat");
    f_out->cd("Cov/BGMCStat");    

    TH2D* h_cov_bg_mc_stat = (TH2D*)f_hist->Get("Reco/Cov/MCStat/Cov_AllBG");
    CrossSectionCov(h_cov_bg_mc_stat,POT);
    mchm.Restore(h_cov_bg_mc_stat);
    h_cov_bg_mc_stat->Write("h_Cov");

    f_out->cd();

    // Calculate background predictions and various non-generator specific objects 

    std::map<std::string,std::vector<TH1D*>> h_bg_m;

    for(int i_s=0;i_s<kSystMAX;i_s++){
      std::string sys = sys_str.at(i_s);
      f_out->cd();
      f_out->mkdir(("Vars/"+sys+"/BG").c_str());
      f_out->cd(("Vars/"+sys+"/BG").c_str());
      std::vector<TH1D*> h_bg_v; 
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        h_bg_v.push_back((TH1D*)f_hist->Get(("Reco/Vars/"+sys+"/h_AllBG_"+std::to_string(i_u)).c_str()));
        CrossSectionH(h_bg_v.back(),POT);
        mchm.Restore(h_bg_v.back());
        h_bg_v.back()->Write(("BG_"+std::to_string(i_u)).c_str());

      }
      h_bg_m[sys] = h_bg_v;

      f_out->cd();
      f_out->mkdir(("Cov/"+sys+"/BG").c_str());
      f_out->cd(("Cov/"+sys+"/BG").c_str());
      TH2D *c,*fc;
      CalcCovMultisim(sys,h_bg_v,c,fc);
      mchm.Restore(c);
      c->Write("Cov_BG");

      f_out->cd();
      f_out->mkdir(("Vars/"+sys+"/Data").c_str());
      f_out->cd(("Vars/"+sys+"/Data").c_str());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++)
        h_reco_data->Write(("Data_"+std::to_string(i_u)).c_str());

      f_out->cd();
      f_out->mkdir(("Cov/"+sys+"/Data").c_str());
      f_out->cd(("Cov/"+sys+"/Data").c_str());
      TH2D* h_cov_data_tmp = (TH2D*)c->Clone("Cov_Data");
      h_cov_data_tmp->Reset();
      h_cov_data_tmp->Write("Cov_Data");

      f_out->cd();
      f_out->mkdir(("Vars/"+sys+"/BGSData").c_str());
      f_out->cd(("Vars/"+sys+"/BGSData").c_str());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        TH1D* h_bgs_data_tmp = (TH1D*)h_reco_data->Clone(("BGSData_"+std::to_string(i_u)).c_str());
        h_bgs_data_tmp->Add(h_bg_v.at(i_u),-1);
        h_bgs_data_tmp->Write(("BGSData_"+std::to_string(i_u)).c_str());   
      }

      f_out->cd();
      f_out->mkdir(("Cov/"+sys+"/BGSData").c_str());
      f_out->cd(("Cov/"+sys+"/BGSData").c_str());
      TH2D* h_cov_bgs_data_tmp = (TH2D*)c->Clone("Cov_BGSData");
      h_cov_bgs_data_tmp->Write("Cov_BGSData");

    }

    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      std::string sys = unisims_str.at(i_s);
      f_out->cd();
      f_out->mkdir(("Vars/"+sys+"/BG").c_str());
      f_out->cd(("Vars/"+sys+"/BG").c_str());
      TH1D* h = (TH1D*)f_hist->Get(("Reco/Vars/"+sys+"/h_AllBG").c_str());
      CrossSectionH(h,POT);
      mchm.Restore(h);
      h->Write("BG");

      f_out->cd();
      f_out->mkdir(("Cov/"+sys+"/BG").c_str());
      f_out->cd(("Cov/"+sys+"/BG").c_str());
      TH2D *c,*fc;
      CalcCovUnisim(sys,h_bg_cv,h,c,fc);
      mchm.Restore(c);
      c->Write("Cov_BG");

      f_out->cd();
      f_out->mkdir(("Vars/"+sys+"/Data").c_str());
      f_out->cd(("Vars/"+sys+"/Data").c_str());
      h_reco_data->Write("Data");

      f_out->cd();
      f_out->mkdir(("Cov/"+sys+"/Data").c_str());
      f_out->cd(("Cov/"+sys+"/Data").c_str());
      TH2D* h_cov_data_tmp = (TH2D*)c->Clone("Cov_Data");
      h_cov_data_tmp->Reset();
      h_cov_data_tmp->Write("Cov_Data");

      f_out->cd();
      f_out->mkdir(("Vars/"+sys+"/BGSData").c_str());
      f_out->cd(("Vars/"+sys+"/BGSData").c_str());
      TH1D* h_bgs_data_tmp = (TH1D*)h_reco_data->Clone("BGSData");
      h_bgs_data_tmp->Add(h,-1);
      h_bgs_data_tmp->Write("BGSData");   
      
      f_out->cd();
      f_out->mkdir(("Cov/"+sys+"/BGSData").c_str());
      f_out->cd(("Cov/"+sys+"/BGSData").c_str());
      TH2D* h_cov_bgs_data_tmp = (TH2D*)c->Clone("Cov_BGSData");
      h_cov_bgs_data_tmp->Write("Cov_BGSData");

    }

    // Calculate the data and BG subtracted data when using different total fluxes
    f_out->cd();
    f_out->mkdir("Vars/Flux/Data");
    f_out->cd("Vars/Flux/Data");
    std::vector<TH1D*> h_data_v; 
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++){

      double ratio = (h_numu_integrals->GetBinContent(i_u+1)*flux_numu + h_numubar_integrals->GetBinContent(i_u+1)*flux_numubar)/(flux_numu+flux_numubar);

      h_data_v.push_back(blinded ? (TH1D*)f_hist->Get("Reco/CV/h_Tot")->Clone(("h_reco_data_"+std::to_string(i_u)).c_str())
                                 : (TH1D*)f_hist->Get("Reco/CV/h_Data")->Clone(("h_reco_data_"+std::to_string(i_u)).c_str()));

      CrossSectionH(h_data_v.back(),POT*ratio);
      mchm.Restore(h_data_v.back());
      h_data_v.back()->Write(("Data_"+std::to_string(i_u)).c_str());
    }

    f_out->cd();
    f_out->mkdir("Cov/Flux/Data");
    f_out->cd("Cov/Flux/Data");
    TH2D *c,*fc;
    CalcCovMultisim("Flux",h_data_v,c,fc);
    mchm.Restore(c);
    c->Write("Cov_Data");

    f_out->cd();
    f_out->mkdir("Vars/Flux/BGSData");
    f_out->cd("Vars/Flux/BGSData");
    std::vector<TH1D*> h_bgs_data_v; 
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++){

      double ratio = (h_numu_integrals->GetBinContent(i_u+1)*flux_numu + h_numubar_integrals->GetBinContent(i_u+1)*flux_numubar)/(flux_numu+flux_numubar);

      h_bgs_data_v.push_back(blinded ? (TH1D*)f_hist->Get("Reco/CV/h_Tot")->Clone(("h_bgs_data_"+std::to_string(i_u)).c_str())
                                     : (TH1D*)f_hist->Get("Reco/CV/h_Data")->Clone(("h_bgs_data_"+std::to_string(i_u)).c_str()));

      TH1D* h_bg = (TH1D*)f_hist->Get(("Reco/Vars/Flux/h_AllBG_"+std::to_string(i_u)).c_str());
      h_bgs_data_v.back()->Add(h_bg,-1);
      CrossSectionH(h_bgs_data_v.back(),POT*ratio);

      mchm.Restore(h_bgs_data_v.back());
      h_bgs_data_v.back()->Write(("BGSData_"+std::to_string(i_u)).c_str());
    }

    f_out->cd();
    f_out->mkdir("Cov/Flux/BGSData");
    f_out->cd("Cov/Flux/BGSData");
    TH2D *c_flux,*fc_flux;
    CalcCovMultisim("Flux",h_data_v,c_flux,fc_flux);
    mchm.Restore(c_flux);
    c_flux->Write("Cov_BGSData");


    // Fold the generator predictions in different universes
    for(size_t i_g=0;i_g<generators.size();i_g++){

      std::string gen = generators.at(i_g);
      std::cout << gen << std::endl;

      const TH1D* h_gen_truth = (TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str());

      // Multisims
      for(int i_s=0;i_s<kSystMAX;i_s++){
        if(i_s == kFlux) continue;
        std::string sys = sys_str.at(i_s);
        f_out->cd();
        f_out->mkdir(("Vars/"+sys+"/"+gen).c_str());
        f_out->cd(("Vars/"+sys+"/"+gen).c_str());
        std::vector<TH1D*> h;
        for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
          TH2D* h_res = (TH2D*)f_hist->Get(("Response/Vars/"+sys+"/h_Signal_"+std::to_string(i_u)).c_str());
          h.push_back(Multiply(h_gen_truth,h_res,"Pred_"+std::to_string(i_u)));
          mchm.Restore(h.back());
          h.back()->Write(("Pred_"+std::to_string(i_u)).c_str());
        }
        f_out->cd();
        f_out->mkdir(("Cov/"+sys+"/"+gen).c_str());
        f_out->cd(("Cov/"+sys+"/"+gen).c_str());
        TH2D *c,*fc;
        CalcCovMultisim(sys,h,c,fc);
        mchm.Restore(c);
        c->Write("Cov_Pred");
      }

      // Unisims
      for(int i_s=0;i_s<kUnisimMAX;i_s++){
        std::string sys = unisims_str.at(i_s);
        f_out->cd();
        f_out->mkdir(("Vars/"+sys+"/"+gen).c_str());
        f_out->cd(("Vars/"+sys+"/"+gen).c_str());
        TH1D* h = Multiply(h_gen_truth,(TH2D*)f_hist->Get(("Response/Vars/"+sys+"/h_Signal").c_str()),"Pred");
        mchm.Restore(h);
        h->Write("Pred");

        f_out->cd();
        f_out->mkdir(("Cov/"+sys+"/"+gen).c_str());
        f_out->cd(("Cov/"+sys+"/"+gen).c_str());
        TH2D *c,*fc;
        CalcCovUnisim(sys,h_gen_ff_cv_v.at(i_g),h,c,fc);
        mchm.Restore(c);
        c->Write("Cov_Pred");
      }



      // Calculate the generator predictions in each flux universe
      const TH2D* h_gen_truth_2d = (TH2D*)f_gen->Get(("h_xsec_2D_"+var+"_"+gen).c_str());
      f_out->cd();
      f_out->mkdir(("Vars/Flux/"+gen).c_str());
      f_out->cd(("Vars/Flux/"+gen).c_str());
      std::vector<TH1D*> h_gen_ff_flux_v;
      for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++){
        TH1D* h_flux_ratio = (TH1D*)f_flux_ratios->Get(Form("ShapeRatios/NuMu/h_NuMu_FluxRatio_%i",i_u)); // Ratio of alt flux to CV in true nu_e
        TH1D* h_gen_truth_flux = Multiply(h_flux_ratio,h_gen_truth_2d,Form("h_gen_truth_f%i",i_u)); // Gen truth in flux universe
        h_gen_ff_flux_v.push_back(Multiply(h_gen_truth_flux,h_res_cv,Form("h_gen_f%i",i_u))); // Gen in reco space in flux universe
        delete h_gen_truth_flux;
        h_gen_ff_flux_v.back()->Write(("Pred_"+std::to_string(i_u)).c_str());
      }

      f_out->cd();
      f_out->mkdir(("Cov/Flux/"+gen).c_str());
      f_out->cd(("Cov/Flux/"+gen).c_str());
      TH2D *c,*fc;
      CalcCovMultisim("Flux",h_gen_ff_flux_v,c,fc);

    }

    f_hist->Close();
    f_gen->Close();

    std::cout << "Finished cleaning up" << std::endl;
  }

}
