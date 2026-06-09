#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "Systematics.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"

using namespace syst;

// Calculate the background subtracted cross section prediction
// in reco space. Covariance is built by subtracting the per-universe
// background from the data, matching the approach of FFTest_CVSpecRes.C.

void CalcCrossSection(){

  bool add_detvars = false;
  bool blinded = true; // if blinded, use CV prediction as data
  bool draw_underflow = false;
  bool draw_overflow = false;
  bool divide_by_bin_width = true;
  bool include_data_stat = true;

  //std::vector<std::string> vars = var_names;
  std::vector<std::string> vars = {"MuonMom","MuonCosTheta"};
  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  for(size_t i_f=0;i_f<vars.size();i_f++){
    std::string label = vars.at(i_f);

    TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
    TFile* f_in_detvar = add_detvars ? TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str()) : nullptr;

    const double POT = ((TH1D*)f_in_hist->Get("Meta/POT"))->GetBinContent(1);
    std::cout << "POT = " << POT << std::endl;

    hist::MultiChannelHistogramManager mchm(label,false);
    mchm.SetTrueChannelList(channels_t);
    mchm.SetRecoChannelList(channels_r);
    mchm.LoadTemplates();

    TH1D* h_CV_AllBG = (TH1D*)f_in_hist->Get("Reco/CV/h_AllBG");
    TH1D* h_CV_Reco  = (TH1D*)f_in_hist->Get("Reco/CV/h_Tot");
    TH1D* h_Data     = (TH1D*)f_in_hist->Get(blinded ? "Reco/CV/h_Tot" : "Reco/CV/h_Data");

    // Background subtracted data in reco space (event counts)
    TH1D* h_Data_BGS = (TH1D*)h_Data->Clone("h_Data_BGS");
    h_Data_BGS->Add(h_CV_AllBG,-1);
    CrossSectionH(h_Data_BGS,POT);
    h_Data_BGS->GetYaxis()->SetTitle("d#sigma (10^{-40} cm^{2}/nucleon)");

    // Build covariance in event-count space before converting to cross section units
    TH2D* h_Cov = (TH2D*)f_in_hist->Get("Reco/Cov/Total/Cov_Tot");
    h_Cov->Reset();

    // Multisims: BGS data in each universe = h_Data - h_AllBG_universe
    for(int i_s=0;i_s<kSystMAX;i_s++){
      if(i_s == kFlux) continue;
      std::vector<TH1D*> h;
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        std::string name = "h_BGS_"+sys_str.at(i_s)+"_"+std::to_string(i_u);
        TH1D* h_u = (TH1D*)h_Data->Clone(name.c_str());
        h_u->Add((TH1D*)f_in_hist->Get(("Reco/Vars/"+sys_str.at(i_s)+"/h_AllBG_"+std::to_string(i_u)).c_str()),-1);
        CrossSectionH(h_u,POT);
        h.push_back(h_u);
      }
      TH2D *c,*fc;
      CalcCovMultisim(sys_str.at(i_s),h,c,fc);
      h_Cov->Add(c);
      for(TH1D* hh : h) delete hh;
      h.clear();
      delete c;
      delete fc;
    }

    // Unisims: BGS data in each universe = h_Data - h_AllBG_unisim
    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      std::string name = "h_BGS_"+unisims_str.at(i_s);
      TH1D* h = (TH1D*)h_Data->Clone(name.c_str());
      h->Add((TH1D*)f_in_hist->Get(("Reco/Vars/"+unisims_str.at(i_s)+"/h_AllBG").c_str()),-1);
      CrossSectionH(h,POT);
      TH2D *c,*fc;
      CalcCovUnisim(unisims_str.at(i_s),h_Data_BGS,h,c,fc);
      h_Cov->Add(c);
      delete h;
      delete c;
      delete fc;
    }

    // Flux - calculate the covariance in the predicted background and 
    // the fractional covariance in the selected signal
    if(true){
      TFile* f_flux = TFile::Open("Flux/FluxRatios.root");
      TH1D* h_numu_flux_int_ratio = (TH1D*)f_flux->Get("IntegralRatios/h_NuMu_IntegralRatio");
      TH1D* h_numubar_flux_int_ratio = (TH1D*)f_flux->Get("IntegralRatios/h_NuMuBar_IntegralRatio");

      std::vector<TH1D*> h;
      for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++){
        std::string name = "h_BGS_"+sys_str.at(kFlux)+"_"+std::to_string(i_u);
        TH1D* h_u = (TH1D*)f_in_hist->Get(("Reco/Vars/"+sys_str.at(kFlux)+"/h_Tot_"+std::to_string(i_u)).c_str());
        double flux_scale = (h_numu_flux_int_ratio->GetBinContent(i_u+1)*flux_numu + h_numubar_flux_int_ratio->GetBinContent(i_u+1)*flux_numubar)/(flux_numu + flux_numubar);
        CrossSectionH(h_u,POT*flux_scale);
        h.push_back(h_u);
      }
      TH2D *c,*fc;
      CalcCovMultisim(sys_str.at(kFlux),h,c,fc);
      h_Cov->Add(c);
      for(TH1D* hh : h) delete hh;
      h.clear();
      delete c;
      delete fc;
    }  

    // Set errors on BGS data from covariance diagonal
    for(int i_b=0;i_b<h_Data_BGS->GetNbinsX()+2;i_b++)
      h_Data_BGS->SetBinError(i_b,sqrt(h_Cov->GetBinContent(i_b,i_b)));

    mchm.Restore(h_Data_BGS);

    if(divide_by_bin_width){
      DivideByBinWidth(h_Data_BGS);
      h_Data_BGS->GetYaxis()->SetTitle("d#sigma/dvar (10^{-40} cm^{2}/nucleon/unit)");
    }

    std::cout << "Total Cross Section = " << IntegralWithOU(h_Data_BGS) << " 1e-40 cm2" << std::endl;  

    std::string plot_dir = "Analysis/"+label+"/Plots/CalcCrossSection/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    std::string leg = blinded ? "BGS Asimov Data" : "BGS Data";
    pfs::DrawUnstacked2({h_Data_BGS},{kBlack},{leg},plot_dir+"XS_BGS_"+label+".png",true);

    delete h_Data_BGS;
    f_in_hist->Close();
    if(f_in_detvar) f_in_detvar->Close();
  }

}
