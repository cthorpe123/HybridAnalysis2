#include "Funcs.h"
#include "BranchList.h"
#include "EnergyEstimatorFuncs.h"
#include "Systematics.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"
#include "WeightFuncs.h"

using namespace syst;

// Forward fold the CV truth through the CV response and the NuWro_0 response,
// then draw both results on the same figure for comparison.
// Also draws the efficiency (ProjectionX of response) for both universes.
// Computes a chi2 between the two FF predictions (+ CV background) using the
// response MC stat uncertainty plus the full systematic covariance.

void FFTest_ResComp(){

  const bool draw_underflow = false;
  const bool draw_overflow = false;
  const bool dbbw = true;
  bool add_detvars = false;
  const bool include_data_stat = true;

  std::vector<std::string> vars = {"MuonMom"};
  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  const std::string spec = "NuWro_0";

  for(size_t i_f=0;i_f<vars.size();i_f++){

    std::string label = vars.at(i_f);
    std::cout << label << std::endl;

    std::string plot_dir = "Analysis/"+label+"/Plots/FFTest/ResComp/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    TFile* f_in = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
    TFile* f_in_detvar = add_detvars ? TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str()) : nullptr;

    hist::MultiChannelHistogramManager mchm(label,true);
    mchm.SetTrueChannelList(channels_t);
    mchm.SetRecoChannelList(channels_r);
    mchm.LoadTemplates();

    TH1D* h_CV_Truth       = (TH1D*)f_in->Get("Truth/CV/h_Signal");
    TH1D* h_CV_Reco        = (TH1D*)f_in->Get("Reco/CV/h_Tot");
    TH1D* h_CV_Reco_AllBG  = (TH1D*)f_in->Get("Reco/CV/h_AllBG");
    TH2D* h_CV_Res         = (TH2D*)f_in->Get("Response/CV/h_Signal");
    TH2D* h_Spec_Res       = (TH2D*)f_in->Get(("Response/Special/"+spec+"/h_Signal").c_str());

    TH1D* h_CVT_CVRes   = Multiply(h_CV_Truth, h_CV_Res,   "h_CVT_CVRes");
    TH1D* h_CVT_SpecRes = Multiply(h_CV_Truth, h_Spec_Res, "h_CVT_SpecRes");

    // Fractional stat error on the FF prediction from response MC stats:
    // σ²(FF[j]) = Σᵢ Truth[i]² × σ²(Response[i,j])
    // Computed before Restore so the raw prediction values serve as the denominator.
    std::vector<double> reco_bins_v;
    for(int i=1;i<h_CV_Res->GetNbinsY()+2;i++) reco_bins_v.push_back(h_CV_Res->GetYaxis()->GetBinLowEdge(i));
    int n_reco_bins = reco_bins_v.size()-1;

    TH1D* h_CV_FracStat   = new TH1D("h_CV_FracStat",   ";Reco Variable;Fractional Stat. Error", n_reco_bins, &reco_bins_v[0]);
    TH1D* h_Spec_FracStat = new TH1D("h_Spec_FracStat", ";Reco Variable;Fractional Stat. Error", n_reco_bins, &reco_bins_v[0]);

    // Diagonal stat covariance for the chi2: σ²[j,j] = σ²_CV(j) + σ²_Spec(j)
    // (combined because chi2 tests the difference between the two FF predictions)
    TH2D* h_Stat_Cov = new TH2D("h_Stat_Cov", ";Reco Bin;Reco Bin", n_reco_bins, &reco_bins_v[0], n_reco_bins, &reco_bins_v[0]);

    for(int j=0;j<n_reco_bins+2;j++){
      double err2_CV = 0.0, err2_Spec = 0.0;
      for(int i=0;i<h_CV_Truth->GetNbinsX()+2;i++){
        double truth = h_CV_Truth->GetBinContent(i);
        err2_CV   += truth*truth * pow(h_CV_Res->GetBinError(i,j),  2);
        err2_Spec += truth*truth * pow(h_Spec_Res->GetBinError(i,j),2);
      }
      double pred_CV   = h_CVT_CVRes->GetBinContent(j);
      double pred_Spec = h_CVT_SpecRes->GetBinContent(j);
      h_CV_FracStat->SetBinContent(j,   pred_CV   > 0 ? sqrt(err2_CV)  /pred_CV   : 0.0);
      h_Spec_FracStat->SetBinContent(j, pred_Spec > 0 ? sqrt(err2_Spec)/pred_Spec : 0.0);
      h_Stat_Cov->SetBinContent(j, j, err2_CV + err2_Spec);
    }

    // Build full systematic covariance (mirrors FFTest_CVSpecRes).
    // Uses the same template TH2D from file (Reset to zero) so binning is guaranteed consistent.
    TH2D* h_Cov = (TH2D*)f_in->Get("Reco/Cov/Total/Cov_Tot");
    h_Cov->Reset();

    for(int i_s=0;i_s<kSystMAX;i_s++){
      if(i_s == kFlux) continue;
      std::vector<TH1D*> h;
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        std::string name = "h_Signal_FF_"+sys_str.at(i_s)+"_"+std::to_string(i_u);
        h.push_back(Multiply(h_CV_Truth,(TH2D*)f_in->Get(("Response/Vars/"+sys_str.at(i_s)+"/h_Signal_"+std::to_string(i_u)).c_str()),name.c_str()));
        ForceAddTH1D(h.back(),(TH1D*)f_in->Get(("Reco/Vars/"+sys_str.at(i_s)+"/h_AllBG_"+std::to_string(i_u)).c_str()));
      }
      TH2D *c,*fc;
      CalcCovMultisim(sys_str.at(i_s),h,c,fc);
      h_Cov->Add(c);
      for(TH1D* hh : h) delete hh;
      h.clear();
      delete c;
      delete fc;
    }

    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      std::string name = "h_Signal_FF_"+unisims_str.at(i_s);
      TH1D* h = Multiply(h_CV_Truth,(TH2D*)f_in->Get(("Response/Vars/"+unisims_str.at(i_s)+"/h_Signal").c_str()),name.c_str());
      ForceAddTH1D(h,(TH1D*)f_in->Get(("Reco/Vars/"+unisims_str.at(i_s)+"/h_AllBG").c_str()));
      TH2D *c,*fc;
      CalcCovUnisim(unisims_str.at(i_s),h_CV_Reco,h,c,fc);
      h_Cov->Add(c);
      delete h;
      delete c;
      delete fc;
    }

    if(add_detvars){
      TH1D* h_CV_Truth_Detvar = (TH1D*)f_in_detvar->Get("Truth/CV/h_Signal");
      TH1D* h_CV_Reco_Detvar  = (TH1D*)f_in_detvar->Get("Reco/CV/h_Tot");
      TH1D* h_CV_tmp = (TH1D*)f_in->Get("Reco/CV/h_Tot")->Clone("h_CV_tmp");
      for(int i_s=0;i_s<kDetvarMAX;i_s++){
        std::string name = "h_Signal_FF_"+detvar_str.at(i_s);
        TH1D* h = Multiply(h_CV_Truth_Detvar,(TH2D*)f_in_detvar->Get(("Response/Vars/"+detvar_str.at(i_s)+"/h_Signal").c_str()),name.c_str());
        ForceAddTH1D(h,(TH1D*)f_in_detvar->Get(("Reco/Vars/"+detvar_str.at(i_s)+"/h_AllBG").c_str()));
        TH2D *c,*fc;
        CalcCovUnisim(detvar_str.at(i_s),h_CV_Reco_Detvar,h,c,fc);
        for(int i=0;i<h_CV_Reco_Detvar->GetNbinsX()+2;i++)
          for(int j=0;j<h_CV_Reco_Detvar->GetNbinsX()+2;j++)
            c->SetBinContent(i,j,fc->GetBinContent(i,j)*h_CV_tmp->GetBinContent(i)*h_CV_tmp->GetBinContent(j));
        h_Cov->Add(c);
        delete h;
        delete c;
        delete fc;
      }
      delete h_CV_tmp;
    }

    TH2D* h_Cov_Flux = (TH2D*)f_in->Get("Reco/Cov/Flux/Cov_Tot");
    for(int i=0;i<h_Cov_Flux->GetNbinsX()+2;i++)
      for(int j=0;j<h_Cov_Flux->GetNbinsY()+2;j++)
        if(i != j) h_Cov_Flux->SetBinContent(i,j,0.0);
    h_Cov->Add(h_Cov_Flux);

    h_Cov->Add((TH2D*)f_in->Get("Reco/Cov/MCStat/Cov_AllBG"));

    if(include_data_stat) h_Cov->Add((TH2D*)f_in->Get("Reco/Cov/EstDataStat/Cov_Tot"));

    // Chi2 between the two FF predictions with full covariance:
    // response MC stat (h_Stat_Cov, diagonal) + systematics (h_Cov).
    // Background is added to both sides so absolute predictions match FFTest_CVSpecRes convention.
    h_Stat_Cov->Add(h_Cov);
    TH1D* h_CVT_CVRes_chi2   = (TH1D*)h_CVT_CVRes->Clone("h_CVT_CVRes_chi2");
    TH1D* h_CVT_SpecRes_chi2 = (TH1D*)h_CVT_SpecRes->Clone("h_CVT_SpecRes_chi2");
    h_CVT_CVRes_chi2->Add(h_CV_Reco_AllBG);
    h_CVT_SpecRes_chi2->Add(h_CV_Reco_AllBG);
    for(int j=0;j<h_CVT_CVRes->GetNbinsX()+2;j++){
      h_CVT_CVRes_chi2->SetBinError(j,   sqrt(h_Stat_Cov->GetBinContent(j,j)));
      h_CVT_SpecRes_chi2->SetBinError(j, 1e-10);
    }
    std::pair<double,int> chi2 = Chi2(h_CVT_CVRes_chi2, h_CVT_SpecRes_chi2, h_Stat_Cov, draw_overflow, draw_underflow);
    std::cout << spec << " chi2 = " << chi2.first
              << "  ndof = " << chi2.second
              << "  chi2/ndof = " << chi2.first/chi2.second << std::endl;

    // Load CV category histograms for the stacked plot (mirrors FFTest_CVSpecRes)
    std::vector<TH1D*> h_v;
    std::vector<int>   fill_colors;
    std::vector<std::string> legs;
    for(int i_c=0;i_c<kData;i_c++){
      h_v.push_back((TH1D*)f_in->Get(("Reco/CV/h_"+categories.at(i_c)).c_str()));
      fill_colors.push_back(cat_colors[i_c]);
      legs.push_back(categories.at(i_c));
      mchm.Restore(h_v.back());
      if(dbbw) DivideByBinWidth(h_v.back());
    }

    mchm.Restore(h_CVT_CVRes);
    mchm.Restore(h_CVT_SpecRes);
    mchm.Restore(h_CV_FracStat);
    mchm.Restore(h_Spec_FracStat);
    mchm.Restore(h_CVT_CVRes_chi2);
    mchm.Restore(h_CVT_SpecRes_chi2);

    if(dbbw) DivideByBinWidth(h_CVT_CVRes);
    if(dbbw) DivideByBinWidth(h_CVT_SpecRes);
    if(dbbw) DivideByBinWidth(h_CVT_CVRes_chi2);
    if(dbbw) DivideByBinWidth(h_CVT_SpecRes_chi2);
    // Fractional errors are dimensionless — no dbbw

    // Stacked comparison: CV categories + total-with-errors vs spec universe prediction
    pfs::DrawStacked(h_v, fill_colors, legs, h_CVT_CVRes_chi2, h_CVT_SpecRes_chi2,
                     draw_overflow, draw_underflow, plot_dir+spec+"_StackedComp.png", chi2);

    delete h_CVT_CVRes_chi2;
    delete h_CVT_SpecRes_chi2;
    delete h_Stat_Cov;

    // Signal-only comparison (no background, no systematics band)
    h_CVT_CVRes->SetLineStyle(1);
    h_CVT_SpecRes->SetLineStyle(2);

    std::vector<TH1D*> h_sig_v   = {h_CVT_CVRes, h_CVT_SpecRes};
    std::vector<int>   sig_colors = {kBlue+1, kRed+1};
    std::vector<std::string> sig_legs = {"CV Truth #times CV Response", "CV Truth #times "+spec+" Response"};

    pfs::DrawUnstacked2(h_sig_v, sig_colors, sig_legs, plot_dir+spec+"_ResComp.png", false);

    delete h_CVT_CVRes;
    delete h_CVT_SpecRes;

    h_CV_FracStat->SetLineStyle(1);
    h_Spec_FracStat->SetLineStyle(2);

    std::vector<TH1D*> h_frac_v    = {h_CV_FracStat, h_Spec_FracStat};
    std::vector<int>   frac_colors = {kBlue+1, kRed+1};
    std::vector<std::string> frac_legs = {"CV Frac. Stat. Error", spec+" Frac. Stat. Error"};

    pfs::DrawUnstacked2(h_frac_v, frac_colors, frac_legs, plot_dir+spec+"_FracStatErr.png", false);

    delete h_CV_FracStat;
    delete h_Spec_FracStat;

    // Efficiency: sum response over all reco bins (ProjectionX) then restore to physical binning
    int reco_bins = h_CV_Res->GetNbinsY() + 1;
    TH1D* h_CV_Eff   = h_CV_Res->ProjectionX("h_CV_Eff",   0, reco_bins);
    TH1D* h_Spec_Eff = h_Spec_Res->ProjectionX("h_Spec_Eff", 0, reco_bins);

    mchm.Restore(h_CV_Eff,   channels_t.at(0), true);
    mchm.Restore(h_Spec_Eff, channels_t.at(0), true);

    h_CV_Eff->SetLineStyle(1);
    h_Spec_Eff->SetLineStyle(2);

    std::vector<TH1D*> h_eff_v   = {h_CV_Eff, h_Spec_Eff};
    std::vector<int>   eff_colors = {kBlue+1, kRed+1};
    std::vector<std::string> eff_legs = {"CV Efficiency", spec+" Efficiency"};

    pfs::DrawUnstacked2(h_eff_v, eff_colors, eff_legs, plot_dir+spec+"_EffComp.png", true);

    delete h_CV_Eff;
    delete h_Spec_Eff;

    f_in->Close();

  }

}
