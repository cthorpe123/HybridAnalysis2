#include "Funcs.h"
#include "BranchList.h"
#include "EnergyEstimatorFuncs.h"
#include "Systematics.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"
#include "WeightFuncs.h"

using namespace syst;

// Forward fold the CV truth through the CV response and each special universe response,
// then draw both results on the same figure for comparison for each universe.
// Also draws the efficiency (ProjectionX of response) and fractional MC-stat error
// for both universe responses.
// Computes a chi2 between the two FF predictions (+ CV background) using the
// response MC stat uncertainty plus the full systematic covariance.
// A chi2 vs universe index curve is drawn at the end of each special universe family.

void FFTest_ResComp(){

  TCanvas* c = new TCanvas("c","c");

  const bool draw_underflow = false;
  const bool draw_overflow = false;
  const bool dbbw = true;
  bool add_detvars = false;
  const bool include_data_stat = true;
  const bool draw_chi2_curve = true;
  const bool add_nuwro_fd = true;


  std::vector<std::string> vars = {"MuonMom"};
  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  weight::SetWeightFuncs();
  std::vector<std::string> special_univs;
  for(const auto &item : weight::r_m)
    special_univs.push_back(item.first);
  if(add_nuwro_fd) special_univs.push_back("NuWro");

  int pts = weight::spline_pts;

  for(size_t i_f=0;i_f<vars.size();i_f++){

    std::string label = vars.at(i_f);
    std::cout << label << std::endl;

    std::string plot_dir = "Analysis/"+label+"/Plots/FFTest/";
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

    // CV truth forward-folded through CV response (raw bin space; kept raw for reuse across universes)
    TH1D* h_CVT_CVRes = Multiply(h_CV_Truth, h_CV_Res, "h_CVT_CVRes");

    // Reco bin edges for stat covariance histograms
    std::vector<double> reco_bins_v;
    for(int i=1;i<h_CV_Res->GetNbinsY()+2;i++) reco_bins_v.push_back(h_CV_Res->GetYaxis()->GetBinLowEdge(i));
    int n_reco_bins = reco_bins_v.size()-1;

    // Precompute CV response MC stat error² per reco bin:
    // σ²_CV(j) = Σᵢ Truth[i]² × σ²(Response_CV[i,j])
    // This is independent of the special universe so we compute it once.
    std::vector<double> err2_CV_v(n_reco_bins+2, 0.0);
    for(int j=0;j<n_reco_bins+2;j++){
      for(int i=0;i<h_CV_Truth->GetNbinsX()+2;i++){
        double truth = h_CV_Truth->GetBinContent(i);
        err2_CV_v[j] += truth*truth * pow(h_CV_Res->GetBinError(i,j), 2);
      }
    }

    // Build full systematic covariance (once per variable, reused across all universe plots)
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
      TH2D *cov,*fc;
      CalcCovMultisim(sys_str.at(i_s),h,cov,fc);
      h_Cov->Add(cov);
      for(TH1D* hh : h) delete hh;
      h.clear();
      delete cov;
      delete fc;
    }

    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      std::string name = "h_Signal_FF_"+unisims_str.at(i_s);
      TH1D* h = Multiply(h_CV_Truth,(TH2D*)f_in->Get(("Response/Vars/"+unisims_str.at(i_s)+"/h_Signal").c_str()),name.c_str());
      ForceAddTH1D(h,(TH1D*)f_in->Get(("Reco/Vars/"+unisims_str.at(i_s)+"/h_AllBG").c_str()));
      TH2D *cov,*fc;
      CalcCovUnisim(unisims_str.at(i_s),h_CV_Reco,h,cov,fc);
      h_Cov->Add(cov);
      delete h;
      delete cov;
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
        TH2D *cov,*fc;
        CalcCovUnisim(detvar_str.at(i_s),h_CV_Reco_Detvar,h,cov,fc);
        for(int i=0;i<h_CV_Reco_Detvar->GetNbinsX()+2;i++)
          for(int j=0;j<h_CV_Reco_Detvar->GetNbinsX()+2;j++)
            cov->SetBinContent(i,j,fc->GetBinContent(i,j)*h_CV_tmp->GetBinContent(i)*h_CV_tmp->GetBinContent(j));
        h_Cov->Add(cov);
        delete h;
        delete cov;
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

    // Load CV category histograms for stacked plots (once per variable, reused across universes)
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

    // CV fractional MC-stat error on the FF prediction (once per variable).
    // Computed in raw bin space, then restored for reuse via clones in each universe.
    TH1D* h_CV_FracStat = new TH1D("h_CV_FracStat", ";Reco Variable;Fractional Stat. Error", n_reco_bins, &reco_bins_v[0]);
    for(int j=0;j<n_reco_bins+2;j++){
      double pred_CV = h_CVT_CVRes->GetBinContent(j);
      h_CV_FracStat->SetBinContent(j, pred_CV > 0 ? sqrt(err2_CV_v[j])/pred_CV : 0.0);
    }
    mchm.Restore(h_CV_FracStat);

    // CV efficiency: ProjectionX of response, restored to physical true binning.
    // Kept for reuse via clones in each universe.
    int reco_bins_count = h_CV_Res->GetNbinsY() + 1;
    TH1D* h_CV_Eff = h_CV_Res->ProjectionX("h_CV_Eff", 0, reco_bins_count);
    mchm.Restore(h_CV_Eff, channels_t.at(0), true);

    // ── Loop over special universe families ──────────────────────────────────────
    for(std::string s : special_univs){

      gSystem->Exec(("mkdir -p "+plot_dir+"/"+s).c_str());
      std::vector<std::pair<double,int>> spec_chi2;

      for(int i_u=0;i_u<pts;i_u++){

        std::string spec = s + "_" + std::to_string(i_u);
        std::cout << spec << std::endl;

        TH2D* h_Spec_Res = (TH2D*)f_in->Get(("Response/Special/"+spec+"/h_Signal").c_str());

        // CV truth forward-folded through special universe response (raw bin space)
        TH1D* h_CVT_SpecRes = Multiply(h_CV_Truth, h_Spec_Res, "h_CVT_SpecRes");

        // ── Diagonal response MC-stat covariance for this universe pair ──────────
        // σ²(j) = σ²_CV(j) + σ²_Spec(j),
        // where σ²_Spec(j) = Σᵢ Truth[i]² × σ²(Response_Spec[i,j])
        TH2D* h_Stat_Cov = new TH2D("h_Stat_Cov", ";Reco Bin;Reco Bin",
                                     n_reco_bins, &reco_bins_v[0],
                                     n_reco_bins, &reco_bins_v[0]);
        TH1D* h_Spec_FracStat = new TH1D("h_Spec_FracStat", ";Reco Variable;Fractional Stat. Error",
                                          n_reco_bins, &reco_bins_v[0]);

        for(int j=0;j<n_reco_bins+2;j++){
          double err2_Spec = 0.0;
          for(int i=0;i<h_CV_Truth->GetNbinsX()+2;i++){
            double truth = h_CV_Truth->GetBinContent(i);
            err2_Spec += truth*truth * pow(h_Spec_Res->GetBinError(i,j), 2);
          }
          double pred_Spec = h_CVT_SpecRes->GetBinContent(j);
          h_Spec_FracStat->SetBinContent(j, pred_Spec > 0 ? sqrt(err2_Spec)/pred_Spec : 0.0);
          h_Stat_Cov->SetBinContent(j, j, err2_CV_v[j] + err2_Spec);
        }

        // Combine response MC-stat covariance with full systematic covariance
        h_Stat_Cov->Add(h_Cov);

        // ── Stacked comparison (with background, chi2) ───────────────────────────
        TH1D* h_CVT_CVRes_chi2   = (TH1D*)h_CVT_CVRes->Clone("h_CVT_CVRes_chi2");
        TH1D* h_CVT_SpecRes_chi2 = (TH1D*)h_CVT_SpecRes->Clone("h_CVT_SpecRes_chi2");
        h_CVT_CVRes_chi2->Add(h_CV_Reco_AllBG);
        h_CVT_SpecRes_chi2->Add(h_CV_Reco_AllBG);
        for(int j=0;j<h_CVT_CVRes->GetNbinsX()+2;j++){
          h_CVT_CVRes_chi2->SetBinError(j,   sqrt(h_Stat_Cov->GetBinContent(j,j)));
          h_CVT_SpecRes_chi2->SetBinError(j, 1e-10);
        }

        std::pair<double,int> chi2 = Chi2(h_CVT_CVRes_chi2, h_CVT_SpecRes_chi2, h_Stat_Cov, draw_overflow, draw_underflow);
        spec_chi2.push_back(chi2);
        std::cout << spec << " chi2 = " << chi2.first
                  << "  ndof = " << chi2.second
                  << "  chi2/ndof = " << chi2.first/chi2.second << std::endl;

        mchm.Restore(h_CVT_CVRes_chi2);
        mchm.Restore(h_CVT_SpecRes_chi2);
        if(dbbw) DivideByBinWidth(h_CVT_CVRes_chi2);
        if(dbbw) DivideByBinWidth(h_CVT_SpecRes_chi2);

        pfs::DrawStacked(h_v, fill_colors, legs, h_CVT_CVRes_chi2, h_CVT_SpecRes_chi2,
                         draw_overflow, draw_underflow,
                         plot_dir+"/"+s+"/"+spec+"_StackedComp.png", chi2);

        delete h_CVT_CVRes_chi2;
        delete h_CVT_SpecRes_chi2;
        delete h_Stat_Cov;

        // ── Signal-only comparison (no background) ───────────────────────────────
        TH1D* h_CVT_CVRes_sig   = (TH1D*)h_CVT_CVRes->Clone("h_CVT_CVRes_sig");
        TH1D* h_CVT_SpecRes_sig = (TH1D*)h_CVT_SpecRes->Clone("h_CVT_SpecRes_sig");
        mchm.Restore(h_CVT_CVRes_sig);
        mchm.Restore(h_CVT_SpecRes_sig);
        if(dbbw) DivideByBinWidth(h_CVT_CVRes_sig);
        if(dbbw) DivideByBinWidth(h_CVT_SpecRes_sig);
        h_CVT_CVRes_sig->SetLineStyle(1);
        h_CVT_SpecRes_sig->SetLineStyle(2);

        std::vector<TH1D*> h_sig_v    = {h_CVT_CVRes_sig, h_CVT_SpecRes_sig};
        std::vector<int>   sig_colors  = {kBlue+1, kRed+1};
        std::vector<std::string> sig_legs = {"CV Truth #times CV Response",
                                             "CV Truth #times "+spec+" Response"};
        pfs::DrawUnstacked2(h_sig_v, sig_colors, sig_legs,
                            plot_dir+"/"+s+"/"+spec+"_ResComp.png", false);
        delete h_CVT_CVRes_sig;
        delete h_CVT_SpecRes_sig;

        // ── Fractional MC-stat error comparison ──────────────────────────────────
        mchm.Restore(h_Spec_FracStat);
        TH1D* h_CV_FracStat_tmp = (TH1D*)h_CV_FracStat->Clone("h_CV_FracStat_tmp");
        h_CV_FracStat_tmp->SetLineStyle(1);
        h_Spec_FracStat->SetLineStyle(2);

        std::vector<TH1D*> h_frac_v    = {h_CV_FracStat_tmp, h_Spec_FracStat};
        std::vector<int>   frac_colors  = {kBlue+1, kRed+1};
        std::vector<std::string> frac_legs = {"CV Frac. Stat. Error", spec+" Frac. Stat. Error"};
        pfs::DrawUnstacked2(h_frac_v, frac_colors, frac_legs,
                            plot_dir+"/"+s+"/"+spec+"_FracStatErr.png", false);
        delete h_CV_FracStat_tmp;
        delete h_Spec_FracStat;

        // ── Efficiency comparison ────────────────────────────────────────────────
        TH1D* h_Spec_Eff = h_Spec_Res->ProjectionX("h_Spec_Eff", 0, reco_bins_count);
        mchm.Restore(h_Spec_Eff, channels_t.at(0), true);
        TH1D* h_CV_Eff_tmp = (TH1D*)h_CV_Eff->Clone("h_CV_Eff_tmp");
        h_CV_Eff_tmp->SetLineStyle(1);
        h_Spec_Eff->SetLineStyle(2);

        std::vector<TH1D*> h_eff_v    = {h_CV_Eff_tmp, h_Spec_Eff};
        std::vector<int>   eff_colors  = {kBlue+1, kRed+1};
        std::vector<std::string> eff_legs = {"CV Efficiency", spec+" Efficiency"};
        pfs::DrawUnstacked2(h_eff_v, eff_colors, eff_legs,
                            plot_dir+"/"+s+"/"+spec+"_EffComp.png", true);
        delete h_CV_Eff_tmp;
        delete h_Spec_Eff;

        delete h_CVT_SpecRes;

        if(s == "NuWro") break;

      } // end universe index loop

      // Draw chi2/ndof vs universe index curve for this special universe family
      if(draw_chi2_curve){
        TH1D* h_chi2 = new TH1D("h_chi2", ";Universe;#chi^{2}/ndof",
                                  spec_chi2.size(), 0.5, spec_chi2.size()+0.5);
        for(size_t i=0;i<spec_chi2.size();i++)
          h_chi2->SetBinContent(i+1, spec_chi2.at(i).first/spec_chi2.at(i).second);
        h_chi2->Draw("HIST");
        h_chi2->SetStats(0);
        c->Print((plot_dir+s+"_chi2.png").c_str());
        c->Clear();
        delete h_chi2;
      }

    } // end special universe family loop

    delete h_CVT_CVRes;
    delete h_CV_FracStat;
    delete h_CV_Eff;

    f_in->Close();

  } // end variable loop

}
