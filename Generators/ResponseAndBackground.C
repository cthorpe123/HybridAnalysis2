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

// Compare the fractional errors in the combined FF generators and background 
// doing the calculation using the full set of response matrices/background
// predictions, and calculating the covariance in the selected background and 
// response separately

void ResponseAndBackground(){

  bool load_asimov = true;
  std::vector<std::string> vars = {"Norm","Enu","MuonMom","MuonCosTheta","LeadProtonKE","LeadPionE"};
  //std::vector<std::string> vars = {"MuonMom"};
  //std::vector<std::string> vars = var_names;
  //vars.push_back("Enu");
  //vars.push_back("Norm");
  
  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro","GiBUU"};
  //std::vector<std::string> generators = {"Untunedv3.0.6"};

  bool blinded = true;
  bool restore = true;
  bool dbbw = true;
  bool add_detvars = false;
  bool draw_o = false;
  bool draw_u = false;
  bool diag_only = false;

  for(const std::string& var : vars){
    
    std::cout << var << std::endl;

    std::string plot_dir = "Analysis/"+var+"/Plots/ResponseAndBackground/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    hist::MultiChannelHistogramManager mchm(var);
    mchm.LoadTemplates();    

    std::vector<std::string> legs_full;
    std::vector<std::string> legs_nobg;
    std::vector<std::string> chi2s_full;
    std::vector<std::string> chi2s_nobg;
    std::vector<int> cols;

    TFile* f_hist = TFile::Open(("Analysis/"+var+"/rootfiles/Histograms.root").c_str());
    TFile* f_gen = TFile::Open(("Analysis/"+var+"/rootfiles/GeneratorXSec.root").c_str());
    const double POT = ((TH1D*)f_hist->Get("Meta/POT"))->GetBinContent(1);

    // Load the reconstructed data with stat errors
    TH1D* h_reco_data = blinded ? (TH1D*)f_hist->Get("Reco/CV/h_Tot")->Clone("h_reco_data") : (TH1D*)f_hist->Get("Reco/CV/h_Data")->Clone("h_reco_data");
    h_reco_data->SetDirectory(0);

    TH2D* h_res_cv = (TH2D*)f_hist->Get("Response/CV/h_Signal");
    TH2D* h_cov_template = (TH2D*)f_hist->Get("Reco/Cov/Total/Cov_Tot");
    h_cov_template->Reset();

    TH2D* h_cov_data_stat = (TH2D*)f_hist->Get("Reco/Cov/EstDataStat/Cov_Tot");
    h_cov_data_stat->Add((TH2D*)f_hist->Get("Reco/Cov/Flux/Cov_Tot"));

    for(int i=0;i<h_reco_data->GetNbinsX()+2;i++) 
      h_reco_data->SetBinError(i,sqrt(h_cov_data_stat->GetBinContent(i,i)));

    CrossSectionH(h_reco_data,POT);
    CrossSectionCov(h_cov_data_stat,POT);

    TH2D* h_cov_bg_mc_stat = (TH2D*)f_hist->Get("Reco/Cov/MCStat/Cov_AllBG");
    CrossSectionCov(h_cov_bg_mc_stat,POT);

    TH1D* h_bg_cv = (TH1D*)f_hist->Get("Reco/CV/h_AllBG");
    CrossSectionH(h_bg_cv,POT);

    // Fold the generators through the response and combine with the background
    std::vector<TH1D*> h_gen_ff_cv_v;
    for(size_t i_g=0;i_g<generators.size();i_g++){
      std::string gen = generators.at(i_g);
      h_gen_ff_cv_v.push_back(Multiply((TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str()),h_res_cv,"h_xsec_ff_"+gen)); 
      ForceAddTH1D(h_gen_ff_cv_v.back(),h_bg_cv);
      cols.push_back(i_g+2);
    }

    // Get the covariances in the background from each source of systematic
    TH2D* h_cov_bg_tot = (TH2D*)f_hist->Get("Reco/Cov/Total/Cov_AllBG");
    h_cov_bg_tot->Reset();

    std::map<std::string,TH2D*> h_cov_bg_m;

    for(int i_s=0;i_s<kSystMAX;i_s++){
      //if(i_s == kFlux) continue;
      h_cov_bg_m[sys_str.at(i_s)] = (TH2D*)f_hist->Get(("Reco/Cov/"+sys_str.at(i_s)+"/Cov_AllBG").c_str());
      CrossSectionCov(h_cov_bg_m.at(sys_str.at(i_s)),POT);
      h_cov_bg_tot->Add(h_cov_bg_m.at(sys_str.at(i_s)));
    }
        
    h_cov_bg_m["Unisims"] = (TH2D*)h_cov_template->Clone("h_cov_bg_unisim");
    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      h_cov_bg_m[unisims_str.at(i_s)] = (TH2D*)f_hist->Get(("Reco/Cov/"+unisims_str.at(i_s)+"/Cov_AllBG").c_str());
      CrossSectionCov(h_cov_bg_m.at(unisims_str.at(i_s)),POT);
      h_cov_bg_m.at("Unisims")->Add(h_cov_bg_m.at(unisims_str.at(i_s)));
      h_cov_bg_tot->Add(h_cov_bg_m.at(unisims_str.at(i_s)));
    }
    
    h_cov_bg_m["Total"] = h_cov_bg_tot;

    // Multiply the generator predictions by the response in each genie, g4 and detvar universe, combine with background
    std::vector<TH2D*> h_cov_full_tot_v;
    std::vector<TH2D*> h_cov_nobg_tot_v;
    for(size_t i_g=0;i_g<generators.size();i_g++){

      std::string gen = generators.at(i_g);
      std::cout << gen << std::endl;

      // load the generator prediction in truth space
      const TH1D* h_gen_truth = (TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str());
      
      TH2D* h_cov_full_tot = (TH2D*)h_cov_template->Clone("h_cov_full_tot");
      TH2D* h_cov_nobg_tot = (TH2D*)h_cov_template->Clone("h_cov_nobg_tot");

      std::map<std::string,TH2D*> h_cov_full_m; // Cov calculated including background vars
      std::map<std::string,TH2D*> h_cov_nobg_m; // Cov calculated neglecting BG contribution

      // Multisims
      
      for(int i_s=0;i_s<kSystMAX;i_s++){
        //if(i_s == kFlux) continue;
        std::vector<TH1D*> h_full;
        std::vector<TH1D*> h_nobg;
        std::vector<TH1D*> h_bg;
        for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
          std::string name = "h_"+gen+"_ff_"+sys_str.at(i_s)+"_"+std::to_string(i_u);
         
          h_full.push_back(Multiply(h_gen_truth,(TH2D*)f_hist->Get(("Response/Vars/"+sys_str.at(i_s)+"/h_Signal_"+std::to_string(i_u)).c_str()),name.c_str()));
          h_nobg.push_back((TH1D*)h_full.back()->Clone((name+"_nobg").c_str()));
          h_bg.push_back((TH1D*)f_hist->Get(("Reco/Vars/"+sys_str.at(i_s)+"/h_AllBG_"+std::to_string(i_u)).c_str()));
          CrossSectionH(h_bg.back(),POT);
          ForceAddTH1D(h_full.back(),h_bg.back());

        }
      
        TH2D *c_full,*fc_full;
        CalcCovMultisim(sys_str.at(i_s),h_full,c_full,fc_full);
        h_cov_full_m[sys_str.at(i_s)] = ((TH2D*)c_full->Clone(("h_cov_full_"+sys_str.at(i_s)+"_"+gen).c_str()));
        h_cov_full_m.at(sys_str.at(i_s))->SetDirectory(0);
        h_cov_full_tot->Add(h_cov_full_m.at(sys_str.at(i_s)));
        delete c_full;
        delete fc_full;

        TH2D *c_nobg,*fc_nobg;
        CalcCovMultisim(sys_str.at(i_s),h_nobg,c_nobg,fc_nobg);
        h_cov_nobg_m[sys_str.at(i_s)] = ((TH2D*)c_nobg->Clone(("h_cov_nobg_"+sys_str.at(i_s)+"_"+gen).c_str()));
        h_cov_nobg_m.at(sys_str.at(i_s))->SetDirectory(0);
        h_cov_nobg_tot->Add(h_cov_nobg_m.at(sys_str.at(i_s)));
        delete c_nobg;
        delete fc_nobg;

        for(TH1D* hh : h_full) delete hh;
        for(TH1D* hh : h_bg) delete hh;
        for(TH1D* hh : h_nobg) delete hh;

      }
      
      // Unisims - group together into one systematic
      TH1D* h_cv_full = (TH1D*)h_gen_ff_cv_v.at(i_g)->Clone("h_cv"); // h_gen_ff_cv_v already has background
      TH1D* h_cv_nobg = Multiply(h_gen_truth,h_res_cv,"h_cv_nobg"); // nobg cv does not include bg

      h_cov_full_m["Unisims"] = (TH2D*)h_cov_template->Clone(("h_cov_full_Unisims_"+gen).c_str());
      h_cov_nobg_m["Unisims"] = (TH2D*)h_cov_template->Clone(("h_cov_nobg_Unisims_"+gen).c_str());

      for(int i_s=0;i_s<kUnisimMAX;i_s++){
        std::string name = "h_"+gen+"_ff_"+unisims_str.at(i_s);
      
        // FF generator pred in this universe
        TH1D* h_full = Multiply(h_gen_truth,(TH2D*)f_hist->Get(("Response/Vars/"+unisims_str.at(i_s)+"/h_Signal").c_str()),name.c_str());
        TH1D* h_nobg = (TH1D*)h_full->Clone("h_nobg");
        TH1D* h_bg = (TH1D*)f_hist->Get(("Reco/Vars/"+unisims_str.at(i_s)+"/h_AllBG").c_str());
        CrossSectionH(h_bg,POT);
        ForceAddTH1D(h_full,h_bg);

        TH2D *c_full,*fc_full;
        CalcCovUnisim(unisims_str.at(i_s),h_cv_full,h_full,c_full,fc_full);
        h_cov_full_m[unisims_str.at(i_s)] = ((TH2D*)c_full->Clone(("h_cov_full_"+unisims_str.at(i_s)+"_"+gen).c_str()));
        h_cov_full_m.at(unisims_str.at(i_s))->SetDirectory(0);
        h_cov_full_m.at("Unisims")->Add(c_full);
        h_cov_full_tot->Add(c_full);

        delete c_full;
        delete fc_full;
        delete h_full;

        TH2D *c_nobg,*fc_nobg;
        CalcCovUnisim(unisims_str.at(i_s),h_cv_nobg,h_nobg,c_nobg,fc_nobg);
        h_cov_nobg_m[unisims_str.at(i_s)] = ((TH2D*)c_nobg->Clone(("h_cov_nobg_"+unisims_str.at(i_s)+"_"+gen).c_str()));
        h_cov_nobg_m.at(unisims_str.at(i_s))->SetDirectory(0);
        h_cov_nobg_m.at("Unisims")->Add(c_nobg);
        h_cov_nobg_tot->Add(c_nobg);  

        delete c_nobg;
        delete fc_nobg;
        delete h_nobg;

        delete h_bg;
        
      }
      
      // Include the data and BG MC stat in the mix, and build total

      h_cov_full_tot->Add(h_cov_bg_mc_stat);
      h_cov_nobg_tot->Add(h_cov_bg_mc_stat);

      h_cov_full_m["BGMCStat"] = (TH2D*)h_cov_bg_mc_stat->Clone(("h_cov_full_BGMCStat_"+gen).c_str());
      h_cov_nobg_m["BGMCStat"] = (TH2D*)h_cov_bg_mc_stat->Clone(("h_cov_nobg_BGMCStat_"+gen).c_str());
      
      TH2D* h_cov_full_tot_wdatastat = (TH2D*)h_cov_full_tot->Clone("h_cov_full_tot_wdatastat");
      TH2D* h_cov_nobg_tot_wdatastat = (TH2D*)h_cov_nobg_tot->Clone("h_cov_nobg_tot_wdatastat");

      h_cov_full_tot_wdatastat->Add(h_cov_data_stat);
      h_cov_nobg_tot_wdatastat->Add(h_cov_data_stat);

      h_cov_full_m["DataStat"] = (TH2D*)h_cov_data_stat->Clone(("h_cov_full_DataStat_"+gen).c_str());
      h_cov_nobg_m["DataStat"] = (TH2D*)h_cov_data_stat->Clone(("h_cov_nobg_DataStat_"+gen).c_str());

      //pfs::Draw2DHist(h_cov_full_tot_wdatastat,plot_dir+"Cov_Full_"+gen+".png");
      //pfs::Draw2DHist(h_cov_nobg_tot_wdatastat,plot_dir+"Cov_NoBG_"+gen+".png");

      h_cov_full_m["Total"] = h_cov_full_tot_wdatastat;
      h_cov_nobg_m["Total"] = h_cov_nobg_tot_wdatastat;

      // Now try the two different ways of assembling the systematics budget

      // Add the background uncertainty from earlier into the nobg covariances, and calc FE
      int color = 1;
      std::vector<TH1D*> h_fe_full_v;
      std::vector<TH1D*> h_fe_nobg_v;
      std::vector<std::string> sys_v;
      std::vector<int> sys_cols;
      std::map<std::string,TH2D*>::iterator it;
      for(std::map<std::string,TH2D*>::iterator it = h_cov_nobg_m.begin();it != h_cov_nobg_m.end();it++){
        std::string sys = it->first;

        if(sys != "DataStat" && sys != "BGMCStat")
          h_cov_nobg_m.at(sys)->Add(h_cov_bg_m.at(sys));
          
        //pfs::Draw2DHist(h_cov_full_m.at(sys),plot_dir+"Cov_"+sys+"_"+gen+"_full.png");
        //pfs::Draw2DHist(h_cov_nobg_m.at(sys),plot_dir+"Cov_"+sys+"_"+gen+"_nobg.png");

        TH1D* h_fe_full = (TH1D*)h_gen_ff_cv_v.back()->Clone(("h_fe_full_"+sys).c_str());
        TH1D* h_fe_nobg = (TH1D*)h_gen_ff_cv_v.back()->Clone(("h_fe_nobg_"+sys).c_str());
        for(int i_b=0;i_b<h_fe_full->GetNbinsX()+2;i_b++){
          double x = h_gen_ff_cv_v.at(i_g)->GetBinContent(i_b);
          h_fe_full->SetBinContent(i_b,sqrt(h_cov_full_m.at(sys)->GetBinContent(i_b,i_b))/x);
          h_fe_nobg->SetBinContent(i_b,sqrt(h_cov_nobg_m.at(sys)->GetBinContent(i_b,i_b))/x);
        }



        mchm.Restore(h_fe_full);
        mchm.Restore(h_fe_nobg);

        h_fe_full->GetYaxis()->SetTitle("Frac. Unc.");

        h_fe_nobg->SetLineStyle(2);
        h_fe_nobg->GetYaxis()->SetTitle("Frac. Unc.");

        pfs::DrawUnstacked({h_fe_full,h_fe_nobg},{color,color},{"Full","Sep. BG"},draw_o,draw_u,false,false,plot_dir+"FE_"+sys+"_"+gen+".png");

        color++;
        if(in_vec(unisims_str,sys)) continue;
        h_fe_nobg->SetLineStyle(1);
        h_fe_full_v.push_back(h_fe_full);
        h_fe_nobg_v.push_back(h_fe_nobg);
        sys_v.push_back(sys);
        sys_cols.push_back(color-1);
         
      }

      pfs::DrawUnstacked(h_fe_full_v,sys_cols,sys_v,draw_o,draw_u,false,false,plot_dir+"FE_Full_"+gen+".png");
      pfs::DrawUnstacked(h_fe_nobg_v,sys_cols,sys_v,draw_o,draw_u,false,false,plot_dir+"FE_NoBG_"+gen+".png");
      
      TH2D* h_corr_full_tot_wdatastat = CalcCorrelationMatrix("Full",h_cov_full_tot_wdatastat);
      TH2D* h_corr_nobg_tot_wdatastat = CalcCorrelationMatrix("NoBG",h_cov_nobg_tot_wdatastat);
      pfs::Draw2DHist(h_corr_full_tot_wdatastat,plot_dir+"Corr_Full_"+gen+".png");
      pfs::Draw2DHist(h_corr_nobg_tot_wdatastat,plot_dir+"Corr_NoBG_"+gen+".png");  

      std::pair<double,int> chi2_full = Chi2(h_gen_ff_cv_v.at(i_g),h_reco_data,h_cov_full_tot_wdatastat,draw_o,draw_u,diag_only);
      chi2s_full.push_back(to_string_with_precision(chi2_full.first/chi2_full.second,2));
      legs_full.push_back(generators.at(i_g)+ ", #chi^{2}/n = " + chi2s_full.at(i_g));

      std::pair<double,int> chi2_nobg = Chi2(h_gen_ff_cv_v.at(i_g),h_reco_data,h_cov_nobg_tot_wdatastat,draw_o,draw_u,diag_only);
      chi2s_nobg.push_back(to_string_with_precision(chi2_nobg.first/chi2_nobg.second,2));
      legs_nobg.push_back(generators.at(i_g)+ ", #chi^{2}/n = " + chi2s_nobg.at(i_g));

      delete h_cov_full_tot_wdatastat;
      delete h_cov_nobg_tot_wdatastat;

      h_cov_full_tot_v.push_back(h_cov_full_tot);
      h_cov_nobg_tot_v.push_back(h_cov_nobg_tot);

    }

    h_gen_ff_cv_v.push_back(h_reco_data);
    cols.push_back(1);
    legs_full.push_back(!blinded ? "uboone Data" : "Asimov Data");
    legs_nobg.push_back(!blinded ? "uboone Data" : "Asimov Data");

    mchm.Restore(h_gen_ff_cv_v.back()); // Restore binning for data as well

    for(size_t i_g=0;i_g<generators.size();i_g++){
      for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++) 
        h_gen_ff_cv_v.at(i_g)->SetBinError(i_b,sqrt(h_cov_full_tot_v.at(i_g)->GetBinContent(i_b,i_b)));
      mchm.Restore(h_gen_ff_cv_v.at(i_g));
    }
    pfs::DrawUnstacked(h_gen_ff_cv_v,cols,legs_full,draw_o,draw_u,true,dbbw,plot_dir+"Test_Full.png");

    for(size_t i_g=0;i_g<generators.size();i_g++){
      for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++) 
        h_gen_ff_cv_v.at(i_g)->SetBinError(i_b,sqrt(h_cov_nobg_tot_v.at(i_g)->GetBinContent(i_b,i_b)));
    mchm.Restore(h_gen_ff_cv_v.at(i_g));
   }



    pfs::DrawUnstacked(h_gen_ff_cv_v,cols,legs_nobg,draw_o,draw_u,true,dbbw,plot_dir+"Test_NoBG.png");  

    f_hist->Close();
    f_gen->Close();

  }

}