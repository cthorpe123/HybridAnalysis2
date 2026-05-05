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

void FFTest(){

  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");

  bool add_detvars = false;
  const bool include_data_stat = true;
  const bool draw_underflow = true;
  const bool draw_overflow = false;
  const bool dbbw = true;
  const bool draw_chi2_curve = true;
  const bool diag_only = false;
  const bool draw_cov = false;

  std::vector<std::string> vars = {"LeadPionE"};
  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  weight::SetWeightFuncs();
  std::vector<std::string> special_univs;
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

    std::vector<TH1D*> h_v;
    std::vector<int> fill_colors;
    std::vector<std::string> legs;

    TH1D* h_CV_Truth = (TH1D*)f_in->Get("Truth/CV/h_Signal");
    TH1D* h_CV_Reco = (TH1D*)f_in->Get("Reco/CV/h_Tot");
    TH2D* h_CV_Res = (TH2D*)f_in->Get("Response/CV/h_Signal");
  
    // Get the CV background prediction
    TH1D* h_CV_Reco_AllBG = (TH1D*)f_in->Get("Reco/CV/h_AllBG");

    for(int i_c=0;i_c<kData;i_c++){
      h_v.push_back((TH1D*)f_in->Get(("Reco/CV/h_"+categories.at(i_c)).c_str()));
      fill_colors.push_back(cat_colors[i_c]);
      legs.push_back(categories.at(i_c));
    }

    // Calculate the covariance encoding systemaics in the CV prediction
    // Build this by adding the total background to the FF signal in each universe
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

   // Do the same with the unisims
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

    // Need to be a bit careful with detvars - calculate fractional covariance 
    // in prediction using only detvar info, then scale this to the normal sample
    if(add_detvars){

      TH1D* h_CV_Truth_Detvar = (TH1D*)f_in_detvar->Get("Truth/CV/h_Signal");
      TH1D* h_CV_Reco_Detvar = (TH1D*)f_in_detvar->Get("Reco/CV/h_Tot");
      TH1D* h_CV_tmp = (TH1D*)f_in->Get("Reco/CV/h_Tot")->Clone("h_CV_tmp");

      for(int i_s=0;i_s<kDetvarMAX;i_s++){
        std::string name = "h_Signal_FF_"+detvar_str.at(i_s);
        TH1D* h = Multiply(h_CV_Truth_Detvar,(TH2D*)f_in_detvar->Get(("Response/Vars/"+detvar_str.at(i_s)+"/h_Signal").c_str()),name.c_str());
        ForceAddTH1D(h,(TH1D*)f_in_detvar->Get(("Reco/Vars/"+detvar_str.at(i_s)+"/h_AllBG").c_str()));
        TH2D *c,*fc; 
        CalcCovUnisim(detvar_str.at(i_s),h_CV_Reco_Detvar,h,c,fc); 
        for(int i=0;i<h_CV_Reco_Detvar->GetNbinsX()+2;i++){
          for(int j=0;j<h_CV_Reco_Detvar->GetNbinsX()+2;j++){
            c->SetBinContent(i,j,fc->GetBinContent(i,j)*h_CV_tmp->GetBinContent(i)*h_CV_tmp->GetBinContent(j));
          }
        }
        h_Cov->Add(c);
        delete h;
        delete c; 
        delete fc;
      }

    }

    // Add the covariance from the flux
    TH2D* h_Cov_Flux = (TH2D*)f_in->Get("Reco/Cov/Flux/Cov_Tot");
    for(int i=0;i<h_Cov_Flux->GetNbinsX()+2;i++)
      for(int j=0;j<h_Cov_Flux->GetNbinsY()+2;j++)
        if(i != j) h_Cov_Flux->SetBinContent(i,j,0.0);
    h_Cov->Add(h_Cov_Flux);

    // Add the MC stat error on the BG to the covariance matrix 
    h_Cov->Add((TH2D*)f_in->Get("Reco/Cov/MCStat/Cov_AllBG"));

    // If requested, also add the estimated stat error on the data
    if(include_data_stat) h_Cov->Add((TH2D*)f_in->Get("Reco/Cov/EstDataStat/Cov_Tot"));

    for(std::string s : special_univs){

      gSystem->Exec(("mkdir -p "+plot_dir+"/"+s).c_str());
      //std::cout << s << std::endl;
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

        // Show the spec truth folded through the CV response, to show how the special universe prediction differs from the CV in reco space
        TH1D* h_CV_Reco_tmp = (TH1D*)h_CV_Reco->Clone("h_CV_Reco_tmp");
        TH1D* h_SpecT_CVRes = Multiply(h_Spec_Truth,h_CV_Res,"h_SpecT_CVRes");
        ForceAddTH1D(h_SpecT_CVRes,h_CV_Reco_AllBG);
        mchm.Restore(h_SpecT_CVRes);
        mchm.Restore(h_CV_Reco_tmp);
        h_SpecT_CVRes->SetLineStyle(2);
        if(dbbw) DivideByBinWidth(h_SpecT_CVRes);
        if(dbbw) DivideByBinWidth(h_CV_Reco_tmp);

        pfs::DrawUnstacked2({h_CV_Reco_tmp,h_SpecT_CVRes},{1,1},{"CV","Spec Folded Through CV Response"},plot_dir+"/"+s+"/"+spec+"_SpecTimesCVRes.png",false);
        delete h_SpecT_CVRes;
        
        // Try folding the CV truth through the spec response, test if the change in model in truth 
        // space impacts the response in a dangerous way

        // Make clones of histograms needed
        TH1D* h_CV_Reco_tmp2 = (TH1D*)h_CV_Reco->Clone("h_CV_Reco_tmp2");
        std::vector<TH1D*> h_v_tmp;
        for(TH1D* h : h_v){
           h_v_tmp.push_back((TH1D*)h->Clone((string(h->GetName())+"_tmp").c_str()));
           mchm.Restore(h_v_tmp.back());
           if(dbbw) DivideByBinWidth(h_v_tmp.back());
        }

        TH2D* h_Res_Spec = (TH2D*)f_in->Get(("Response/Special/"+spec+"/h_Signal").c_str());
        TH1D* h_CVT_SpecRes = Multiply(h_CV_Truth,h_Res_Spec,"h_CVT_SpecRes");
        h_CVT_SpecRes->Add(h_CV_Reco_AllBG);

        // Stat error in the difference between the FF signal in the CV universe and 
        // the FF signal in the alternative universe
        TH2D* h_Stat_Cov = (TH2D*)f_in->Get(("Reco/Special/"+spec+"/SpecialStatCov/Cov_SpecialStat").c_str());
        h_Stat_Cov->Add(h_Cov);

        for(int i=0;i<h_CV_Reco->GetNbinsX()+2;i++){
          h_CVT_SpecRes->SetBinError(i,1e-10);
          h_CV_Reco_tmp2->SetBinError(i,sqrt(h_Stat_Cov->GetBinContent(i,i)));
        }

        std::pair<double,int> chi2 = Chi2(h_CV_Reco_tmp2,h_CVT_SpecRes,h_Stat_Cov,draw_overflow,draw_underflow);
        spec_chi2.push_back(chi2);
        std::cout << "chi2 = " << chi2.first << " ndof = " << chi2.second << " chi2/ndof = " << chi2.first/chi2.second << std::endl;
        
        mchm.Restore(h_CV_Reco_tmp2);
        mchm.Restore(h_CVT_SpecRes);
        if(dbbw) DivideByBinWidth(h_CV_Reco_tmp2);
        if(dbbw) DivideByBinWidth(h_CVT_SpecRes);

        pfs::DrawStacked(h_v_tmp,fill_colors,legs,h_CV_Reco_tmp2,h_CVT_SpecRes,draw_overflow,draw_underflow,plot_dir+"/"+s+"/"+spec+"_CVTimesSpecRes.png",chi2); 

        delete h_CVT_SpecRes;
        delete h_CV_Reco_tmp2;
        for(TH1D* hh : h_v_tmp) delete hh;
      }

      if(draw_chi2_curve){

        TH1D* h_chi2 = new TH1D("h_chi2",";Universe;#chi^{2}/ndof",spec_chi2.size(),0.5,spec_chi2.size()+0.5);
        std::map<std::string,std::pair<double,int>>::iterator it;
        for(size_t i=0;i<spec_chi2.size();i++)
          h_chi2->SetBinContent(i+1,spec_chi2.at(i).first/spec_chi2.at(i).second);

        h_chi2->Draw("HIST");
        h_chi2->SetStats(0);
        c->Print((plot_dir+s+"_chi2.png").c_str());   
        c->Clear();
        
        delete h_chi2;

      }

    }
    
    f_in->Close();

  }

}

