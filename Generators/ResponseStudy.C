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

void ResponseStudy(){

  bool load_asimov = true;
  std::vector<std::string> vars = {"Norm","Enu","MuonMom","MuonCosTheta","LeadProtonKE","LeadPionE"};
  //std::vector<std::string> vars = var_names;
  //vars.push_back("Enu");
  //vars.push_back("Norm");
  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro","GiBUU"};
  
  bool blinded = true;
  bool restore = true;
  bool dbbw = true;
  bool add_detvars = false;
  bool draw_o = false;
  bool draw_u = false;

  for(const std::string& var : vars){
    
    std::cout << var << std::endl;

    std::string plot_dir = "Analysis/"+var+"/Plots/ResponseStudy/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    hist::MultiChannelHistogramManager mchm(var);
    mchm.LoadTemplates();    

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
    h_reco_data->SetDirectory(0);

    // Subtract the predicted background from the data
    h_reco_data->Add((TH1D*)f_hist->Get("Reco/CV/h_AllBG"),-1);

    // Convert BG subtracted data to cross section space
    CrossSectionH(h_reco_data,POT);
    h_reco_data->GetYaxis()->SetTitle("d#sigma (10^{-38} cm^{2}/unit)");

    // Data just gets stat error for now, plus MC stat error from BG subtracted
    TH2D* h_cov_data_stat = blinded ? (TH2D*)f_hist->Get("Reco/Cov/EstDataStat/Cov_Tot") : (TH2D*)f_hist->Get("Reco/Cov/MCStat/Cov_Data");
    CrossSectionCov(h_cov_data_stat,POT);
    h_cov_breakdown_v["DataStat"] = h_cov_data_stat;   
    h_cov_data->Add(h_cov_data_stat); 

    // Add errors to the background subtracted data
    for(int i=0;i<h_reco_data->GetNbinsX()+2;i++) 
      h_reco_data->SetBinError(i,sqrt(h_cov_data->GetBinContent(i,i)));

    // Now take the generator prediction and fold through the CV detector resposne,
    // and create the total covariance matrix for each
    TH2D* h_res_cv = (TH2D*)f_hist->Get("Response/CV/h_Signal");
    TFile* f_gen = TFile::Open(("Analysis/"+var+"/rootfiles/GeneratorXSec.root").c_str());
    std::vector<TH1D*> h_gen_ff_cv_v;
    std::vector<TH2D*> h_gen_cov_v;

    if(restore) mchm.Restore(h_reco_data);

    // Use the covariance from the first generator as a template for the others
    TH2D* h_gen_fcov_ref = nullptr;

    for(size_t i_g=0;i_g<generators.size();i_g++){
      std::string gen = generators.at(i_g);
      h_gen_ff_cv_v.push_back(Multiply((TH1D*)f_gen->Get(("h_xsec_"+var+"_"+gen).c_str()),h_res_cv,"h_xsec_ff_"+gen)); 
      h_gen_cov_v.push_back((TH2D*)h_cov_data->Clone(("h_gen_cov_"+gen).c_str()));
      h_gen_cov_v.back()->Reset();
      legs.push_back(gen);
      cols.push_back(i_g+2);
    }

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
        pfs::Draw2DHist(c,plot_dir+"Cov_"+sys_str.at(i_s)+"_"+gen+".png");
        pfs::Draw2DHist(fc,plot_dir+"FCov_"+sys_str.at(i_s)+"_"+gen+".png");
        h_gen_cov_v.at(i_g)->Add(c);
        h_cov_breakdown_v[sys_str.at(i_s)] = (TH2D*)c->Clone(("h_cov_"+sys_str.at(i_s)).c_str());
        h_cov_breakdown_v.at(sys_str.at(i_s))->SetDirectory(0);
        for(TH1D* hh : h) delete hh;
        h.clear();
        delete c;
        delete fc;
      }

      
      // Unisims - group together into one systematic
      TH1D* h_cv = (TH1D*)h_gen_ff_cv_v.at(i_g)->Clone("h_cv");
      
      TH2D* h_cov_unisim = (TH2D*)h_cov_data->Clone("h_cov_unisim");
      TH2D* h_fcov_unisim = (TH2D*)h_cov_data->Clone("h_cov_unisim");
      h_cov_unisim->Reset();
      h_fcov_unisim->Reset();

      for(int i_s=0;i_s<kUnisimMAX;i_s++){

        std::string name = "h_"+gen+"_ff_"+unisims_str.at(i_s);
      
        // FF generator pred in this universe
        TH1D* h = Multiply(h_gen_truth,(TH2D*)f_hist->Get(("Response/Vars/"+unisims_str.at(i_s)+"/h_Signal").c_str()),name.c_str());

        TH2D *c,*fc;
        CalcCovUnisim(unisims_str.at(i_s),h_cv,h,c,fc);

        pfs::DrawUnstacked({h_cv,h},{1,2},{"CV","Alt"},draw_o,draw_u,false,dbbw,plot_dir+"Univs_"+unisims_str.at(i_s)+".png");

        pfs::Draw2DHist(c,plot_dir+"Cov_"+unisims_str.at(i_s)+"_"+gen+".png");
        pfs::Draw2DHist(fc,plot_dir+"FCov_"+unisims_str.at(i_s)+"_"+gen+".png");        
        h_gen_cov_v.at(i_g)->Add(c);
        h_cov_unisim->Add(c);
        h_fcov_unisim->Add(fc);
        delete h;
        delete c;
        delete fc;
      }

      pfs::Draw2DHist(h_cov_unisim,plot_dir+"Cov_Unisim_"+gen+".png");
      pfs::Draw2DHist(h_fcov_unisim,plot_dir+"FCov_Unisim_"+gen+".png");

      delete h_cv; 
      h_cov_breakdown_v["Unisims"] = (TH2D*)h_cov_unisim->Clone("h_cov_Unisims");
      h_cov_breakdown_v.at("Unisims")->SetDirectory(0);
      delete h_cov_unisim;
      delete h_fcov_unisim;

      for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++)
        h_gen_ff_cv_v.at(i_g)->SetBinError(i_b,sqrt(h_gen_cov_v.at(i_g)->GetBinContent(i_b,i_b)));

      // For the chi2, include the data statistics
      TH2D* h_cov = (TH2D*)h_gen_cov_v.at(i_g)->Clone("h_cov"); 
      h_cov->Add(h_cov_data);

      // If the first generator in the stack, copy and store the fractional covariance for comparisons later
      if(i_g == 0){
        h_gen_fcov_ref = (TH2D*)h_gen_cov_v.at(i_g)->Clone("h_gen_fcov_ref");
        for(int i_b=0;i_b<h_gen_fcov_ref->GetNbinsX()+2;i_b++)
          for(int j_b=0;j_b<h_gen_fcov_ref->GetNbinsX()+2;j_b++)
            h_gen_fcov_ref->SetBinContent(i_b,j_b,h_gen_fcov_ref->GetBinContent(i_b,j_b)/h_gen_ff_cv_v.at(i_g)->GetBinContent(i_b)/h_gen_ff_cv_v.at(i_g)->GetBinContent(j_b));
      }

      // Calculate the chi2 of each generator with the data
      std::pair<double,int> chi2 = Chi2(h_gen_ff_cv_v.at(i_g),h_reco_data,h_cov,false,false);
      chi2s.push_back(to_string_with_precision(chi2.first/chi2.second,2));
      legs.at(i_g) += ", #chi^{2}/n = " + chi2s.at(i_g);

      if(restore) mchm.Restore(h_gen_ff_cv_v.at(i_g));
      h_gen_ff_cv_v.at(i_g)->GetYaxis()->SetTitle("d#sigma (10^{-38} cm^{2}/unit)");
      
      pfs::DrawUnstacked({h_gen_ff_cv_v.at(i_g),h_reco_data},{cols.at(i_g),1},{legs.at(i_g),"Data"},draw_o,draw_u,true,dbbw,plot_dir+"Test_"+gen+".png");

      // Draw breakdown of systematics budget
      std::vector<TH1D*> h_sys;
      std::vector<std::string> sys_legs;
      std::vector<int> sys_cols;
      int i_s=0;
      for(auto const& item : h_cov_breakdown_v){
        if(item.first == "DataStat") continue;
        TH1D* h = (TH1D*)h_reco_data->Clone(("h_sys_"+item.first).c_str()); 
        for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++){
          h->SetBinContent(i_b,sqrt(item.second->GetBinContent(i_b,i_b))/h_gen_ff_cv_v.at(i_g)->GetBinContent(i_b));
        } 
        if(restore) mchm.Restore(h);
        h_sys.push_back(h);
        sys_legs.push_back(item.first);
        sys_cols.push_back(i_s+2);
        i_s++;
      } 
      
      TH1D* h = (TH1D*)h_reco_data->Clone("h_sys_tot"); 
      for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++) 
        h->SetBinContent(i_b,sqrt(h_cov->GetBinContent(i_b,i_b))/h_gen_ff_cv_v.at(i_g)->GetBinContent(i_b));
      if(restore) mchm.Restore(h);
      h_sys.push_back(h); 
      sys_legs.push_back("Total");
      sys_cols.push_back(1);
    
      // Dummy histogram to force equal axis scale
      TH1D* h_dummy = (TH1D*)h_sys.back()->Clone("h_dummy");
      for(int i_b=1;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+1;i_b++) h_dummy->SetBinContent(i_b,0.07);
      h_sys.push_back(h_dummy); 
      sys_legs.push_back("");
      sys_cols.push_back(0);

      pfs::DrawUnstacked(h_sys,sys_cols,sys_legs,draw_o,draw_u,false,false,plot_dir+"SysBreakdown_"+gen+".png");

      delete h_dummy;

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
 
    for(TH1D*& h : h_gen_ff_cv_v){
      if(restore) mchm.Restore(h);
      h->GetYaxis()->SetTitle("d#sigma (10^{-38} cm^{2}/unit)");
    }

    pfs::DrawUnstacked(h_gen_ff_cv_v,cols,legs,draw_o,draw_u,true,dbbw,plot_dir+"Test.png");

    // Now instead draw the data and generators, using the reference for the fractional covariance
    // on each generator
    chi2s.clear();
    legs.clear();
    for(size_t i_g=0;i_g<generators.size();i_g++){
      TH2D* h_cov = (TH2D*)h_gen_fcov_ref->Clone("h_cov"); 
      
      for(int i_b=0;i_b<h_gen_fcov_ref->GetNbinsX()+2;i_b++)
          for(int j_b=0;j_b<h_gen_fcov_ref->GetNbinsX()+2;j_b++)
            h_cov->SetBinContent(i_b,j_b,h_gen_fcov_ref->GetBinContent(i_b,j_b)*h_gen_ff_cv_v.at(i_g)->GetBinContent(i_b)*h_gen_ff_cv_v.at(i_g)->GetBinContent(j_b));
     
      for(int i_b=0;i_b<h_gen_ff_cv_v.at(i_g)->GetNbinsX()+2;i_b++) 
        h_gen_ff_cv_v.at(i_g)->SetBinError(i_b,sqrt(h_cov->GetBinContent(i_b,i_b))); 

      h_gen_ff_cv_v.at(i_g)->SetTitle("d#sigma (10^{-38} cm^{2}/unit)");
      
      h_cov->Add(h_cov_data);
      std::cout << h_gen_ff_cv_v.at(i_g)->GetNbinsX() << std::endl;
      std::cout << h_reco_data->GetNbinsX() << std::endl;
      std::cout << h_cov->GetNbinsX() << std::endl;
      std::pair<double,int> chi2 = Chi2(h_gen_ff_cv_v.at(i_g),h_reco_data,h_cov,false,false);
      chi2s.push_back(to_string_with_precision(chi2.first/chi2.second,2));
      legs.push_back(generators.at(i_g)+ ", #chi^{2}/n = " + chi2s.at(i_g));
      
    }
    
    legs.push_back(!blinded ? "uboone Data" : "Asimov Data");
    pfs::DrawUnstacked(h_gen_ff_cv_v,cols,legs,draw_o,draw_u,true,dbbw,plot_dir+"Test_Ref.png");

    f_hist->Close();
    f_gen->Close();

  }

}
