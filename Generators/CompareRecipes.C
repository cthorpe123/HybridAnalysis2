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

void CompareRecipes(){

  std::vector<std::string> vars = {"MuonMom","MuonCosTheta","LeadProtonKE","ProtonKE"};
  //std::vector<std::string> vars = var_names;
  vars.push_back("Enu");
  vars.push_back("Norm");
  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro","GiBUU"};
  std::vector<std::string> recipes = {"Recipe1","Recipe2","Recipe3","Recipe1_5"};
  std::vector<int> styles = {2,3,4,6};
  int n_r=recipes.size();
  bool add_detvars = false;
  bool draw_o = false;
  bool draw_u = false;
  bool shape_only = false;

  std::vector<std::string> sys_v = {"Total","DataStat","BGMCStat"};
  for(int i_s=0;i_s<kSystMAX;i_s++) sys_v.push_back(sys_str.at(i_s));
  for(int i_s=0;i_s<kUnisimMAX;i_s++) sys_v.push_back(unisims_str.at(i_s));

  for(const std::string& var : vars){
    std::cout << var << std::endl;

    std::string plot_dir = !shape_only ? "Analysis/"+var+"/Plots/CompareRecipes/" 
                                       : "Analysis/"+var+"/Plots/CompareRecipesShape/";
    
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    // Open the files contianing the different recipes 
    std::vector<TFile*> f_in_v;
    for(std::string r : recipes){
      if(!shape_only) f_in_v.push_back(TFile::Open(("Analysis/"+var+"/rootfiles/"+r+".root").c_str()));
      else f_in_v.push_back(TFile::Open(("Analysis/"+var+"/rootfiles/"+r+"_ShapeOnly.root").c_str()));
    }

    // Big map storing the covariance calculated for each systematic, generator and recipe
    // h_cov_m[sys][gen][rec]
    std::map<std::string,std::vector<std::vector<TH2D*>>> h_cov_m;
    std::map<std::string,std::vector<std::vector<TH1D*>>> h_fe_m;
    for(std::string sys : sys_v){
      h_cov_m[sys] = std::vector<std::vector<TH2D*>>(generators.size(),std::vector<TH2D*>(recipes.size(),nullptr));
      h_fe_m[sys] = std::vector<std::vector<TH1D*>>(generators.size(),std::vector<TH1D*>(recipes.size(),nullptr));
      for(size_t i_g=0;i_g<generators.size();i_g++){
        for(size_t i_r=0;i_r<recipes.size();i_r++){
          h_cov_m[sys].at(i_g).at(i_r) = (TH2D*)f_in_v.at(i_r)->Get((generators.at(i_g)+"/"+"Cov_"+sys).c_str());
          h_fe_m[sys].at(i_g).at(i_r) = (TH1D*)f_in_v.at(i_r)->Get((generators.at(i_g)+"/"+"FE_"+sys).c_str());
          //h_fe_m[sys].at(i_g).at(i_r)->SetLineStyle(styles.at(i_r));
        }
      }
    }

    TH1D* h_data = (TH1D*)f_in_v.at(0)->Get((generators.at(0)+"/BGSData").c_str()); 

    std::vector<std::vector<std::pair<double,int>>> chi2_v;
    for(size_t i_g=0;i_g<generators.size();i_g++){
      std::string gen = generators.at(i_g);
      std::string plot_dir_gen = !shape_only ? "Analysis/"+var+"/Plots/CompareRecipes/"+gen+"/" 
                                             : "Analysis/"+var+"/Plots/CompareRecipesShape/"+gen+"/";
      gSystem->Exec(("mkdir -p " + plot_dir_gen).c_str());

      TH1D* h_pred = (TH1D*)f_in_v.at(0)->Get((gen+"/Pred").c_str()); 

      // Compare the fractional errors by category between methods
      for(std::string sys : sys_v){
        std::vector<TH1D*> fe_v = h_fe_m[sys].at(i_g);
        fe_v.at(0)->GetYaxis()->SetTitle("Frac. Unc.");
        pfs::DrawUnstacked(fe_v,styles,recipes,draw_o,draw_u,false,false,plot_dir_gen+"FE_"+sys+".png");
      }

      // Compute the total FE from unisims
      std::vector<TH2D*> h_cov_unisim_v;
      std::vector<TH1D*> h_fe_unisim_v;
      for(size_t i_r=0;i_r<recipes.size();i_r++){
        h_cov_unisim_v.push_back((TH2D*)h_cov_m["Total"].at(i_g).at(i_r)->Clone(("Unisim_"+recipes.at(i_r)).c_str()));
        h_cov_unisim_v.back()->Reset();
        for(std::string u : unisims_str) h_cov_unisim_v.back()->Add(h_cov_m[u].at(i_g).at(i_r));
        h_fe_unisim_v.push_back((TH1D*)h_pred->Clone(("h_fe_unisim_"+recipes.at(i_r)).c_str()));
        MakeFEHist(h_fe_unisim_v.back(),h_pred,h_cov_unisim_v.back()); 
      }
      pfs::DrawUnstacked(h_fe_unisim_v,styles,recipes,draw_o,draw_u,false,false,plot_dir_gen+"FE_Unisim.png");


      // Calculate the chi2 with the asimov data/beam data using each recipe and compare
      chi2_v.push_back(std::vector<std::pair<double,int>>());
      for(size_t i_r=0;i_r<recipes.size();i_r++){
        
        h_pred = (TH1D*)f_in_v.at(i_r)->Get((gen+"/Pred").c_str()); 
        h_data = (TH1D*)f_in_v.at(i_r)->Get((generators.at(i_g)+"/BGSData").c_str()); 
        const TH2D* h_cov =  h_cov_m["Total"].at(i_g).at(i_r);
        const TH2D* h_cov_data = h_cov_m["DataStat"].at(i_g).at(i_r);
        TH2D* h_cov_pred = (TH2D*)h_cov->Clone("h_cov_pred"); 
        h_cov_pred->Add(h_cov_m["DataStat"].at(i_g).at(i_r),-1);

        std::pair<double,int> chi2 = Chi2(h_pred,h_data,h_cov,draw_o,draw_u,false);
        chi2_v.back().push_back(chi2);
        std::cout << "Recipe " << i_r+1 << " "  << var << " " << gen << ", chi2/ndof = " << chi2.first << "/" << chi2.second << " = " << chi2.first/chi2.second << std::endl;
        
        std::pair<double,int> chi2_diag = Chi2(h_pred,h_data,h_cov,draw_o,draw_u,true);
        std::cout << "Recipe " << i_r+1 << " "  << var << " " << gen << ", chi2/ndof = " << chi2_diag.first << "/" << chi2_diag.second << " = " << chi2_diag.first/chi2_diag.second << std::endl;

        for(int i=0;i<h_data->GetNbinsX()+2;i++){
          h_data->SetBinError(i,sqrt(h_cov_data->GetBinContent(i,i)));
          h_pred->SetBinError(i,sqrt(h_cov_pred->GetBinContent(i,i)));
        }

        std::string leg = gen + ", chi2/ndof = " + to_string_with_precision(chi2.first/chi2.second,2);
        pfs::DrawUnstacked({h_data,h_pred},{1,(int)i_g+2},{"Data",leg},draw_o,draw_u,true,true,plot_dir_gen+"Test_"+recipes.at(i_r)+".png");

        for(int i=0;i<h_data->GetNbinsX()+2;i++){
          h_data->SetBinError(i,0);
          h_pred->SetBinError(i,sqrt(h_cov->GetBinContent(i,i)));
        }

        pfs::DrawUnstacked({h_data,h_pred},{1,(int)i_g+2},{"Data",leg},draw_o,draw_u,true,true,plot_dir_gen+"Test_Alt_"+recipes.at(i_r)+".png");

        leg = gen + ", chi2/ndof = " + to_string_with_precision(chi2_diag.first/chi2_diag.second,2);
        pfs::DrawUnstacked({h_data,h_pred},{1,(int)i_g+2},{"Data",leg},draw_o,draw_u,true,true,plot_dir_gen+"Test_DiagOnly_"+recipes.at(i_r)+".png");
        

        delete h_cov_pred;

      }

    }

    // Draw all generators together with their chi2s
    for(size_t i_r=0;i_r<recipes.size();i_r++){
      std::vector<TH1D*> h_pred_v;
      std::vector<std::string> legs;
      std::vector<int> cols;
      for(size_t i_g=0;i_g<generators.size();i_g++){
        std::string gen =  generators.at(i_g);
        cols.push_back(i_g+2);
        h_pred_v.push_back((TH1D*)f_in_v.at(0)->Get((gen+"/Pred").c_str())->Clone(("h_Pred_"+gen).c_str()));
        legs.push_back(gen+" #chi^{2}/n="+to_string_with_precision(chi2_v.at(i_g).at(i_r).first/chi2_v.at(i_g).at(i_r).second,2));
        const TH2D* h_cov =  h_cov_m["Total"].at(i_g).at(i_r);
        for(int i_b=0;i_b<h_pred_v.back()->GetNbinsX()+2;i_b++) h_pred_v.back()->SetBinError(i_b,sqrt(h_cov->GetBinContent(i_b,i_b)));
      }

      h_pred_v.push_back((TH1D*)f_in_v.at(i_r)->Get((generators.at(0)+"/BGSData").c_str())->Clone("h_data"));
      for(int i_b=0;i_b<h_pred_v.back()->GetNbinsX()+2;i_b++) h_pred_v.back()->SetBinError(i_b,0);
      cols.push_back(1);
      legs.push_back("Asimov Data");

      pfs::DrawUnstacked(h_pred_v,cols,legs,draw_o,draw_u,true,true,plot_dir+"Test_"+recipes.at(i_r)+".png");

      for(TH1D* hh : h_pred_v) delete hh;
    }

    for(TFile* f : f_in_v) f->Close();

  }

}