#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "EnergyEstimatorFuncs.h"
#include "Systematics.h"
#include "PlotFuncs.h"

using namespace syst;

// Try forward folding the CV truth through the response
// calulated in the CV and special universes, calculate chi2
// between the CV and each special prediction

void FFTest(){

  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");

  const bool include_data_stat = false;
  const bool draw_underflow = false;
  const bool draw_overflow = false;
  const bool dbbw = false;
  const bool draw_chi2_curve = true;

  std::vector<std::string> label_v = {"MuonMom","ProtonE","W","PeLEELike0Pi"};

  // Special universe setup
  const int spline_pts = 100;
  std::vector<std::string> special_univ = {"ExtraPi","ExtraP","ExtraG"};

  for(size_t i_f=0;i_f<label_v.size();i_f++){

    std::string label = label_v.at(i_f);
    std::cout << "label = " << label << std::endl;
    std::string plot_dir = "Analysis/"+label+"/Plots/FFTest/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    TFile* f_in = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());

    // Get the CV reco prediction with signal only
    TH1D* h_CV = (TH1D*)f_in->Get("Reco/CV/h_Tot");
    TH1D* h_CV_Truth = (TH1D*)f_in->Get("Truth/CV/h_Signal");

    std::cout << h_CV->GetNbinsX() << "  " << h_CV_Truth->GetNbinsX() <<  std::endl;

    // Get the CV background prediction
    TH1D* h_CV_BG = (TH1D*)h_CV->Clone("h_CV_BG");
    std::vector<TH1D*> h_v;
    std::vector<int> fill_colors;
    std::vector<std::string> legs;
    h_CV_BG->Reset();
    for(int i_c=0;i_c<kData;i_c++){
      h_v.push_back((TH1D*)f_in->Get(("Reco/CV/h_"+categories.at(i_c)).c_str()));
      fill_colors.push_back(cat_colors[i_c]); 
      legs.push_back(categories.at(i_c));
      if(i_c == kSignal) continue;
      h_CV_BG->Add(h_v.back());
    } 

    // Calculate the covariance encoding systemaics in the CV prediction
    // Build this by adding the total background to the FF signal in each
    // universe
    TH2D* h_Cov = (TH2D*)f_in->Get("Reco/Cov/Total/Cov_Tot");   
    h_Cov->Reset();
    for(int i_s=0;i_s<kSystMAX;i_s++){
      if(i_s == kFlux) continue;
      std::vector<TH1D*> h;
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        h.push_back(Multiply(h_CV_Truth,(TH2D*)f_in->Get(("Response/Vars/"+sys_str.at(i_s)+"/h_Signal_"+std::to_string(i_u)).c_str()),("h_Signal_FF_"+sys_str.at(i_s)+"_"+std::to_string(i_u)).c_str()));
        for(int i_c=0;i_c<kData;i_c++){
          if(i_c == kSignal) continue;
          h.back()->Add((TH1D*)f_in->Get(("Reco/Vars/"+sys_str.at(i_s)+"/h_"+categories.at(i_c)+"_"+std::to_string(i_u)).c_str()));
        } 
      }
      TH2D *c,*fc; 
      CalcCovMultisim(sys_str.at(i_s),h,c,fc); 
      h_Cov->Add(c);
    }

    // Add the covariance from the flux (for signal and background) to the covariance matrix
    h_Cov->Add((TH2D*)f_in->Get("Reco/Cov/Flux/Cov_Tot"));

    // Add the MC stat error on the BG to the covariance matrix 
    for(int i_c=0;i_c<kData;i_c++){
      if(i_c == kSignal) continue;
      h_Cov->Add((TH2D*)f_in->Get(("Reco/Cov/MCStat/Cov_"+categories.at(i_c)).c_str()));
    }

    // If requested, also add the estimated stat error on the data
    if(include_data_stat) h_Cov->Add((TH2D*)f_in->Get("Reco/Cov/EstDataStat/Cov_Tot"));

    if(dbbw) for(TH1D* h : h_v) DivideByBinWidth(h);

    for(std::string s : special_univ){

      gSystem->Exec(("mkdir -p "+plot_dir+"/"+s).c_str());
      std::cout << s << std::endl;

      std::vector<std::pair<double,int>> spec_chi2;

      for(int i=0;i<spline_pts;i++){

        std::string spec = s + "_" + std::to_string(i); 

        // CV folded through response calculated in special universe
        TH1D* h_FF_Spec = (TH1D*)f_in->Get(("Reco/Special/"+spec+"/FoldedCV/h_Signal").c_str());
        h_FF_Spec->Add(h_CV_BG);

        // Stat error in the difference between the FF signal in the CV universe and 
        // the FF signal in the alternative universe
        TH2D* h_Stat_Cov = (TH2D*)f_in->Get(("Reco/Special/"+spec+"/SpecialStatCov/Cov_SpecialStat").c_str());
        h_Stat_Cov->Add(h_Cov);

        for(int i=0;i<h_FF_Spec->GetNbinsX()+2;i++){
          h_FF_Spec->SetBinError(i,1e-10);
          h_CV->SetBinError(i,sqrt(h_Stat_Cov->GetBinContent(i,i)));
        }

        std::pair<double,int> chi2 = Chi2(h_CV,h_FF_Spec,h_Stat_Cov,draw_overflow,draw_underflow);
        spec_chi2.push_back(chi2);
        std::cout << "chi2 = " << chi2.first << " ndof = " << chi2.second << " chi2/ndof = " << chi2.first/chi2.second << std::endl;

        TH1D* h_CV_tmp = (TH1D*)h_CV->Clone("h_CV_tmp");    

        if(dbbw){
          DivideByBinWidth(h_CV_tmp);
          DivideByBinWidth(h_FF_Spec);
        }
       
        pfs::DrawStacked(h_v,fill_colors,legs,h_CV_tmp,h_FF_Spec,draw_overflow,draw_underflow,plot_dir+"/"+s+"/"+"FF_Signal_"+spec+".png",chi2); 
        delete h_CV_tmp; 

      }

      if(draw_chi2_curve){

        TH1D* h_chi2 = new TH1D("h_chi2",";Universe;#chi^{2}/ndof",spec_chi2.size(),0.5,spec_chi2.size()+0.5);
        std::map<std::string,std::pair<double,int>>::iterator it;
        for(size_t i=0;i<spec_chi2.size();i++)
          h_chi2->SetBinContent(i,spec_chi2.at(i).first/spec_chi2.at(i).second);

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

