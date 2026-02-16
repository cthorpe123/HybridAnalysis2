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
  bool draw_underflow = false;
  bool draw_overflow = true;

  std::vector<std::string> label_v = {"MuonMom"};
  std::vector<std::string> special_univ_v = {"ExtraPi","ExtraPi2","ExtraP","ExtraP2","ExtraGamma","ExtraGamma2"};

  for(size_t i_f=0;i_f<label_v.size();i_f++){

    std::string label = label_v.at(i_f);
    std::string plot_dir = "Analysis/"+label+"/Plots/ResponseCheck/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    TFile* f_in = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());

    // Get the CV reco prediction with signal only
    TH1D* h_CV = (TH1D*)f_in->Get("Reco/CV/h_Tot");
    TH1D* h_CV_Truth = (TH1D*)f_in->Get("Truth/CV/h_Signal");

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


    // Add the MC stat error on the BG to the covariance matrix 
    for(int i_c=0;i_c<kData;i_c++){
      if(i_c == kSignal) continue;
      h_Cov->Add((TH2D*)f_in->Get(("Reco/Cov/MCStat/Cov_"+categories.at(i_c)).c_str()));
    }
 
    // If requested, also add the estimated stat error on the data
    if(include_data_stat) h_Cov->Add((TH2D*)f_in->Get("Reco/Cov/EstDataStat/Cov_Tot"));


    for(std::string spec : special_univ_v){

      std::cout << spec << std::endl;

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
      std::cout << "chi2 = " << chi2.first << " ndof = " << chi2.second << " chi2/ndof = " << chi2.first/chi2.second << std::endl;
    
      TH1D* h_CV_tmp = (TH1D*)h_CV->Clone("h_CV_tmp");    
      pfs::DrawStacked(h_v,fill_colors,legs,h_CV_tmp,h_FF_Spec,draw_overflow,draw_underflow,plot_dir+"FF_Signal_"+spec+".png"); 
      delete h_CV_tmp; 

    }


  }

}

