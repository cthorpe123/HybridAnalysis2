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
#include "DetvarHistograms.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"

using namespace syst;

void CheckEff(){

  // Label and set the branches defining the selection and systematics
  bool draw_truth = true;
  bool draw_hist = true; // Grab the CV from the non-detvar file
  bool draw_o=false,draw_u=false;

  std::vector<std::string> vars = {"MuonMom","MuonCosTheta","NProt","NPi","NSh","ProtonKE","PionE","PiZeroE","W"};
  for(int i_e=0;i_e<ee::kMAX;i_e++)
    vars.push_back(ee::estimators_str.at(i_e));

  std::vector<std::string> channels_t = {"1p0pi0g","2p0pi0g","1p1pi0g","2p1pi0g","1p0pi1g","2p0pi1g","1p1pi1g","2p1pi1g","1p0pi2g","2p0pi2g","1p1pi2g","2p1pi2g"};
  std::vector<std::string> channels_r = {"All"};

  for(std::string label : vars){

    std::string plot_dir = "Analysis/"+label+"/Plots/CheckEff/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    TFile* f_in = draw_hist ? TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str()) : nullptr;

    hist::MultiChannelHistogramManager mchm(label,true);
    mchm.SetTrueChannelList(channels_t);
    mchm.SetRecoChannelList(channels_r);
    mchm.LoadTemplates();

    TH2D* h_CV_Response_Signal = (TH2D*)f_in->Get("Response/CV/h_Signal");
    int bins = h_CV_Response_Signal->GetNbinsY()+1;
    for(int i=0;i<h_CV_Response_Signal->GetNbinsX()+2;i++){
      double tot = 0.0;  
      for(int j=0;j<h_CV_Response_Signal->GetNbinsY()+2;j++){
        tot += h_CV_Response_Signal->GetBinContent(i,j);
        //std::cout << i << " " << j << " " << h_CV_Response_Signal->GetBinContent(i,j) << std::endl;
      }
        std::cout << i << " "  << tot << std::endl;
         }


    TH1D* h_CV_Eff_Signal = h_CV_Response_Signal->ProjectionX("h_CV_Eff",0,bins); 

    //for(int i=0;i<h_CV_Eff_Signal->GetNbinsX()+2;i++) std::cout << i << " " <<  h_CV_Eff_Signal->GetBinContent(i) << std::endl;

    std::vector<TH1D*> h_eff_v;
    std::vector<int> fill_colors;
    std::vector<std::string> legs;
    for(size_t i_ch=0;i_ch<channels_t.size();i_ch++){
      std::string ch = channels_t.at(i_ch);
      //std::cout << ch << std::endl;
      h_eff_v.push_back((TH1D*)h_CV_Eff_Signal->Clone(("Eff_"+ch).c_str()));
      mchm.Restore(h_eff_v.back(),ch,true);
      //for(int i=0;i<h_eff_v.back()->GetNbinsX()+1;i++) std::cout << i <<  " " << h_eff_v.back()->GetBinContent(i) << std::endl;
      fill_colors.push_back(i_ch+1); 
      legs.push_back(ch);
    } 

    pfs::DrawUnstacked(h_eff_v,fill_colors,legs,draw_o,draw_u,plot_dir+"Efficiency.png"); 

  }


}
