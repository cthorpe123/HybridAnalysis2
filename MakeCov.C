#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "Systematics.h"

using namespace syst;

// Try using different combinations of cuts and frameworks to calculate W

void MakeCov(){

  std::string label = "RecoW";
  double* var = &W_lt; 
  bool* sel = &sel_h8;

  const double scale = 1.0;

  const double data_POT = 1.332E+20;
  const double data_Trig = 31582916.0;

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";

  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root"
  };

  std::vector<double> weights = {
    scale*data_POT/7.88166e+20,
    scale*data_POT/3.06E+20,
    scale*data_Trig/88445969.0  
  };

  std::vector<bool> load_syst_v = {
    true,
    false,
    false
  };

  std::vector<bool> is_data_v = {
    false,
    false,
    true
  };

  // Giant contained for all of the systematics
  std::vector<std::vector<std::vector<TH1D*>>> h_Vars;
  std::vector<TH1D*> h_CV;

  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_CV.push_back(new TH1D(("h_CV_"+categories.at(i_c)).c_str(),";Reco W (GeV);Events",50,0.95,5.0));
    h_Vars.push_back(std::vector<std::vector<TH1D*>>());
    for(int i_s=0;i_s<kSystMAX;i_s++){
      h_Vars.back().push_back(std::vector<TH1D*>());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        std::string label = "h_Vars_"+categories.at(i_c)+"_"+sys_str.at(i_s)+"_"+std::to_string(i_u);
        h_Vars.back().back().push_back(new TH1D(label.c_str(),";Reco W (GeV);Events",50,0.95,5.0));
      }
    }
  }

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    LoadTreeFiltered(file,f_in,t_in,is_data_v.at(i_f),load_syst_v.at(i_f));

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 10000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);

      if(!*sel) continue;

      if(category == -1) std::cout << "Bad event" << std::endl;

      h_CV.at(category)->Fill(*var,weights.at(i_f));

      if(load_syst_v.at(i_f)){
        for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) h_Vars.at(category).at(kGenie).at(i_u)->Fill(*var,weights.at(i_f)*weightsGenie->at(i_u)/1000);
        for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) h_Vars.at(category).at(kReint).at(i_u)->Fill(*var,weights.at(i_f)*weightsReint->at(i_u)/1000);
        for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) h_Vars.at(category).at(kFlux).at(i_u)->Fill(*var,weights.at(i_f)*weightsFlux->at(i_u)/1000);
      }
      else {
        for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) h_Vars.at(category).at(kGenie).at(i_u)->Fill(W_lt,weights.at(i_f));
        for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) h_Vars.at(category).at(kReint).at(i_u)->Fill(W_lt,weights.at(i_f));
        for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) h_Vars.at(category).at(kFlux).at(i_u)->Fill(W_lt,weights.at(i_f));
      } 

    }

  }


  gSystem->Exec(("mkdir -p rootfiles/"+label).c_str());
  TFile* f_out = TFile::Open(("rootfiles/"+label+"/Histograms.root").c_str(),"RECREATE");

  for(size_t i_c=0;i_c<categories.size();i_c++)
    h_CV.at(i_c)->Write();

  for(size_t i_c=0;i_c<categories.size();i_c++)
    for(int i_s=0;i_s<kSystMAX;i_s++)
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++)
        h_Vars.at(i_c).at(i_s).at(i_u)->Write();        

  std::vector<std::vector<TH2D*>> h_Cov; 
  std::vector<std::vector<TH2D*>> h_FCov; 

  for(size_t i_c=0;i_c<categories.size();i_c++){
     h_Cov.push_back(std::vector<TH2D*>());
     h_FCov.push_back(std::vector<TH2D*>());
    for(int i_s=0;i_s<kSystMAX;i_s++){
      h_Cov.back().push_back(nullptr);
      h_FCov.back().push_back(nullptr);
      CalcCovMultisim(categories.at(i_c)+"_"+sys_str.at(i_s),h_CV.at(i_c),h_Vars.at(i_c).at(i_s),h_Cov.back().back(),h_FCov.back().back());
      h_Cov.back().back()->Write();
      h_FCov.back().back()->Write();
    }
  }       

  f_out->Close();

}
