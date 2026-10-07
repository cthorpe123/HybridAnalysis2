#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "BranchList.h"
#include "Systematics.h"
#include "EnergyEstimatorFuncs.h"
#include "MultiChannelHistograms.h"
#include "WeightFuncs.h"

void MakeCovNuWro(){

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/retupled/";

  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  std::vector<std::string> vars = {"MuonMom"};
  //std::vector<std::string> vars = var_names;
  std::vector<std::string> int_vars = {"NProt","NPi","NSh","NPi0"};

  std::map<std::string,hist::MultiChannelHistogramManager> h_m;
  for(std::string var : vars){
    h_m.emplace(var,hist::MultiChannelHistogramManager(var,true));
    h_m.at(var).SetTrueChannelList(channels_t);
    h_m.at(var).SetRecoChannelList(channels_r);
    h_m.at(var).KeepAll();
    if(!in_vec(int_vars,var)) h_m.at(var).LoadTemplates();
    else h_m.at(var).SetTemplates("",4,-0.5,3.5,4,-0.5,3.5); 
    h_m.at(var).MakeHM();
  } 
  
  for(std::string var : vars)
    h_m.at(var).AddSpecialUniv("NuWro_0");

  // then analyse the nuwro files
  std::vector<std::string> files_v = {
    "run4c/Filtered_Merged_checkout_MCC9.10_Run45_v10_04_07_23_BNB_nuwro_overlay_surprise_reco2_hist_4c.root",
    "run5/Filtered_Merged_checkout_MCC9.10_Run45_v10_04_07_23_BNB_nuwro_overlay_surprise_reco2_hist_5.root"
  };

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 20000) break;
      if (ievent % 1000 == 0) std::cout << "  " << ievent << " / " << t_in->GetEntries() << "\r" << std::flush;
      t_in->GetEntry(ievent);
      
      std::string channel_t = "All";
      std::string channel_h8 = "All";

      weightSplineTimesTune = 1.0;

      // Fill the histograms
      for(const auto &item : h_m){
        std::string var = item.first;
        if(vars_t->find(var) == vars_t->end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
        if(vars_h8->find(var) == vars_h8->end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
        const double& t = vars_t->at(var);
        const double& r = vars_h8->at(var);
        h_m.at(var).FillSpecialHistograms2D("NuWro_0",is_signal_t,sel_h8,t,r,1.0,channel_t,channel_h8);
      }
    }

  }

  for(const auto &item : h_m){
    std::string var = item.first;
    h_m.at(var).Write("NuWroFD.root");
  }

  // Normalise the NuWro plots to the same as the truth CV from the 
  // main histograms file
  // Append the result of the NuWro FD to the main Histograms.root file
  for(const auto &item : h_m){
    std::string var = item.first;
    TFile* f_nuwro = TFile::Open((AnalysisDir()+"/"+var+"/rootfiles/NuWroFD.root").c_str());
    TFile* f_hist = TFile::Open((AnalysisDir()+"/"+var+"/rootfiles/Histograms.root").c_str(),"UPDATE");
 
    // Make new dirs in the main file if needed
    f_hist->cd();
    if(f_hist->GetDirectory("Truth/Special") == nullptr){ 
      f_hist->mkdir("Truth/Special");
      f_hist->mkdir("Reco/Special");
      f_hist->mkdir("Response/Special");
      f_hist->mkdir("Joint/Special");

    }

    if(f_hist->GetDirectory("Truth/Special/NuWro_0") == nullptr){ 
      f_hist->mkdir("Truth/Special/NuWro_0");
      f_hist->mkdir("Reco/Special/NuWro_0");
      f_hist->mkdir("Joint/Special/NuWro_0");
      f_hist->mkdir("Response/Special/NuWro_0");
    }
    
    TH1D* h_Truth = (TH1D*)f_nuwro->Get("Truth/Special/NuWro_0/h_Signal");
    TH1D* h_Reco = (TH1D*)f_nuwro->Get("Reco/Special/NuWro_0/h_Signal");
    TH2D* h_Joint = (TH2D*)f_nuwro->Get("Joint/Special/NuWro_0/h_Signal");
    TH2D* h_Response = (TH2D*)f_nuwro->Get("Response/Special/NuWro_0/h_Signal");

    // Make clones of them to normalise to the CV truth
    TH1D* h_Truth_Norm = (TH1D*)h_Truth->Clone("h_Truth_Norm");
    TH1D* h_Reco_Norm = (TH1D*)h_Reco->Clone("h_Reco_Norm");
    TH2D* h_Joint_Norm = (TH2D*)h_Joint->Clone("h_Reco_Norm");

    TH1D* h_CV_Truth = (TH1D*)f_hist->Get("Truth/CV/h_Signal");
    double scale = IntegralWithOU(h_CV_Truth)/IntegralWithOU(h_Truth_Norm);
    h_Truth_Norm->Scale(scale);
    h_Reco_Norm->Scale(scale);
    h_Joint_Norm->Scale(scale);

    // Write everything
    f_hist->cd();
    f_hist->cd("Truth/Special/NuWro_0");
    h_Truth_Norm->Write("h_Signal",TObject::kOverwrite);
    h_Truth->Write("h_Signal_NoNorm",TObject::kOverwrite);

    f_hist->cd();
    f_hist->cd("Reco/Special/NuWro_0");
    h_Reco_Norm->Write("h_Signal",TObject::kOverwrite);
    h_Reco->Write("h_Signal_NoNorm",TObject::kOverwrite);

    f_hist->cd();
    f_hist->cd("Joint/Special/NuWro_0");
    h_Joint_Norm->Write("h_Signal",TObject::kOverwrite);
    h_Joint->Write("h_Signal_NoNorm",TObject::kOverwrite);

    f_hist->cd();
    f_hist->cd("Response/Special/NuWro_0");
    h_Response->Write("h_Signal",TObject::kOverwrite);
    h_Response->Write("h_Signal_NoNorm",TObject::kOverwrite);

    f_hist->Close();
    f_nuwro->Close();

  }


}