#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "Systematics.h"
#include "EnergyEstimatorFuncs.h"

using namespace syst;

// Try using different combinations of cuts and frameworks to calculate W

void MakeCov(){

  std::string label = "RecoE_Test";
  std::string axis_title = ";Reco E (GeV);Events/GeV";
  TFile* f_tp = TFile::Open(("Analysis/"+label+"/rootfiles/BinningTemplate.root").c_str());
  TH1D* h_tp = (TH1D*)f_tp->Get("h_template");
  h_tp->SetDirectory(0);
  f_tp->Close();

  h_tp->GetXaxis()->SetTitle("Reco E (GeV)");
  h_tp->GetYaxis()->SetTitle("Events/GeV");

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

  //TH1D* h_CV_Tot = new TH1D("h_CV_Tot",axis_title.c_str(),bins.size()-1,&bins[0]);
  TH1D* h_CV_Tot = (TH1D*)h_tp->Clone("h_CV_Tot");
  h_CV_Tot->Sumw2();
  std::vector<std::vector<TH1D*>> h_Vars_Tot;
  for(int i_s=0;i_s<kSystMAX;i_s++){
    h_Vars_Tot.push_back(std::vector<TH1D*>());
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      std::string label = "h_Vars_Tot_"+sys_str.at(i_s)+"_"+std::to_string(i_u);
      //h_Vars_Tot.back().push_back(new TH1D(label.c_str(),axis_title.c_str(),bins.size()-1,&bins[0]));
      h_Vars_Tot.back().push_back((TH1D*)h_tp->Clone(label.c_str()));
    }
  }
  
  for(size_t i_c=0;i_c<categories.size();i_c++){
    //h_CV.push_back(new TH1D(("h_CV_"+categories.at(i_c)).c_str(),axis_title.c_str(),bins.size()-1,&bins[0]));
    h_CV.push_back((TH1D*)h_tp->Clone(("h_CV_"+categories.at(i_c)).c_str()));
    h_CV.back()->Sumw2();
    h_Vars.push_back(std::vector<std::vector<TH1D*>>());
    for(int i_s=0;i_s<kSystMAX;i_s++){
      h_Vars.back().push_back(std::vector<TH1D*>());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        std::string label = "h_Vars_"+categories.at(i_c)+"_"+sys_str.at(i_s)+"_"+std::to_string(i_u);
        //h_Vars.back().back().push_back(new TH1D(label.c_str(),axis_title.c_str(),bins.size()-1,&bins[0]));
        h_Vars.back().back().push_back((TH1D*)h_tp->Clone(label.c_str()));
      }
    }
  }

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    LoadTreeFiltered(file,f_in,t_in,is_data_v.at(i_f),load_syst_v.at(i_f));

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 100000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);

      if(!sel_h8) continue;
      double var = est_nu_e_h8->at(ee::kMuonKin);

      if(category == -1) std::cout << "Bad event" << std::endl;

      if(std::isnan(weightSpline) || std::isinf(weightSpline)) continue;
           
      h_CV.at(category)->Fill(var,weights.at(i_f)*weightSpline);
      if(!is_data) h_CV_Tot->Fill(var,weights.at(i_f)*weightSpline);

      if(load_syst_v.at(i_f)){
        for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) h_Vars_Tot.at(kGenie).at(i_u)->Fill(var,(double)weights.at(i_f)*weightsGenie->at(i_u)/1000);
        for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) h_Vars_Tot.at(kReint).at(i_u)->Fill(var,(double)weights.at(i_f)*weightsReint->at(i_u)*weightSpline/1000);
        for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) h_Vars_Tot.at(kFlux).at(i_u)->Fill(var,(double)weights.at(i_f)*weightsFlux->at(i_u)*weightSpline/1000);
      }
      else {
        for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) h_Vars_Tot.at(kGenie).at(i_u)->Fill(var,(double)weights.at(i_f)*weightSpline);
        for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) h_Vars_Tot.at(kReint).at(i_u)->Fill(var,(double)weights.at(i_f)*weightSpline);
        for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) h_Vars_Tot.at(kFlux).at(i_u)->Fill(var,(double)weights.at(i_f)*weightSpline);
      } 

      if(load_syst_v.at(i_f)){
        for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) h_Vars.at(category).at(kGenie).at(i_u)->Fill(var,(double)weights.at(i_f)*weightsGenie->at(i_u)/1000);
        for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) h_Vars.at(category).at(kReint).at(i_u)->Fill(var,(double)weights.at(i_f)*weightsReint->at(i_u)*weightSpline/1000);
        for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) h_Vars.at(category).at(kFlux).at(i_u)->Fill(var,(double)weights.at(i_f)*weightsFlux->at(i_u)*weightSpline/1000);
      }
      else {
        for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) h_Vars.at(category).at(kGenie).at(i_u)->Fill(var,(double)weights.at(i_f)*weightSpline);
        for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) h_Vars.at(category).at(kReint).at(i_u)->Fill(var,(double)weights.at(i_f)*weightSpline);
        for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) h_Vars.at(category).at(kFlux).at(i_u)->Fill(var,(double)weights.at(i_f)*weightSpline);
      } 

    }

  }

  gSystem->Exec(("mkdir -p Analysis/"+label+"/rootfiles/").c_str());
  TFile* f_out = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str(),"RECREATE");

  DivideByBinWidth(h_CV_Tot);
  h_CV_Tot->Write();

  for(size_t i_c=0;i_c<categories.size();i_c++){
    DivideByBinWidth(h_CV.at(i_c));
    h_CV.at(i_c)->Write();
   }

  for(int i_s=0;i_s<kSystMAX;i_s++){
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      DivideByBinWidth(h_Vars_Tot.at(i_s).at(i_u));
      h_Vars_Tot.at(i_s).at(i_u)->Write();     
    }
  }   

  for(size_t i_c=0;i_c<categories.size();i_c++){
    for(int i_s=0;i_s<kSystMAX;i_s++){
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        DivideByBinWidth(h_Vars.at(i_c).at(i_s).at(i_u)); 
        h_Vars.at(i_c).at(i_s).at(i_u)->Write();        
      }
    }
  }

  // Systematics covariance matrices
  std::vector<TH2D*> h_Cov_Tot; 
  std::vector<TH2D*> h_FCov_Tot; 
  for(int i_s=0;i_s<kSystMAX;i_s++){
    h_Cov_Tot.push_back(nullptr);
    h_FCov_Tot.push_back(nullptr);
    CalcCovMultisim(sys_str.at(i_s),h_CV_Tot,h_Vars_Tot.at(i_s),h_Cov_Tot.back(),h_FCov_Tot.back());
    h_Cov_Tot.back()->Write();
    h_FCov_Tot.back()->Write();
  }
  
  // Statistics covariance matrices
  TH2D* h_Cov_MCStat = (TH2D*)h_Cov_Tot.at(0)->Clone("Cov_MCStat");
  TH2D* h_FCov_MCStat = (TH2D*)h_Cov_Tot.at(0)->Clone("FCov_MCStat");
  TH2D* h_Cov_EstDataStat = (TH2D*)h_Cov_Tot.at(0)->Clone("Cov_EstDataStat");
  TH2D* h_FCov_EstDataStat = (TH2D*)h_Cov_Tot.at(0)->Clone("FCov_EstDataStat");
  h_Cov_MCStat->Reset();
  h_FCov_MCStat->Reset();
  h_Cov_EstDataStat->Reset();
  h_FCov_EstDataStat->Reset();
  for(int i=1;i<h_CV_Tot->GetNbinsX()+1;i++){
    h_Cov_MCStat->SetBinContent(i,i,h_CV_Tot->GetBinError(i)*h_CV_Tot->GetBinError(i));
    h_FCov_MCStat->SetBinContent(i,i,h_CV_Tot->GetBinError(i)*h_CV_Tot->GetBinError(i)/h_CV_Tot->GetBinContent(i)/h_CV_Tot->GetBinContent(i));
    h_Cov_EstDataStat->SetBinContent(i,i,h_CV_Tot->GetBinContent(i)/h_CV_Tot->GetBinWidth(i));
    h_FCov_EstDataStat->SetBinContent(i,i,h_CV_Tot->GetBinContent(i)*h_CV_Tot->GetBinContent(i)/h_Cov_EstDataStat->GetBinContent(i,i));
  } 
  h_Cov_MCStat->Write();
  h_FCov_MCStat->Write();
  h_Cov_EstDataStat->Write();
  h_FCov_EstDataStat->Write();

  std::vector<TH2D*> h_Cov_MCStat_Cat; 
  std::vector<TH2D*> h_FCov_MCStat_Cat; 
  std::vector<TH2D*> h_Cov_EstDataStat_Cat; 
  std::vector<TH2D*> h_FCov_EstDataStat_Cat; 
  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_Cov_MCStat_Cat.push_back((TH2D*)h_Cov_Tot.at(0)->Clone(("Cov_"+categories.at(i_c)+"_MCStat").c_str()));
    h_FCov_MCStat_Cat.push_back((TH2D*)h_Cov_Tot.at(0)->Clone(("FCov_"+categories.at(i_c)+"_MCStat").c_str()));
    h_Cov_EstDataStat_Cat.push_back((TH2D*)h_Cov_Tot.at(0)->Clone(("Cov_"+categories.at(i_c)+"_EstDataStat").c_str()));
    h_FCov_EstDataStat_Cat.push_back((TH2D*)h_Cov_Tot.at(0)->Clone(("FCov_"+categories.at(i_c)+"_EstDataStat").c_str()));
    h_Cov_MCStat_Cat.back()->Reset();  
    h_FCov_MCStat_Cat.back()->Reset();  
    h_Cov_EstDataStat_Cat.back()->Reset();  
    h_FCov_EstDataStat_Cat.back()->Reset();  
    for(int i=1;i<h_CV_Tot->GetNbinsX()+1;i++){
      h_Cov_MCStat_Cat.back()->SetBinContent(i,i,h_CV.at(i_c)->GetBinError(i)*h_CV.at(i_c)->GetBinError(i));
      h_FCov_MCStat_Cat.back()->SetBinContent(i,i,h_CV.at(i_c)->GetBinError(i)*h_CV.at(i_c)->GetBinError(i)/h_CV.at(i_c)->GetBinContent(i)/h_CV.at(i_c)->GetBinContent(i));
      h_Cov_EstDataStat_Cat.back()->SetBinContent(i,i,h_CV.at(i_c)->GetBinContent(i));
      h_FCov_EstDataStat_Cat.back()->SetBinContent(i,i,1.0/h_CV.at(i_c)->GetBinContent(i));
    } 
    h_Cov_MCStat_Cat.back()->Write();
    h_FCov_MCStat_Cat.back()->Write();
    h_Cov_EstDataStat_Cat.back()->Write();
    h_FCov_EstDataStat_Cat.back()->Write();
  }

  TH2D* h_Cov = static_cast<TH2D*>(h_Cov_Tot.at(0)->Clone("Cov"));
  TH2D* h_FCov = static_cast<TH2D*>(h_FCov_Tot.at(0)->Clone("FCov"));
  h_Cov->Reset();
  h_FCov->Reset();
 
  for(int i_s=0;i_s<kSystMAX;i_s++){
    h_Cov->Add(h_Cov_Tot.at(i_s));    
    h_FCov->Add(h_FCov_Tot.at(i_s));    
  }
  h_Cov->Add(h_Cov_MCStat);
  h_FCov->Add(h_FCov_MCStat);
  h_Cov->Write();
  h_FCov->Write();

  std::vector<std::vector<TH2D*>> h_Cov_Cat; 
  std::vector<std::vector<TH2D*>> h_FCov_Cat; 
  for(size_t i_c=0;i_c<categories.size();i_c++){
     h_Cov_Cat.push_back(std::vector<TH2D*>());
     h_FCov_Cat.push_back(std::vector<TH2D*>());
    for(int i_s=0;i_s<kSystMAX;i_s++){
      h_Cov_Cat.back().push_back(nullptr);
      h_FCov_Cat.back().push_back(nullptr);
      CalcCovMultisim(categories.at(i_c)+"_"+sys_str.at(i_s),h_CV.at(i_c),h_Vars.at(i_c).at(i_s),h_Cov_Cat.back().back(),h_FCov_Cat.back().back());
      h_Cov_Cat.back().back()->Write();
      h_FCov_Cat.back().back()->Write();
    }
  }       

  f_out->Close();

}
