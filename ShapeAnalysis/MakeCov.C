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
#include "Histograms.h"

void MakeCov(){

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root"
  };

  hist::HistogramManager h_MuonMom("MuonMom");
  h_MuonMom.LoadTemplate();
  h_MuonMom.DBBW();
  h_MuonMom.ShapeOnly();

  hist::HistogramManager h_MuonCosTheta("MuonCosTheta");
  h_MuonCosTheta.LoadTemplate();
  h_MuonCosTheta.DBBW();
  h_MuonCosTheta.ShapeOnly();

  hist::HistogramManager h_ProtonE("ProtonE");
  h_ProtonE.LoadTemplate();
  h_ProtonE.DBBW();
  h_ProtonE.ShapeOnly();

  hist::HistogramManager h_PionE("PionE");
  h_PionE.LoadTemplate();
  h_PionE.DBBW();
  h_PionE.ShapeOnly();

  hist::HistogramManager h_PiZeroE("PiZeroE");
  h_PiZeroE.LoadTemplate();
  h_PiZeroE.DBBW();
  h_PiZeroE.ShapeOnly();

  hist::HistogramManager h_NProt("NProt");
  h_NProt.SetTemplate(";Num of FS Protons;",6,-0.5,5.5);
  h_NProt.ShapeOnly(); 

  hist::HistogramManager h_NPi("NPi");
  h_NPi.SetTemplate(";Num of FS Pions;",6,-0.5,5.5);
  h_NPi.ShapeOnly(); 

  hist::HistogramManager h_NShr("NShr");
  h_NShr.SetTemplate(";Num of Showers;",6,-0.5,5.5);
  h_NShr.ShapeOnly(); 

  hist::HistogramManager h_W("W");
  h_W.LoadTemplate();
  h_W.DBBW();
  h_W.ShapeOnly();  
 
  std::vector<hist::HistogramManager> h_NuE_v;
  for(int i_e=0;i_e<ee::kMAX;i_e++){
    h_NuE_v.push_back(hist::HistogramManager(ee::estimators_str.at(i_e)));
    h_NuE_v.back().LoadTemplate();
    h_NuE_v.back().DBBW();
    h_NuE_v.back().ShapeOnly();
  }

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //std::cout << ievent << std::endl;

      //if(ievent > 10000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);

      if(!sel_h8) continue;

      h_MuonMom.FillHistograms(muon_mom_h8->Mag(),load_syst);
      h_MuonCosTheta.FillHistograms(muon_mom_h8->CosTheta(),load_syst);
      h_ProtonE.FillHistograms(proton_p4_h8->E() - nprot_h8*Mp,load_syst);
      if(npi_h8 > 0) h_PionE.FillHistograms(pion_p4_h8->E(),load_syst);
      if(nsh_h8 > 0) h_PiZeroE.FillHistograms(shower_p4_h8->E(),load_syst);
  
      h_NProt.FillHistograms(nprot_h8,load_syst);
      h_NPi.FillHistograms(npi_h8,load_syst);
      h_NShr.FillHistograms(nsh_h8,load_syst);

      h_W.FillHistograms(W_h8,load_syst);

      for(int i_e=0;i_e<ee::kMAX;i_e++)
        h_NuE_v.at(i_e).FillHistograms(est_nu_e_h8->at(i_e),load_syst);

    }

  }

  h_MuonMom.Write();
  h_MuonCosTheta.Write();
  h_ProtonE.Write();
  h_PionE.Write();
  h_PiZeroE.Write();

  h_NProt.Write();
  h_NPi.Write();
  h_NShr.Write();

  h_W.Write();

  for(int i_e=0;i_e<ee::kMAX;i_e++)
    h_NuE_v.at(i_e).Write();
}
