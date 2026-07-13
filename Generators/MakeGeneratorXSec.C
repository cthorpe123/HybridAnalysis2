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

void MakeGeneratorXSec(){

  bool load_asimov = true;
  std::vector<std::string> vars = {"Enu","MuonMom","MuonCosTheta","Norm"};
  bool dbbw = true;
  
  std::vector<std::string> generators = {"v3.0.6","NuWro"};

  std::map<std::string, std::map<std::string, TH1D*>> h_m;
  std::map<std::string,std::map<std::string,TH2D*>> h_m_2d;

  for(const std::string& var : vars){
    for(const std::string& gen : generators){
      TFile* f_tp_truth = TFile::Open(("Analysis/"+var+"/rootfiles/TruthBinningTemplate.root").c_str());
      h_m[var][gen] = (TH1D*)f_tp_truth->Get("h_template_All")->Clone(("h_xsec_"+var+"_"+gen).c_str());
      auto* xbins = h_m[var][gen]->GetXaxis()->GetXbins()->GetArray();
      h_m_2d[var][gen] = new TH2D(("h_xsec_2D_"+var+"_"+gen).c_str(),";;True Neutrino Energy (GeV)",h_m[var][gen]->GetNbinsX(),xbins,100,0.0,3.0);
      h_m[var][gen]->SetDirectory(0);
      h_m_2d[var][gen]->SetDirectory(0);
      f_tp_truth->Close();
    }
  }

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Generators/";
  std::vector<std::string> files_v = {
    "14_1000180400_CC_v3_0_6_G18_10a_02_11aEvents.root",
    "NuWroEvents.root"
  };

  for(int i_f=0;i_f<(int)files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);
    std::string gen = generators.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    LoadGeneratorTree(file,f_in,t_in);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      if (ievent % 1000 == 0) std::cout << "  " << ievent << " / " << t_in->GetEntries() << "\r" << std::flush;
      t_in->GetEntry(ievent);

      vars_t->emplace("Norm",0.5);
      vars_t->emplace("Enu",nu_e);

      if(is_signal_t){
        for(const std::string& var : vars){
          if(vars_t->find(var) == vars_t->end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
          h_m.at(var).at(gen)->Fill(vars_t->at(var), gen_weight*scale*1e38*40);
          h_m_2d.at(var).at(gen)->Fill(vars_t->at(var),nu_e,gen_weight*scale*1e38*40);
        }
      }

    }

    f_in->Close();
  }

  for(const std::string& var : vars){
    TFile* f_out = TFile::Open(("Analysis/"+var+"/rootfiles/GeneratorXSec.root").c_str(), "RECREATE");
    for(const std::string& gen : generators){
      h_m.at(var).at(gen)->GetYaxis()->SetTitle("d#sigma (10^{-38} cm^{2})");
      f_out->cd();
      h_m.at(var).at(gen)->Write();
      h_m_2d.at(var).at(gen)->Write();
    }
    f_out->Close();
  }

  for(const std::string& var : vars){
    gSystem->Exec(("mkdir -p Analysis/"+var+"/Plots/MakeGeneratorXSec").c_str());
    std::vector<TH1D*> h_v;
    std::vector<std::string> legs = generators;
    std::vector<int> cols;
    int ctr = 0;
    for(const std::string& gen : generators){
      if(dbbw) DivideByBinWidth(h_m.at(var).at(gen));
      h_v.push_back(h_m.at(var).at(gen));
      cols.push_back(ctr+2);
      ctr++;
    }

    TH1D* h_asimov = nullptr;
    if(load_asimov){
      TFile* f_in_hist = TFile::Open(("Analysis/"+var+"/rootfiles/Histograms.root").c_str());
      const double POT = ((TH1D*)f_in_hist->Get("Meta/POT"))->GetBinContent(1);
      hist::MultiChannelHistogramManager mchm(var,true);
      mchm.LoadTemplates();
      h_asimov = (TH1D*)f_in_hist->Get("Truth/CV/h_Signal");
      h_asimov->SetDirectory(0);
      CrossSectionH(h_asimov,POT);
      mchm.Restore(h_asimov,"All",true);
      f_in_hist->Close();
      if(dbbw) DivideByBinWidth(h_asimov);
    }
    if(h_asimov != nullptr){
      h_v.push_back(h_asimov);
      legs.push_back("MicroBooNE");
      cols.push_back(kBlack);
    }
    pfs::DrawUnstacked(h_v, cols, legs,true,true,false,"Analysis/"+var+"/Plots/MakeGeneratorXSec/GeneratorXSec.png");

    // Shape comparison: normalise each histogram to unit area
    std::vector<TH1D*> h_shape_v;
    for(TH1D* h : h_v){
      TH1D* h_norm = (TH1D*)h->Clone((std::string(h->GetName())+"_shape").c_str());
      h_norm->SetDirectory(0);
      const double integral = h_norm->Integral();
      if(integral > 0) h_norm->Scale(1.0/integral);
      h_norm->GetYaxis()->SetTitle("Normalised");
      h_shape_v.push_back(h_norm);
    }
    pfs::DrawUnstacked(h_shape_v, cols, legs,true,true,false,"Analysis/"+var+"/Plots/MakeGeneratorXSec/GeneratorXSec_Shape.png");
    for(TH1D* h : h_shape_v) delete h;
  }

}
