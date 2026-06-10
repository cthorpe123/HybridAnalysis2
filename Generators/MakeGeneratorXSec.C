#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "BranchList.h"
#include "Systematics.h"
#include "EnergyEstimatorFuncs.h"
#include "BinningFuncs.h"

using namespace syst;
using namespace binning;

// Tune binning so all bins have data FE of this value

void MakeGeneratorXSec(){

  std::map<std::string,TH1D*> h_m;

  h_m["MuonMom"]                = new TH1D("h_MuonMom",                ";Muon Momentum [GeV];Events/GeV",              100, 0,   2);
  h_m["MuonCosTheta"]           = new TH1D("h_MuonCosTheta",           ";Muon cos#theta;Events/Unit",                   100,-1,   1);

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Generators/";
  std::vector<std::string> files_v = {
    "GENIEEvents.root"
  };

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    LoadGeneratorTree(file,f_in,t_in);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 10000) break;
      if (ievent % 1000 == 0) std::cout << "  " << ievent << " / " << t_in->GetEntries() << "\r" << std::flush;
      t_in->GetEntry(ievent);

      if(is_signal_t){
        for(const auto &item : h_m){
          std::string var = item.first;
          if(vars_t->find(var) == vars_t->end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
          h_m.at(var)->Fill(vars_t->at(var),scale*1e38*40);
        }
      }

    }

  }


  TCanvas* c = new TCanvas("c","c");
  for(const auto &item : h_m){
    std::string var = item.first;
    h_m.at(var)->Draw("HIST");
    h_m.at(var)->GetYaxis()->SetTitle("d#sigma/dvar (10^{-38} cm^{2}/unit)");
    c->Print((var+".png").c_str());
  }



}
