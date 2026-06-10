#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;
#include "Funcs.h"
#include "WeightFuncs.h"
#include "EnergyEstimatorFuncs.h"

const std::string generator = "GENIE";

// GiBUU ntuple is made of many files - scaling branch needs to be corrected
// for this, use this variable to do so
const int num_files = 1;

void NuisanceToFlat(){

  gInterpreter->GenerateDictionary("std::vector<TLorentzVector>","vector;TLorentzVector.h");

  TFile* f_in = TFile::Open("/exp/uboone/data/users/cthorpe/DIS/Generators/Unprocessed/GENIE.flat.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("FlatTree_VARS")) ;

  Char_t          cc;
  Int_t           PDGnu;
  Float_t         Enu_true;
  Int_t           nfsp;
  Float_t         px[100];   //[nfsp]
  Float_t         py[100];   //[nfsp]
  Float_t         pz[100];   //[nfsp]
  Float_t         E[100];   //[nfsp]
  Int_t           pdg_in[100];   //[nfsp]
  Float_t         Weight;
  Double_t        fScaleFactor;

  t_in->SetBranchAddress("cc", &cc);
  t_in->SetBranchAddress("PDGnu", &PDGnu);
  t_in->SetBranchAddress("Enu_true", &Enu_true);
  t_in->SetBranchAddress("nfsp", &nfsp);
  t_in->SetBranchAddress("px", px);
  t_in->SetBranchAddress("py", py);
  t_in->SetBranchAddress("pz", pz);
  t_in->SetBranchAddress("E", E);
  t_in->SetBranchAddress("pdg", pdg_in);
  t_in->SetBranchAddress("Weight", &Weight);
  t_in->SetBranchAddress("fScaleFactor", &fScaleFactor);

  TFile* f_out = new TFile((generator+"Events.root").c_str(),"RECREATE");
  TTree* t_out = new TTree("eventtree","eventtree");

  Double_t scale;
  Double_t gen_weight;
  Float_t nu_e;
  Int_t ccnc;
  Int_t nu_pdg;

  Int_t lepton_pdg;
  TLorentzVector lepton_p4;

  std::vector<int> pdg;
  std::vector<TLorentzVector> p4;

  bool is_signal_t;
  bool has_muon_t;
  std::map<std::string,double> vars_t_map;

  TVector3 muon_mom_t;
  int nprot_t, npi_t, npi0_t, nsh_t;
  TLorentzVector proton_p4_t, pion_p4_t, pi0_p4_t, gamma_p4_t;
  std::vector<TLorentzVector> protons_t, pions_t, pi0s_t, gammas_t;
  double W_t;
  int ch_t;

  t_out->Branch("scale",&scale);
  t_out->Branch("gen_weight",&gen_weight);
  t_out->Branch("nu_e",&nu_e);
  t_out->Branch("nu_pdg",&nu_pdg);
  t_out->Branch("ccnc",&ccnc);
  t_out->Branch("lepton_pdg",&lepton_pdg);
  t_out->Branch("lepton_p4",&lepton_p4);
  t_out->Branch("pdg",&pdg);
  t_out->Branch("p4",&p4,32000,0);
  t_out->Branch("is_signal_t",&is_signal_t);
  t_out->Branch("has_muon_t",&has_muon_t);
  t_out->Branch("vars_t",&vars_t_map);
  t_out->Branch("muon_mom_t",&muon_mom_t);
  t_out->Branch("nprot_t",&nprot_t);
  t_out->Branch("npi_t",&npi_t);
  t_out->Branch("npi0_t",&npi0_t);
  t_out->Branch("nsh_t",&nsh_t);
  t_out->Branch("proton_p4_t",&proton_p4_t);
  t_out->Branch("pion_p4_t",&pion_p4_t);
  t_out->Branch("pi0_p4_t",&pi0_p4_t);
  t_out->Branch("gamma_p4_t",&gamma_p4_t);
  t_out->Branch("protons_t",&protons_t,32000,0);
  t_out->Branch("pions_t",&pions_t,32000,0);
  t_out->Branch("pi0s_t",&pi0s_t,32000,0);
  t_out->Branch("gammas_t",&gammas_t,32000,0);
  t_out->Branch("W_t",&W_t);
  t_out->Branch("ch_t",&ch_t);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){
    //if(ievent > 10000) break;
    if(ievent % 10000 == 0) std::cout << "Event " << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    p4.clear();
    pdg.clear();

    scale = fScaleFactor / num_files;
    gen_weight = Weight;
    nu_e = Enu_true;
    ccnc = cc;
    nu_pdg = PDGnu;

    lepton_pdg = 0;
    lepton_p4 = TLorentzVector(0,0,0,0);

    int nlep = 0;
    for(size_t i_p=0;i_p<nfsp;i_p++){
      if(abs(pdg_in[i_p]) == 13 || abs(pdg_in[i_p]) == 11){
        lepton_pdg = pdg_in[i_p];
        lepton_p4 = TLorentzVector(px[i_p],py[i_p],pz[i_p],E[i_p]);
        nlep++;
      }
      else {
        pdg.push_back(pdg_in[i_p]);
        p4.push_back(TLorentzVector(px[i_p],py[i_p],pz[i_p],E[i_p]));
      }
    }

    if(nlep > 1){
      std::cout << "Event has more than one FS lepton " << nlep << std::endl;
      continue;
    }

    // Truth-level final state kinematics
    muon_mom_t = TVector3(0,0,0);
    nprot_t = 0; npi_t = 0; npi0_t = 0; nsh_t = 0;
    proton_p4_t = TLorentzVector(0,0,0,0);
    pion_p4_t   = TLorentzVector(0,0,0,0);
    pi0_p4_t    = TLorentzVector(0,0,0,0);
    gamma_p4_t  = TLorentzVector(0,0,0,0);
    protons_t.clear(); pions_t.clear(); pi0s_t.clear(); gammas_t.clear();

    if(abs(lepton_pdg) == 13){
      muon_mom_t = lepton_p4.Vect();
    }

    for(size_t i_p=0;i_p<p4.size();i_p++){
      TVector3 mom = p4.at(i_p).Vect();
      if(pdg.at(i_p) == 2212 && mom.Mag() > thresholds.at(2212).first){
        nprot_t++;
        proton_p4_t += p4.at(i_p);
        protons_t.push_back(p4.at(i_p));
      }
      if(abs(pdg.at(i_p)) == 211 && mom.Mag() > thresholds.at(211).first){
        npi_t++;
        pion_p4_t += p4.at(i_p);
        pions_t.push_back(p4.at(i_p));
      }
      // Pi0s used as single shower proxies since NUISANCE flat tree does not track photon daughters
      if(pdg.at(i_p) == 111){
        npi0_t++;
        pi0_p4_t += p4.at(i_p);
        pi0s_t.push_back(p4.at(i_p));
        if(mom.Mag() > thresholds.at(22).first){
          nsh_t++;
          gamma_p4_t += p4.at(i_p);
          gammas_t.push_back(p4.at(i_p));
        }
      }
    }

    SortTLorentzVector(protons_t);
    SortTLorentzVector(pions_t);
    SortTLorentzVector(pi0s_t);
    SortTLorentzVector(gammas_t);

    W_t = (proton_p4_t + pion_p4_t + gamma_p4_t).M();
    has_muon_t = abs(lepton_pdg) == 13 && muon_mom_t.Mag() > thresholds.at(13).first;
    is_signal_t = abs(nu_pdg) == 14 && ccnc == 1 && has_muon_t && nprot_t > 0;

    ch_t = -1;
    std::string ch_str_t = channel_str(nprot_t,npi_t,nsh_t);
    for(size_t i_ch=0;i_ch<channels.size();i_ch++){
      if(ch_str_t == channels.at(i_ch)){ ch_t = i_ch; break; }
    }

    std::vector<double> est_nu_e_t(ee::kMAX,-1);
    if(is_signal_t){
      TLorentzVector plepton_t(muon_mom_t.X(),muon_mom_t.Y(),muon_mom_t.Z(),sqrt(muon_mom_t.Mag2()+ml*ml));
      est_nu_e_t = ee::GetEnergyEst(plepton_t,W_t,proton_p4_t,nprot_t,pion_p4_t,npi_t,gamma_p4_t,npi0_t);
    }

    vars_t_map = {
      {"MuonMom",muon_mom_t.Mag()},
      {"MuonCosTheta",muon_mom_t.CosTheta()},
      {"LeadProtonKE",-1000},
      {"LeadPionE",-1000},
      {"1p1piOpeningAngle",-1000},
      {"1p1piAsym",-1000},
      {"MuonProtonOpeningAngle",-1000},
      {"2pOpeningAngle",-1000},
      {"2pAsym",-1000},
      {"2shwOpenAngle",-1000},
      {"2shwAsym",-1000},
      {"NProt",(double)nprot_t},
      {"NPi",(double)npi_t},
      {"NSh",(double)nsh_t},
      {"NPi0",(double)npi0_t},
      {"ProtonKE",proton_p4_t.E()-nprot_t*Mp},
      {"PionE",pion_p4_t.E()},
      {"PiZeroE",gamma_p4_t.E()},
      {"W",W_t},
      {"Channel",(double)ch_t}
    };
    for(int i_e=0;i_e<ee::kMAX;i_e++) vars_t_map[ee::estimators_str.at(i_e)] = est_nu_e_t.at(i_e);

    if(is_signal_t){
      vars_t_map.at("LeadProtonKE") = protons_t.at(0).E() - Mp;
      vars_t_map.at("MuonProtonOpeningAngle") = 180/3.142*muon_mom_t.Angle(protons_t.at(0).Vect());
      if(nprot_t == 2){
        vars_t_map.at("2pOpeningAngle") = 180/3.142*protons_t.at(0).Vect().Angle(protons_t.at(1).Vect());
        vars_t_map.at("2pAsym") = weight::Asymmetry3({protons_t.at(0)},{protons_t.at(1)});
      }
      if(nsh_t == 2){
        vars_t_map.at("2shwOpenAngle") = 180/3.142*gammas_t.at(0).Vect().Angle(gammas_t.at(1).Vect());
        vars_t_map.at("2shwAsym") = weight::Asymmetry3({gammas_t.at(0)},{gammas_t.at(1)});
      }
      if(npi_t > 0) vars_t_map.at("LeadPionE") = pions_t.at(0).E();
      if(nprot_t == 1 && npi_t == 1){
        vars_t_map.at("1p1piOpeningAngle") = 180/3.142*protons_t.at(0).Vect().Angle(pions_t.at(0).Vect());
        vars_t_map.at("1p1piAsym") = weight::Asymmetry3({protons_t.at(0)},{pions_t.at(0)});
      }
    }

    t_out->Fill();

  }

  t_out->Write();
  f_out->Close();

  f_in->Close();

}
