#include "TLorentzVector.h"
#pragma link C++ class std::vector<TLorentzVector>+;
#include "Funcs.h"
#include "WeightFuncs.h"
#include "EnergyEstimatorFuncs.h"

void GENIEToFlat(){

  gInterpreter->GenerateDictionary("std::vector<TLorentzVector>","vector;TLorentzVector.h");

  TFile* f_in = TFile::Open("/exp/uboone/data/users/cthorpe/DIS/Generators/Unprocessed/genie_events_numu_Ar40.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("gst")) ;

  Double_t        wght;
  Int_t           neu;
  Bool_t          cc;
  Double_t        Ev;
  Int_t           fspl;
  Double_t        El;
  Double_t        pxl;
  Double_t        pyl;
  Double_t        pzl;
  Int_t           nf;
  Int_t           pdgf[100];   //[nf]
  Double_t        Ef[100];   //[nf]
  Double_t        pxf[100];   //[nf]
  Double_t        pyf[100];   //[nf]
  Double_t        pzf[100];   //[nf]
  Int_t           ni;
  Int_t           pdgi[100];   //[ni]
  Double_t        Ei[100];   //[ni]
  Double_t        pxi[100];   //[ni]
  Double_t        pyi[100];   //[ni]
  Double_t        pzi[100];   //[ni]

  t_in->SetBranchAddress("wght",&wght);
  t_in->SetBranchAddress("neu",&neu);
  t_in->SetBranchAddress("cc",&cc);
  t_in->SetBranchAddress("Ev",&Ev);
  t_in->SetBranchAddress("fspl",&fspl);
  t_in->SetBranchAddress("El",&El);
  t_in->SetBranchAddress("pxl",&pxl);
  t_in->SetBranchAddress("pyl",&pyl);
  t_in->SetBranchAddress("pzl",&pzl);
  t_in->SetBranchAddress("nf",&nf);
  t_in->SetBranchAddress("pdgf",pdgf);
  t_in->SetBranchAddress("Ef",Ef);
  t_in->SetBranchAddress("pxf",pxf);
  t_in->SetBranchAddress("pyf",pyf);
  t_in->SetBranchAddress("pzf",pzf);
  t_in->SetBranchAddress("ni",&ni);
  t_in->SetBranchAddress("pdgi",pdgi);
  t_in->SetBranchAddress("Ei",Ei);
  t_in->SetBranchAddress("pxi",pxi);
  t_in->SetBranchAddress("pyi",pyi);
  t_in->SetBranchAddress("pzi",pzi);

  TFile* f_out = new TFile("GENIEEvents.root","RECREATE");
  TTree* t_out = new TTree("eventtree","eventtree");

  Double_t scale = 1;
  Double_t weight;
  Double_t nu_e;
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
  t_out->Branch("weight",&weight);
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
    weight = wght;
    nu_e = Ev;
    ccnc = cc;
    nu_pdg = neu;

    lepton_pdg = fspl;
    lepton_p4 = TLorentzVector(pxl,pyl,pzl,El);

    for(size_t i_p=0;i_p<nf;i_p++){
      pdg.push_back(pdgf[i_p]);
      p4.push_back(TLorentzVector(pxf[i_p],pyf[i_p],pzf[i_p],Ef[i_p]));
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
      // Pi0s are used as single shower proxies since GENIE GST does not track photon daughters
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
