#ifndef _WC_Funcs_h_
#define _WC_Funcs_h_

#include "Funcs.h"

namespace wc {

// Simple selection for finding muon 
int SimpleMuonSelection(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4], const Float_t endxyzt[][4]){
  //double min_E = sqrt(thresholds.at(13).first*thresholds.at(13).first+ml*ml);
  double min_E = 0;
  double max_E = sqrt(thresholds.at(13).second*thresholds.at(13).second+ml*ml);
  int longest_idx = -1;
  int n = 0;
  for(size_t i_tr=0;i_tr<ntracks;i_tr++){
    if(mother[i_tr] == 0 && abs(pdg[i_tr]) == 13 && mom[i_tr][3] > min_E){
      min_E = mom[i_tr][3];
      longest_idx = i_tr;
      n++;
    }
  }
  //if(n > 1) return -1;
  return longest_idx;
}

int SimpleNueSelection(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4], const Float_t endxyzt[][4]){
  double min_E = 0;
  int longest_idx = -1;
  int n = 0;
  for(size_t i_tr=0;i_tr<ntracks;i_tr++){
    if(mother[i_tr] == 0 && abs(pdg[i_tr]) == 11 && mom[i_tr][3] > min_E){
      min_E = mom[i_tr][3];
      longest_idx = i_tr;
      n++;
    }
  }
  //if(n > 1) return -1;
  return longest_idx;
}

// Simple selection for finding leading charged pion 
int SimplePionSelection(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4]){
  double min_E = sqrt(thresholds.at(211).first*thresholds.at(211).first+mpi*mpi);
  int longest_idx = -1;
  for(size_t i_tr=0;i_tr<ntracks;i_tr++){
    if(mother[i_tr] == 0 && abs(pdg[i_tr]) == 211 && mom[i_tr][3] > min_E){
      min_E = mom[i_tr][3];
      longest_idx = i_tr;
    }
  }
  return longest_idx;
}

// Simple selection for finding proton 
int SimpleProtonSelection(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4]){
  double min_E = sqrt(thresholds.at(2212).first*thresholds.at(2212).first+Mp*Mp);
  int longest_idx = -1;
  for(size_t i_tr=0;i_tr<ntracks;i_tr++){
    if(mother[i_tr] == 0 && abs(pdg[i_tr]) == 2212 && mom[i_tr][3] > min_E){
      min_E = mom[i_tr][3];
      longest_idx = i_tr;
    }
  }
  return longest_idx;
}

double CalcTrueW(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4],const Int_t id[]){

  // Calculating the true invariant mass
  TLorentzVector p4_tot(0,0,0,0); 
  std::vector<Int_t> pi0_ids;  
  for(size_t i_p=0;i_p<ntracks;i_p++){
    if(mother[i_p] != 0) continue;
    if(abs(pdg[i_p]) == 13) continue;
    if(pdg[i_p] == 111){
      pi0_ids.push_back(id[i_p]);
      continue;
    }

    double  p = TVector3(mom[i_p][0],mom[i_p][1],mom[i_p][2]).Mag();

    if(thresholds.find(abs(pdg[i_p])) != thresholds.end() && p > thresholds.at(abs(pdg[i_p])).first && p < thresholds.at(abs(pdg[i_p])).second)
      p4_tot += TLorentzVector(mom[i_p][0],mom[i_p][1],mom[i_p][2],mom[i_p][3]);
  }

  for(size_t i_pi0=0;i_pi0<pi0_ids.size();i_pi0++){
    for(size_t i_p=0;i_p<ntracks;i_p++){
      if(mother[i_p] != pi0_ids.at(i_pi0) || pdg[i_p] != 22) continue;
      p4_tot += TLorentzVector(mom[i_p][0],mom[i_p][1],mom[i_p][2],mom[i_p][3]);
    }
  } 

  return p4_tot.M();

}


std::vector<TLorentzVector> RecoProton4MomV(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4],const Float_t endpos[][4],const Int_t id[],int muon_idx){
  std::vector<TLorentzVector> p4; 
  for(size_t i_p=0;i_p<ntracks;i_p++){
    if(mother[i_p] != 0) continue;
    if(i_p == muon_idx || pdg[i_p] != 2212) continue;
    p4.push_back(TLorentzVector(mom[i_p][0],mom[i_p][1],mom[i_p][2],mom[i_p][3]));
  }
  SortTLorentzVector(p4); 
  return p4;
}

TLorentzVector RecoProton4Mom(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4],const Float_t endpos[][4],const Int_t id[],int muon_idx){
  return SumTLorentzVector(RecoProton4MomV(ntracks,pdg,mother,mom,endpos,id,muon_idx));
}

std::vector<TLorentzVector> RecoPion4MomV(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4],const Float_t endpos[][4],const Int_t id[],int muon_idx){
  std::vector<TLorentzVector> p4; 
  for(size_t i_p=0;i_p<ntracks;i_p++){
    if(mother[i_p] != 0) continue;
    if(i_p == muon_idx || abs(pdg[i_p]) != 211) continue;
    p4.push_back(TLorentzVector(mom[i_p][0],mom[i_p][1],mom[i_p][2],mom[i_p][3]));
  }
  SortTLorentzVector(p4); 
  return p4;
}

TLorentzVector RecoPion4Mom(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4],const Float_t endpos[][4],const Int_t id[],int muon_idx){
  return SumTLorentzVector(RecoPion4MomV(ntracks,pdg,mother,mom,endpos,id,muon_idx));
}

std::vector<TLorentzVector> RecoShower4MomV(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4],const Float_t endpos[][4],const Int_t id[],int muon_idx){

  std::map<int,int> part_by_id;
  for(size_t i_p=0;i_p<ntracks;i_p++)
    part_by_id[id[i_p]] = pdg[i_p];

  std::vector<TLorentzVector> p4;
  for(size_t i_tr=0;i_tr<ntracks;i_tr++){
    if (part_by_id.find(mother[i_tr]) == part_by_id.end()) continue;
    if(pdg[i_tr] == 22 && mother[i_tr] != 0 && (part_by_id.at(mother[i_tr]) != 2212 && part_by_id.at(mother[i_tr]) != 211 && part_by_id.at(mother[i_tr]) != 11 && part_by_id.at(mother[i_tr]) != 13)){
      p4.push_back(TLorentzVector(mom[i_tr][0],mom[i_tr][1],mom[i_tr][2],mom[i_tr][3]));
    }
  }
  SortTLorentzVector(p4); 
  return p4;
}

TLorentzVector RecoShower4Mom(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4],const Float_t endpos[][4],const Int_t id[],int muon_idx){
  return SumTLorentzVector(RecoShower4MomV(ntracks,pdg,mother,mom,endpos,id,muon_idx));
}

double CalcRecoW(const Int_t& ntracks,const Int_t pdg[],const Int_t mother[],const Float_t mom[][4],const Float_t endpos[][4],const Int_t id[],int muon_idx){

  TLorentzVector p4_tot = RecoProton4Mom(ntracks,pdg,mother,mom,endpos,id,muon_idx)
                        + RecoPion4Mom(ntracks,pdg,mother,mom,endpos,id,muon_idx)
                        + RecoShower4Mom(ntracks,pdg,mother,mom,endpos,id,muon_idx);

  return p4_tot.M();

}

}

#endif
