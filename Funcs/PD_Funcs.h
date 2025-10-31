#ifndef _PD_Funcs_h_
#define _PD_Funcs_h_

#include "Funcs.h"

namespace pd {

double ProtonMom(double range){
  double Ep = 29.9*pow(range,0.586)/1e3;
  return sqrt(Ep*Ep + 2*Mp*Ep); 
}

double PionMom(double range){
   return 0.25798 + 0.0024088*range - 0.18828*pow(range,-0.11687);
}

double CalcTrueW(const std::vector<int>* pdg,const std::vector<double>* px,const std::vector<double>* py,const std::vector<double>* pz,std::vector<double>* E){
  TLorentzVector p4_tot(0,0,0,0); 
  for(size_t i_p=0;i_p<pdg->size();i_p++){
    if(abs(pdg->at(i_p)) == 13 || abs(pdg->at(i_p)) == 11) continue;
    double  p = TVector3(px->at(i_p),py->at(i_p),pz->at(i_p)).Mag();
    if(thresholds.find(abs(pdg->at(i_p))) != thresholds.end() && p > thresholds.at(abs(pdg->at(i_p))).first && p < thresholds.at(abs(pdg->at(i_p))).second)
      p4_tot += TLorentzVector(px->at(i_p),py->at(i_p),pz->at(i_p),E->at(i_p));
  }
  return p4_tot.M();
}

// Simple selection for finding muon 
int SimpleMuonSelection(const std::vector<float>* pid_v,const std::vector<float>* len_v){
  double longest_len = 0.0; // enforce minimum length of 10cm 
  int longest_idx = -1;
  for(size_t i_tr=0;i_tr<pid_v->size();i_tr++){
    if(pid_v->at(i_tr) > 0.2 && len_v->at(i_tr) > longest_len){
      longest_len = len_v->at(i_tr);
      longest_idx = i_tr;
    }
  }
  return longest_idx;
}

// Simple selection for finding leading reconstructed proton - PID requirement enforces orthogonality with muon selection
int SimpleProtonSelection(const std::vector<float>* pid_v,const std::vector<float>* len_v){
  double longest_len = 2.0; // enforce minimum length of 2cm 
  int longest_idx = -1;
  for(size_t i_tr=0;i_tr<pid_v->size();i_tr++){
    if(pid_v->at(i_tr) < 0.2 && len_v->at(i_tr) > longest_len){
      longest_len = len_v->at(i_tr);
      longest_idx = i_tr;
    }
  }
  return longest_idx;
}

// Simple selection for leading reconstructed charged pion
int SimpleChargedPionSelection(const std::vector<float>* pid_v,const std::vector<float>* len_v,int muon_idx){
  double longest_len = -1;  
  int longest_idx = -1;
  for(size_t i_tr=0;i_tr<pid_v->size();i_tr++){
    if(i_tr != muon_idx && pid_v->at(i_tr) > 0.2 && len_v->at(i_tr) > longest_len){
      longest_len = len_v->at(i_tr);
      longest_idx = i_tr;
    }
  }
  return longest_idx;
}


std::vector<TLorentzVector> RecoProton4MomV(const std::vector<float>* pid_v,const std::vector<float>* len_v,const std::vector<float>* trk_dir_x_v,const std::vector<float>* trk_dir_y_v,const std::vector<float>* trk_dir_z_v,int muon_idx){
  std::vector<TLorentzVector> p4;
  for(size_t i_tr=0;i_tr<pid_v->size();i_tr++){
    if(i_tr == muon_idx) continue;
    if(pid_v->at(i_tr) < 0.2){
      double p = ProtonMom(len_v->at(i_tr));
      p4.push_back(TLorentzVector(p*trk_dir_x_v->at(i_tr),p*trk_dir_y_v->at(i_tr),p*trk_dir_z_v->at(i_tr),sqrt(Mp*Mp+p*p)));
    }
  }
  SortTLorentzVector(p4); 
  return p4;
}

TLorentzVector RecoProton4Mom(const std::vector<float>* pid_v,const std::vector<float>* len_v,const std::vector<float>* trk_dir_x_v,const std::vector<float>* trk_dir_y_v,const std::vector<float>* trk_dir_z_v,int muon_idx){
  return SumTLorentzVector(RecoProton4MomV(pid_v,len_v,trk_dir_x_v,trk_dir_y_v,trk_dir_z_v,muon_idx));
}

std::vector<TLorentzVector> RecoPion4MomV(const std::vector<float>* pid_v,const std::vector<float>* len_v,const std::vector<float>* trk_dir_x_v,const std::vector<float>* trk_dir_y_v,const std::vector<float>* trk_dir_z_v,int muon_idx){
  std::vector<TLorentzVector> p4;
  for(size_t i_tr=0;i_tr<pid_v->size();i_tr++){
    if(i_tr == muon_idx) continue;
    if(pid_v->at(i_tr) > 0.2){
      double p = PionMom(len_v->at(i_tr));
      p4.push_back(TLorentzVector(p*trk_dir_x_v->at(i_tr),p*trk_dir_y_v->at(i_tr),p*trk_dir_z_v->at(i_tr),sqrt(mpi*mpi+p*p)));
    }
  }
  SortTLorentzVector(p4); 
  return p4;
}

TLorentzVector RecoPion4Mom(const std::vector<float>* pid_v,const std::vector<float>* len_v,const std::vector<float>* trk_dir_x_v,const std::vector<float>* trk_dir_y_v,const std::vector<float>* trk_dir_z_v,int muon_idx){
  return SumTLorentzVector(RecoPion4MomV(pid_v,len_v,trk_dir_x_v,trk_dir_y_v,trk_dir_z_v,muon_idx));
}

std::vector<TLorentzVector> RecoShower4MomV(const std::vector<float>* shr_x,const std::vector<float>* shr_y,const std::vector<float>* shr_z,const std::vector<float>* shr_e){
  std::vector<TLorentzVector> p4;
  for(size_t i_shr=0;i_shr<shr_x->size();i_shr++){
    double E = (shr_e->at(i_shr))/1e3/0.83; 
    if(E < 0) continue;
    p4.push_back(TLorentzVector(shr_x->at(i_shr)*E,shr_y->at(i_shr)*E,shr_z->at(i_shr)*E,E));
  } 
  SortTLorentzVector(p4); 
  return p4;
}

TLorentzVector RecoShower4Mom(const std::vector<float>* shr_x,const std::vector<float>* shr_y,const std::vector<float>* shr_z,const std::vector<float>* shr_e){
  return SumTLorentzVector(RecoShower4MomV(shr_x,shr_y,shr_z,shr_e));
}


double ReconstructedW(const std::vector<float>* pid_v,const std::vector<float>* len_v,const std::vector<float>* trk_dir_x_v,const std::vector<float>* trk_dir_y_v,const std::vector<float>* trk_dir_z_v,int muon_idx,const std::vector<float>* shr_x,const std::vector<float>* shr_y,const std::vector<float>* shr_z,const std::vector<float>* shr_e){

  TLorentzVector p4 = RecoProton4Mom(pid_v,len_v,trk_dir_x_v,trk_dir_y_v,trk_dir_z_v,muon_idx) 
                    + RecoPion4Mom(pid_v,len_v,trk_dir_x_v,trk_dir_y_v,trk_dir_z_v,muon_idx)
                    + RecoShower4Mom(shr_x,shr_y,shr_z,shr_e);

  return p4.M();

}

}

#endif
