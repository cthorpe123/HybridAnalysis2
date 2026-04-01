#ifndef _LT_Funcs_h_
#define _LT_Funcs_h_

#include "Funcs.h"

namespace lt {

int FindLeadingTruePart(Int_t nTruePrimParts,Int_t truePrimPartPDG[],Float_t truePrimPartPx[],Float_t truePrimPartPy[],Float_t truePrimPartPz[],Int_t pdg){

  int leading = -1;
  double leading_p = -1;
  for(int i=0;i<nTruePrimParts;i++){
    double p = TVector3(truePrimPartPx[i],truePrimPartPy[i],truePrimPartPz[i]).Mag();
    if(abs(truePrimPartPDG[i]) == pdg && p > leading_p){
      leading = i;
      leading_p = p;
    } 
  } 
  if(leading_p > thresholds.at(abs(pdg)).first)  return leading;
  else return -1;
}

int SimpleMuonSelection(Int_t nTracks,Int_t trackIsSecondary[],Int_t trackPID[]){
  int reco_muon = -1; 
  for(int i=0;i<nTracks;i++){
    if(trackIsSecondary[i]) continue;
    if(abs(trackPID[i]) == 13){
      reco_muon = i;
      break;
    }
  } 
  return reco_muon;
}

int SimpleNueSelection(Int_t nShowers,Int_t showerIsSecondary[],Int_t showerPID[]){
  int reco_electron = -1; 
  for(int i=0;i<nShowers;i++){
    if(showerIsSecondary[i]) continue;
    if(abs(showerPID[i]) == 11){
      reco_electron = i;
      break;
    }
  } 
  return reco_electron;
}

  
// Search for reconstruction artefacts
const double cut_c = 0.2;
const double cut_f_low = 0.00, cut_f_high = 0.4;

std::vector<TLorentzVector> RecoProton4MomV(Int_t nTracks,Int_t trackIsSecondary[],Int_t trackPID[],Float_t trackRecoE[],Float_t trackStartDirX[],Float_t trackStartDirY[],Float_t trackStartDirZ[]){
  std::vector<TLorentzVector> p4;
  for(int i=0;i<nTracks;i++){
    if(trackIsSecondary[i]) continue;
    if(abs(trackPID[i]) == 2212){
      double mom = sqrt(trackRecoE[i]*trackRecoE[i]/1e6 + 2*0.938*trackRecoE[i]/1e3);
      p4.push_back(TLorentzVector(mom*trackStartDirX[i],mom*trackStartDirY[i],mom*trackStartDirZ[i],sqrt(mom*mom+0.938*0.938)));
    }
  }
  SortTLorentzVector(p4);

  bool bad = true;
  while(bad){
    bad = false;
    for(size_t i=0;i<p4.size();i++){
      for(size_t j=i+1;j<p4.size();j++){
        double c = p4.at(i).Vect().Angle(p4.at(j).Vect());
        double f = abs(p4.at(i).Vect().Mag() - p4.at(j).Vect().Mag())/(p4.at(i).Vect().Mag() + p4.at(j).Vect().Mag());
        if(c < cut_c && f > cut_f_low && f < cut_f_high){
          p4.erase(p4.begin()+j);
          bad = true;
          break;
        }
      }
      if(bad) break;
    }
  }

  return p4;
}

TLorentzVector RecoProton4Mom(Int_t nTracks,Int_t trackIsSecondary[],Int_t trackPID[],Float_t trackRecoE[],Float_t trackStartDirX[],Float_t trackStartDirY[],Float_t trackStartDirZ[]){
  return SumTLorentzVector(RecoProton4MomV(nTracks,trackIsSecondary,trackPID,trackRecoE,trackStartDirX,trackStartDirY,trackStartDirZ));
}



std::vector<TLorentzVector> RecoPion4MomV(Int_t nTracks,Int_t trackIsSecondary[],Int_t trackPID[],Float_t trackRecoE[],Float_t trackStartDirX[],Float_t trackStartDirY[],Float_t trackStartDirZ[]){
  std::vector<TLorentzVector> p4;
  for(int i=0;i<nTracks;i++){
    if(trackIsSecondary[i]) continue;
    if(abs(trackPID[i]) == 211){
      double mom = sqrt(trackRecoE[i]*trackRecoE[i]/1e6 + 2*mpi*trackRecoE[i]/1e3);
      p4.push_back(TLorentzVector(mom*trackStartDirX[i],mom*trackStartDirY[i],mom*trackStartDirZ[i],sqrt(mom*mom+mpi*mpi)));
    }
  }
  SortTLorentzVector(p4);
  return p4;
}

TLorentzVector RecoPion4Mom(Int_t nTracks,Int_t trackIsSecondary[],Int_t trackPID[],Float_t trackRecoE[],Float_t trackStartDirX[],Float_t trackStartDirY[],Float_t trackStartDirZ[]){
  return SumTLorentzVector(RecoPion4MomV(nTracks,trackIsSecondary,trackPID,trackRecoE,trackStartDirX,trackStartDirY,trackStartDirZ));
}


std::vector<TLorentzVector> RecoShower4MomV(Int_t nShowers,Int_t showerIsSecondary[],Int_t showerPID[],Float_t showerRecoE[],Float_t showerStartDirX[],Float_t showerStartDirY[],Float_t showerStartDirZ[]){
  std::vector<TLorentzVector> p4;
  for(int i=0;i<nShowers;i++){
    if(showerPID[i] == 22 && !showerIsSecondary[i]){
      p4.push_back(TLorentzVector(showerRecoE[i]*showerStartDirX[i]/1e3,showerRecoE[i]*showerStartDirY[i]/1e3,showerRecoE[i]*showerStartDirZ[i]/1e3,showerRecoE[i]/1e3));
    } 
  }    
  SortTLorentzVector(p4);
  return p4;
}

TLorentzVector RecoShower4Mom(Int_t nShowers,Int_t showerIsSecondary[],Int_t showerPID[],Float_t showerRecoE[],Float_t showerStartDirX[],Float_t showerStartDirY[],Float_t showerStartDirZ[]){
  return SumTLorentzVector(RecoShower4MomV(nShowers,showerIsSecondary,showerPID,showerRecoE,showerStartDirX,showerStartDirY,showerStartDirZ));
}

double CalcRecoW(Int_t nTracks,Int_t trackIsSecondary[],Int_t trackPID[],Float_t trackRecoE[],Float_t trackStartDirX[],Float_t trackStartDirY[],Float_t trackStartDirZ[],Int_t nShowers,Int_t showerIsSecondary[],Int_t showerPID[],Float_t showerRecoE[],Float_t showerStartDirX[],Float_t showerStartDirY[],Float_t showerStartDirZ[]){ 

  TLorentzVector reco_p4 = RecoProton4Mom(nTracks,trackIsSecondary,trackPID,trackRecoE,trackStartDirX,trackStartDirY,trackStartDirZ)
                         + RecoPion4Mom(nTracks,trackIsSecondary,trackPID,trackRecoE,trackStartDirX,trackStartDirY,trackStartDirZ)
                         + RecoShower4Mom(nShowers,showerIsSecondary,showerPID,showerRecoE,showerStartDirX,showerStartDirY,showerStartDirZ);

  return reco_p4.M();

}


}

#endif
