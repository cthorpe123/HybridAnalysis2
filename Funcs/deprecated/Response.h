#ifndef _Response_h_
#define _Response_h_

#include "Systematics.h"
#include "BranchList.h"

using namespace syst;

namespace hist {

class Response {

  public: 

    Response(std::string label);
    void LoadTemplate();
    void SetTemplate(std::string axis_title,int nbins_truth,double low_truth,double high_truth,int nbins_reco,double low_reco,double high_reco);
    void AddSpecialUniv(std::string name);
    void Fill(double var_t,double var_r,bool sig,bool sel,double weight=1.0);
    void FillSpecialUniv(std::string name,double var_t,double var_r,bool sig,bool sel,double weight=1.0);
    void DBBW() { _divide_by_bin_width = true; }
    void ShapeOnly() { _shape_only = true; }
    void KeepOU(){ _keep_overflow_underflow_ = true; }
    void DiscardEmpty(){ _discard_empty = true; }
    void KeepAll(){ _keep_all = true; }
    void CalcCov() { _calc_cov = true; }
    void Write();

  private:

    const std::string _label;
    const bool _get_res;
    bool _divide_by_bin_width = false;
    bool _shape_only = false; 
    bool _keep_overflow_underflow_ = false;
    bool _discard_empty = false;
    bool _keep_all = false;
    bool _calc_cov = false;

    TH1D* _h_tp_truth = nullptr;
    TH1D* _h_tp_reco = nullptr;
    TH1D* _h_CV_Truth;
    TH1D* _h_CV_Reco;
    TH1D* _h_CV_BG;
    TH2D* _h_CV_Joint;

    std::vector<std::vector<TH1D*>> _h_Vars_Truth;
    std::vector<std::vector<TH1D*>> _h_Vars_Reco;
    std::vector<std::vector<TH1D*>> _h_Vars_BG;
    std::vector<std::vector<TH2D*>> _h_Vars_Joint;

    std::map<std::string,TH1D*> _h_Special_Truth;
    std::map<std::string,TH1D*> _h_Special_Reco;
    std::map<std::string,TH1D*> _h_Special_BG;
    std::map<std::string,TH2D*> _h_Special_Joint;

    std::vector<double> _bin_edges_truth;
    std::vector<double> _bin_edges_reco;
    void SetupHistograms();

};

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

Response::Response(std::string label) : _label(label)
{
  std::cout << "Setting up Response for " << _label << std::endl;
}

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void Response::LoadTemplate()
{
  TFile* f_tp_truth = TFile::Open(("Analysis/"+_label+"/rootfiles/TruthBinningTemplate.root").c_str());
  _h_tp_truth = (TH1D*)f_tp_truth->Get(("h_template_"+_label).c_str());
  _h_tp_truth->SetDirectory(0);
  f_tp_truth->Close();
  _h_tp_truth->SetName(("h_template_truth_"+_label).c_str()); 

  TFile* f_tp_reco = TFile::Open(("Analysis/"+_label+"/rootfiles/BinningTemplate.root").c_str());
  _h_tp_reco = (TH1D*)f_tp_reco->Get(("h_template_"+_label).c_str());
  _h_tp_reco->SetDirectory(0);
  f_tp_reco->Close();

  SetupHistograms();
}

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void Response::SetTemplate(std::string axis_title,int nbins_truth,double low_truth,double high_truth,int nbins_reco,double low_reco,double high_reco)
{
  _h_tp_truth = new TH1D(("h_template_truth_"+_label).c_str(),axis_title.c_str(),nbins_truth,low_truth,high_truth);
  _h_tp_reco = new TH1D(("h_template_"+_label).c_str(),axis_title.c_str(),nbins_reco,low_reco,high_reco);

  SetupHistograms();
}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void Response::SetupHistograms()
{

  _bin_edges_truth.clear();
  _bin_edges_reco.clear();
  for(int i=1;i<_h_tp_truth->GetNbinsX()+2;i++) _bin_edges_truth.push_back(_h_tp_truth->GetXaxis()->GetBinLowEdge(i));
  for(int i=1;i<_h_tp_reco->GetNbinsX()+2;i++) _bin_edges_reco.push_back(_h_tp_reco->GetXaxis()->GetBinLowEdge(i));

  // Storing totals
  _h_CV_Truth = (TH1D*)_h_tp_truth->Clone(("h_CV_Truth_"+_label).c_str());
  _h_CV_Truth->Sumw2();

  _h_CV_Reco = (TH1D*)_h_tp_reco->Clone(("h_CV_Reco_"+_label).c_str());
  _h_CV_Reco->Sumw2();

  _h_CV_BG = (TH1D*)_h_tp_reco->Clone(("h_CV_BG_"+_label).c_str());
  _h_CV_BG->Sumw2();

  _h_CV_Joint = new TH2D(("h_CV_Joint_"+_label).c_str(),"",_bin_edges_truth.size()-1,&_bin_edges_truth[0],_bin_edges_reco.size()-1,&_bin_edges_reco[0]);
  _h_CV_Joint->Sumw2();

  for(int i_s=0;i_s<kSystMAX;i_s++){
    _h_Vars_Truth.push_back(std::vector<TH1D*>());
    _h_Vars_Reco.push_back(std::vector<TH1D*>());
    _h_Vars_BG.push_back(std::vector<TH1D*>());
    _h_Vars_Joint.push_back(std::vector<TH2D*>());
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      _h_Vars_Truth.back().push_back((TH1D*)_h_tp_truth->Clone(("h_Vars_Truth_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
      _h_Vars_Reco.back().push_back((TH1D*)_h_tp_reco->Clone(("h_Vars_Reco_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
      _h_Vars_BG.back().push_back((TH1D*)_h_tp_reco->Clone(("h_Vars_BG_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
      _h_Vars_Joint.back().push_back(new TH2D(("h_Vars_Joint_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str(),"",_bin_edges_truth.size()-1,&_bin_edges_truth[0],_bin_edges_reco.size()-1,&_bin_edges_reco[0]));
    }
  }

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void Response::AddSpecialUniv(std::string name){

  if(_h_tp_truth == nullptr || _h_tp_reco == nullptr) 
    throw std::invalid_argument("Trying to call AddSpecialUniv before binning has been set, call LoadTemplate or SetTemplate first!");

  _h_Special_Truth[name] = (TH1D*)_h_tp_truth->Clone(("h_Special_Truth_"+name).c_str()); 
  _h_Special_Reco[name] = (TH1D*)_h_tp_reco->Clone(("h_Special_Reco_"+name).c_str()); 
  _h_Special_BG[name] = (TH1D*)_h_tp_reco->Clone(("h_Special_BG_"+name).c_str()); 
  _h_Special_Joint[name] = (TH2D*)_h_CV_Joint->Clone(("h_Special_Joint_"+name).c_str());

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void Response::Fill(double var_t,double var_r,bool sig,bool sel,double weight)
{

  if(_h_tp_truth == nullptr || _h_tp_reco == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(std::isnan(weightSpline) || std::isinf(weightSpline)) return;


  if(sig){
    _h_CV_Truth->Fill(var_t,POT_weight*weightSpline);
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Truth.at(kGenie).at(i_u)->Fill(var_t,(double)POT_weight*weightsGenie->at(i_u)/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Truth.at(kReint).at(i_u)->Fill(var_t,(double)POT_weight*weightsReint->at(i_u)*weightSpline/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Truth.at(kFlux).at(i_u)->Fill(var_t,(double)POT_weight*weightsFlux->at(i_u)*weightSpline/1000);
  }

  if(!sel) return;

  if(sig){
    _h_CV_Reco->Fill(var_r,POT_weight*weightSpline);
    _h_CV_Joint->Fill(var_t,var_r,POT_weight*weightSpline);
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Reco.at(kGenie).at(i_u)->Fill(var_r,(double)POT_weight*weightsGenie->at(i_u)/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Reco.at(kReint).at(i_u)->Fill(var_r,(double)POT_weight*weightsReint->at(i_u)*weightSpline/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Reco.at(kFlux).at(i_u)->Fill(var_r,(double)POT_weight*weightsFlux->at(i_u)*weightSpline/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Joint.at(kGenie).at(i_u)->Fill(var_t,var_r,(double)POT_weight*weightsGenie->at(i_u)/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Joint.at(kReint).at(i_u)->Fill(var_t,var_r,(double)POT_weight*weightsReint->at(i_u)*weightSpline/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Joint.at(kFlux).at(i_u)->Fill(var_t,var_r,(double)POT_weight*weightsFlux->at(i_u)*weightSpline/1000);
  }
  else {
    _h_CV_BG->Fill(var_r,POT_weight*weightSpline);
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_BG.at(kGenie).at(i_u)->Fill(var_r,(double)POT_weight*weightsGenie->at(i_u)/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_BG.at(kReint).at(i_u)->Fill(var_r,(double)POT_weight*weightsReint->at(i_u)*weightSpline/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_BG.at(kFlux).at(i_u)->Fill(var_r,(double)POT_weight*weightsFlux->at(i_u)*weightSpline/1000);
  }

}

///////////////////////////////////////////////////////////////////////
// Fill special universe histogram 

void Response::FillSpecialUniv(std::string name,double var_t,double var_r,bool sig,bool sel,double weight=1.0)
{

  if(_h_tp_truth == nullptr || _h_tp_reco == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(std::isnan(weightSpline) || std::isinf(weightSpline)) return;


  if(sig) _h_Special_Truth.at(name)->Fill(var_t,POT_weight*weightSpline*weight);  

  if(!sel) return;

  if(sig){
    _h_Special_Reco.at(name)->Fill(var_r,POT_weight*weightSpline*weight);  
    _h_Special_Joint.at(name)->Fill(var_t,var_r,POT_weight*weightSpline*weight);  
  }
  else 
    _h_Special_BG.at(name)->Fill(var_r,POT_weight*weightSpline*weight); 

}

///////////////////////////////////////////////////////////////////////
// Write the histograms to file

void Response::Write()
{

  std::cout << "Writing response histograms for " << _label << std::endl;

  gSystem->Exec(("mkdir -p Analysis/"+_label+"/rootfiles/").c_str());
  TFile* f_out = TFile::Open(("Analysis/"+_label+"/rootfiles/Response.root").c_str(),"RECREATE");

  NormaliseResponse(_h_CV_Truth,_h_CV_Joint);
  _h_CV_Truth->Write("h_CV_Truth");
  _h_CV_Reco->Write("h_CV_Reco");
  _h_CV_BG->Write("h_CV_BG");
  _h_CV_Joint->Write("h_CV_Joint");

  std::vector<std::pair<int,int>> empty_bins;
  for(int i=1;i<_h_CV_Joint->GetNbinsX()+1;i++)
    for(int j=1;j<_h_CV_Joint->GetNbinsY()+1;j++)
      if(!(_h_CV_Joint->GetBinContent(i,j) > 0)) empty_bins.push_back(std::make_pair(i,j));

  for(int i_s=0;i_s<kSystMAX;i_s++){
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      NormaliseResponse(_h_Vars_Truth.at(i_s).at(i_u),_h_Vars_Joint.at(i_s).at(i_u));
      _h_Vars_Truth.at(i_s).at(i_u)->Write();
      _h_Vars_Reco.at(i_s).at(i_u)->Write();
      _h_Vars_BG.at(i_s).at(i_u)->Write();
      _h_Vars_Joint.at(i_s).at(i_u)->Write();
    }
  }

  if(_calc_cov){

    TH1D* h_CV_Joint_Unrolled = _discard_empty ? Unroll2DDist(_h_CV_Joint,empty_bins) : Unroll2DDist(_h_CV_Joint);
    h_CV_Joint_Unrolled->Write("h_CV_Joint_Unrolled");

    std::vector<TH2D*> h_Cov(kSystMAX);
    std::vector<TH2D*> h_FCov(kSystMAX);
    std::vector<std::vector<TH1D*>> h_Vars_Joint_Unrolled(kSystMAX); 

    for(int i_s=0;i_s<kSystMAX;i_s++){
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        h_Vars_Joint_Unrolled.at(i_s).push_back(_discard_empty ? Unroll2DDist(_h_Vars_Joint.at(i_s).at(i_u),empty_bins) : Unroll2DDist(_h_Vars_Joint.at(i_s).at(i_u)));
        if(_keep_all) h_Vars_Joint_Unrolled.at(i_s).back()->Write((string(_h_Vars_Truth.at(i_s).at(i_u)->GetName())+"_Unrolled").c_str());
      }
      std::cout << "Calculating covariance for systematic " << sys_str.at(i_s) << std::endl;
      CalcCovMultisim(sys_str.at(i_s),h_CV_Joint_Unrolled,h_Vars_Joint_Unrolled.at(i_s),h_Cov.at(i_s),h_FCov.at(i_s));
      h_Cov.at(i_s)->Write(("Cov_"+sys_str.at(i_s)).c_str());
      h_FCov.at(i_s)->Write(("FCov_"+sys_str.at(i_s)).c_str());
    }

    TH2D* h_Cov_MCStat = (TH2D*)h_Cov.back()->Clone("h_Cov_MCStat");
    TH2D* h_FCov_MCStat = (TH2D*)h_FCov.back()->Clone("h_FCov_MCStat");

    h_Cov_MCStat->Reset();
    h_FCov_MCStat->Reset();
    for(int i=1;i<h_CV_Joint_Unrolled->GetNbinsX()+1;i++){
      double content = h_CV_Joint_Unrolled->GetBinContent(i);
      double error = h_CV_Joint_Unrolled->GetBinError(i);
      h_Cov_MCStat->SetBinContent(i,i,error*error);
      h_FCov_MCStat->SetBinContent(i,i,error*error/content/content);
    }

    h_Cov_MCStat->Write("Cov_MCStat");
    h_FCov_MCStat->Write("FCov_MCStat");

    TH2D* h_Cov_Tot = (TH2D*)h_Cov_MCStat->Clone("h_Cov_Tot");
    TH2D* h_FCov_Tot = (TH2D*)h_FCov_MCStat->Clone("h_FCov_Tot");

    for(int i_s=0;i_s<kSystMAX;i_s++){
      h_Cov_Tot->Add(h_Cov.at(i_s));
      h_FCov_Tot->Add(h_FCov.at(i_s));
    } 

    h_Cov_Tot->Write("Cov_Tot");
    h_FCov_Tot->Write("FCov_Tot");
  }

  std::map<std::string,TH1D*>::iterator it;
  for(it = _h_Special_Truth.begin();it != _h_Special_Truth.end();it++){
    std::string name = it->first;
    TH1D* truth = _h_Special_Truth.at(name);
    TH1D* reco = _h_Special_Reco.at(name);
    TH1D* bg = _h_Special_BG.at(name);
    TH2D* joint = _h_Special_Joint.at(name);
    NormaliseResponse(truth,joint);

    truth->Write();
    reco->Write();
    bg->Write();
    joint->Write();

    if(_calc_cov){
      TH1D* unrolled = _discard_empty ? Unroll2DDist(joint,empty_bins) : Unroll2DDist(joint); 
      unrolled->Write((string(joint->GetName())+"_Unrolled").c_str());
    }

  }

  f_out->Close();

}

}

#endif
