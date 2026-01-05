#ifndef _Histograms_h_
#define _Histograms_h_

#include "Systematics.h"
#include "BranchList.h"

using namespace syst;

namespace hist {

class HistogramManager {

  public: 

    HistogramManager(std::string label,bool get_res=false);
    void LoadTemplate();
    void SetTemplate(std::string axis_title,int nbins,double low,double high);
    void SetTemplate(std::string axis_title,int nbins_t,double low_t,double high_t,int nbins_r,double low_r,double high_r);
    void FillHistograms(bool sel,double var,bool load_syst,double weight=1.0);   
    void FillHistograms2D(bool sel,double var_t,double var_r,bool load_syst,double weight=1.0);
    void DBBW() { _divide_by_bin_width = true; }
    void ShapeOnly() { _shape_only = true; }
    void KeepOU(){ _keep_overflow_underflow_ = true; }
    void KeepAll(){ _keep_all = true; }
    void Write();
    
  private:

    const std::string _label;
    const bool _get_res;
    bool _divide_by_bin_width = false;
    bool _shape_only = false; 
    bool _keep_overflow_underflow_ = false;
    bool _keep_all = false;

    // Reco histograms 
    TH1D* _h_tp = nullptr;
    TH1D* _h_tp_truth = nullptr;
    TH1D* _h_CV_Tot;
    std::vector<std::vector<TH1D*>> _h_Vars_Tot;
    std::vector<TH1D*> _h_CV;
    std::vector<std::vector<std::vector<TH1D*>>> _h_Vars;
  
    // Response histograms
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
    
    std::vector<double> _bin_edges_t;
    std::vector<double> _bin_edges_r;

    void SetupHistograms();

    TFile* _f_out = nullptr;

};

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

HistogramManager::HistogramManager(std::string label,bool get_res) : _label(label), _get_res(get_res)
{
  std::cout << "Setting up HistogramManager for " << _label << std::endl;
}

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void HistogramManager::LoadTemplate()
{
  TFile* f_tp = TFile::Open(("Analysis/"+_label+"/rootfiles/BinningTemplate.root").c_str());
  _h_tp = (TH1D*)f_tp->Get(("h_template_"+_label).c_str());
  _h_tp->SetDirectory(0);
  f_tp->Close();

  if(_get_res){
    TFile* f_tp_truth = TFile::Open(("Analysis/"+_label+"/rootfiles/TruthBinningTemplate.root").c_str());
    _h_tp_truth = (TH1D*)f_tp_truth->Get(("h_template_"+_label).c_str());
    _h_tp_truth->SetDirectory(0);
    f_tp_truth->Close();
    _h_tp_truth->SetName(("h_template_truth_"+_label).c_str()); 
  }

  SetupHistograms();
}

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void HistogramManager::SetTemplate(std::string axis_title,int nbins,double low,double high)
{
  if(_get_res) 
    throw std::invalid_argument("HistogramManager::SetTemplate: Calling wrong function when response matrices are being set, need to specify truth binning as well");   

  _h_tp = new TH1D(("h_template_"+_label).c_str(),axis_title.c_str(),nbins,low,high);

  SetupHistograms();
}

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void HistogramManager::SetTemplate(std::string axis_title,int nbins_t,double low_t,double high_t,int nbins_r,double low_r,double high_r)

{
  if(!_get_res) 
    throw std::invalid_argument("HistogramManager::SetTemplate: Calling wrong function response matrices not enabled");   

  _h_tp_truth = new TH1D(("h_template_truth_"+_label).c_str(),axis_title.c_str(),nbins_t,low_t,high_t);
  _h_tp = new TH1D(("h_template_"+_label).c_str(),axis_title.c_str(),nbins_r,low_r,high_r);

  SetupHistograms();
}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void HistogramManager::SetupHistograms()
{

  _bin_edges_r.clear();
  for(int i=1;i<_h_tp->GetNbinsX()+2;i++) _bin_edges_r.push_back(_h_tp->GetXaxis()->GetBinLowEdge(i));

  // Storing totals
  _h_CV_Tot = (TH1D*)_h_tp->Clone(("h_CV_Tot_"+_label).c_str());
  _h_CV_Tot->Sumw2();
  for(int i_s=0;i_s<kSystMAX;i_s++){
    _h_Vars_Tot.push_back(std::vector<TH1D*>());
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      _h_Vars_Tot.back().push_back((TH1D*)_h_tp->Clone(("h_Vars_Tot_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
    }
  }
  
  // Storing events broken down by category
  for(size_t i_c=0;i_c<categories.size();i_c++){
    _h_CV.push_back((TH1D*)_h_tp->Clone(("h_CV_"+categories.at(i_c)+"_"+_label).c_str()));
    _h_CV.back()->Sumw2();
    _h_Vars.push_back(std::vector<std::vector<TH1D*>>());
    for(int i_s=0;i_s<kSystMAX;i_s++){
      _h_Vars.back().push_back(std::vector<TH1D*>());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        _h_Vars.back().back().push_back((TH1D*)_h_tp->Clone(("h_Vars_"+categories.at(i_c)+"_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
      }
    }
  }


  if(_get_res){

    _bin_edges_t.clear();
    for(int i=1;i<_h_tp_truth->GetNbinsX()+2;i++) _bin_edges_t.push_back(_h_tp_truth->GetXaxis()->GetBinLowEdge(i));

    // Storing totals
    _h_CV_Truth = (TH1D*)_h_tp_truth->Clone(("h_CV_Truth_"+_label).c_str());
    _h_CV_Truth->Sumw2();

    _h_CV_Reco = (TH1D*)_h_tp->Clone(("h_CV_Reco_"+_label).c_str());
    _h_CV_Reco->Sumw2();

    _h_CV_BG = (TH1D*)_h_tp->Clone(("h_CV_BG_"+_label).c_str());
    _h_CV_BG->Sumw2();

    _h_CV_Joint = new TH2D(("h_CV_Joint_"+_label).c_str(),"",_bin_edges_t.size()-1,&_bin_edges_t[0],_bin_edges_r.size()-1,&_bin_edges_r[0]);
    _h_CV_Joint->Sumw2();

    for(int i_s=0;i_s<kSystMAX;i_s++){
      _h_Vars_Truth.push_back(std::vector<TH1D*>());
      _h_Vars_Reco.push_back(std::vector<TH1D*>());
      _h_Vars_BG.push_back(std::vector<TH1D*>());
      _h_Vars_Joint.push_back(std::vector<TH2D*>());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        _h_Vars_Truth.back().push_back((TH1D*)_h_tp_truth->Clone(("h_Vars_Truth_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
        _h_Vars_Reco.back().push_back((TH1D*)_h_tp->Clone(("h_Vars_Reco_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
        _h_Vars_BG.back().push_back((TH1D*)_h_tp->Clone(("h_Vars_BG_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
        _h_Vars_Joint.back().push_back(new TH2D(("h_Vars_Joint_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str(),"",_bin_edges_t.size()-1,&_bin_edges_t[0],_bin_edges_r.size()-1,&_bin_edges_r[0]));
      }
    }

  }

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void HistogramManager::FillHistograms(bool sel,double var,bool load_syst,double weight)
{

  if(!sel) return;

  if(_h_tp == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(category == -1) std::cout << "Bad event" << std::endl;

  if(std::isnan(weightSpline) || std::isinf(weightSpline)) return;

  _h_CV.at(category)->Fill(var,POT_weight*weightSpline);
  if(!is_data) _h_CV_Tot->Fill(var,POT_weight*weightSpline);

  if(!is_data){
    if(load_syst){
      for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Tot.at(kGenie).at(i_u)->Fill(var,(double)POT_weight*weightsGenie->at(i_u)/1000);
      for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Tot.at(kReint).at(i_u)->Fill(var,(double)POT_weight*weightsReint->at(i_u)*weightSpline/1000);
      for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Tot.at(kFlux).at(i_u)->Fill(var,(double)POT_weight*weightsFlux->at(i_u)*weightSpline/1000);
    }
    else {
      for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Tot.at(kGenie).at(i_u)->Fill(var,(double)POT_weight*weightSpline);
      for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Tot.at(kReint).at(i_u)->Fill(var,(double)POT_weight*weightSpline);
      for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Tot.at(kFlux).at(i_u)->Fill(var,(double)POT_weight*weightSpline);
    } 
  }

  if(load_syst){
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars.at(category).at(kGenie).at(i_u)->Fill(var,(double)POT_weight*weightsGenie->at(i_u)/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars.at(category).at(kReint).at(i_u)->Fill(var,(double)POT_weight*weightsReint->at(i_u)*weightSpline/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars.at(category).at(kFlux).at(i_u)->Fill(var,(double)POT_weight*weightsFlux->at(i_u)*weightSpline/1000);
  }
  else {
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars.at(category).at(kGenie).at(i_u)->Fill(var,(double)POT_weight*weightSpline);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars.at(category).at(kReint).at(i_u)->Fill(var,(double)POT_weight*weightSpline);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars.at(category).at(kFlux).at(i_u)->Fill(var,(double)POT_weight*weightSpline);
  } 

}


///////////////////////////////////////////////////////////////////////
// Fill the histograms

void HistogramManager::FillHistograms2D(bool sel,double var_t,double var_r,bool load_syst,double weight)
{
  if(!_get_res) 
    throw std::invalid_argument("HistogramManager::FillHistograms2D: Trying to fill histograms before they enabled");   
  
  FillHistograms(sel,var_r,load_syst,weight);

  if(category == kSignal){
    _h_CV_Truth->Fill(var_t,POT_weight*weightSpline);
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Truth.at(kGenie).at(i_u)->Fill(var_t,(double)POT_weight*weightsGenie->at(i_u)/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Truth.at(kReint).at(i_u)->Fill(var_t,(double)POT_weight*weightsReint->at(i_u)*weightSpline/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Truth.at(kFlux).at(i_u)->Fill(var_t,(double)POT_weight*weightsFlux->at(i_u)*weightSpline/1000);
  }

  if(!sel) return;

  if(category == kSignal){
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
// Write the histograms to file

void HistogramManager::Write()
{

  std::cout << "Writing histograms for " << _label << std::endl;

  gSystem->Exec(("mkdir -p Analysis/"+_label+"/rootfiles/").c_str());
  TFile* f_out = TFile::Open(("Analysis/"+_label+"/rootfiles/Histograms.root").c_str(),"RECREATE");

  // Store the CV before any scaling so we can calculate the errors correctly
  TH1D* h_CV_NoScale = (TH1D*)_h_CV_Tot->Clone("h_CV_NoScale");
  std::vector<TH1D*> h_CV_Cat_NoScale;
  for(size_t i_c=0;i_c<categories.size();i_c++)
    h_CV_Cat_NoScale.push_back((TH1D*)_h_CV.at(i_c)->Clone(("h_CV_Cat_NoScale_"+categories.at(i_c)).c_str()));

  if(_divide_by_bin_width){
    DivideByBinWidth(_h_CV_Tot);
  }  

  double integral = _h_CV_Tot->Integral("width");
  if(_keep_overflow_underflow_){
    integral = _h_CV_Tot->GetBinContent(0) + _h_CV_Tot->GetBinContent(_h_CV_Tot->GetNbinsX()+1);
    if(_divide_by_bin_width ) for(int i=1;i<_h_CV_Tot->GetNbinsX()+1;i++) integral += _h_CV_Tot->GetBinContent(i)*_h_CV_Tot->GetBinWidth(i); 
    else for(int i=1;i<_h_CV_Tot->GetNbinsX()+1;i++) integral += _h_CV_Tot->GetBinContent(i);
  } 

  if(_shape_only){
    _h_CV_Tot->Scale(1.0/integral);
  }
  _h_CV_Tot->Write("h_CV_Tot");

  for(size_t i_c=0;i_c<categories.size();i_c++){
    if(_divide_by_bin_width) DivideByBinWidth(_h_CV.at(i_c));
    if(_shape_only) _h_CV.at(i_c)->Scale(1.0/integral);
    _h_CV.at(i_c)->Write(("h_CV_"+categories.at(i_c)).c_str());
   }

  std::vector<std::vector<double>> integral_v(kSystMAX);
  for(int i_s=0;i_s<kSystMAX;i_s++){
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      if(_divide_by_bin_width) DivideByBinWidth(_h_Vars_Tot.at(i_s).at(i_u));
      integral_v.at(i_s).push_back(_h_Vars_Tot.at(i_s).at(i_u)->Integral("width"));       

      if(_keep_overflow_underflow_){
        integral_v.at(i_s).at(i_u) = _h_Vars_Tot.at(i_s).at(i_u)->GetBinContent(0) + _h_Vars_Tot.at(i_s).at(i_u)->GetBinContent(_h_CV_Tot->GetNbinsX()+1);
        if(_divide_by_bin_width ) for(int i=1;i<_h_CV_Tot->GetNbinsX()+1;i++) integral_v.at(i_s).at(i_u) += _h_Vars_Tot.at(i_s).at(i_u)->GetBinContent(i)*_h_Vars_Tot.at(i_s).at(i_u)->GetBinWidth(i); 
        else for(int i=1;i<_h_CV_Tot->GetNbinsX()+1;i++) integral_v.at(i_s).at(i_u) += _h_Vars_Tot.at(i_s).at(i_u)->GetBinContent(i);
      } 

      if(_shape_only) _h_Vars_Tot.at(i_s).at(i_u)->Scale(1.0/integral_v.at(i_s).back());
      if(_keep_all) _h_Vars_Tot.at(i_s).at(i_u)->Write(("h_Vars_Tot_"+sys_str.at(i_s)+"_"+std::to_string(i_u)).c_str());     
    }
  }   

  for(size_t i_c=0;i_c<categories.size();i_c++){
    for(int i_s=0;i_s<kSystMAX;i_s++){
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        if(_divide_by_bin_width) DivideByBinWidth(_h_Vars.at(i_c).at(i_s).at(i_u)); 
        if(_shape_only) _h_Vars.at(i_c).at(i_s).at(i_u)->Scale(1.0/integral_v.at(i_s).back());
        if(_keep_all) _h_Vars.at(i_c).at(i_s).at(i_u)->Write(("h_Vars_"+categories.at(i_c)+"_"+sys_str.at(i_s)+"_"+std::to_string(i_u)).c_str());
      }
    }
  }

  // Systematics covariance matrices
  std::vector<TH2D*> h_Cov_Tot; 
  std::vector<TH2D*> h_FCov_Tot; 
  for(int i_s=0;i_s<kSystMAX;i_s++){
    h_Cov_Tot.push_back(nullptr);
    h_FCov_Tot.push_back(nullptr);
    CalcCovMultisim(sys_str.at(i_s),_h_CV_Tot,_h_Vars_Tot.at(i_s),h_Cov_Tot.back(),h_FCov_Tot.back());
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
  for(int i=0;i<_h_CV_Tot->GetNbinsX()+2;i++){
    double bin_error = _h_CV_Tot->GetBinError(i);
    double bin_content = _h_CV_Tot->GetBinContent(i);
    double orig_bin_content = h_CV_NoScale->GetBinContent(i);
    double scale = bin_content/orig_bin_content;
    h_Cov_MCStat->SetBinContent(i,i,bin_error*bin_error);
    h_FCov_MCStat->SetBinContent(i,i,bin_error*bin_error/bin_content/bin_content);
    h_Cov_EstDataStat->SetBinContent(i,i,orig_bin_content*scale*scale);
    h_FCov_EstDataStat->SetBinContent(i,i,1.0/orig_bin_content);
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
    for(int i=0;i<_h_CV_Tot->GetNbinsX()+2;i++){
      double bin_error = _h_CV.at(i_c)->GetBinError(i);
      double bin_content = _h_CV.at(i_c)->GetBinContent(i);
      double orig_bin_content = h_CV_Cat_NoScale.at(i_c)->GetBinContent(i);
      double scale = bin_content/orig_bin_content;
      h_Cov_MCStat_Cat.back()->SetBinContent(i,i,bin_error*bin_error);
      h_FCov_MCStat_Cat.back()->SetBinContent(i,i,bin_error*bin_error/bin_content/bin_content);
      h_Cov_EstDataStat_Cat.back()->SetBinContent(i,i,orig_bin_content*scale*scale);
      h_FCov_EstDataStat_Cat.back()->SetBinContent(i,i,1.0/orig_bin_content);
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
      CalcCovMultisim(categories.at(i_c)+"_"+sys_str.at(i_s),_h_CV.at(i_c),_h_Vars.at(i_c).at(i_s),h_Cov_Cat.back().back(),h_FCov_Cat.back().back());
      h_Cov_Cat.back().back()->Write();
      h_FCov_Cat.back().back()->Write();
    }
  }       

  f_out->Close();

  if(_get_res){

    TFile* f_out_res = TFile::Open(("Analysis/"+_label+"/rootfiles/Response.root").c_str(),"RECREATE");

    NormaliseResponse(_h_CV_Truth,_h_CV_Joint);
    _h_CV_Truth->Write("h_CV_Truth");
    _h_CV_Reco->Write("h_CV_Reco");
    _h_CV_BG->Write("h_CV_BG");
    _h_CV_Joint->Write("h_CV_Joint");
    /*
       std::vector<std::pair<int,int>> empty_bins;
       for(int i=1;i<_h_CV_Joint->GetNbinsX()+1;i++)
       for(int j=1;j<_h_CV_Joint->GetNbinsY()+1;j++)
       if(!(_h_CV_Joint->GetBinContent(i,j) > 0)) empty_bins.push_back(std::make_pair(i,j));
       */
    for(int i_s=0;i_s<kSystMAX;i_s++){
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        NormaliseResponse(_h_Vars_Truth.at(i_s).at(i_u),_h_Vars_Joint.at(i_s).at(i_u));
        _h_Vars_Truth.at(i_s).at(i_u)->Write();
        _h_Vars_Reco.at(i_s).at(i_u)->Write();
        _h_Vars_BG.at(i_s).at(i_u)->Write();
        _h_Vars_Joint.at(i_s).at(i_u)->Write();
      }
    }
    /*
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
       */

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
      /*
         if(_calc_cov){
         TH1D* unrolled = _discard_empty ? Unroll2DDist(joint,empty_bins) : Unroll2DDist(joint); 
         unrolled->Write((string(joint->GetName())+"_Unrolled").c_str());
         }
         */

    }

    f_out_res->Close();

  }

}

}

#endif
