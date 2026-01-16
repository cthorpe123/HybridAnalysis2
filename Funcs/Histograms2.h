#ifndef _Histograms_h_
#define _Histograms_h_

#include "Systematics.h"
#include "BranchList.h"

using namespace syst;

namespace hist {

class HistogramManager {

  public: 

    HistogramManager(std::string label,bool save_truth=false);
    void LoadTemplate();
    void SetTemplate(std::string axis_title,int nbins,double low,double high);
    void SetTemplate(std::string axis_title,int nbins_t,double low_t,double high_t,int nbins_r,double low_r,double high_r);

    void FillTruthHistograms(bool sig,double var_t,bool load_syst,double weight=1.0);
    void FillRecoHistograms(bool sel,double var_r,bool load_syst,double weight=1.0);
    void FillHistograms2D(bool sig,bool sel,double var_t,double var_r,bool load_syst,double weight=1.0);

    void AddSpecialUniv(std::string name);
    void FillSpecialTruthHistograms(std::string name,bool sig,double var_t,double weight);
    void FillSpecialRecoHistograms(std::string name,bool sel,double var_r,double weight);
    void FillSpecialHistograms2D(std::string name,bool sig,bool sel,double var_t,double var_r,double weight);

    void DBBW() { _divide_by_bin_width = true; }
    void ShapeOnly() { _shape_only = true; }
    void KeepOU(){ _keep_overflow_underflow_ = true; }
    void KeepAll(){ _keep_all = true; }

    void Write();
    
  private:

    const std::string _label;
    const bool _save_truth;
    bool _divide_by_bin_width = false;
    bool _shape_only = false; 
    bool _keep_overflow_underflow_ = false;
    bool _keep_all = false;

    // histogram templates
    TH1D* _h_tp = nullptr;
    TH1D* _h_tp_truth = nullptr;

    // Reco histograms 
    TH1D* _h_CV_Reco_Tot;
    std::vector<std::vector<TH1D*>> _h_Vars_Reco_Tot;
    std::vector<TH1D*> _h_Unisim_Vars_Reco_Tot;
    std::map<std::string,TH1D*> _h_Special_Reco_Tot;

    // Reco histograms by category
    std::vector<TH1D*> _h_CV_Reco_Cat;
    std::vector<std::vector<std::vector<TH1D*>>> _h_Vars_Reco_Cat;
    std::vector<std::vector<TH1D*>> _h_Unisim_Vars_Reco_Cat;
    std::map<std::string,std::vector<TH1D*>> _h_Special_Reco_Cat;
 
    // Truth histograms, only filled for the signal
    TH1D* _h_CV_Truth_Signal;
    std::vector<std::vector<TH1D*>> _h_Vars_Truth_Signal;
    std::vector<TH1D*> _h_Unisim_Vars_Truth_Signal;
    std::map<std::string,TH1D*> _h_Special_Truth_Signal;

    // Joint histograms, only filled for the signal
    TH2D* _h_CV_Joint_Signal;
    std::vector<std::vector<TH2D*>> _h_Vars_Joint_Signal;
    std::vector<TH2D*> _h_Unisim_Vars_Joint_Signal;
    std::map<std::string,TH2D*> _h_Special_Joint_Signal;

    int _nbins_t;    
    std::vector<double> _bin_edges_t;

    int _nbins_r;
    std::vector<double> _bin_edges_r;

    void _SetupHistograms();
    void _SetupRecoHistograms();
    void _SetupTruthHistograms();
    void _SetupJointHistograms();
    void _WriteReco();
    void _WriteTruth();
    void _WriteJoint();

    TFile* _f_out = nullptr;

};

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

HistogramManager::HistogramManager(std::string label,bool save_truth) : _label(label), _save_truth(save_truth)
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

  if(_save_truth){ 
    TFile* f_tp_truth = TFile::Open(("Analysis/"+_label+"/rootfiles/TruthBinningTemplate.root").c_str());
    _h_tp_truth = (TH1D*)f_tp_truth->Get(("h_template_"+_label).c_str());
    _h_tp_truth->SetDirectory(0);
    f_tp_truth->Close();
    _h_tp_truth->SetName(("h_template_truth_"+_label).c_str()); 
  }

  _SetupHistograms();
}

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void HistogramManager::SetTemplate(std::string axis_title,int nbins_r,double low_r,double high_r)
{
  if(_save_truth) 
    throw std::invalid_argument("HistogramManager::SetTemplate: Calling wrong function when response matrices are being set, need to specify truth binning as well");   
  _h_tp = new TH1D(("h_template_"+_label).c_str(),axis_title.c_str(),nbins_r,low_r,high_r);

  _SetupHistograms();

}


///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void HistogramManager::SetTemplate(std::string axis_title,int nbins_t,double low_t,double high_t,int nbins_r,double low_r,double high_r)
{
  if(!_save_truth) 
    throw std::invalid_argument("HistogramManager::SetTemplate: Calling wrong function when response matrices are being set, truth binning not enabled");  

  _h_tp = new TH1D(("h_template_"+_label).c_str(),axis_title.c_str(),nbins_r,low_r,high_r);
  _h_tp_truth = new TH1D(("h_template_truth_"+_label).c_str(),axis_title.c_str(),nbins_t,low_t,high_t);

  _SetupHistograms();

}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void HistogramManager::_SetupHistograms()
{
  _SetupRecoHistograms();
  if(_save_truth){
   _SetupTruthHistograms();
   _SetupJointHistograms();
  }
}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void HistogramManager::_SetupRecoHistograms()
{

  _bin_edges_r.clear();
  for(int i=1;i<_h_tp->GetNbinsX()+2;i++) _bin_edges_r.push_back(_h_tp->GetXaxis()->GetBinLowEdge(i));
  _nbins_r = _bin_edges_r.size()-1;

  // Storing totals
  _h_CV_Reco_Tot = (TH1D*)_h_tp->Clone(("h_CV_Reco_Tot_"+_label).c_str());
  _h_CV_Reco_Tot->Sumw2();

  for(size_t i_c=0;i_c<categories.size();i_c++)
    _h_CV_Reco_Cat.push_back((TH1D*)_h_tp->Clone(("h_CV_Reco_"+categories.at(i_c)+"_"+_label).c_str()));

  for(int i_s=0;i_s<kSystMAX;i_s++){
    _h_Vars_Reco_Tot.push_back(std::vector<TH1D*>());
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      _h_Vars_Reco_Tot.back().push_back((TH1D*)_h_tp->Clone(("h_Vars_Reco_Tot_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
    }
  }

  for(int i_s=0;i_s<kUnisimMAX;i_s++)
    _h_Unisim_Vars_Reco_Tot.push_back((TH1D*)_h_tp->Clone(("h_Unisim_Vars_Reco_Tot_"+unisims_str.at(i_s)+"_"+_label).c_str()));

  for(size_t i_c=0;i_c<categories.size();i_c++){
    _h_Vars_Reco_Cat.push_back(std::vector<std::vector<TH1D*>>());
    for(int i_s=0;i_s<kSystMAX;i_s++){
      _h_Vars_Reco_Cat.back().push_back(std::vector<TH1D*>());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        _h_Vars_Reco_Cat.back().back().push_back((TH1D*)_h_tp->Clone(("h_Vars_Reco_"+categories.at(i_c)+"_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
      }
    }
  }

  for(size_t i_c=0;i_c<categories.size();i_c++){
    _h_Unisim_Vars_Reco_Cat.push_back(std::vector<TH1D*>());
    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      _h_Unisim_Vars_Reco_Cat.at(i_c).push_back((TH1D*)_h_tp->Clone(("h_Unisim_Vars_Reco_"+categories.at(i_c)+"_"+unisims_str.at(i_s)+"_"+_label).c_str()));
    }
  }

}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void HistogramManager::_SetupTruthHistograms()
{

  _bin_edges_t.clear();
  for(int i=1;i<_h_tp_truth->GetNbinsX()+2;i++) _bin_edges_t.push_back(_h_tp_truth->GetXaxis()->GetBinLowEdge(i));
  _nbins_t = _bin_edges_t.size()-1;

  // Storing totals
  _h_CV_Truth_Signal = (TH1D*)_h_tp_truth->Clone(("h_CV_Truth_Signal_"+_label).c_str());
  _h_CV_Truth_Signal->Sumw2();

  for(int i_s=0;i_s<kSystMAX;i_s++){
    _h_Vars_Truth_Signal.push_back(std::vector<TH1D*>());
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      _h_Vars_Truth_Signal.back().push_back((TH1D*)_h_tp_truth->Clone(("h_Vars_Truth_Signal_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
    }
  }

  for(int i_s=0;i_s<kUnisimMAX;i_s++)
    _h_Unisim_Vars_Truth_Signal.push_back((TH1D*)_h_tp_truth->Clone(("h_Unisim_Vars_Truth_Signal_"+unisims_str.at(i_s)+"_"+_label).c_str()));
 
}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void HistogramManager::_SetupJointHistograms()
{

  // Storing totals
  std::string t_title = _h_tp_truth->GetXaxis()->GetTitle();
  std::string r_title = _h_tp->GetXaxis()->GetTitle();
  std::string title = ";"+t_title+";"+r_title+";";
  _h_CV_Joint_Signal = new TH2D(("h_CV_Joint_Signal_"+_label).c_str(),title.c_str(),_nbins_t,&_bin_edges_t[0],_nbins_r,&_bin_edges_r[0]);
  _h_CV_Joint_Signal->Sumw2();

  for(int i_s=0;i_s<kSystMAX;i_s++){
    _h_Vars_Joint_Signal.push_back(std::vector<TH2D*>());
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      _h_Vars_Joint_Signal.back().push_back((TH2D*)_h_CV_Joint_Signal->Clone(("h_Vars_Joint_Signal_"+sys_str.at(i_s)+"_"+std::to_string(i_u)+"_"+_label).c_str()));
    }
  }

  for(int i_s=0;i_s<kUnisimMAX;i_s++)
    _h_Unisim_Vars_Joint_Signal.push_back((TH2D*)_h_CV_Joint_Signal->Clone(("h_Unisim_Vars_Joint_Signal_"+unisims_str.at(i_s)+"_"+_label).c_str()));

}

///////////////////////////////////////////////////////////////////////
// Setup histograms for special alternative universe 

void HistogramManager::AddSpecialUniv(std::string name)
{

  _h_Special_Reco_Tot[name] = (TH1D*)_h_tp->Clone(("h_Special_Reco_Tot_"+name+"_"+_label).c_str()); 
  _h_Special_Reco_Cat[name] = std::vector<TH1D*>();
  for(size_t i_c=0;i_c<categories.size();i_c++)
    _h_Special_Reco_Cat.at(name).push_back((TH1D*)_h_tp->Clone(("h_Special_Reco_"+categories.at(i_c)+"_"+name+"_"+_label).c_str()));
 
  if(_save_truth){
    _h_Special_Truth_Signal[name] = (TH1D*)_h_tp->Clone(("h_Special_Truth_Signal_"+name+"_"+_label).c_str()); 
    _h_Special_Joint_Signal[name] = (TH2D*)_h_CV_Joint_Signal->Clone(("h_Special_Truth_Signal_"+name+"_"+_label).c_str()); 
  }
 
}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void HistogramManager::FillRecoHistograms(bool sel,double var_r,bool load_syst,double weight)
{

  if(!sel) return;

  if(_h_tp == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(category == -1) std::cout << "Bad event" << std::endl;

  if(std::isnan(weightSplineTimesTune) || std::isinf(weightSplineTimesTune)) return;

  if(!is_data) _h_CV_Reco_Tot->Fill(var_r,POT_weight*weightSplineTimesTune);
  _h_CV_Reco_Cat.at(category)->Fill(var_r,POT_weight*weightSplineTimesTune);

  if(!is_data){
    if(load_syst){
      for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Reco_Tot.at(kGenie).at(i_u)->Fill(var_r,(double)POT_weight*weightsGenie->at(i_u)/1000);
      for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Reco_Tot.at(kReint).at(i_u)->Fill(var_r,(double)POT_weight*weightsReint->at(i_u)*weightSplineTimesTune/1000);
      for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Reco_Tot.at(kFlux).at(i_u)->Fill(var_r,(double)POT_weight*weightsFlux->at(i_u)*weightSplineTimesTune/1000);
      for(int i_s=0;i_s<kUnisimMAX;i_s++) _h_Unisim_Vars_Reco_Tot.at(i_s)->Fill(var_r,(double)POT_weight*ChooseUnisimWeight(i_s,weightsUnisim));
    }
    else {
      for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Reco_Tot.at(kGenie).at(i_u)->Fill(var_r,(double)POT_weight*weightSplineTimesTune);
      for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Reco_Tot.at(kReint).at(i_u)->Fill(var_r,(double)POT_weight*weightSplineTimesTune);
      for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Reco_Tot.at(kFlux).at(i_u)->Fill(var_r,(double)POT_weight*weightSplineTimesTune);
      for(int i_s=0;i_s<kUnisimMAX;i_s++) _h_Unisim_Vars_Reco_Tot.at(i_s)->Fill(var_r,(double)POT_weight*weightSplineTimesTune);
    } 
  }

  if(load_syst){
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Reco_Cat.at(category).at(kGenie).at(i_u)->Fill(var_r,(double)POT_weight*weightsGenie->at(i_u)/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Reco_Cat.at(category).at(kReint).at(i_u)->Fill(var_r,(double)POT_weight*weightsReint->at(i_u)*weightSplineTimesTune/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Reco_Cat.at(category).at(kFlux).at(i_u)->Fill(var_r,(double)POT_weight*weightsFlux->at(i_u)*weightSplineTimesTune/1000);
    for(int i_s=0;i_s<kUnisimMAX;i_s++) _h_Unisim_Vars_Reco_Cat.at(category).at(i_s)->Fill(var_r,(double)POT_weight*ChooseUnisimWeight(i_s,weightsUnisim));
  }
  else {
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Reco_Cat.at(category).at(kGenie).at(i_u)->Fill(var_r,(double)POT_weight*weightSplineTimesTune);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Reco_Cat.at(category).at(kReint).at(i_u)->Fill(var_r,(double)POT_weight*weightSplineTimesTune);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Reco_Cat.at(category).at(kFlux).at(i_u)->Fill(var_r,(double)POT_weight*weightSplineTimesTune);
    for(int i_s=0;i_s<kUnisimMAX;i_s++) _h_Unisim_Vars_Reco_Cat.at(category).at(i_s)->Fill(var_r,(double)POT_weight*weightSplineTimesTune);
  } 

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void HistogramManager::FillTruthHistograms(bool sig,double var_t,bool load_syst,double weight=1.0)
{

  if(!sig) return;

  if(_h_tp_truth == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(category == -1) std::cout << "Bad event" << std::endl;

  if(std::isnan(weightSplineTimesTune) || std::isinf(weightSplineTimesTune)) return;

  _h_CV_Truth_Signal->Fill(var_t,POT_weight*weightSplineTimesTune);

  if(load_syst){
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Truth_Signal.at(kGenie).at(i_u)->Fill(var_t,(double)POT_weight*weightsGenie->at(i_u)/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Truth_Signal.at(kReint).at(i_u)->Fill(var_t,(double)POT_weight*weightsReint->at(i_u)*weightSplineTimesTune/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Truth_Signal.at(kFlux).at(i_u)->Fill(var_t,(double)POT_weight*weightsFlux->at(i_u)*weightSplineTimesTune/1000);
    for(int i_s=0;i_s<kUnisimMAX;i_s++) _h_Unisim_Vars_Truth_Signal.at(i_s)->Fill(var_t,(double)POT_weight*ChooseUnisimWeight(i_s,weightsUnisim));
  }
  else {
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Truth_Signal.at(kGenie).at(i_u)->Fill(var_t,(double)POT_weight*weightSplineTimesTune);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Truth_Signal.at(kReint).at(i_u)->Fill(var_t,(double)POT_weight*weightSplineTimesTune);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Truth_Signal.at(kFlux).at(i_u)->Fill(var_t,(double)POT_weight*weightSplineTimesTune);
    for(int i_s=0;i_s<kUnisimMAX;i_s++) _h_Unisim_Vars_Truth_Signal.at(i_s)->Fill(var_t,(double)POT_weight*weightSplineTimesTune);
  } 

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void HistogramManager::FillHistograms2D(bool sig,bool sel,double var_t,double var_r,bool load_syst,double weight)
{

 FillTruthHistograms(sig,var_t,load_syst,weight);
 FillRecoHistograms(sel,var_r,load_syst,weight);

 if(!sig || !sel) return;

 if(_h_tp == nullptr || _h_tp_truth == nullptr) 
   throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

 if(category == -1) std::cout << "Bad event" << std::endl;

 if(std::isnan(weightSplineTimesTune) || std::isinf(weightSplineTimesTune)) return;

 _h_CV_Joint_Signal->Fill(var_t,var_r,POT_weight*weightSplineTimesTune);

  if(load_syst){
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Joint_Signal.at(kGenie).at(i_u)->Fill(var_t,var_r,(double)POT_weight*weightsGenie->at(i_u)/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Joint_Signal.at(kReint).at(i_u)->Fill(var_t,var_r,(double)POT_weight*weightsReint->at(i_u)*weightSplineTimesTune/1000);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Joint_Signal.at(kFlux).at(i_u)->Fill(var_t,var_r,(double)POT_weight*weightsFlux->at(i_u)*weightSplineTimesTune/1000);
    for(int i_s=0;i_s<kUnisimMAX;i_s++) _h_Unisim_Vars_Joint_Signal.at(i_s)->Fill(var_t,var_r,(double)POT_weight*ChooseUnisimWeight(i_s,weightsUnisim));
  }
  else {
    for(int i_u=0;i_u<sys_nuniv.at(kGenie);i_u++) _h_Vars_Joint_Signal.at(kGenie).at(i_u)->Fill(var_t,var_r,(double)POT_weight*weightSplineTimesTune);
    for(int i_u=0;i_u<sys_nuniv.at(kReint);i_u++) _h_Vars_Joint_Signal.at(kReint).at(i_u)->Fill(var_t,var_r,(double)POT_weight*weightSplineTimesTune);
    for(int i_u=0;i_u<sys_nuniv.at(kFlux);i_u++) _h_Vars_Joint_Signal.at(kFlux).at(i_u)->Fill(var_t,var_r,(double)POT_weight*weightSplineTimesTune);
    for(int i_s=0;i_s<kUnisimMAX;i_s++) _h_Unisim_Vars_Joint_Signal.at(i_s)->Fill(var_t,var_r,(double)POT_weight*weightSplineTimesTune);
  } 

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void HistogramManager::FillSpecialTruthHistograms(std::string name,bool sig,double var_t,double weight)
{

  if(!sig) return;

  if(_h_tp_truth == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(category == -1) std::cout << "Bad event" << std::endl;

  _h_Special_Truth_Signal.at(name)->Fill(var_t,POT_weight*weightSplineTimesTune*weight);

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void HistogramManager::FillSpecialRecoHistograms(std::string name,bool sel,double var_r,double weight)
{

  if(!sel) return;

  if(_h_tp == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(category == -1) std::cout << "Bad event" << std::endl;

  if(std::isnan(weightSplineTimesTune) || std::isinf(weightSplineTimesTune)) return;

  if(!is_data) _h_Special_Reco_Tot.at(name)->Fill(var_r,POT_weight*weightSplineTimesTune*weight);
  _h_Special_Reco_Cat.at(name).at(category)->Fill(var_r,POT_weight*weightSplineTimesTune*weight);

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void HistogramManager::FillSpecialHistograms2D(std::string name,bool sig,bool sel,double var_t,double var_r,double weight)
{

  FillSpecialTruthHistograms(name,sig,var_t,weight);
  FillSpecialRecoHistograms(name,sel,var_r,weight);

  if(!sig || !sel) return;

  if(_h_tp == nullptr || _h_tp_truth == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(category == -1) std::cout << "Bad event" << std::endl;

  if(std::isnan(weightSplineTimesTune) || std::isinf(weightSplineTimesTune)) return;

  _h_Special_Joint_Signal.at(name)->Fill(var_t,var_r,POT_weight*weightSplineTimesTune*weight);

}

///////////////////////////////////////////////////////////////////////
// Write the histograms to file

void HistogramManager::Write()
{
  std::cout << "Writing histograms for " << _label << std::endl;
  gSystem->Exec(("mkdir -p Analysis/"+_label+"/rootfiles/").c_str());
  _f_out = TFile::Open(("Analysis/"+_label+"/rootfiles/Histograms.root").c_str(),"RECREATE");
  _WriteReco();

  if(_save_truth){
    _WriteTruth();
    _WriteJoint();
  }

  _f_out->Close();
}

///////////////////////////////////////////////////////////////////////
// Write the reco histograms to file

void HistogramManager::_WriteReco()
{
  _f_out->mkdir("Reco");
  _f_out->cd("Reco");

  // Make copies of histograms with error equal to bin count instead of sqrt(sum(weights^2)) for data stat error
  TH1D* h_CV_Reco_Tot_E = (TH1D*)_h_CV_Reco_Tot->Clone("_h_CV_Reco_Tot_E");
  for(int i=0;i<_h_CV_Reco_Tot->GetNbinsX()+2;i++) h_CV_Reco_Tot_E->SetBinError(i,sqrt(_h_CV_Reco_Tot->GetBinContent(i)));
  std::vector<TH1D*> h_CV_Reco_Cat_E;
  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_CV_Reco_Cat_E.push_back((TH1D*)_h_CV_Reco_Cat.at(i_c)->Clone(("h_CV_Reco_"+categories.at(i_c)+"_E").c_str()));
    for(int i=0;i<_h_CV_Reco_Tot->GetNbinsX()+2;i++) h_CV_Reco_Cat_E.back()->SetBinError(i,sqrt(h_CV_Reco_Cat_E.back()->GetBinContent(i)));
  }

  // Histograms to store totals
  TH2D* h_Cov = Make2DHist("h_Cov",_h_tp); 
  TH2D* h_FCov = Make2DHist("h_FCov",_h_tp); 
  std::vector<TH2D*> h_Cov_Cat,h_FCov_Cat;
 
  TH2D* h_Cov_Sys = Make2DHist("h_Cov_Sys",_h_tp); 
  TH2D* h_FCov_Sys = Make2DHist("h_FCov_Sys",_h_tp); 
  std::vector<TH2D*> h_Cov_Sys_Cat,h_FCov_Sys_Cat;

  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_Cov_Cat.push_back(Make2DHist("h_Cov_"+categories.at(i_c),_h_tp));
    h_FCov_Cat.push_back(Make2DHist("h_FCov_"+categories.at(i_c),_h_tp));
    h_Cov_Sys_Cat.push_back(Make2DHist("h_Cov_Sys_"+categories.at(i_c),_h_tp));
    h_FCov_Sys_Cat.push_back(Make2DHist("h_FCov_Sys_"+categories.at(i_c),_h_tp));
  }

  // Apply scaling to the reco if shape only
  if(_shape_only){
    double integral = _h_CV_Reco_Tot->Integral() + _h_CV_Reco_Tot->GetBinContent(0) + _h_CV_Reco_Tot->GetBinContent(_nbins_r+1);
    _h_CV_Reco_Tot->Scale(1.0/integral);
    h_CV_Reco_Tot_E->Scale(1.0/integral);
    for(size_t i_c=0;i_c<categories.size();i_c++) _h_CV_Reco_Cat.at(i_c)->Scale(1.0/integral);
    for(size_t i_c=0;i_c<categories.size();i_c++) h_CV_Reco_Cat_E.at(i_c)->Scale(1.0/integral);
    for(size_t i_s=0;i_s<kSystMAX;i_s++){
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        TH1D* h =_h_Vars_Reco_Tot.at(i_s).at(i_u); 
        integral = h->Integral() + h->GetBinContent(0) + h->GetBinContent(_nbins_r+1);
        h->Scale(1.0/integral);
        for(size_t i_c=0;i_c<categories.size();i_c++){
          _h_Vars_Reco_Cat.at(i_c).at(i_s).at(i_u)->Scale(1.0/integral);
        }
      } 
    }
    for(size_t i_s=0;i_s<kUnisimMAX;i_s++){
      TH1D* h =_h_Unisim_Vars_Reco_Tot.at(i_s); 
      integral = h->Integral() + h->GetBinContent(0) + h->GetBinContent(_nbins_r+1);
      h->Scale(1.0/integral);
      for(size_t i_c=0;i_c<categories.size();i_c++){
        _h_Unisim_Vars_Reco_Cat.at(i_c).at(i_s)->Scale(1.0/integral);
      }
    }
  } 

  if(_divide_by_bin_width){
    DivideByBinWidth(_h_CV_Reco_Tot);
    DivideByBinWidth(h_CV_Reco_Tot_E);
    for(size_t i_c=0;i_c<categories.size();i_c++) DivideByBinWidth(_h_CV_Reco_Cat.at(i_c));
    for(size_t i_c=0;i_c<categories.size();i_c++) DivideByBinWidth(h_CV_Reco_Cat_E.at(i_c));
    for(size_t i_s=0;i_s<kSystMAX;i_s++){
      for(size_t i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        DivideByBinWidth(_h_Vars_Reco_Tot.at(i_s).at(i_u));
        for(size_t i_c=0;i_c<categories.size();i_c++) DivideByBinWidth(_h_Vars_Reco_Cat.at(i_c).at(i_s).at(i_u));
      }
    } 
    for(size_t i_s=0;i_s<kUnisimMAX;i_s++){
      DivideByBinWidth(_h_Unisim_Vars_Reco_Tot.at(i_s));
      for(size_t i_c=0;i_c<categories.size();i_c++) DivideByBinWidth(_h_Unisim_Vars_Reco_Cat.at(i_c).at(i_s));
    } 
  }

 
  // Central value reco histograms
  _f_out->mkdir("Reco/CV");
  _f_out->cd("Reco/CV");
  //_f_out->mkdir("CV/Reco");
  //_f_out->cd("CV/Reco"); 
  _h_CV_Reco_Tot->Write("h_Tot");  
  for(size_t i_c=0;i_c<categories.size();i_c++) _h_CV_Reco_Cat.at(i_c)->Write(("h_"+categories.at(i_c)).c_str());
  _f_out->cd();

  // Variation reco histograms
  if(_keep_all){
    _f_out->mkdir("Reco/Vars");
    for(int i_s=0;i_s<kSystMAX;i_s++){
      _f_out->mkdir(("Reco/Vars/"+sys_str.at(i_s)).c_str());
      _f_out->cd(("Reco/Vars/"+sys_str.at(i_s)).c_str());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        _h_Vars_Reco_Tot.at(i_s).at(i_u)->Write(Form("h_Tot_%i",i_u));
        for(size_t i_c=0;i_c<categories.size();i_c++) 
          _h_Vars_Reco_Cat.at(i_c).at(i_s).at(i_u)->Write(("h_"+categories.at(i_c)+"_"+std::to_string(i_u)).c_str());
      }
      _f_out->cd();
    }
    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      _f_out->mkdir(("Reco/Vars/"+unisims_str.at(i_s)).c_str());
      _f_out->cd(("Reco/Vars/"+unisims_str.at(i_s)).c_str());
      _h_Unisim_Vars_Reco_Tot.at(i_s)->Write("h_Tot");
      for(size_t i_c=0;i_c<categories.size();i_c++) 
        _h_Unisim_Vars_Reco_Cat.at(i_c).at(i_s)->Write(("h_"+categories.at(i_c)).c_str());
      _f_out->cd();
    }
  }

  // Reco covariance matrices - multisims 
  _f_out->mkdir("Reco/Cov");
  for(int i_s=0;i_s<kSystMAX;i_s++){
    _f_out->mkdir(("Reco/Cov/"+sys_str.at(i_s)).c_str());
    _f_out->cd(("Reco/Cov/"+sys_str.at(i_s)).c_str());
    TH2D *C,*FC;
    CalcCovMultisim(sys_str.at(i_s),_h_Vars_Reco_Tot.at(i_s),C,FC);
    C->Write("Cov_Tot");
    FC->Write("FCov_Tot");
    h_Cov->Add(C);
    h_FCov->Add(FC);
    h_Cov_Sys->Add(C);
    h_FCov_Sys->Add(FC);
    for(size_t i_c=0;i_c<categories.size();i_c++){
      TH2D *C_Cat,*FC_Cat;
      CalcCovMultisim(sys_str.at(i_s)+"_"+categories.at(i_c),_h_Vars_Reco_Cat.at(i_c).at(i_s),C_Cat,FC_Cat);
      C_Cat->Write(("Cov_"+categories.at(i_c)).c_str());
      FC_Cat->Write(("FCov_"+categories.at(i_c)).c_str());
      h_Cov_Cat.at(i_c)->Add(C_Cat);
      h_FCov_Cat.at(i_c)->Add(FC_Cat);
      h_Cov_Sys_Cat.at(i_c)->Add(C_Cat);
      h_FCov_Sys_Cat.at(i_c)->Add(FC_Cat);
    } 
    _f_out->cd();
  } 

  // Reco covariance matrices - unisims 
  for(int i_s=0;i_s<kUnisimMAX;i_s++){
    _f_out->mkdir(("Reco/Cov/"+unisims_str.at(i_s)).c_str());
    _f_out->cd(("Reco/Cov/"+unisims_str.at(i_s)).c_str());
    TH2D *C,*FC;
    CalcCovUnisim(unisims_str.at(i_s),_h_CV_Reco_Tot,_h_Unisim_Vars_Reco_Tot.at(i_s),C,FC);
    C->Write("Cov_Tot");
    FC->Write("FCov_Tot");
    h_Cov->Add(C);
    h_FCov->Add(FC);
    h_Cov_Sys->Add(C);
    h_FCov_Sys->Add(FC);
    for(size_t i_c=0;i_c<categories.size();i_c++){
      TH2D *C_Cat,*FC_Cat;
      CalcCovUnisim(unisims_str.at(i_s)+"_"+categories.at(i_c),_h_CV_Reco_Cat.at(i_c),_h_Unisim_Vars_Reco_Cat.at(i_c).at(i_s),C_Cat,FC_Cat);
      C_Cat->Write(("Cov_"+categories.at(i_c)).c_str());
      FC_Cat->Write(("FCov_"+categories.at(i_c)).c_str());
      h_Cov_Cat.at(i_c)->Add(C_Cat);
      h_FCov_Cat.at(i_c)->Add(FC_Cat);
      h_Cov_Sys_Cat.at(i_c)->Add(C_Cat);
      h_FCov_Sys_Cat.at(i_c)->Add(FC_Cat);
    } 
    _f_out->cd();
  } 

  // Reco stat covariance matrices
  _f_out->mkdir("Reco/Cov/MCStat");
  _f_out->cd("Reco/Cov/MCStat");
  TH2D* h_Cov_MCStat_Reco = Make2DHist("Cov_MCStat",_h_tp);
  TH2D* h_FCov_MCStat_Reco = Make2DHist("FCov_MCStat",_h_tp);
  for(int i=0;i<h_Cov_MCStat_Reco->GetNbinsX()+2;i++){
    h_Cov_MCStat_Reco->SetBinContent(i,i,pow(_h_CV_Reco_Tot->GetBinError(i),2));
    h_FCov_MCStat_Reco->SetBinContent(i,i,pow(_h_CV_Reco_Tot->GetBinError(i)/_h_CV_Reco_Tot->GetBinContent(i),2));
  }
  h_Cov->Add(h_Cov_MCStat_Reco);
  h_FCov->Add(h_FCov_MCStat_Reco);
  h_Cov_MCStat_Reco->Write("Cov_Tot"); 
  h_FCov_MCStat_Reco->Write("FCov_Tot");

  for(size_t i_c=0;i_c<categories.size();i_c++){
    TH2D* h_Cov_MCStat_Reco_Cat = Make2DHist("Cov_MCStat_"+categories.at(i_c),_h_tp);
    TH2D* h_FCov_MCStat_Reco_Cat = Make2DHist("FCov_MCStat_"+categories.at(i_c),_h_tp);
    for(int i=0;i<h_Cov_MCStat_Reco->GetNbinsX()+2;i++){
      h_Cov_MCStat_Reco_Cat->SetBinContent(i,i,pow(_h_CV_Reco_Cat.at(i_c)->GetBinError(i),2));
      h_FCov_MCStat_Reco_Cat->SetBinContent(i,i,pow(_h_CV_Reco_Cat.at(i_c)->GetBinError(i)/_h_CV_Reco_Cat.at(i_c)->GetBinContent(i),2));
    }
    h_Cov_MCStat_Reco_Cat->Write(("Cov_"+categories.at(i_c)).c_str());
    h_FCov_MCStat_Reco_Cat->Write(("FCov_"+categories.at(i_c)).c_str());
    h_Cov_Cat.at(i_c)->Add(h_Cov_MCStat_Reco_Cat);
    h_FCov_Cat.at(i_c)->Add(h_FCov_MCStat_Reco_Cat);
  }  
  _f_out->cd();

  _f_out->mkdir("Reco/Cov/EstDataStat");
  _f_out->cd("Reco/Cov/EstDataStat");
  TH2D* h_Cov_EstDataStat_Reco = Make2DHist("Cov_EstDataStat",_h_tp);
  TH2D* h_FCov_EstDataStat_Reco = Make2DHist("FCov_EstDataStat",_h_tp);
  for(int i=0;i<h_Cov_EstDataStat_Reco->GetNbinsX()+2;i++){
    h_Cov_EstDataStat_Reco->SetBinContent(i,i,pow(h_CV_Reco_Tot_E->GetBinError(i),2));
    h_FCov_EstDataStat_Reco->SetBinContent(i,i,pow(h_CV_Reco_Tot_E->GetBinError(i)/h_CV_Reco_Tot_E->GetBinContent(i),2));
  }
  h_Cov_EstDataStat_Reco->Write("Cov_Tot"); 
  h_FCov_EstDataStat_Reco->Write("FCov_Tot");

  for(size_t i_c=0;i_c<categories.size();i_c++){
    TH2D* h_Cov_EstDataStat_Reco_Cat = Make2DHist("Cov_EstDataStat_"+categories.at(i_c),_h_tp);
    TH2D* h_FCov_EstDataStat_Reco_Cat = Make2DHist("FCov_EstDataStat_"+categories.at(i_c),_h_tp);
    for(int i=0;i<h_Cov_EstDataStat_Reco->GetNbinsX()+2;i++){
      h_Cov_EstDataStat_Reco_Cat->SetBinContent(i,i,pow(h_CV_Reco_Cat_E.at(i_c)->GetBinError(i),2));
      h_FCov_EstDataStat_Reco_Cat->SetBinContent(i,i,pow(h_CV_Reco_Cat_E.at(i_c)->GetBinError(i)/h_CV_Reco_Cat_E.at(i_c)->GetBinContent(i),2));
    }
    h_Cov_EstDataStat_Reco_Cat->Write(("Cov_"+categories.at(i_c)).c_str());
    h_FCov_EstDataStat_Reco_Cat->Write(("FCov_"+categories.at(i_c)).c_str());
  }  
  _f_out->cd();

  _f_out->mkdir("Reco/Cov/Total");
  _f_out->cd("Reco/Cov/Total");
  h_Cov->Write("Cov_Tot");
  h_FCov->Write("FCov_Tot");
  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_Cov_Cat.at(i_c)->Write(("Cov_"+categories.at(i_c)).c_str());
    h_FCov_Cat.at(i_c)->Write(("FCov_"+categories.at(i_c)).c_str());
  }
  
  _f_out->mkdir("Reco/Cov/Sys");
  _f_out->cd("Reco/Cov/Sys");
  h_Cov_Sys->Write("Cov_Tot");
  h_FCov_Sys->Write("FCov_Tot");
  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_Cov_Sys_Cat.at(i_c)->Write(("Cov_"+categories.at(i_c)).c_str());
    h_FCov_Sys_Cat.at(i_c)->Write(("FCov_"+categories.at(i_c)).c_str());
  }

  _f_out->cd();

  if(_h_Special_Reco_Tot.size()){
    _f_out->mkdir("Reco/Special");
    std::map<std::string,TH1D*>::iterator it;
    for(it=_h_Special_Reco_Tot.begin();it!=_h_Special_Reco_Tot.end();it++){
      std::string name = it->first;
      _f_out->mkdir(("Reco/Special/"+name).c_str());
      _f_out->cd(("Reco/Special/"+name).c_str());
      _h_Special_Reco_Tot.at(name)->Write("h_Tot");
      for(size_t i_c=0;i_c<categories.size();i_c++)
        _h_Special_Reco_Cat.at(name).at(i_c)->Write(("h_"+categories.at(i_c)).c_str());
    }
    _f_out->cd();
  } 

}

///////////////////////////////////////////////////////////////////////
// Write the truth histograms to file

void HistogramManager::_WriteTruth()
{

  _f_out->cd();
  _f_out->mkdir("Truth");
  _f_out->cd("Truth");

  // Make copies of histograms with error equal to bin count instead of sqrt(sum(weights^2)) for data stat error
  TH1D* h_CV_Truth_Signal_E = (TH1D*)_h_CV_Truth_Signal->Clone("_h_CV_Truth_Signal_E");
  for(int i=0;i<_h_CV_Truth_Signal->GetNbinsX()+2;i++) h_CV_Truth_Signal_E->SetBinError(i,sqrt(_h_CV_Truth_Signal->GetBinContent(i)));

  // Histograms to store totals
  TH2D* h_Cov = Make2DHist("h_Cov_Truth",_h_tp_truth); 
  TH2D* h_FCov = Make2DHist("h_FCov_Truth",_h_tp_truth); 
 
  TH2D* h_Cov_Sys = Make2DHist("h_Cov_Truth_Sys",_h_tp_truth); 
  TH2D* h_FCov_Sys = Make2DHist("h_FCov_Truth_Sys",_h_tp_truth); 

  // Apply scaling to the reco if shape only
  if(_shape_only){
    double integral = _h_CV_Truth_Signal->Integral() + _h_CV_Truth_Signal->GetBinContent(0) + _h_CV_Truth_Signal->GetBinContent(_nbins_t+1);
    _h_CV_Truth_Signal->Scale(1.0/integral);
    h_CV_Truth_Signal_E->Scale(1.0/integral);
    for(int i_s=0;i_s<kSystMAX;i_s++){
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        TH1D* h =_h_Vars_Truth_Signal.at(i_s).at(i_u); 
        integral = h->Integral() + h->GetBinContent(0) + h->GetBinContent(_nbins_t+1);
        h->Scale(1.0/integral);
      } 
    }
    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      TH1D* h =_h_Unisim_Vars_Truth_Signal.at(i_s); 
      integral = h->Integral() + h->GetBinContent(0) + h->GetBinContent(_nbins_t+1);
      h->Scale(1.0/integral);
    }
  } 

  if(_divide_by_bin_width){
    DivideByBinWidth(_h_CV_Truth_Signal);
    DivideByBinWidth(h_CV_Truth_Signal_E);
    for(int i_s=0;i_s<kSystMAX;i_s++){
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        DivideByBinWidth(_h_Vars_Truth_Signal.at(i_s).at(i_u));
      }
    } 
    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      DivideByBinWidth(_h_Unisim_Vars_Truth_Signal.at(i_s));
    } 
  }

  _f_out->mkdir("Truth/CV"); 
  _f_out->cd("Truth/CV");
  _h_CV_Truth_Signal->Write("h_Signal");  
  _f_out->cd();

  // Variation truth histograms
  if(_keep_all){
    _f_out->mkdir("Truth/Vars");
    for(int i_s=0;i_s<kSystMAX;i_s++){
      _f_out->mkdir(("Truth/Vars/"+sys_str.at(i_s)).c_str());
      _f_out->cd(("Truth/Vars/"+sys_str.at(i_s)).c_str());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        _h_Vars_Truth_Signal.at(i_s).at(i_u)->Write(Form("h_Signal_%i",i_u));
      }
      _f_out->cd();
    }
    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      _f_out->mkdir(("Truth/Vars/"+unisims_str.at(i_s)).c_str());
      _f_out->cd(("Truth/Vars/"+unisims_str.at(i_s)).c_str());
      _h_Unisim_Vars_Truth_Signal.at(i_s)->Write("h_Signal");
      _f_out->cd();
    }
  }

  // Truth covariance matrices - multisim 
  _f_out->mkdir("Truth/Cov");
  for(int i_s=0;i_s<kSystMAX;i_s++){
    _f_out->mkdir(("Truth/Cov/"+sys_str.at(i_s)).c_str());
    _f_out->cd(("Truth/Cov/"+sys_str.at(i_s)).c_str());
    TH2D *C,*FC;
    CalcCovMultisim(sys_str.at(i_s),_h_Vars_Truth_Signal.at(i_s),C,FC);
    C->Write("Cov_Signal");
    FC->Write("FCov_Signal");
    h_Cov->Add(C);
    h_FCov->Add(FC);
    h_Cov_Sys->Add(C);
    h_FCov_Sys->Add(FC);
    _f_out->cd();
  } 

  // Truth covariance matrices - unisim 
  for(int i_s=0;i_s<kUnisimMAX;i_s++){
    _f_out->mkdir(("Truth/Cov/"+unisims_str.at(i_s)).c_str());
    _f_out->cd(("Truth/Cov/"+unisims_str.at(i_s)).c_str());
    TH2D *C,*FC;
    CalcCovUnisim(unisims_str.at(i_s),_h_CV_Truth_Signal,_h_Unisim_Vars_Truth_Signal.at(i_s),C,FC);
    C->Write("Cov_Signal");
    FC->Write("FCov_Signal");
    h_Cov->Add(C);
    h_FCov->Add(FC);
    h_Cov_Sys->Add(C);
    h_FCov_Sys->Add(FC);
    _f_out->cd();
  } 

  // Truth stat covariance matrices
  _f_out->mkdir("Truth/Cov/MCStat");
  _f_out->cd("Truth/Cov/MCStat");
  TH2D* h_Cov_MCStat_Truth = Make2DHist("Cov_MCStat",_h_tp_truth);
  TH2D* h_FCov_MCStat_Truth = Make2DHist("FCov_MCStat",_h_tp_truth);
  for(int i=0;i<h_Cov_MCStat_Truth->GetNbinsX()+2;i++){
    h_Cov_MCStat_Truth->SetBinContent(i,i,pow(_h_CV_Truth_Signal->GetBinError(i),2));
    h_FCov_MCStat_Truth->SetBinContent(i,i,pow(_h_CV_Truth_Signal->GetBinError(i)/_h_CV_Truth_Signal->GetBinContent(i),2));
  }
  h_Cov->Add(h_Cov_MCStat_Truth);
  h_FCov->Add(h_FCov_MCStat_Truth);
  h_Cov_MCStat_Truth->Write("Cov_Signal"); 
  h_FCov_MCStat_Truth->Write("FCov_Signal");

  _f_out->cd();
  _f_out->mkdir("Truth/Cov/EstDataStat");
  _f_out->cd("Truth/Cov/EstDataStat");
  TH2D* h_Cov_EstDataStat_Truth = Make2DHist("Cov_EstDataStat",_h_tp);
  TH2D* h_FCov_EstDataStat_Truth = Make2DHist("FCov_EstDataStat",_h_tp);
  for(int i=0;i<h_Cov_EstDataStat_Truth->GetNbinsX()+2;i++){
    h_Cov_EstDataStat_Truth->SetBinContent(i,i,pow(h_CV_Truth_Signal_E->GetBinError(i),2));
    h_FCov_EstDataStat_Truth->SetBinContent(i,i,pow(h_CV_Truth_Signal_E->GetBinError(i)/h_CV_Truth_Signal_E->GetBinContent(i),2));
  }
  h_Cov_EstDataStat_Truth->Write("Cov_Signal"); 
  h_FCov_EstDataStat_Truth->Write("FCov_Signal");

  _f_out->mkdir("Truth/Cov/Total");
  _f_out->cd("Truth/Cov/Total");
  h_Cov->Write("Cov_Signal");
  h_FCov->Write("FCov_Signal");

  _f_out->mkdir("Truth/Cov/Sys");
  _f_out->cd("Truth/Cov/Sys");
  h_Cov_Sys->Write("Cov_Signal");
  h_FCov_Sys->Write("FCov_Signal");

  _f_out->cd();

  if(_h_Special_Truth_Signal.size()){
    _f_out->mkdir("Truth/Special");
    _f_out->cd("Truth/Special");
    std::map<std::string,TH1D*>::iterator it;
    for(it=_h_Special_Truth_Signal.begin();it!=_h_Special_Truth_Signal.end();it++){
      std::string name = it->first;
      _f_out->mkdir(("Truth/Special/"+name).c_str());
      _f_out->cd(("Truth/Special/"+name).c_str());
      _h_Special_Truth_Signal.at(name)->Write("h_Signal");
    }
  } 

}

///////////////////////////////////////////////////////////////////////
// Write the truth histograms to file

void HistogramManager::_WriteJoint()
{

  _f_out->cd();
  _f_out->mkdir("Joint");

  if(_divide_by_bin_width){
    for(int i_t=0;i_t<_nbins_t;i_t++){
      for(int i_r=0;i_r<_nbins_t;i_r++){
        double a = _h_tp->GetBinWidth(i_r)*_h_tp_truth->GetBinWidth(i_t);
        _h_CV_Joint_Signal->SetBinContent(i_t,i_r,_h_CV_Joint_Signal->GetBinContent(i_t,i_r)/a);
        for(int i_s=0;i_s<kSystMAX;i_s++){
          for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
            _h_Vars_Joint_Signal.at(i_s).at(i_u)->SetBinContent(i_t,i_r,_h_Vars_Joint_Signal.at(i_s).at(i_u)->GetBinContent(i_t,i_r)/a);
          }
        }
        for(int i_s=0;i_s<kUnisimMAX;i_s++){
          _h_Unisim_Vars_Joint_Signal.at(i_s)->SetBinContent(i_t,i_r,_h_Unisim_Vars_Joint_Signal.at(i_s)->GetBinContent(i_t,i_r)/a);
        }
      }
    }
  }

  _f_out->mkdir("Joint/CV");
  _f_out->cd("Joint/CV");
  _h_CV_Joint_Signal->Write("h_Signal");  
  _f_out->cd();

  // Variation truth histograms
  if(_keep_all){
    _f_out->mkdir("Joint/Vars");
    for(int i_s=0;i_s<kSystMAX;i_s++){
      _f_out->mkdir(("Joint/Vars/"+sys_str.at(i_s)).c_str());
      _f_out->cd(("Joint/Vars/"+sys_str.at(i_s)).c_str());
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        _h_Vars_Joint_Signal.at(i_s).at(i_u)->Write(Form("h_Signal_%i",i_u));
      }
      _f_out->cd();
    }
    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      _f_out->mkdir(("Joint/Vars/"+unisims_str.at(i_s)).c_str());
      _f_out->cd(("Joint/Vars/"+unisims_str.at(i_s)).c_str());
      _h_Unisim_Vars_Joint_Signal.at(i_s)->Write("h_Signal");
      _f_out->cd();
    }
  }

  if(_h_Special_Joint_Signal.size()){
    _f_out->mkdir("Joint/Special");
    _f_out->cd("Joint/Special");
    std::map<std::string,TH2D*>::iterator it;
    for(it=_h_Special_Joint_Signal.begin();it!=_h_Special_Joint_Signal.end();it++){
      std::string name = it->first;
      _f_out->mkdir(("Joint/Special/"+name).c_str());
      _f_out->cd(("Joint/Special/"+name).c_str());
      _h_Special_Joint_Signal.at(name)->Write("h_Signal");
    }
  } 
  
}

///////////////////////////////////////////////////////////////////////

}


#endif
