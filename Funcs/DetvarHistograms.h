#ifndef _DetvarHistograms_h_
#define _DetvarHistograms_h_

#include "Systematics.h"
#include "BranchList.h"

using namespace syst;

namespace hist {

///////////////////////////////////////////////////////////////////////
// Class for managing detvar histograms

class DetvarHistogramManager {

  public: 

    DetvarHistogramManager(std::string label);
    void LoadTemplate();
    void SetTemplate(std::string axis_title,int nbins,double low,double high);
    void FillHistograms(double var,double weight=1.0);
    void DBBW() { _divide_by_bin_width = true; }
    void ShapeOnly() { _shape_only = true; }
    void Write();

 private:

    const std::string _label;
    bool _divide_by_bin_width = false;
    bool _shape_only = false; 
 
    TH1D* _h_tp = nullptr;
    TH1D* _h_CV_Tot;
    std::vector<TH1D*> _h_CV;
    std::vector<TH1D*> _h_Vars_Tot;
    std::vector<std::vector<TH1D*>> _h_Vars;

    void SetupHistograms();

};

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

DetvarHistogramManager::DetvarHistogramManager(std::string label) : _label(label)
{
}

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void DetvarHistogramManager::LoadTemplate()
{
  TFile* f_tp = TFile::Open(("Analysis/"+_label+"/rootfiles/BinningTemplate.root").c_str());
  _h_tp = (TH1D*)f_tp->Get("h_template");
  _h_tp->SetDirectory(0);
  f_tp->Close();

  SetupHistograms();
}

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void DetvarHistogramManager::SetTemplate(std::string axis_title,int nbins,double low,double high)
{
  _h_tp = new TH1D("h_template",axis_title.c_str(),nbins,low,high);

  SetupHistograms();
}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void DetvarHistogramManager::SetupHistograms()
{

  _h_CV_Tot = (TH1D*)_h_tp->Clone(("h_Detvar_CV_Tot"+_label).c_str());
  
  for(int i_s=0;i_s<kDetvarMAX;i_s++)
    _h_Vars_Tot.push_back((TH1D*)_h_tp->Clone(("h_Detvar_Vars_Tot_"+detvar_str.at(i_s)+"_"+_label).c_str()));
 
  for(size_t i_c=0;i_c<categories.size();i_c++){
    _h_CV.push_back((TH1D*)_h_tp->Clone(("h_Detvar_CV_"+categories.at(i_c)+"_"+_label).c_str()));
    _h_Vars.push_back(std::vector<TH1D*>());
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      _h_Vars.back().push_back((TH1D*)_h_tp->Clone(("h_Detvar_Vars_"+categories.at(i_c)+"_"+detvar_str.at(i_s)+"_"+_label).c_str()));
    }
  }

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void DetvarHistogramManager::FillHistograms(double var,double weight)
{

  if(_h_tp == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(detvar_univ == -1){
    if(!is_data) _h_CV_Tot->Fill(var,POT_weight);
    _h_CV.at(category)->Fill(var,POT_weight);
  }
  else {
    _h_Vars_Tot.at(detvar_univ)->Fill(var,POT_weight);
    _h_Vars.at(category).at(detvar_univ)->Fill(var,POT_weight);
  }

}

///////////////////////////////////////////////////////////////////////
// Write the histograms 

void DetvarHistogramManager::Write()
{

  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    _h_Vars_Tot.at(i_s)->Add(_h_CV.at(kEXT));
    _h_Vars_Tot.at(i_s)->Add(_h_CV.at(kDirt));
    _h_Vars.at(kEXT).at(i_s)->Add(_h_CV.at(kEXT));
    _h_Vars.at(kDirt).at(i_s)->Add(_h_CV.at(kDirt));
  }

  gSystem->Exec(("mkdir -p Analysis/"+_label+"/rootfiles/").c_str());
  TFile* f_out = TFile::Open(("Analysis/"+_label+"/rootfiles/Detvars.root").c_str(),"RECREATE");

  if(_divide_by_bin_width) DivideByBinWidth(_h_CV_Tot);
  double integral = _h_CV_Tot->Integral("width");
  if(_shape_only) _h_CV_Tot->Scale(1.0/integral);
  _h_CV_Tot->Write("h_Detvar_CV_Tot");

  for(size_t i_c=0;i_c<categories.size();i_c++){
    if(_divide_by_bin_width) DivideByBinWidth(_h_CV.at(i_c));
    if(_shape_only) _h_CV.at(i_c)->Scale(1.0/integral);
    _h_CV.at(i_c)->Write(("h_Detvar_CV_"+categories.at(i_c)).c_str());
  }

  std::vector<double> integral_v;
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    if(_divide_by_bin_width) DivideByBinWidth(_h_Vars_Tot.at(i_s));
    integral_v.push_back(_h_Vars_Tot.at(i_s)->Integral("width"));
    if(_shape_only) _h_Vars_Tot.at(i_s)->Scale(1.0/integral_v.back());
    _h_Vars_Tot.at(i_s)->Write(("h_Detvar_Vars_Tot_"+detvar_str.at(i_s)).c_str());
  }  

  for(size_t i_c=0;i_c<categories.size();i_c++){
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      if(_divide_by_bin_width) DivideByBinWidth(_h_Vars.at(i_c).at(i_s));
      if(_shape_only) _h_Vars.at(i_c).at(i_s)->Scale(1.0/integral_v.at(i_s));
      _h_Vars.at(i_c).at(i_s)->Write(("h_Detvar_Vars_"+categories.at(i_c)+"_"+detvar_str.at(i_s)).c_str());  
    }
  }      

  std::vector<TH2D*> h_Cov_Tot; 
  std::vector<TH2D*> h_FCov_Tot; 
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    h_Cov_Tot.push_back(nullptr);
    h_FCov_Tot.push_back(nullptr);
    CalcCovUnisim(detvar_str.at(i_s),_h_CV_Tot,_h_Vars_Tot.at(i_s),h_Cov_Tot.back(),h_FCov_Tot.back());
    h_Cov_Tot.back()->Write();
    h_FCov_Tot.back()->Write();
  }

  TH2D* h_Cov = static_cast<TH2D*>(h_Cov_Tot.at(0)->Clone("Cov"));
  TH2D* h_FCov = static_cast<TH2D*>(h_FCov_Tot.at(0)->Clone("FCov"));
  h_Cov->Reset();
  h_FCov->Reset();
  for(int i_s=0;i_s<kSystMAX;i_s++){
    h_Cov->Add(h_Cov_Tot.at(i_s));    
    h_FCov->Add(h_FCov_Tot.at(i_s));    
  }
  h_Cov->Write();
  h_FCov->Write();

  std::vector<std::vector<TH2D*>> h_Cov_Cat; 
  std::vector<std::vector<TH2D*>> h_FCov_Cat; 

  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_Cov_Cat.push_back(std::vector<TH2D*>());
    h_FCov_Cat.push_back(std::vector<TH2D*>());
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      h_Cov_Cat.back().push_back(nullptr);
      h_FCov_Cat.back().push_back(nullptr);
      CalcCovUnisim(categories.at(i_c)+"_"+detvar_str.at(i_s),_h_CV.at(i_c),_h_Vars.at(i_c).at(i_s),h_Cov_Cat.back().back(),h_FCov_Cat.back().back());
      h_Cov_Cat.back().back()->Write();
      h_FCov_Cat.back().back()->Write();
    }
  }       

  f_out->Close();

}

}

#endif
