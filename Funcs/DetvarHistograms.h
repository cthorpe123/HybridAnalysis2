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

    DetvarHistogramManager(std::string label,bool save_truth=false);

    void LoadTemplate(std::string t = "");
    void SetTemplate(std::string axis_title,int nbins,double low,double high);
    void SetTemplate(std::string axis_title,int nbins_t,double low_t,double high_t,int nbins_r,double low_r,double high_r);

    void FillTruthHistograms(bool sig,double var_t,bool load_syst,double weight=1.0);
    void FillRecoHistograms(bool sel,double var_r,bool load_syst,double weight=1.0);
    void FillHistograms2D(bool sig,bool sel,double var_t,double var_r,bool load_syst,double weight=1.0);

    void KeepAll(){ _keep_all = true; }
    void Write();

  private:

    const std::string _label;
    const bool _save_truth;
    bool _keep_all = false;

    // histogram templates
    TH1D* _h_tp = nullptr;
    TH1D* _h_tp_truth = nullptr;
    TH2D* _h_tp_joint = nullptr;

    // Reco histograms 
    TH1D* _h_CV_Reco_Tot;
    std::vector<TH1D*> _h_Vars_Reco_Tot;

    // Reco histograms by category
    std::vector<TH1D*> _h_CV_Reco_Cat;
    std::vector<std::vector<TH1D*>> _h_Vars_Reco_Cat;
 
    // Truth histograms, only filled for the signal
    TH1D* _h_CV_Truth_Signal;
    std::vector<TH1D*> _h_Vars_Truth_Signal;

    // Joint histograms, only filled for the signal
    TH2D* _h_CV_Joint_Signal;
    std::vector<TH2D*> _h_Vars_Joint_Signal;

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

DetvarHistogramManager::DetvarHistogramManager(std::string label,bool save_truth) : _label(label), _save_truth(save_truth)
{
  std::cout << "Setting up DetvarHistogramManager for " << _label << std::endl;
}

///////////////////////////////////////////////////////////////////////
// initial setup - load template histogram from file

void DetvarHistogramManager::LoadTemplate(std::string t)
{
  if(t == "") t = _label;

  TFile* f_tp = TFile::Open(("Analysis/"+t+"/rootfiles/BinningTemplate.root").c_str());
  _h_tp = (TH1D*)f_tp->Get("h_template");
  _h_tp->SetDirectory(0);
  f_tp->Close();
  _h_tp->SetName(("h_template_"+_label).c_str()); 

  if(_save_truth){ 
    TFile* f_tp_truth = TFile::Open(("Analysis/"+t+"/rootfiles/TruthBinningTemplate.root").c_str());
    _h_tp_truth = (TH1D*)f_tp_truth->Get("h_template");
    _h_tp_truth->SetDirectory(0);
    f_tp_truth->Close();
    _h_tp_truth->SetName(("h_template_truth_"+_label).c_str()); 
  }

  _SetupHistograms();
}

///////////////////////////////////////////////////////////////////////
// initial setup - Set binning scheme manually

void DetvarHistogramManager::SetTemplate(std::string axis_title,int nbins_r,double low_r,double high_r)
{
  if(_save_truth) 
    throw std::invalid_argument("HistogramManager::SetTemplate: Calling wrong function when response matrices are being set, need to specify truth binning as well");   
  _h_tp = new TH1D(("h_template_"+_label).c_str(),axis_title.c_str(),nbins_r,low_r,high_r);

  _SetupHistograms();

}


///////////////////////////////////////////////////////////////////////
// initial setup - Set binning scheme manually 

void DetvarHistogramManager::SetTemplate(std::string axis_title,int nbins_t,double low_t,double high_t,int nbins_r,double low_r,double high_r)
{
  if(!_save_truth) 
    throw std::invalid_argument("HistogramManager::SetTemplate: Calling wrong function when response matrices are being set, truth binning not enabled");  

  _h_tp = new TH1D(("h_template_"+_label).c_str(),axis_title.c_str(),nbins_r,low_r,high_r);
  _h_tp_truth = new TH1D(("h_template_truth_"+_label).c_str(),axis_title.c_str(),nbins_t,low_t,high_t);

  _SetupHistograms();

}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void DetvarHistogramManager::_SetupHistograms()
{
  _SetupRecoHistograms();
  if(_save_truth){
   _SetupTruthHistograms();
   _SetupJointHistograms();
  }
}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void DetvarHistogramManager::_SetupRecoHistograms()
{

  _bin_edges_r.clear();
  for(int i=1;i<_h_tp->GetNbinsX()+2;i++) _bin_edges_r.push_back(_h_tp->GetXaxis()->GetBinLowEdge(i));
  _nbins_r = _bin_edges_r.size()-1;

  _h_CV_Reco_Tot = (TH1D*)_h_tp->Clone(("h_CV_Reco_Tot_"+_label).c_str());
  
  for(int i_s=0;i_s<kDetvarMAX;i_s++)
    _h_Vars_Reco_Tot.push_back((TH1D*)_h_tp->Clone(("h_Vars_Reco_Tot_"+detvar_str.at(i_s)+"_"+_label).c_str()));
 
  for(size_t i_c=0;i_c<categories.size();i_c++){
    _h_CV_Reco_Cat.push_back((TH1D*)_h_tp->Clone(("h_CV_Reco_"+categories.at(i_c)+"_"+_label).c_str()));
    _h_Vars_Reco_Cat.push_back(std::vector<TH1D*>());
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      _h_Vars_Reco_Cat.back().push_back((TH1D*)_h_tp->Clone(("h_Vars_Reco_"+categories.at(i_c)+"_"+detvar_str.at(i_s)+"_"+_label).c_str()));
    }
  }

}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void DetvarHistogramManager::_SetupTruthHistograms()
{

  _bin_edges_t.clear();
  for(int i=1;i<_h_tp_truth->GetNbinsX()+2;i++) _bin_edges_t.push_back(_h_tp_truth->GetXaxis()->GetBinLowEdge(i));
  _nbins_t = _bin_edges_t.size()-1;

  _h_CV_Truth_Signal = (TH1D*)_h_tp_truth->Clone(("h_CV_Truth_Signal_"+_label).c_str());
  
  for(int i_s=0;i_s<kDetvarMAX;i_s++)
    _h_Vars_Truth_Signal.push_back((TH1D*)_h_tp_truth->Clone(("h_Vars_Truth_Signal_"+detvar_str.at(i_s)+"_"+_label).c_str()));
 
}

///////////////////////////////////////////////////////////////////////
// Create all of the vectors of histogram pointers

void DetvarHistogramManager::_SetupJointHistograms()
{
  // Storing totals
  std::string t_title = _h_tp_truth->GetXaxis()->GetTitle();
  std::string r_title = _h_tp->GetXaxis()->GetTitle();
  std::string title = ";"+t_title+";"+r_title+";";
  _h_CV_Joint_Signal = new TH2D(("h_CV_Joint_Signal_"+_label).c_str(),title.c_str(),_nbins_t,&_bin_edges_t[0],_nbins_r,&_bin_edges_r[0]);

  _h_tp_joint = new TH2D(("h_template_joint_"+_label).c_str(),title.c_str(),_nbins_t,&_bin_edges_t[0],_nbins_r,&_bin_edges_r[0]);

  for(int i_s=0;i_s<kDetvarMAX;i_s++)
    _h_Vars_Joint_Signal.push_back((TH2D*)_h_CV_Joint_Signal->Clone(("h_Vars_Joint_Signal_"+unisims_str.at(i_s)+"_"+_label).c_str()));
 
}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void DetvarHistogramManager::FillTruthHistograms(bool sig,double var_t,bool load_syst,double weight=1.0)
{
  if(!sig) return;

  if(_h_tp_truth == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(category == -1) std::cout << "Bad event" << std::endl;

  if(std::isnan(weightSplineTimesTune) || std::isinf(weightSplineTimesTune)) return;

  if(detvar_univ == -1) _h_CV_Truth_Signal->Fill(var_t,POT_weight*weightSplineTimesTune);
  else _h_Vars_Truth_Signal.at(detvar_univ)->Fill(var_t,POT_weight*weightSplineTimesTune);

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void DetvarHistogramManager::FillRecoHistograms(bool sel,double var_r,bool load_syst,double weight=1.0)
{
  if(!sel) return;
  
  if(_h_tp == nullptr) 
    throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

  if(category == -1) std::cout << "Bad event" << std::endl;

  if(std::isnan(weightSplineTimesTune) || std::isinf(weightSplineTimesTune)) return;

  if(detvar_univ == -1){
    if(!is_data) _h_CV_Reco_Tot->Fill(var_r,POT_weight*weightSplineTimesTune);
    _h_CV_Reco_Cat.at(category)->Fill(var_r,POT_weight*weightSplineTimesTune);
  }
  else {
    _h_Vars_Reco_Tot.at(detvar_univ)->Fill(var_r,POT_weight*weightSplineTimesTune);
    _h_Vars_Reco_Cat.at(category).at(detvar_univ)->Fill(var_r,POT_weight*weightSplineTimesTune);
  }

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void DetvarHistogramManager::FillHistograms2D(bool sig,bool sel,double var_t,double var_r,bool load_syst,double weight=1.0)
{

 FillTruthHistograms(sig,var_t,weight);
 FillRecoHistograms(sel,var_r,weight);

 if(!sig || !sel) return;

 if(_h_tp == nullptr || _h_tp_truth == nullptr) 
   throw std::invalid_argument("Trying to fill histograms before binning has been set, call LoadTemplate or SetTemplate first!");

 if(category == -1) std::cout << "Bad event" << std::endl;

 if(std::isnan(weightSplineTimesTune) || std::isinf(weightSplineTimesTune)) return;

  if(detvar_univ == -1) _h_CV_Joint_Signal->Fill(var_t,var_r,POT_weight*weightSplineTimesTune);
  else _h_Vars_Joint_Signal.at(detvar_univ)->Fill(var_t,var_r,POT_weight*weightSplineTimesTune);

}

///////////////////////////////////////////////////////////////////////
// Write the histograms to file

void DetvarHistogramManager::Write()
{
  std::cout << "Writing histograms for " << _label << std::endl;
  gSystem->Exec(("mkdir -p Analysis/"+_label+"/rootfiles/").c_str());
  _f_out = TFile::Open(("Analysis/"+_label+"/rootfiles/Detvars.root").c_str(),"RECREATE");

  _WriteReco();
  if(_save_truth){
    _WriteTruth();
    _WriteJoint();
  }

  _f_out->Close();
}

///////////////////////////////////////////////////////////////////////
// Write the reco histograms to file

void DetvarHistogramManager::_WriteReco()
{
  _f_out->mkdir("Reco");

  // Add the EXT and dirt backgrounds to the variations
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    _h_Vars_Reco_Tot.at(i_s)->Add(_h_CV_Reco_Cat.at(kEXT));
    _h_Vars_Reco_Tot.at(i_s)->Add(_h_CV_Reco_Cat.at(kDirt));
    _h_Vars_Reco_Cat.at(kEXT).at(i_s)->Add(_h_CV_Reco_Cat.at(kEXT));
    _h_Vars_Reco_Cat.at(kDirt).at(i_s)->Add(_h_CV_Reco_Cat.at(kDirt));
  }

  // Store the total covariance as well 
  TH2D* h_Cov = Make2DHist("h_Cov",_h_tp); 
  TH2D* h_FCov = Make2DHist("h_FCov",_h_tp); 
  std::vector<TH2D*> h_Cov_Cat,h_FCov_Cat;
  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_Cov_Cat.push_back(Make2DHist("h_Cov_"+categories.at(i_c),_h_tp));
    h_FCov_Cat.push_back(Make2DHist("h_FCov_"+categories.at(i_c),_h_tp));
  }

  _f_out->mkdir("Reco/CV");
  _f_out->cd("Reco/CV");

  _h_CV_Reco_Tot->Write("h_Tot");
  for(size_t i_c=0;i_c<categories.size();i_c++)
    _h_CV_Reco_Cat.at(i_c)->Write(("h_"+categories.at(i_c)).c_str());
   
  _f_out->cd();

  if(_keep_all){
    _f_out->mkdir("Reco/Vars");
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      _f_out->mkdir(("Reco/Vars/"+detvar_str.at(i_s)).c_str());
      _f_out->cd(("Reco/Vars/"+detvar_str.at(i_s)).c_str());
      _h_Vars_Reco_Tot.at(i_s)->Write("h_Tot");
      for(size_t i_c=0;i_c<categories.size();i_c++)
        _h_Vars_Reco_Cat.at(i_c).at(i_s)->Write(("h_"+categories.at(i_c)).c_str());
      _f_out->cd();
    }
  }

  _f_out->mkdir("Reco/Cov");
  _f_out->cd("Reco/Cov");
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    _f_out->mkdir(("Reco/Cov/"+detvar_str.at(i_s)).c_str());
    _f_out->cd(("Reco/Cov/"+detvar_str.at(i_s)).c_str());
    TH2D *C,*FC;
    CalcCovUnisim(detvar_str.at(i_s),_h_CV_Reco_Tot,_h_Vars_Reco_Tot.at(i_s),C,FC);
    C->Write("Cov_Tot");
    FC->Write("FCov_Tot");
    h_Cov->Add(C);
    h_FCov->Add(FC);
    for(size_t i_c=0;i_c<categories.size();i_c++){
      TH2D *C_Cat,*FC_Cat;
      CalcCovUnisim(detvar_str.at(i_s)+"_"+categories.at(i_c),_h_CV_Reco_Cat.at(i_c),_h_Vars_Reco_Cat.at(i_c).at(i_s),C_Cat,FC_Cat);
      C_Cat->Write(("Cov_"+categories.at(i_c)).c_str());
      FC_Cat->Write(("FCov_"+categories.at(i_c)).c_str());
      h_Cov_Cat.at(i_c)->Add(C_Cat);
      h_FCov_Cat.at(i_c)->Add(FC_Cat);
    }
    _f_out->cd();
  }

  _f_out->mkdir("Reco/Cov/Total");
  _f_out->cd("Reco/Cov/Total");
  h_Cov->Write("Cov_Tot"); 
  h_FCov->Write("FCov_Tot"); 
  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_Cov_Cat.at(i_c)->Write(("Cov_"+categories.at(i_c)).c_str());
    h_FCov_Cat.at(i_c)->Write(("FCov_"+categories.at(i_c)).c_str());
  }

  _f_out->cd();
 
}

///////////////////////////////////////////////////////////////////////
// Write the truth histograms to file

void DetvarHistogramManager::_WriteTruth()
{
  _f_out->mkdir("Truth");

  // Histograms to store totals
  TH2D* h_Cov = Make2DHist("h_Cov_Truth",_h_tp_truth); 
  TH2D* h_FCov = Make2DHist("h_FCov_Truth",_h_tp_truth); 

  _f_out->mkdir("Truth/CV");
  _f_out->cd("Truth/CV");

  _h_CV_Truth_Signal->Write("h_Signal");
  _f_out->cd();

  if(_keep_all){
    _f_out->mkdir("Truth/Vars");
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      _f_out->mkdir(("Truth/Vars/"+detvar_str.at(i_s)).c_str());
      _f_out->cd(("Truth/Vars/"+detvar_str.at(i_s)).c_str());
      _h_Vars_Truth_Signal.at(i_s)->Write("h_Signal");
      _f_out->cd();
    }
  }

  _f_out->mkdir("Truth/Cov");
  _f_out->cd("Truth/Cov");
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    _f_out->mkdir(("Truth/Cov/"+detvar_str.at(i_s)).c_str());
    _f_out->cd(("Truth/Cov/"+detvar_str.at(i_s)).c_str());
    TH2D *C,*FC;
    CalcCovUnisim(detvar_str.at(i_s),_h_CV_Truth_Signal,_h_Vars_Truth_Signal.at(i_s),C,FC);
    C->Write("Cov_Signal");
    FC->Write("FCov_Signal");
    h_Cov->Add(C);
    h_FCov->Add(FC);
    _f_out->cd();
  }

  _f_out->mkdir("Truth/Cov/Total");
  _f_out->cd("Truth/Cov/Total");
  h_Cov->Write("Cov_Tot"); 
  h_FCov->Write("FCov_Tot"); 
  _f_out->cd();

}

///////////////////////////////////////////////////////////////////////
// Write the truth histograms to file

void DetvarHistogramManager::_WriteJoint()
{
  _f_out->mkdir("Response");

  _f_out->mkdir("Response/CV");
  _f_out->cd("Response/CV");
  TH2D* h_CV_Response_Signal = (TH2D*)_h_CV_Joint_Signal->Clone("h_CV_Response_Signal");  
  NormaliseResponse(_h_CV_Truth_Signal,h_CV_Response_Signal); 
  h_CV_Response_Signal->Write("h_Signal");

  if(_keep_all){
    _f_out->mkdir("Response/Vars");
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      _f_out->mkdir(("Response/Vars/"+detvar_str.at(i_s)).c_str());
      _f_out->cd(("Response/Vars/"+detvar_str.at(i_s)).c_str());
      TH2D* h = (TH2D*)_h_Vars_Joint_Signal.at(i_s)->Clone("h_CV_Response_Signal");  
      NormaliseResponse(_h_Vars_Truth_Signal.at(i_s),h);
      h->Write("h_Signal");
      _f_out->cd();
    }
  }

  _f_out->mkdir("Joint/CV");
  _f_out->cd("Joint/CV");
  _h_CV_Joint_Signal->Write("h_Signal");

  _f_out->cd();

  if(_keep_all){
    _f_out->mkdir("Joint/Vars");
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      _f_out->mkdir(("Joint/Vars/"+detvar_str.at(i_s)).c_str());
      _f_out->cd(("Joint/Vars/"+detvar_str.at(i_s)).c_str());
      _h_Vars_Joint_Signal.at(i_s)->Write("h_Signal");
      _f_out->cd();
    }
  }

}

///////////////////////////////////////////////////////////////////////

}

#endif
