#ifndef _MultiChannelHistograms_h_
#define _MultiChannelHistograms_h_

#include "Histograms2.h"

// Wrapper class that acts as an interface between functions trying to 
// performa analysis with multiple channels and the histogram manager 

namespace hist {

class MultiChannelHistogramManager {

  public:

    MultiChannelHistogramManager(std::string label,bool save_truth=false);
    void SetRecoChannelList(std::vector<std::string> ch_v={});
    void SetTrueChannelList(std::vector<std::string> ch_v={});
    void MakeHM();    

    void FillTruthHistograms(bool sig,double var_t,bool load_syst,std::string ch="",double weight=1.0);
    void FillRecoHistograms(bool sel,double var_r,bool load_syst,std::string ch="",double weight=1.0);
    void FillHistograms2D(bool sig,bool sel,double var_t,double var_r,bool load_syst,std::string ch_t="",std::string ch_r="",double weight=1.0);
    void KeepAll(){ _keep_all = true; _hm.KeepAll(); }
    void Write(){ _hm.Write(); }

    void AddSpecialUniv(std::string name){ _hm.AddSpecialUniv(name); }
    void FillSpecialTruthHistograms(std::string name,bool sig,double var_t,double weight,std::string ch="");
    void FillSpecialRecoHistograms(std::string name,bool sel,double var_r,double weight,std::string ch="");
    void FillSpecialHistograms2D(std::string name,bool sig,bool sel,double var_t,double var_r,double weight,std::string ch_t="",std::string ch_r="");

  private:

    const std::string _label;
    const bool _save_truth;
    bool _keep_all;
   
    std::vector<std::string> _ch_list_r;
    std::vector<std::string> _ch_list_t;
    int _n_ch_r;
    int _n_ch_t;
    std::vector<TH1D*> _h_tp_v;
    std::vector<TH1D*> _h_tp_truth_v;
    std::vector<int> _offset_r;
    std::vector<int> _offset_t;
    
    HistogramManager _hm;

};

///////////////////////////////////////////////////////////////////////
// initial setup 

MultiChannelHistogramManager::MultiChannelHistogramManager(std::string label,bool save_truth) : _label(label), _save_truth(save_truth), _hm(label,save_truth)
{
  std::cout << "Setting up MultiChannelHistogramManager for " << _label << std::endl;
}

///////////////////////////////////////////////////////////////////////
// Set reco channels 

void MultiChannelHistogramManager::SetRecoChannelList(std::vector<std::string> ch_v)
{
  _ch_list_r = ch_v;
  _n_ch_r = _ch_list_r.size();
}

///////////////////////////////////////////////////////////////////////
// Set true channels 

void MultiChannelHistogramManager::SetTrueChannelList(std::vector<std::string> ch_v)
{
  _ch_list_t = ch_v;
  _n_ch_t = _ch_list_t.size();
}

///////////////////////////////////////////////////////////////////////
// Create the binning scheme 

void MultiChannelHistogramManager::MakeHM()
{

  TFile* f_tp = TFile::Open(("Analysis/"+_label+"/rootfiles/BinningTemplate.root").c_str());

  if(!_ch_list_r.size()) _h_tp_v.push_back((TH1D*)f_tp->Get("h_template"));
  else for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++) _h_tp_v.push_back((TH1D*)f_tp->Get(("h_template_"+_ch_list_r.at(i_ch)).c_str()));

  for(TH1D* h : _h_tp_v) h->SetDirectory(0);

  f_tp->Close();

  for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++)
    _h_tp_v.at(i_ch)->SetName(("h_template_"+_ch_list_r.at(i_ch)+"_"+_label).c_str()); 

  if(_save_truth){

    TFile* f_tp_truth = TFile::Open(("Analysis/"+_label+"/rootfiles/TruthBinningTemplate.root").c_str());

    if(!_ch_list_t.size()) _h_tp_truth_v.push_back((TH1D*)f_tp_truth->Get("h_template"));
    else for(size_t i_ch=0;i_ch<_ch_list_t.size();i_ch++) _h_tp_truth_v.push_back((TH1D*)f_tp_truth->Get(("h_template_"+_ch_list_t.at(i_ch)).c_str()));

    for(TH1D* h : _h_tp_truth_v) h->SetDirectory(0);

    f_tp_truth->Close();

    for(size_t i_ch=0;i_ch<_ch_list_t.size();i_ch++)
      _h_tp_truth_v.at(i_ch)->SetName(("h_template_"+_ch_list_t.at(i_ch)+"_"+_label).c_str()); 

  }

  int tot_bins_r = _h_tp_v.at(0)->GetNbinsX()+2;
  _offset_r.push_back(0);
  for(size_t i_ch=1;i_ch<_ch_list_r.size();i_ch++){
    _offset_r.push_back(tot_bins_r); 
    tot_bins_r += _h_tp_v.at(i_ch)->GetNbinsX()+2; 
  }
 
  int tot_bins_t = 0;
  if(_save_truth){
    tot_bins_t = _h_tp_truth_v.at(0)->GetNbinsX()+2;
    _offset_t.push_back(0);
    for(size_t i_ch=1;i_ch<_ch_list_t.size();i_ch++){
      _offset_t.push_back(tot_bins_t); 
      tot_bins_t += _h_tp_truth_v.at(i_ch)->GetNbinsX()+2; 
    }
  }

  if(_save_truth) _hm.SetTemplate("",tot_bins_t,-0.5,tot_bins_t-0.5,tot_bins_r,-0.5,tot_bins_r-0.5);
  else _hm.SetTemplate("",tot_bins_r,-0.5,tot_bins_r-0.5);

}    

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void MultiChannelHistogramManager::FillRecoHistograms(bool sel,double var_r,bool load_syst,std::string ch,double weight)
{
  int ch_idx = 0;
  for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++)
    if(ch == _ch_list_r.at(i_ch)) ch_idx = i_ch;  
  int bin_r = _offset_r.at(ch_idx)+_h_tp_v.at(ch_idx)->FindBin(var_r);

  _hm.FillRecoHistograms(sel,bin_r,load_syst,weight);
}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void MultiChannelHistogramManager::FillTruthHistograms(bool sig,double var_t,bool load_syst,std::string ch,double weight)
{
  int ch_idx = 0;
  for(size_t i_ch=0;i_ch<_ch_list_t.size();i_ch++)
    if(ch == _ch_list_t.at(i_ch)) ch_idx = i_ch;  
  int bin_t = _offset_t.at(ch_idx)+_h_tp_truth_v.at(ch_idx)->FindBin(var_t);

  _hm.FillTruthHistograms(sig,bin_t,load_syst,weight);

}

///////////////////////////////////////////////////////////////////////
// Fill the histograms

void MultiChannelHistogramManager::FillHistograms2D(bool sig,bool sel,double var_t,double var_r,bool load_syst,std::string ch_t,std::string ch_r,double weight)
{
  int ch_idx_t = 0;
  for(size_t i_ch=0;i_ch<_ch_list_t.size();i_ch++)
    if(ch_t == _ch_list_t.at(i_ch)) ch_idx_t = i_ch;  
  int bin_t = _offset_t.at(ch_idx_t)+_h_tp_truth_v.at(ch_idx_t)->FindBin(var_t);

  int ch_idx_r = 0;
  for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++)
    if(ch_r == _ch_list_r.at(i_ch)) ch_idx_r = i_ch;  
  int bin_r = _offset_r.at(ch_idx_r)+_h_tp_v.at(ch_idx_r)->FindBin(var_r);
   
  _hm.FillHistograms2D(sig,sel,bin_t,bin_r,load_syst,weight);

}

///////////////////////////////////////////////////////////////////////
// Fill the special truth histograms

void MultiChannelHistogramManager::FillSpecialTruthHistograms(std::string name,bool sig,double var_t,double weight,std::string ch)
{
  int ch_idx = 0;
  for(size_t i_ch=0;i_ch<_ch_list_t.size();i_ch++)
    if(ch == _ch_list_t.at(i_ch)) ch_idx = i_ch;  
  int bin_t = _offset_t.at(ch_idx)+_h_tp_truth_v.at(ch_idx)->FindBin(var_t);

  _hm.FillSpecialTruthHistograms(name,sig,var_t,weight);
}

///////////////////////////////////////////////////////////////////////
// Fill the special reco histograms

void MultiChannelHistogramManager::FillSpecialRecoHistograms(std::string name,bool sel,double var_r,double weight,std::string ch)
{
  int ch_idx = 0;
  for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++)
    if(ch == _ch_list_r.at(i_ch)) ch_idx = i_ch;  
  int bin_r = _offset_r.at(ch_idx)+_h_tp_v.at(ch_idx)->FindBin(var_r);

  _hm.FillSpecialRecoHistograms(name,sel,var_r,weight);
}

///////////////////////////////////////////////////////////////////////
// Fill the special histograms

void MultiChannelHistogramManager::FillSpecialHistograms2D(std::string name,bool sig,bool sel,double var_t,double var_r,double weight,std::string ch_t,std::string ch_r)
{
  int ch_idx_t = 0;
  for(size_t i_ch=0;i_ch<_ch_list_t.size();i_ch++)
    if(ch_t == _ch_list_t.at(i_ch)) ch_idx_t = i_ch;  
  int bin_t = _offset_t.at(ch_idx_t)+_h_tp_truth_v.at(ch_idx_t)->FindBin(var_t);

  int ch_idx_r = 0;
  for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++)
    if(ch_r == _ch_list_r.at(i_ch)) ch_idx_r = i_ch;  
  int bin_r = _offset_r.at(ch_idx_r)+_h_tp_v.at(ch_idx_r)->FindBin(var_r);

  _hm.FillSpecialHistograms2D(name,sig,sel,var_t,var_r,weight);
}

///////////////////////////////////////////////////////////////////////

}

#endif
