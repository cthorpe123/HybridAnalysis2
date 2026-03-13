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
    void LoadTemplates();
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

    TH1D* RestoreRecoBinning(const TH1D* h,std::string name) const;
    TH2D* RestoreRecoBinning(const TH2D* h,std::string name) const;

  private:

    const std::string _label;
    const bool _save_truth;
    bool _keep_all;
    bool _hm_loaded = false;  
 
    std::vector<std::string> _ch_list_r;
    std::vector<std::string> _ch_list_t;
    int _n_ch_r;
    int _n_ch_t;
    std::vector<TH1D*> _h_tp_v;
    std::vector<TH1D*> _h_tp_truth_v;
    std::vector<int> _offset_r;
    std::vector<int> _offset_t;
    int _tot_bins_r = 0;
    int _tot_bins_t = 0; 

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

void MultiChannelHistogramManager::LoadTemplates()
{

  TFile* f_tp = TFile::Open(("Analysis/"+_label+"/rootfiles/BinningTemplate.root").c_str());

  if(!_ch_list_r.size()){
    _h_tp_v.push_back((TH1D*)f_tp->Get("h_template"));
    _ch_list_r.push_back("All");
  }
  else for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++) _h_tp_v.push_back((TH1D*)f_tp->Get(("h_template_"+_ch_list_r.at(i_ch)).c_str()));

  for(TH1D* h : _h_tp_v) h->SetDirectory(0);

  f_tp->Close();

  for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++)
    _h_tp_v.at(i_ch)->SetName(("h_template_"+_ch_list_r.at(i_ch)+"_"+_label).c_str()); 

  if(_save_truth){

    TFile* f_tp_truth = TFile::Open(("Analysis/"+_label+"/rootfiles/TruthBinningTemplate.root").c_str());

    if(!_ch_list_t.size()){
      _h_tp_truth_v.push_back((TH1D*)f_tp_truth->Get("h_template"));
      _ch_list_t.push_back("All");
    }
    else for(size_t i_ch=0;i_ch<_ch_list_t.size();i_ch++) _h_tp_truth_v.push_back((TH1D*)f_tp_truth->Get(("h_template_"+_ch_list_t.at(i_ch)).c_str()));

    for(TH1D* h : _h_tp_truth_v) h->SetDirectory(0);

    f_tp_truth->Close();

    for(size_t i_ch=0;i_ch<_ch_list_t.size();i_ch++)
      _h_tp_truth_v.at(i_ch)->SetName(("h_template_"+_ch_list_t.at(i_ch)+"_"+_label).c_str()); 

  }

  _tot_bins_r = 0;
  for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++){
    _offset_r.push_back(_tot_bins_r); 
    _tot_bins_r += _h_tp_v.at(i_ch)->GetNbinsX()+2; 
  }
 
  if(_save_truth){
    _tot_bins_t = 0;
    for(size_t i_ch=0;i_ch<_ch_list_t.size();i_ch++){
      _offset_t.push_back(_tot_bins_t); 
      _tot_bins_t += _h_tp_truth_v.at(i_ch)->GetNbinsX()+2; 
    }
  }

}    

///////////////////////////////////////////////////////////////////////
// Create the histogram manager object 

void MultiChannelHistogramManager::MakeHM()
{
  // Arrange the bins such that 0 corresponds to underflow of first channel, 1 the first bin of the first channel
  // and the last bin is the overflow of the last channel, and the second to last bin the hhighest bing of the last channel 
  if(_save_truth) _hm.SetTemplate("",_tot_bins_t-2,0.5,_tot_bins_t-1.5,_tot_bins_r-2,0.5,_tot_bins_r-1.5);
  else _hm.SetTemplate("",_tot_bins_r-2,0.5,_tot_bins_r-1.5);

  bool _hm_loaded = true;  
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

  _hm.FillSpecialTruthHistograms(name,sig,bin_t,weight);
}

///////////////////////////////////////////////////////////////////////
// Fill the special reco histograms

void MultiChannelHistogramManager::FillSpecialRecoHistograms(std::string name,bool sel,double var_r,double weight,std::string ch)
{
  int ch_idx = 0;
  for(size_t i_ch=0;i_ch<_ch_list_r.size();i_ch++)
    if(ch == _ch_list_r.at(i_ch)) ch_idx = i_ch;  
  int bin_r = _offset_r.at(ch_idx)+_h_tp_v.at(ch_idx)->FindBin(var_r);

  _hm.FillSpecialRecoHistograms(name,sel,bin_r,weight);
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

  _hm.FillSpecialHistograms2D(name,sig,sel,bin_t,bin_r,weight);
}

///////////////////////////////////////////////////////////////////////
// Take a histogram with the numeric binning scheme in reco space and 
// convert it back to the physics variable binning scheme - currently
// only implemented for single reco channel

TH1D* MultiChannelHistogramManager::RestoreRecoBinning(const TH1D* h,std::string name) const
{
  if(_ch_list_r.size() > 1)
    throw std::invalid_argument("MultiChannelHistogramManager::RestoreRecoBinning only implemented for one reco channel at the moment");


  TH1D* h_out = (TH1D*)_h_tp_v.at(0)->Clone(name.c_str());  
  std::cout << "1D: h_out->GetNbinsX() = " << h_out->GetNbinsX() << std::endl;

  // bin 1 of numeric binning scheme is underflow physical binning scheme 
  for(int i=0;i<h_out->GetNbinsX()+2;i++){
    h_out->SetBinContent(i,h->GetBinContent(i));
    h_out->SetBinError(i,h->GetBinError(i));
  }
  
  return h_out; 

}

///////////////////////////////////////////////////////////////////////

TH2D* MultiChannelHistogramManager::RestoreRecoBinning(const TH2D* h,std::string name) const
{
  if(_ch_list_r.size() > 1)
    throw std::invalid_argument("MultiChannelHistogramManager::RestoreRecoBinning only implemented for one reco channel at the moment");

  TH2D* h_out = Make2DHist(name,_h_tp_v.at(0));
  std::cout << "2D: h_out->GetNbinsX() = " << h_out->GetNbinsX() << std::endl;

  for(int i=0;i<h->GetNbinsX()+2;i++){
    for(int j=0;j<h->GetNbinsY()+2;j++){
      h_out->SetBinContent(i,j,h->GetBinContent(i,j));
      h_out->SetBinError(i,j,h->GetBinError(i,j));
    }
    std::cout << i << " " << sqrt(h_out->GetBinContent(i,i)) << std::endl;
  } 

  return h_out; 

}

///////////////////////////////////////////////////////////////////////

}

#endif
