#ifndef _BinningFuncs_h_
#define _BinningFuncs_h_

namespace binning {

const double _EPSILON_ = 1e-10;

std::vector<double> MakeEdges(const TH1D* h_data,double target_fe){

  if(h_data->Integral() < _EPSILON_ || 1.0/sqrt(h_data->Integral()) >= target_fe) return {};

  std::vector<double> bin_edges;

  double events = 0;
  for(int i=1;i<h_data->GetNbinsX()+1;i++){

    if(h_data->GetBinContent(i) < _EPSILON_) continue;

    if(!bin_edges.size()) bin_edges.push_back(h_data->GetBinLowEdge(i));

    events += h_data->GetBinContent(i);
    if(1.0/sqrt(events) < target_fe){
      bin_edges.push_back(h_data->GetBinLowEdge(i+1));
      events = 0;
    } 

  }

  if(!(bin_edges.size()-1)){
    std::cout << "No bin edges generated" << std::endl;
    return {};
  }

  std::cout << "Number of bins generated " << bin_edges.size()-1 << std::endl;
  std::cout << "Low edge = " << bin_edges.front() << " Hgh edge = " << bin_edges.back() << std::endl;

  return bin_edges;

} 

///////////////////////////////////////////////////////////////////////
// Generate variably binned histogram  

bool MakeBinningTemplate(std::string label,TH1D* h_data,bool truth=false,double target_fe=0.05){

  std::cout << "Generating binning template for " << label << std::endl;

  gSystem->Exec(("mkdir -p Analysis/"+label+"/rootfiles/").c_str());
  TFile* f_out = !truth ? TFile::Open(("Analysis/"+label+"/rootfiles/BinningTemplate.root").c_str(),"RECREATE")
                        : TFile::Open(("Analysis/"+label+"/rootfiles/TruthBinningTemplate.root").c_str(),"RECREATE");
  
  // If stats in channel are too low to have more than 1 bin
  if(h_data->Integral() < _EPSILON_ || 1.0/sqrt(h_data->Integral()) >= target_fe){
    std::cout << "Histogram is empty or doesn't contain enough data to generate at least one bin with target FE, saving single bin" << std::endl;
    TH1D* h_template = new TH1D(("h_template_"+label).c_str(),"",1,h_data->GetBinLowEdge(1),h_data->GetBinLowEdge(h_data->GetNbinsX()+1)); 
    h_template->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
    h_template->GetYaxis()->SetTitle(h_data->GetYaxis()->GetTitle());
    h_template->Write("h_template");
    f_out->Close();
    if(h_data->Integral() > _EPSILON_) return true;   
    else return false;
  }

  std::vector<double> bin_edges = MakeEdges(h_data,target_fe);

  if(bin_edges.size() < 2) 
    throw std::invalid_argument("BinningFuncs::MakeBinningTemplate: Too few edges");

  TH1D* h_template = new TH1D(("h_template_"+label).c_str(),"",bin_edges.size()-1,&bin_edges[0]); 
  h_template->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
  h_template->GetYaxis()->SetTitle(h_data->GetYaxis()->GetTitle());
  h_template->Write("h_template_All");
  f_out->Close();

  return true;

}

///////////////////////////////////////////////////////////////////////
// Generate several templates in a single file for multiple channels 

bool MakeMultiChannelTemplate(std::string label,std::map<std::string,TH1D*> h_data_m,bool truth=false,double target_fe=0.05){

  std::cout << "Generating binning template for " << label << std::endl;

  gSystem->Exec(("mkdir -p Analysis/"+label+"/rootfiles/").c_str());
  TFile* f_out = !truth ? TFile::Open(("Analysis/"+label+"/rootfiles/BinningTemplate.root").c_str(),"RECREATE")
                        : TFile::Open(("Analysis/"+label+"/rootfiles/TruthBinningTemplate.root").c_str(),"RECREATE");
 
  std::map<std::string,TH1D*>::iterator it;
  for(it = h_data_m.begin();it != h_data_m.end();it++){
       
    std::string ch = it->first;
    std::cout << "Channel " << ch << std::endl;
    TH1D* h_data = it->second;
        
    // If stats in channel are too low to have more than 1 bin
    if(h_data->Integral() < _EPSILON_ || 1.0/sqrt(h_data->Integral()) >= target_fe){
      std::cout << "Histogram is empty or doesn't contain enough data to generate at least one bin with target FE, saving single bin" << std::endl;
      TH1D* h_template = new TH1D(("h_template_"+ch+"_"+label).c_str(),"",1,h_data->GetBinLowEdge(1),h_data->GetBinLowEdge(h_data->GetNbinsX()+1)); 
      h_template->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
      h_template->GetYaxis()->SetTitle(h_data->GetYaxis()->GetTitle());
      h_template->Write(("h_template_"+ch).c_str());
      continue;
    }

    std::vector<double> bin_edges = MakeEdges(h_data,target_fe);

    if(bin_edges.size() < 2) 
      throw std::invalid_argument("BinningFuncs::MakeMultichannelTemplate: Too few edges");

    TH1D* h_template = new TH1D(("h_template_"+ch+"_"+label).c_str(),"",bin_edges.size()-1,&bin_edges[0]); 
    h_template->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
    h_template->GetYaxis()->SetTitle(h_data->GetYaxis()->GetTitle());
    h_template->Write(("h_template_"+ch).c_str());

  }

  f_out->Close();

  return true;

}

///////////////////////////////////////////////////////////////////////
// Root TH1::Add occasionally fails to add two histograms with same 
// limits - seems to be a bug.  

void ForceAdd(TH1D* h1, TH1D* h2){
  if(h1->GetNbinsX() != h2->GetNbinsX())
    throw std::invalid_argument("binning::ForceAdd: histograms have different numbers of bins");
  for(int i=0;i<h1->GetNbinsX()+2;i++){
    double c1 = h1->GetBinContent(i);
    double c2 = h2->GetBinContent(i);
    double e1 = h1->GetBinError(i);
    double e2 = h2->GetBinError(i);
    h1->SetBinContent(i,c1+c2);
    h1->SetBinError(i,sqrt(e1*e1+e2*e2));
  }
}

///////////////////////////////////////////////////////////////////////

}

#endif
