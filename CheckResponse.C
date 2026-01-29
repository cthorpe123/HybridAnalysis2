#include "Systematics.h"
#include "EnergyEstimatorFuncs.h"
#include "Histograms2.h"
#include "Funcs.h"
#include "PlotFuncs.h"

// Validate the multiplication of truth and response gives the reco in some universes

void CheckUniverse(std::string label,TFile* fh,std::string histstr,bool draw_o,bool draw_u,std::string univ=""){

  TH1D* h_truth = (TH1D*)fh->Get(("Truth/"+histstr+"/h_Signal"+univ).c_str()); 
  TH1D* h_reco = (TH1D*)fh->Get(("Reco/"+histstr+"/h_Signal"+univ).c_str()); 
  TH2D* h_joint = (TH2D*)fh->Get(("Joint/"+histstr+"/h_Signal"+univ).c_str()); 
  TH2D* h_res = (TH2D*)fh->Get(("Response/"+histstr+"/h_Signal"+univ).c_str()); 
  
  std::string plot_dir = "Analysis/"+label+"/Plots/CheckResponse/" + histstr + univ + "/";
  gSystem->Exec(("mkdir -p "+plot_dir).c_str());

  pfs::DrawUnstacked({h_truth},{1},{"Truth"},draw_o,draw_u,plot_dir+"Truth.png");
  pfs::DrawUnstacked({h_reco},{1},{"Reco"},draw_o,draw_u,plot_dir+"Reco.png");
  pfs::Draw2DHist(h_joint,plot_dir+"Joint.png");
  pfs::Draw2DHist(h_res,plot_dir+"Response.png");
 
  TH1D* h_folded_res = Multiply(h_truth,h_res); 
  pfs::DrawUnstacked({h_reco,h_folded_res},{1,2},{"Reco","Folded Truth"},draw_o,draw_u,plot_dir+"FoldedReco.png");

}

void CheckResponse(){

  std::vector<std::string> label_v = {"MuonMom"};
  std::vector<std::string> special_univ = {"ExtraPi","ExtraPi2"};

  bool draw_underflow = false;
  bool draw_overflow = false;

  for(std::string label : label_v){

    TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());

    CheckUniverse(label,f_in_hist,"CV",draw_underflow,draw_overflow);

    // Check 10 universes of multisim systematics
    for(int i_s=0;i_s<syst::kSystMAX;i_s++)
      for(int i=0;i<10;i++)
        CheckUniverse(label,f_in_hist,"Vars/"+syst::sys_str.at(i_s),draw_underflow,draw_overflow,"_"+std::to_string(i));

    // Check the unisim universes
    for(int i_s=0;i_s<syst::kUnisimMAX;i_s++)
      CheckUniverse(label,f_in_hist,"Vars/"+syst::unisims_str.at(i_s),draw_underflow,draw_overflow);

    if(!special_univ.size()) continue;

    for(std::string su : special_univ)
      CheckUniverse(label,f_in_hist,"Special/"+su,draw_underflow,draw_overflow); 


    // build the prediction by folding the CV through response matrices  
    TH1D* h_reco = (TH1D*)f_in_hist->Get("Reco/CV/h_Signal");
    TH1D* h_truth = (TH1D*)f_in_hist->Get("Truth/CV/h_Signal");
    TH2D* h_res = (TH2D*)f_in_hist->Get("Response/CV/h_Signal");   
 
    TH2D* h_Cov = (TH2D*)f_in_hist->Get("Reco/Cov/Total/Cov_Signal");   
    h_Cov->Reset();
    for(int i_s=0;i_s<syst::kSystMAX;i_s++){
      std::vector<TH1D*> h;
      for(int i_u=0;i_u<syst::sys_nuniv.at(i_s);i_u++)
        h.push_back(Multiply(h_truth,(TH2D*)f_in_hist->Get(("Response/Vars/"+syst::sys_str.at(i_s)+"/h_Signal_"+std::to_string(i_u)).c_str())));
       TH2D *c,*fc; 
       CalcCovMultisim(syst::sys_str.at(i_s),h,c,fc); 
       h_Cov->Add(c);
    }
  
    TH1D* h_reco_ff = Multiply(h_truth,h_res);
    for(int i=0;i<h_reco->GetNbinsX()+2;i++) h_reco_ff->SetBinError(i,sqrt(h_Cov->GetBinContent(i,i)));

    std::string plot_dir = "Analysis/"+label+"/Plots/CheckResponse/Folding/"; 
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());
  
    for(std::string su : special_univ){ 
      TH1D* h_reco_test = Multiply(h_truth,(TH2D*)f_in_hist->Get(("Response/Special/"+su+"/h_Signal").c_str()));
      pfs::DrawStacked({h_reco},{cat_colors[kSignal]},{"Folded Signal"},h_reco_ff,h_reco_test,draw_underflow,draw_overflow,plot_dir+"Reco_"+su+".png"); 
    }
    
    f_in_hist->Close();

  }

}
