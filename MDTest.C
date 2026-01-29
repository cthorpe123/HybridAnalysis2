#include "Systematics.h"
#include "EnergyEstimatorFuncs.h"
#include "Histograms2.h"
#include "Funcs.h"
#include "PlotFuncs.h"

// Split MC into two sets and use first to calculate Covariance matrix 
// between elements of response matrix. The calculate mahalonobis distance
// of any special universes using this to see if response is roust 

void MDTest(){

  std::vector<std::string> label_v = {"MuonMom"};
  std::vector<std::string> special_univ = {"ExtraPi","ExtraPi2","ExtraP","ExtraP2"};

  bool draw_underflow = false;
  bool draw_overflow = false;

  std::vector<std::string> to_check = {"CV"};
  for(std::string u : special_univ) to_check.push_back("Special/"+u);

  for(std::string label : label_v){

    std::cout << "Doing MD check with " << label << std::endl;

    TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
    TFile* f_in_hist_2 = TFile::Open(("Analysis/"+label+"2/rootfiles/Histograms.root").c_str());

    TH2D* h_res = (TH2D*)f_in_hist->Get("Response/CV/h_Signal"); 
    std::vector<std::pair<int,int>> empty;
    for(int i=0;i<h_res->GetNbinsX()+2;i++)
      for(int j=0;j<h_res->GetNbinsY()+2;j++)
        if(h_res->GetBinContent(i,j) < 1e-10) empty.push_back(std::make_pair(i,j)); 

    std::cout << "Ignoring " << empty.size() << " bins out of " << (h_res->GetNbinsX()+2)*(h_res->GetNbinsY()+2) << " in response matrix" << std::endl;

    TH1D* h_res_u = Unroll2DDist(h_res,"h_res_u",empty);

    TH2D* h_Cov = Make2DHist("h_Cov_Res_Unroll",h_res_u);
    for(int i=0;i<h_res_u->GetNbinsX()+2;i++)
      h_Cov->SetBinContent(i,i,h_res_u->GetBinError(i)*h_res_u->GetBinError(i));

    for(int i_s=0;i_s<syst::kSystMAX;i_s++){
      std::cout << "Calculating cov for " << syst::sys_str.at(i_s) << std::endl;
      std::vector<TH1D*> h;
      for(int i_u=0;i_u<syst::sys_nuniv.at(i_s);i_u++)
        h.push_back(Unroll2DDist((TH2D*)f_in_hist->Get(("Response/Vars/"+syst::sys_str.at(i_s)+"/h_Signal_"+std::to_string(i_u)).c_str()),Form("h_%i",i_u),empty));
       TH2D *c,*fc; 
       CalcCovMultisim(syst::sys_str.at(i_s),h,c,fc); 
       h_Cov->Add(c);
       for(TH1D* x : h) delete x;       
       h.clear();
    }


    for(int i_s=0;i_s<syst::kUnisimMAX;i_s++){
      TH1D* h = Unroll2DDist((TH2D*)f_in_hist->Get(("Response/Vars/"+syst::unisims_str.at(i_s)+"/h_Signal").c_str()),"h",empty);
      TH2D *c,*fc; 
      CalcCovUnisim(syst::unisims_str.at(i_s),h_res_u,h,c,fc); 
      h_Cov->Add(c);
      delete h;
    }

    for(std::string u : to_check){ 

      std::cout << "Checking universe " << u << std::endl; 

      TH2D* h_Cov_tmp = (TH2D*)h_Cov->Clone("h_Cov_tmp");

      // Add the stat error from the CV of the second subset of events
      TH2D* h_res2 = (TH2D*)f_in_hist_2->Get(("Response/"+u+"/h_Signal").c_str()); 
      TH1D* h_res_u2 = Unroll2DDist(h_res2,"h_res_u2",empty);
      TH2D* h_Cov_Stat2 = Make2DHist("h_Cov_Res_Unroll_2",h_res_u2);
      for(int i=0;i<h_res_u2->GetNbinsX()+2;i++)
        h_Cov_Stat2->SetBinContent(i,i,h_res_u2->GetBinError(i)*h_res_u2->GetBinError(i));

      h_Cov_tmp->Add(h_Cov_Stat2);

      TMatrixDSym m_Cov = MakeCovMat(h_Cov_tmp);
      m_Cov.Invert();

      // Calculate squared Mahalanobis distance of special universe
      double md = 0.0;
      for(int i=1;i<h_res_u2->GetNbinsX()+1;i++)
        for(int j=1;j<h_res_u2->GetNbinsX()+1;j++)
          md += (h_res_u2->GetBinContent(i) - h_res_u->GetBinContent(i))*m_Cov[i-1][j-1]*(h_res_u2->GetBinContent(j) - h_res_u->GetBinContent(j));

      std::cout << md << "/" << h_res_u2->GetNbinsX() << " = " << md/h_res_u2->GetNbinsX() << std::endl;

      delete h_Cov_tmp;
      delete h_res_u2;
      delete h_Cov_Stat2;

    } 
 

  }

}
