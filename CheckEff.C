#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16 + ;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "EnergyEstimatorFuncs.h"
#include "Systematics.h"
#include "DetvarHistograms.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"

using namespace syst;

void CheckEff()
{

  // Label and set the branches defining the selection and systematics
  bool draw_truth = true;
  bool draw_hist = true; // Grab the CV from the non-detvar file
  bool draw_o = false, draw_u = false;

  std::vector<std::string> channels_t = {"1p", "2p"};
  std::vector<std::string> channels_r = {"1p", "2p"};

  std::vector<std::string> vars = {"NProt", "ProtonKE", "Channel"};

  for (std::string label : vars)
  {

    std::string plot_dir = "Analysis/" + label + "/Plots/CheckEff/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    TFile *f_in = draw_hist ? TFile::Open(("Analysis/" + label + "/rootfiles/Histograms.root").c_str()) : nullptr;

    hist::MultiChannelHistogramManager mchm(label, true);
    mchm.SetTrueChannelList(channels_t);
    mchm.SetRecoChannelList(channels_r);
    mchm.LoadTemplates();

    TH2D *h_CV_Response_Signal = (TH2D *)f_in->Get("Response/CV/h_Signal");
    int bins = h_CV_Response_Signal->GetNbinsY() + 1;

    TH1D *h_CV_Eff_Signal = h_CV_Response_Signal->ProjectionX("h_CV_Eff", 0, bins);

    std::vector<TH1D *> h_v;
    std::vector<int> fill_colors;
    std::vector<std::string> legs;
    for (size_t i_ch = 0; i_ch < channels_t.size(); i_ch++)
    {
      std::string ch = channels_t.at(i_ch);
      // std::cout << ch << std::endl;
      h_v.push_back((TH1D *)h_CV_Eff_Signal->Clone(("Eff_" + ch).c_str()));
      mchm.Restore(h_v.back(), ch, true);
      // for(int i=0;i<h_v.back()->GetNbinsX()+1;i++) std::cout << i <<  " " << h_v.back()->GetBinContent(i) << std::endl;
      fill_colors.push_back(i_ch + 1);
      legs.push_back(ch);
    }

    pfs::DrawUnstacked2(h_v, fill_colors, legs, plot_dir + "Efficiency.png", true);

    for (TH1D *h : h_v)
      delete h;
    h_v.clear();
    fill_colors.clear();
    legs.clear();

    // Calculate the fraction of events within X% of the true value (assumes true and reco are the same variable)
    int ctr = 1;
    for (size_t i_ch_t = 0; i_ch_t < channels_t.size(); i_ch_t++)
    {
      for (size_t i_ch_r = 0; i_ch_r < channels_r.size(); i_ch_r++)
      {
        // if(i_ch_t != i_ch_r) continue;
        std::string ch_t = channels_t.at(i_ch_t);
        std::string ch_r = channels_r.at(i_ch_r);
        TH2D *h_res = (TH2D *)h_CV_Response_Signal->Clone("h_res");
        mchm.RestoreRes(h_res, ch_t, ch_r);

        TH1D *h = (TH1D *)f_in->Get("Truth/CV/h_Signal");
        mchm.Restore(h, ch_t, true);
        h_v.push_back(h);
        fill_colors.push_back(ctr);
        legs.push_back(ch_t + " " + ch_r);

        pfs::Draw2DHist(h_res, plot_dir + "Response_" + ch_t + "_" + ch_r + ".png");

        for (int i = 1; i < h_res->GetNbinsX() + 1; i++)
        {
          double good = 0.0;
          double center = h_res->GetXaxis()->GetBinCenter(i);
          for (int j = 1; j < h_res->GetNbinsY() + 1; j++)
          {
            if (abs(h_res->GetYaxis()->GetBinCenter(j) - center) / center < 0.2)
              good += h_res->GetBinContent(i, j);
          }
          h->SetBinContent(i, good);
        }
        delete h_res;
        ctr++;
      }
    }

    pfs::DrawUnstacked2(h_v, fill_colors, legs, plot_dir + "FOM.png", false);

    h_v.clear();
    fill_colors.clear();
    legs.clear();

    // Try calculating the bias and variance in the observable for each channel
    ctr = 1;
    std::vector<TH1D *> h_bias, h_variance, h_bias_weighted;
    for (size_t i_ch_t = 0; i_ch_t < channels_t.size(); i_ch_t++)
    {
      for (size_t i_ch_r = 0; i_ch_r < channels_r.size(); i_ch_r++)
      {
        // if(i_ch_t != i_ch_r) continue;
        std::string ch_t = channels_t.at(i_ch_t);
        std::string ch_r = channels_r.at(i_ch_r);
        TH2D *h_res = (TH2D *)h_CV_Response_Signal->Clone("h_res");
        mchm.RestoreRes(h_res, ch_t, ch_r);

        h_bias.push_back((TH1D *)f_in->Get("Truth/CV/h_Signal")->Clone(("h_bias_" + ch_t + "_" + ch_r).c_str()));
        h_variance.push_back((TH1D *)f_in->Get("Truth/CV/h_Signal")->Clone(("h_variance_" + ch_t + "_" + ch_r).c_str()));
        h_bias_weighted.push_back((TH1D *)f_in->Get("Truth/CV/h_Signal")->Clone(("h_bias_weighted_" + ch_t + "_" + ch_r).c_str()));
        mchm.Restore(h_bias.back(), ch_t, true);
        mchm.Restore(h_variance.back(), ch_t, true);
        mchm.Restore(h_bias_weighted.back(), ch_t, true);

        GetBiasVariance(h_res, h_bias.back(), h_variance.back());
        GetWeightedBias(h_res, h_bias_weighted.back());

       // for(int i=1;i<h_bias_weighted.back()->GetNbinsX()+1;i++)
       //  std::cout << ch_t << " " << ch_r << " " << h_bias_weighted.back()->GetBinContent(i) << std::endl;


        fill_colors.push_back(ctr);
        legs.push_back(ch_t + " " + ch_r);
        ctr++;
      }
    }

    pfs::DrawUnstacked2(h_bias, fill_colors, legs, plot_dir + "Bias.png", false);
    pfs::DrawUnstacked2(h_variance, fill_colors, legs, plot_dir + "Variance.png", false);
    pfs::DrawUnstacked2(h_bias_weighted, fill_colors, legs, plot_dir + "WeightedBias.png", true);
  }
}
