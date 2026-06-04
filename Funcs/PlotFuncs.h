#ifndef _PlotFuncs_h_
#define _PlotFuncs_h_

#include "Funcs.h"

namespace pfs {

///////////////////////////////////////////////////////////////////////
// Draw 2D histogram

void Draw2DHist(TH2D* h,std::string name){
  TCanvas* c = new TCanvas("c2","c2");
  h->Draw("colz");
  h->SetStats(0);
  c->Print(name.c_str());
  c->Close();
}

///////////////////////////////////////////////////////////////////////
// Draw a set of histograms on one canvas, not stacked together,
// as lines

void DrawUnstacked(std::vector<TH1D*> h_v,std::vector<int> colors,std::vector<std::string> legs,bool draw_o,bool draw_u,std::string name){

  if(!h_v.size()) return;

  const std::string title = ";"+string(h_v.at(0)->GetXaxis()->GetTitle())+";"+string(h_v.at(0)->GetYaxis()->GetTitle());

  THStack* hs_middle = new THStack("hs_middle",title.c_str());
  THStack* hs_U = new THStack("hs_U",";;FE");
  THStack* hs_O = new THStack("hs_O",";;FE");
  //TLegend* l2 = new TLegend(0.75,0.75,0.98,0.98);
  TLegend* l2 = new TLegend(0.0,0.9,1.0,1.0);
  l2->SetNColumns(legs.size());
  l2->SetBorderSize(0);

  std::vector<TH1D*> h_O(h_v.size());
  std::vector<TH1D*> h_U(h_v.size());
  for(size_t i_s=0;i_s<h_v.size();i_s++){
    MakeOU(h_v.at(i_s)->GetName(),h_v.at(i_s),h_U.at(i_s),h_O.at(i_s));
    h_v.at(i_s)->SetLineColor(colors.at(i_s));
    h_U.at(i_s)->SetLineColor(colors.at(i_s));
    h_O.at(i_s)->SetLineColor(colors.at(i_s));
    h_v.at(i_s)->SetLineWidth(2);
    h_U.at(i_s)->SetLineWidth(2);
    h_O.at(i_s)->SetLineWidth(2);
    hs_middle->Add(h_v.at(i_s));
    hs_U->Add(h_U.at(i_s));
    hs_O->Add(h_O.at(i_s));
    l2->AddEntry(h_v.at(i_s),legs.at(i_s).c_str(),"L");
  }

  TCanvas* c2 = (draw_u || draw_o) ? new TCanvas("c2","c2",1000,600) : new TCanvas("c2","c2");
  double  split_low = draw_u ? 0.18 : 0.0;
  double  split_high = draw_o ? 0.82 : 1.0;
  TPad* p_U = draw_u ? new TPad("p_U","p_U",0.0,0.0,split_low,1.0) : nullptr;
  TPad* p_middle = new TPad("p_middle","p_middle",split_low,0.0,split_high,1.0);
  TPad* p_O = draw_o ? new TPad("p_O","p_O",split_high,0.0,1.0,1.0) : nullptr;

  c2->cd();
  p_middle->Draw();
  p_middle->SetRightMargin(0.03);
  p_middle->SetLeftMargin(0.11);
  if(draw_u && draw_o) p_middle->SetLeftMargin(0.13);
  if(draw_u){
    p_U->Draw();
    p_U->SetLeftMargin(0.44);
    p_U->SetRightMargin(0.1);
    p_U->cd();
    hs_U->Draw("nostack HIST");       
    //hs_U->SetMaximum(GetMax(h_v_Tot_U)*1.1);
    hs_U->GetXaxis()->SetLabelSize(0.16);
    hs_U->GetYaxis()->SetTitleSize(0.1);
    hs_U->GetYaxis()->SetLabelSize(0.11);
    hs_U->GetYaxis()->SetLabelOffset(0.04);
    //hs_U->GetYaxis()->SetTitleOffset(1.8);
    c2->cd();
  }
  if(draw_o){
    p_O->Draw(); 
    p_O->SetRightMargin(0.44);
    p_O->SetLeftMargin(0.1);
    p_O->cd();
    hs_O->Draw("nostack HIST Y+");       
    //hs_O->SetMaximum(GetMax(h_v_Tot_O)*1.1);
    hs_O->GetYaxis()->SetTitleSize(0.1);
    hs_O->GetXaxis()->SetLabelSize(0.16);
    hs_O->GetYaxis()->SetTitleOffset(1.8);
    hs_O->GetYaxis()->SetLabelSize(0.11);
    hs_O->GetYaxis()->SetLabelOffset(0.04);
    c2->cd();
  }    

  p_middle->cd(); 
  hs_middle->Draw("nostack HIST");
  //hs_middle->SetMaximum(GetMax(h_v_Tot)*1.1);
  l2->Draw();

  c2->cd();
  c2->Print(name.c_str());
  delete c2;

}

///////////////////////////////////////////////////////////////////////
// Draw a set of histograms on one canvas, not stacked together,
// as lines

void DrawStacked(std::vector<TH1D*> h_v,std::vector<int> colors,std::vector<std::string> legs,TH1D* h_tot,TH1D* h_data,bool draw_o,bool draw_u,std::string name,std::pair<double,int> chi2={0,-1}){

  if(!h_v.size()) return;
  
  // Sometimes the first hist and h_tot are the same pointer that we now need to behave separately, make a clone of the 
  // tot and work on that
  TH1D* h_tot_tmp = (TH1D*)h_tot->Clone("h_tot_tmp");

  THStack* hs_middle = new THStack("hs_middle",(";"+string(h_tot_tmp->GetXaxis()->GetTitle())+";"+string(h_tot_tmp->GetYaxis()->GetTitle())).c_str());
  THStack* hs_U = new THStack("hs_U",";;Events");
  THStack* hs_O = new THStack("hs_O",";;Events");
  TLegend* l2 = draw_o && draw_u ? new TLegend(0.13,0.91,0.97,1.0): new TLegend(0.11,0.91,0.97,1.0);
  l2->SetNColumns(legs.size());
  l2->SetBorderSize(0);

  std::vector<TH1D*> h_O(h_v.size());
  std::vector<TH1D*> h_U(h_v.size());
  for(size_t i_s=0;i_s<h_v.size();i_s++){
    MakeOU(h_v.at(i_s)->GetName(),h_v.at(i_s),h_U.at(i_s),h_O.at(i_s));
    h_v.at(i_s)->SetFillColor(colors.at(i_s));
    h_U.at(i_s)->SetFillColor(colors.at(i_s));
    h_O.at(i_s)->SetFillColor(colors.at(i_s));
    hs_middle->Add(h_v.at(i_s));
    hs_U->Add(h_U.at(i_s));
    hs_O->Add(h_O.at(i_s));
    l2->AddEntry(h_v.at(i_s),legs.at(i_s).c_str(),"F");
  }

  TH1D *h_tot_O,*h_tot_U;
  MakeOU(h_tot_tmp->GetName(),h_tot_tmp,h_tot_U,h_tot_O,"","");
  h_tot_U->SetFillStyle(3253);
  h_tot_U->SetFillColor(1);
  h_tot_O->SetFillStyle(3253);
  h_tot_O->SetFillColor(1);
  h_tot_tmp->SetFillStyle(3253);
  h_tot_tmp->SetFillColor(1);

  TH1D *h_data_U=nullptr,*h_data_O=nullptr;
  if(h_data != nullptr){
    MakeOU(h_data->GetName(),h_data,h_data_U,h_data_O);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(0.8);
    h_data->SetLineColor(1);
    h_data_U->SetMarkerStyle(20);
    h_data_U->SetMarkerSize(0.8);
    h_data_U->SetLineColor(1);
    h_data_O->SetMarkerStyle(20);
    h_data_O->SetMarkerSize(0.8);
    h_data_O->SetLineColor(1);
  }

  TCanvas* c2 = (draw_u || draw_o) ? new TCanvas("c2","c2",1000,600) : new TCanvas("c2","c2");
  double  split_low = draw_u ? 0.18 : 0.0;
  double  split_high = draw_o ? 0.82 : 1.0;
  TPad* p_U = draw_u ? new TPad("p_U","p_U",0.0,0.0,split_low,1.0) : nullptr;
  TPad* p_middle = new TPad("p_middle","p_middle",split_low,0.0,split_high,1.0);
  TPad* p_O = draw_o ? new TPad("p_O","p_O",split_high,0.0,1.0,1.0) : nullptr;

  c2->cd();
  p_middle->Draw();
  p_middle->SetRightMargin(0.03);
  p_middle->SetLeftMargin(0.11);
  if(draw_u && draw_o) p_middle->SetLeftMargin(0.13);
  if(draw_u){
    p_U->Draw();
    p_U->SetLeftMargin(0.44);
    p_U->SetRightMargin(0.1);
    p_U->cd();
    hs_U->Draw("HIST");       
    h_tot_U->Draw("same e2");
    if(h_data != nullptr) h_data_U->Draw("same e1");
    hs_U->SetMaximum(GetMax(h_tot_U)*1.1);
    hs_U->GetXaxis()->SetLabelSize(0.16);
    hs_U->GetYaxis()->SetTitleSize(0.11);
    hs_U->GetYaxis()->SetLabelSize(0.11);
    hs_U->GetYaxis()->SetLabelOffset(0.04);
    //hs_U->GetYaxis()->SetTitleOffset(1.8);
    c2->cd();
  }
  if(draw_o){
    p_O->Draw(); 
    p_O->SetRightMargin(0.44);
    p_O->SetLeftMargin(0.1);
    p_O->cd();
    hs_O->Draw("HIST Y+");       
    h_tot_O->Draw("same e2");
    if(h_data != nullptr) h_data_O->Draw("same e1");
    hs_O->SetMaximum(GetMax(h_tot_O)*1.1);
    hs_O->GetYaxis()->SetTitleSize(0.11);
    hs_O->GetXaxis()->SetLabelSize(0.16);
    hs_O->GetYaxis()->SetTitleOffset(1.8);
    hs_O->GetYaxis()->SetLabelSize(0.11);
    hs_O->GetYaxis()->SetLabelOffset(0.04);
    c2->cd();
  }    

  p_middle->cd(); 

  TLegend *l_Chi2 = draw_o && draw_u ? new TLegend(0.15,0.83,0.37,0.895) : new TLegend(0.13,0.83,0.37,0.895);
  l_Chi2->SetBorderSize(0);
  l_Chi2->SetMargin(0.005);
  l_Chi2->SetTextAlign(12);
  l_Chi2->SetTextSize(0.05);
  l_Chi2->SetHeader(("#chi^{2}/ndof = " + to_string_with_precision(chi2.first,1) + "/" + std::to_string(chi2.second)).c_str());

  hs_middle->Draw("HIST");
  h_tot_tmp->Draw("same e2");
  if(h_data != nullptr) h_data->Draw("same e1");
  hs_middle->SetMaximum(GetMax(h_tot_tmp)*1.15);
  l2->Draw();
  if(chi2.second > 0) l_Chi2->Draw();
  c2->cd();
  c2->Print(name.c_str());
  delete c2;

  delete hs_U;
  delete hs_O;
  delete h_tot_U;
  delete h_tot_O;

  if(h_data != nullptr){
    delete h_data_U;
    delete h_data_O;
  }

  for(size_t i_s=0;i_s<h_v.size();i_s++){
    delete h_U.at(i_s);
    delete h_O.at(i_s);
  }

   h_U.clear();
   h_O.clear();
   delete h_tot_tmp;

}

///////////////////////////////////////////////////////////////////////
// Like DrawStacked but adds a row of ratio panels below showing
// data/MC, with a hatched band for the MC uncertainty (same style as
// the upper plot) and a reference line at 1.

void DrawStackedRatio(std::vector<TH1D*> h_v,std::vector<int> colors,std::vector<std::string> legs,TH1D* h_tot,TH1D* h_data,bool draw_o,bool draw_u,std::string name,std::pair<double,int> chi2={0,-1}){

  if(!h_v.size()) return;

  TH1D* h_tot_tmp = (TH1D*)h_tot->Clone("h_tot_tmp_r");

  THStack* hs_middle = new THStack("hs_middle_r",(";"+string(h_tot_tmp->GetXaxis()->GetTitle())+";"+string(h_tot_tmp->GetYaxis()->GetTitle())).c_str());
  THStack* hs_U      = new THStack("hs_U_r",";;Events");
  THStack* hs_O      = new THStack("hs_O_r",";;Events");
  TLegend* l2 = draw_o && draw_u ? new TLegend(0.13,0.91,0.97,1.0) : new TLegend(0.11,0.91,0.97,1.0);
  l2->SetNColumns(legs.size());
  l2->SetBorderSize(0);

  std::vector<TH1D*> h_O(h_v.size()), h_U(h_v.size());
  for(size_t i_s=0;i_s<h_v.size();i_s++){
    MakeOU(h_v.at(i_s)->GetName(),h_v.at(i_s),h_U.at(i_s),h_O.at(i_s));
    h_v.at(i_s)->SetFillColor(colors.at(i_s));
    h_U.at(i_s)->SetFillColor(colors.at(i_s));
    h_O.at(i_s)->SetFillColor(colors.at(i_s));
    hs_middle->Add(h_v.at(i_s));
    hs_U->Add(h_U.at(i_s));
    hs_O->Add(h_O.at(i_s));
    l2->AddEntry(h_v.at(i_s),legs.at(i_s).c_str(),"F");
  }

  TH1D *h_tot_O,*h_tot_U;
  MakeOU(h_tot_tmp->GetName(),h_tot_tmp,h_tot_U,h_tot_O,"","");
  h_tot_U->SetFillStyle(3253); h_tot_U->SetFillColor(1);
  h_tot_O->SetFillStyle(3253); h_tot_O->SetFillColor(1);
  h_tot_tmp->SetFillStyle(3253); h_tot_tmp->SetFillColor(1);

  TH1D *h_data_U=nullptr, *h_data_O=nullptr;
  if(h_data != nullptr){
    MakeOU(h_data->GetName(),h_data,h_data_U,h_data_O);
    h_data->SetMarkerStyle(20);   h_data->SetMarkerSize(0.8);   h_data->SetLineColor(1);
    h_data_U->SetMarkerStyle(20); h_data_U->SetMarkerSize(0.8); h_data_U->SetLineColor(1);
    h_data_O->SetMarkerStyle(20); h_data_O->SetMarkerSize(0.8); h_data_O->SetLineColor(1);
  }

  // Build data/MC ratio histogram; error = data_err / tot
  auto MakeRatio = [](const char* hname, TH1D* data_h, TH1D* tot_h) -> TH1D* {
    TH1D* h = (TH1D*)tot_h->Clone(hname);
    h->Reset();
    h->SetFillStyle(0);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->SetLineColor(1);
    h->SetStats(0);
    for(int i=0; i<=h->GetNbinsX()+1; ++i){
      double tot = tot_h->GetBinContent(i);
      if(tot > 0){
        h->SetBinContent(i, data_h->GetBinContent(i) / tot);
        h->SetBinError(i,   data_h->GetBinError(i)   / tot);
      }
    }
    return h;
  };

  // Build MC uncertainty band for the ratio panel: centred at 1, err = tot_err / tot
  auto MakeBand = [](const char* hname, TH1D* tot_h) -> TH1D* {
    TH1D* h = (TH1D*)tot_h->Clone(hname);
    h->Reset();
    h->SetStats(0);
    h->SetFillStyle(3253);
    h->SetFillColor(1);
    h->SetLineColor(0);
    for(int i=0; i<=h->GetNbinsX()+1; ++i){
      double tot = tot_h->GetBinContent(i);
      if(tot > 0){
        h->SetBinContent(i, 1.0);
        h->SetBinError(i,   tot_h->GetBinError(i) / tot);
      }
    }
    return h;
  };

  TH1D* h_ratio   = h_data ? MakeRatio("h_ratio_mid_r", h_data,   h_tot_tmp) : nullptr;
  TH1D* h_ratio_U = h_data ? MakeRatio("h_ratio_U_r",   h_data_U, h_tot_U)   : nullptr;
  TH1D* h_ratio_O = h_data ? MakeRatio("h_ratio_O_r",   h_data_O, h_tot_O)   : nullptr;
  TH1D* h_band    = MakeBand("h_band_mid_r", h_tot_tmp);
  TH1D* h_band_U  = MakeBand("h_band_U_r",  h_tot_U);
  TH1D* h_band_O  = MakeBand("h_band_O_r",  h_tot_O);

  // Canvas: split 70 % top (main) / 30 % bottom (ratio)
  const double rs         = 0.3;
  const double split_low  = draw_u ? 0.18 : 0.0;
  const double split_high = draw_o ? 0.82 : 1.0;

  TCanvas* c2 = (draw_u || draw_o) ? new TCanvas("c2","c2",1000,800) : new TCanvas("c2","c2",700,800);

  TPad* p_U_top = draw_u ? new TPad("p_U_top_r",  "",0.0,       rs,split_low, 1.0) : nullptr;
  TPad* p_U_bot = draw_u ? new TPad("p_U_bot_r",  "",0.0,      0.0,split_low,  rs) : nullptr;
  TPad* p_mid_top =         new TPad("p_mid_top_r","",split_low, rs,split_high,1.0);
  TPad* p_mid_bot =         new TPad("p_mid_bot_r","",split_low,0.0,split_high, rs);
  TPad* p_O_top = draw_o ? new TPad("p_O_top_r",  "",split_high, rs,1.0,      1.0) : nullptr;
  TPad* p_O_bot = draw_o ? new TPad("p_O_bot_r",  "",split_high,0.0,1.0,       rs) : nullptr;

  c2->cd();
  p_mid_top->Draw(); p_mid_bot->Draw();
  p_mid_top->SetRightMargin(0.03); p_mid_top->SetBottomMargin(0.06);
  p_mid_bot->SetRightMargin(0.03); p_mid_bot->SetTopMargin(0.01); p_mid_bot->SetBottomMargin(0.2);
  if(draw_u && draw_o){ p_mid_top->SetLeftMargin(0.13); p_mid_bot->SetLeftMargin(0.13); }
  else                { p_mid_top->SetLeftMargin(0.11); p_mid_bot->SetLeftMargin(0.11); }

  if(draw_u){
    p_U_top->Draw(); p_U_bot->Draw();
    p_U_top->SetLeftMargin(0.44); p_U_top->SetRightMargin(0.1); p_U_top->SetBottomMargin(0.06);
    p_U_bot->SetLeftMargin(0.44); p_U_bot->SetRightMargin(0.1); p_U_bot->SetTopMargin(0.01); p_U_bot->SetBottomMargin(0.2);
  }
  if(draw_o){
    p_O_top->Draw(); p_O_bot->Draw();
    p_O_top->SetRightMargin(0.44); p_O_top->SetLeftMargin(0.1); p_O_top->SetBottomMargin(0.06);
    p_O_bot->SetRightMargin(0.44); p_O_bot->SetLeftMargin(0.1); p_O_bot->SetTopMargin(0.01); p_O_bot->SetBottomMargin(0.2);
  }

  
  // Helper: draw reference line at 1 on the current pad
  auto DrawRatioLine = [](TH1D* h){
    double xlo = h->GetXaxis()->GetXmin(), xhi = h->GetXaxis()->GetXmax();
    TLine* l1 = new TLine(xlo, 1, xhi, 1); l1->SetLineWidth(2); l1->Draw();
  };

  // Underflow column
  if(draw_u){
    p_U_top->cd();
    hs_U->Draw("HIST");
    h_tot_U->Draw("same e2");
    if(h_data) h_data_U->Draw("same e1");
    hs_U->SetMaximum(GetMax(h_tot_U)*1.1);
    hs_U->GetXaxis()->SetLabelSize(0);
    hs_U->GetYaxis()->SetTitleSize(0.11);
    hs_U->GetYaxis()->SetLabelSize(0.11);
    hs_U->GetYaxis()->SetLabelOffset(0.04);

    if(h_ratio_U){
      p_U_bot->cd();
      double dev_U = std::abs(h_ratio_U->GetBinContent(1) - 1.0) + h_ratio_U->GetBinError(1);
      dev_U = std::max(dev_U, h_band_U->GetBinError(1));
      double rmax_U = 1.0 + dev_U * 1.3;
      double rmin_U = 1.0 - dev_U * 1.3;
      if(rmax_U == rmin_U){ rmax_U = 1.5; rmin_U = 0.5; }
      h_ratio_U->SetMaximum(rmax_U); h_ratio_U->SetMinimum(rmin_U);
      h_ratio_U->Draw("e1");
      h_band_U->Draw("e2 same");
      h_ratio_U->Draw("e1 same");
      h_ratio_U->GetXaxis()->SetLabelSize(0.16);
      h_ratio_U->GetYaxis()->SetLabelSize(0.11);
      h_ratio_U->GetYaxis()->SetLabelOffset(0.04);
      gPad->RedrawAxis();
      DrawRatioLine(h_ratio_U);
    }
    c2->cd();
  }

  // Overflow column
  if(draw_o){
    p_O_top->cd();
    hs_O->Draw("HIST Y+");
    h_tot_O->Draw("same e2");
    if(h_data) h_data_O->Draw("same e1");
    hs_O->SetMaximum(GetMax(h_tot_O)*1.1);
    hs_O->GetXaxis()->SetLabelSize(0);
    hs_O->GetYaxis()->SetTitleSize(0.11);
    hs_O->GetYaxis()->SetTitleOffset(1.8);
    hs_O->GetYaxis()->SetLabelSize(0.11);
    hs_O->GetYaxis()->SetLabelOffset(0.04);

    if(h_ratio_O){
      p_O_bot->cd();
      double dev_O = std::abs(h_ratio_O->GetBinContent(1) - 1.0) + h_ratio_O->GetBinError(1);
      dev_O = std::max(dev_O, h_band_O->GetBinError(1));
      double rmax_O = 1.0 + dev_O * 1.3;
      double rmin_O = 1.0 - dev_O * 1.3;
      if(rmax_O == rmin_O){ rmax_O = 1.5; rmin_O = 0.5; }
      h_ratio_O->SetMaximum(rmax_O); h_ratio_O->SetMinimum(rmin_O);
      h_ratio_O->Draw("e1 Y+");
      h_band_O->Draw("e2 same");
      h_ratio_O->Draw("e1 Y+ same");
      h_ratio_O->GetXaxis()->SetLabelSize(0.16);
      h_ratio_O->GetYaxis()->SetLabelSize(0.11);
      h_ratio_O->GetYaxis()->SetLabelOffset(0.04);
      gPad->RedrawAxis();
      DrawRatioLine(h_ratio_O);
    }
    c2->cd();
  }

  // Main column — top (main plot)
  p_mid_top->cd();

  TLegend* l_Chi2 = draw_o && draw_u ? new TLegend(0.15,0.83,0.37,0.895) : new TLegend(0.13,0.83,0.37,0.895);
  l_Chi2->SetBorderSize(0);
  l_Chi2->SetMargin(0.005);
  l_Chi2->SetTextAlign(12);
  l_Chi2->SetTextSize(0.05);
  l_Chi2->SetHeader(("#chi^{2}/ndof = "+to_string_with_precision(chi2.first,1)+"/"+std::to_string(chi2.second)).c_str());

  hs_middle->Draw("HIST");
  h_tot_tmp->Draw("same e2");
  if(h_data) h_data->Draw("same e1");
  hs_middle->SetMaximum(GetMax(h_tot_tmp)*1.15);
  hs_middle->GetXaxis()->SetLabelSize(0);
  hs_middle->GetXaxis()->SetTitleSize(0);
  l2->Draw();
  if(chi2.second > 0) l_Chi2->Draw();

  // Main column — bottom (ratio)
  if(h_ratio){
    p_mid_bot->cd();

    double dev = 0;
    for(int i=1; i<=h_ratio->GetNbinsX(); ++i){
      double v = std::abs(h_ratio->GetBinContent(i) - 1.0) + h_ratio->GetBinError(i);
      double b = h_band->GetBinError(i);
      if(v > dev) dev = v;
      if(b > dev) dev = b;
    }
    dev = dev * 1.1;
    if(dev == 0) dev = 0.5;
    h_ratio->SetMaximum(1.0 + dev); h_ratio->SetMinimum(1.0 - dev);
    h_ratio->Draw("e1");
    h_band->Draw("e2 same");
    h_ratio->Draw("e1 same");
    h_ratio->GetYaxis()->SetTitle("Data/MC");
    h_ratio->GetYaxis()->SetTitleSize(0.075);
    h_ratio->GetYaxis()->SetLabelSize(0.08);
    h_ratio->GetYaxis()->SetLabelOffset(0.01);
    h_ratio->GetYaxis()->SetTitleOffset(0.8);
    h_ratio->GetXaxis()->SetTitle(h_tot_tmp->GetXaxis()->GetTitle());
    h_ratio->GetXaxis()->SetTitleSize(0.075);
    h_ratio->GetXaxis()->SetTitleOffset(0.95);
    h_ratio->GetXaxis()->SetLabelSize(0.08);
    gPad->RedrawAxis();
    DrawRatioLine(h_ratio);
  }

  c2->cd();
  c2->Print(name.c_str());
  delete c2;

  delete hs_U; delete hs_O;
  delete h_tot_U; delete h_tot_O;
  if(h_data){ delete h_data_U; delete h_data_O; }
  if(h_ratio)   delete h_ratio;
  if(h_ratio_U) delete h_ratio_U;
  if(h_ratio_O) delete h_ratio_O;
  delete h_band; delete h_band_U; delete h_band_O;
  for(size_t i_s=0;i_s<h_v.size();i_s++){ delete h_U.at(i_s); delete h_O.at(i_s); }
  h_U.clear(); h_O.clear();
  delete h_tot_tmp;

}

///////////////////////////////////////////////////////////////////////
// Draw a set of histograms on one canvas with different binning
// ranges, as lines, optionally with error bars

void DrawUnstacked2(std::vector<TH1D*> h_v,std::vector<int> colors,std::vector<std::string> legs,std::string name,bool draw_errors){

  if(!h_v.size()) return;

  const std::string title = ";"+string(h_v.at(0)->GetXaxis()->GetTitle())+";"+string(h_v.at(0)->GetYaxis()->GetTitle());

  // Find the ranges needed and generate a template with the right axis ranges
  double xmin = h_v.at(0)->GetBinLowEdge(1), xmax = h_v.at(0)->GetBinLowEdge(h_v.at(0)->GetNbinsX()+1);
  double ymin = h_v.at(0)->GetMinimum(), ymax = h_v.at(0)->GetMaximum();

  for(TH1D* h : h_v){
    xmin = std::min(xmin,h->GetBinLowEdge(1));
    xmax = std::max(xmax,h->GetBinLowEdge(h->GetNbinsX()+1));
    for(int i=1;i<=h->GetNbinsX();i++){
      double val = h->GetBinContent(i);
      double err = draw_errors ? h->GetBinError(i) : 0.0;
      ymin = std::min(ymin, val - err);
      ymax = std::max(ymax, val + err);
    }
  }
  double xrange = xmax - xmin;
  double yrange = ymax - ymin;
 
  THStack* hs_middle = new THStack("hs_middle",title.c_str());
  //TLegend* l2 = new TLegend(0.75,0.75,0.98,0.98);
  TLegend* l2 = new TLegend(0.0,0.9,1.0,1.0);
  l2->SetNColumns(legs.size());
  l2->SetBorderSize(0);

  for(size_t i_s=0;i_s<h_v.size();i_s++){
    h_v.at(i_s)->SetLineColor(colors.at(i_s));
    h_v.at(i_s)->SetLineWidth(2);
    hs_middle->Add(h_v.at(i_s));
    l2->AddEntry(h_v.at(i_s),legs.at(i_s).c_str(),"L");
  }

  TCanvas* c2 = new TCanvas("c2","c2");

  TH1D* h = new TH1D("h",title.c_str(),1,xmin-0.05*xrange,xmax+0.05*xrange);
  h->Draw();
  h->SetStats(0);
  if(!draw_errors)hs_middle->Draw("nostack HIST same");
  else hs_middle->Draw("nostack E same");
  h->SetMaximum(ymax+0.05*yrange);
  h->SetMinimum(ymin-0.05*yrange);

  l2->Draw();

  c2->cd();
  c2->Print(name.c_str());
  delete c2;
  delete h;

}

}

#endif
