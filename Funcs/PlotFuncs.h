#ifndef _PlotFuncs_h_
#define _PlotFuncs_h_

namespace pfs {

///////////////////////////////////////////////////////////////////////
// Draw 2D histogram

void Draw2DHist(TH2D* h,std::string name){
  TCanvas* c = new TCanvas("c","c");
  h->Draw("colz");
  h->SetStats(0);
  c->Print(name.c_str());
  c->Close();
}

///////////////////////////////////////////////////////////////////////
// Draw a set of histograms on one canvas, not stacked together,
// as lines

void DrawUnstacked(std::vector<TH1D*> h_v,std::vector<int> colors,std::vector<std::string> legs,bool draw_o,bool draw_u,std::string name){

  THStack* hs_middle = new THStack("hs_middle",(";"+string(h_v.at(0)->GetXaxis()->GetTitle())+";FE").c_str());
  THStack* hs_U = new THStack("hs_U",";;FE");
  THStack* hs_O = new THStack("hs_O",";;FE");
  TLegend* l2 = new TLegend(0.75,0.75,0.98,0.98);
  l2->SetNColumns(2);

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
  if(draw_u){
    p_U->Draw();
    p_U->SetLeftMargin(0.37);
    p_U->SetRightMargin(0.04);
    p_U->cd();
    hs_U->Draw("nostack HIST");       
    //hs_U->SetMaximum(GetMax(h_v_Tot_U)*1.1);
    hs_U->GetXaxis()->SetLabelSize(0.16);
    hs_U->GetYaxis()->SetTitleSize(0.1);
    hs_U->GetYaxis()->SetLabelSize(0.11);
    hs_U->GetYaxis()->SetLabelOffset(0.04);
    hs_U->GetYaxis()->SetTitleOffset(1.8);
    c2->cd();
  }
  if(draw_o){
    p_O->Draw(); 
    p_O->SetRightMargin(0.37);
    p_O->SetLeftMargin(0.04);
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

void DrawStacked(std::vector<TH1D*> h_v,std::vector<int> colors,std::vector<std::string> legs,TH1D* h_tot,TH1D* h_data,bool draw_o,bool draw_u,std::string name){

  THStack* hs_middle = new THStack("hs_middle",(";"+string(h_v.at(0)->GetXaxis()->GetTitle())+";Events").c_str());
  THStack* hs_U = new THStack("hs_U",";;Events");
  THStack* hs_O = new THStack("hs_O",";;Events");
  TLegend* l2 = new TLegend(0.75,0.75,0.98,0.98);
  l2->SetNColumns(2);

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
  MakeOU(h_tot->GetName(),h_tot,h_tot_U,h_tot_O,"","");
  h_tot_U->SetFillStyle(3253);
  h_tot_U->SetFillColor(1);
  h_tot_O->SetFillStyle(3253);
  h_tot_O->SetFillColor(1);
  h_tot->SetFillStyle(3253);
  h_tot->SetFillColor(1);

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
  if(draw_u){
    p_U->Draw();
    p_U->SetLeftMargin(0.37);
    p_U->SetRightMargin(0.04);
    p_U->cd();
    hs_U->Draw("HIST");       
    h_tot_U->Draw("same e2");
    if(h_data != nullptr) h_data_U->Draw("same e1");
    hs_U->SetMaximum(GetMax(h_tot_U)*1.1);
    hs_U->GetXaxis()->SetLabelSize(0.16);
    hs_U->GetYaxis()->SetTitleSize(0.1);
    hs_U->GetYaxis()->SetLabelSize(0.11);
    hs_U->GetYaxis()->SetLabelOffset(0.04);
    hs_U->GetYaxis()->SetTitleOffset(1.8);
    c2->cd();
  }
  if(draw_o){
    p_O->Draw(); 
    p_O->SetRightMargin(0.37);
    p_O->SetLeftMargin(0.04);
    p_O->cd();
    hs_O->Draw("HIST Y+");       
    h_tot_O->Draw("same e2");
    if(h_data != nullptr) h_data_O->Draw("same e1");
    hs_O->SetMaximum(GetMax(h_tot_O)*1.1);
    hs_O->GetYaxis()->SetTitleSize(0.1);
    hs_O->GetXaxis()->SetLabelSize(0.16);
    hs_O->GetYaxis()->SetTitleOffset(1.8);
    hs_O->GetYaxis()->SetLabelSize(0.11);
    hs_O->GetYaxis()->SetLabelOffset(0.04);
    c2->cd();
  }    

  p_middle->cd(); 
  hs_middle->Draw("HIST");
  h_tot->Draw("same e2");
  if(h_data != nullptr) h_data->Draw("same e1");
  hs_middle->SetMaximum(GetMax(h_tot)*1.1);
  l2->Draw();

  c2->cd();
  c2->Print(name.c_str());
  delete c2;

}

}

#endif
