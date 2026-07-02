#include "Funcs.h"

void print_uboone_flux(){

  TFile* f_in = TFile::Open("/pnfs/uboone/persistent/uboonebeam/bnb_gsimple/bnb_gsimple_fluxes_01.09.2019_463_hist/MCC9_FluxHist_volTPCActive.root");
  TH1D* h_numu = static_cast<TH1D*>(f_in->Get("hEnumubar_cv"));
  h_numu->Scale(1/(4997.*5e8)/(256.35*233.)); // Accodring to readme this gives flux in nu/POT/bin/cm2

  TCanvas* c = new TCanvas("c","c");

  h_numu->Draw("HIST");
  DivideByBinWidth(h_numu);
  //h_numu->Scale(460.0*460.0/541.0/541.0); // To scale to Miniboone baseline
  h_numu->GetYaxis()->SetTitle("#nu/POT/GeV/cm^{2}");
  c->Print("uboone_numu_flux_scaled.png");

  std::cout << "Total flux: " << h_numu->Integral("width") << " nu/POT/cm2" << std::endl; 

  // Sample plot with different units - nu/50MeV/m2/10^6 POT
  h_numu->Scale(1e6*100*100/20);
  h_numu->GetYaxis()->SetTitle("#nu/50 MeV/m^{2}/10^{6} POT");

  c->SetLogy();
  h_numu->GetXaxis()->SetRangeUser(0.0,3.0);
  h_numu->Draw("HIST");
  c->Print("uboone_numu_flux_scaled_other_units.png");

  f_in->Close();
 

}
