#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs/BranchList.h"
#include "Funcs/PlotFuncs.h"

// Counts events (total, signal, selected per reco) in a filtered ntuple and
// normalises by the sample exposure (POT for overlay/dirt, triggers for beam-off).
// Produces one figure per sample type (overlay, dirt, beam-off) comparing
// event rates across runs.

struct RateSummary {
  double total  = 0;
  double signal = 0; // only valid for overlay
  double sel_pd = 0;
  double sel_wc = 0;
  double sel_lt = 0;
  double sel_h8 = 0;
  // Poisson statistical errors: sqrt(N) / exposure
  double err_total  = 0;
  double err_signal = 0;
  double err_sel_pd = 0;
  double err_sel_wc = 0;
  double err_sel_lt = 0;
  double err_sel_h8 = 0;
};

RateSummary ProcessFile(const std::string& filepath) {

  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  bool is_overlay = false, load_syst = false;
  LoadTreeFiltered(filepath, f_in, t_in, is_overlay, load_syst);

  // POT is constant per file — read it from the first entry
  t_in->GetEntry(0);
  double sample_exposure = POT;

  RateSummary r;
  long long n = t_in->GetEntries();
  double sw_total  = 0, sw_signal  = 0, sw_pd  = 0, sw_wc  = 0, sw_lt  = 0, sw_h8  = 0;
  double sw2_total = 0, sw2_signal = 0, sw2_pd = 0, sw2_wc = 0, sw2_lt = 0, sw2_h8 = 0;
  for (long long i = 0; i < n; ++i) {
    if (i % 1000 == 0)
      std::cout << "  " << i << " / " << n << "\r" << std::flush;
    t_in->GetEntry(i);
    double w  = weightSplineTimesTune;
    if (std::isnan(w) || std::isinf(w)) continue;
    double w2 = w * w;
    sw_total  += w; sw2_total  += w2;
    if (is_overlay && is_signal_t) { sw_signal += w; sw2_signal += w2; }
    if (sel_pd) { sw_pd += w; sw2_pd += w2; }
    if (sel_wc) { sw_wc += w; sw2_wc += w2; }
    if (sel_lt) { sw_lt += w; sw2_lt += w2; }
    if (sel_h8) { sw_h8 += w; sw2_h8 += w2; }
  }
  std::cout << std::endl;

  r.total   = sw_total   / sample_exposure;
  r.signal  = sw_signal  / sample_exposure;
  r.sel_pd  = sw_pd      / sample_exposure;
  r.sel_wc  = sw_wc      / sample_exposure;
  r.sel_lt  = sw_lt      / sample_exposure;
  r.sel_h8  = sw_h8      / sample_exposure;

  r.err_total   = std::sqrt(sw2_total)   / sample_exposure;
  r.err_signal  = std::sqrt(sw2_signal)  / sample_exposure;
  r.err_sel_pd  = std::sqrt(sw2_pd)      / sample_exposure;
  r.err_sel_wc  = std::sqrt(sw2_wc)      / sample_exposure;
  r.err_sel_lt  = std::sqrt(sw2_lt)      / sample_exposure;
  r.err_sel_h8  = std::sqrt(sw2_h8)      / sample_exposure;

  f_in->Close();
  return r;
}

// Build a styled TH1D with run bin labels
TH1D* MakeLabelledHist(const char* name, const std::vector<std::string>& labels) {
  int n = labels.size();
  TH1D* h = new TH1D(name, "", n, 0, n);
  for (int i = 0; i < n; ++i)
    h->GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  h->GetXaxis()->SetLabelSize(0.07);
  h->GetYaxis()->SetTitleOffset(1.3);
  return h;
}


void FillHist(TH1D* h, const std::vector<double>& vals, const std::vector<double>& errs) {
  for (int i = 0; i < (int)vals.size(); ++i) {
    h->SetBinContent(i + 1, vals[i]);
    h->SetBinError(i + 1, errs[i]);
  }
}

// Build, fill, and save one canvas per series for a given sample type
void DrawFigure(const std::string& title,
                const std::string& y_label,
                const std::string& out_path,
                const std::vector<std::string>& run_labels,
                const std::vector<RateSummary>& rates,
                bool show_signal) {

  TH1D* h_total  = MakeLabelledHist("h_total",  run_labels);
  TH1D* h_signal = MakeLabelledHist("h_signal", run_labels);
  TH1D* h_sel_pd = MakeLabelledHist("h_sel_pd", run_labels);
  TH1D* h_sel_wc = MakeLabelledHist("h_sel_wc", run_labels);
  TH1D* h_sel_lt = MakeLabelledHist("h_sel_lt", run_labels);
  TH1D* h_sel_h8 = MakeLabelledHist("h_sel_h8", run_labels);

  for (TH1D* h : {h_total, h_signal, h_sel_pd, h_sel_wc, h_sel_lt, h_sel_h8})
    h->GetYaxis()->SetTitle(y_label.c_str());

  std::vector<double> vt, vs, vpd, vwc, vlt, vh8;
  std::vector<double> et, es, epd, ewc, elt, eh8;
  for (const auto& r : rates) {
    vt.push_back(r.total);   vs.push_back(r.signal);
    vpd.push_back(r.sel_pd); vwc.push_back(r.sel_wc);
    vlt.push_back(r.sel_lt); vh8.push_back(r.sel_h8);
    et.push_back(r.err_total);   es.push_back(r.err_signal);
    epd.push_back(r.err_sel_pd); ewc.push_back(r.err_sel_wc);
    elt.push_back(r.err_sel_lt); eh8.push_back(r.err_sel_h8);
  }

  FillHist(h_total,  vt,  et);
  FillHist(h_signal, vs,  es);
  FillHist(h_sel_pd, vpd, epd);
  FillHist(h_sel_wc, vwc, ewc);
  FillHist(h_sel_lt, vlt, elt);
  FillHist(h_sel_h8, vh8, eh8);

  auto Save = [&](TH1D* h, int color, const std::string& label, const std::string& suffix) {
    pfs::DrawUnstacked2({h}, {color}, {label}, out_path + suffix + ".png", true);
  };

  Save(h_total,  kBlack,    title + " - Total",               "_total");
  if (show_signal)
    Save(h_signal, kAzure+1, title + " - Signal",             "_signal");
  Save(h_sel_pd, kRed+1,    title + " - Selected (Pandora)",  "_sel_pd");
  Save(h_sel_wc, kGreen+2,  title + " - Selected (WireCell)", "_sel_wc");
  Save(h_sel_lt, kOrange+7, title + " - Selected (Lantern)",  "_sel_lt");
  Save(h_sel_h8, kViolet+1, title + " - Selected (H8)",       "_sel_h8");

  delete h_total; delete h_signal;
  delete h_sel_pd; delete h_sel_wc; delete h_sel_lt; delete h_sel_h8;
}

// Holds per-variable histograms for one file (normalised by POT)
struct FileVarHists {
  std::map<std::string, TH1D*> h8;    // from vars_h8, filled for all sample types
  std::map<std::string, TH1D*> truth; // from vars_t,  filled for overlay only
  bool is_overlay = false;
};

// Process one file: fill per-variable histograms (weighted, normalised by POT).
// Histograms are cloned from pre-loaded templates (same pattern as HistogramManager::LoadTemplate).
// name_prefix must be unique per file to avoid ROOT name collisions.
FileVarHists ProcessFileVars(const std::string& filepath,
                              const std::map<std::string, TH1D*>& h8_templates,
                              const std::map<std::string, TH1D*>& truth_templates,
                              const std::string& name_prefix) {
  TFile* f_in = nullptr; TTree* t_in = nullptr;
  bool is_overlay = false, load_syst = false;
  LoadTreeFiltered(filepath, f_in, t_in, is_overlay, load_syst);

  t_in->GetEntry(0);
  double exposure = POT;

  FileVarHists result;
  result.is_overlay = is_overlay;

  for (const auto& kv : h8_templates) {
    const std::string& var = kv.first;
    auto* hh = (TH1D*)kv.second->Clone((name_prefix + "_h8_" + var).c_str());
    hh->Reset();
    hh->Sumw2();
    hh->SetDirectory(nullptr);
    result.h8[var] = hh;
  }
  if (is_overlay) {
    for (const auto& kv : truth_templates) {
      const std::string& var = kv.first;
      auto* ht = (TH1D*)kv.second->Clone((name_prefix + "_t_" + var).c_str());
      ht->Reset();
      ht->Sumw2();
      ht->SetDirectory(nullptr);
      result.truth[var] = ht;
    }
  }

  long long n = t_in->GetEntries();
  for (long long i = 0; i < n; ++i) {
    if (i % 1000 == 0)
      std::cout << "  " << i << " / " << n << "\r" << std::flush;
    t_in->GetEntry(i);
    double w = weightSplineTimesTune;
    if (std::isnan(w) || std::isinf(w)) continue;

    if (vars_h8) {
      for (const auto& kv : *vars_h8) {
        auto it = result.h8.find(kv.first);
        if (it != result.h8.end() && std::isfinite(kv.second) && kv.second > -900)
          it->second->Fill(kv.second, w);
      }
    }
    if (is_overlay && vars_t) {
      for (const auto& kv : *vars_t) {
        auto it = result.truth.find(kv.first);
        if (it != result.truth.end() && std::isfinite(kv.second) && kv.second > -900)
          it->second->Fill(kv.second, w);
      }
    }
  }
  std::cout << std::endl;

  for (auto& kv : result.h8)    kv.second->Scale(1.0 / exposure);
  for (auto& kv : result.truth) kv.second->Scale(1.0 / exposure);

  f_in->Close();
  return result;
}

// Draw one variable for one sample type: each run is a separate line
void DrawVarFigure(const std::string& out_path,
                   const std::vector<std::string>& run_labels,
                   const std::vector<TH1D*>& hists) {
  if (hists.empty()) return;
  const std::vector<int> colors = {kBlue+1, kRed+1, kGreen+2, kOrange+7};
  std::vector<int> cols;
  for (int i = 0; i < (int)hists.size(); ++i)
    cols.push_back(colors[i % (int)colors.size()]);
  pfs::DrawUnstacked2(hists, cols, run_labels, out_path + ".png", true);
}

void EventRates() {

  gStyle->SetOptStat(0);

  std::vector<std::string> vars = {"MuonMom","MuonCosTheta"};
  //std::vector<std::string> vars = var_names;

  
  struct RunDef {
    std::string label;
    std::string overlay;
    std::string dirt;
    std::string beam_off;
  };
  
  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/retupled/";
  const std::vector<RunDef> runs = {
    { "Run 4b",
      "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root",
      "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
      "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_beam_off_metapatch_retuple_retuple_hist.root" },
    { "Run 4c",
      "run4c/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4c.root",
      "run4c/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4c.root",
      "run4c/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4c.root" },
    { "Run 4d",
      "run4d/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4d.root",
      "run4d/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4d.root",
      "run4d/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4d.root" },
    { "Run 5",
      "run5/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_5.root",
      "run5/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_5.root",
      "run5/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_5.root" }
  };
  

  // Alternative: pre-retuple files produced directly by Filters/deprecated/Filter.C
  /*
  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/old/";
  const std::vector<RunDef> runs = {
    { "Run 4b",
      "run4b/Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
      "run4b/Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
      "run4b/Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root" },
    { "Run 4c",
      "run4c/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4c.root",
      "run4c/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4c.root",
      "run4c/Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4c.root" },
    { "Run 4d",
      "run4d/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4d.root",
      "run4d/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4d.root",
      "run4d/Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4d.root" },
    { "Run 5",
      "run5/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_5.root",
      "run5/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_5.root",
      "run5/Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_5.root" }
  };
  */

  std::vector<std::string>  run_labels;
  std::vector<RateSummary>  overlay_rates, dirt_rates, ext_rates;

  for (const auto& r : runs) {
    run_labels.push_back(r.label);

    std::cout << "\n=== " << r.label << " overlay ===" << std::endl;
    overlay_rates.push_back(ProcessFile(in_dir + r.overlay));

    std::cout << "=== " << r.label << " dirt ===" << std::endl;
    dirt_rates.push_back(ProcessFile(in_dir + r.dirt));

    std::cout << "=== " << r.label << " beam-off ===" << std::endl;
    ext_rates.push_back(ProcessFile(in_dir + r.beam_off));
  }

  // Print summary table
  std::cout << "\n--- Event rates per POT / trigger ---\n";
  std::cout << std::left << std::setw(10) << "Run"
            << std::setw(14) << "Type"
            << std::setw(14) << "Total"
            << std::setw(14) << "Signal"
            << std::setw(14) << "Sel(PD)"
            << std::setw(14) << "Sel(WC)"
            << std::setw(14) << "Sel(LT)"
            << std::setw(14) << "Sel(H8)" << "\n";
  for (int i = 0; i < (int)runs.size(); ++i) {
    auto Print = [&](const std::string& type, const RateSummary& rs) {
      std::cout << std::left << std::setw(10) << run_labels[i]
                << std::setw(14) << type
                << std::setw(14) << rs.total
                << std::setw(14) << rs.signal
                << std::setw(14) << rs.sel_pd
                << std::setw(14) << rs.sel_wc
                << std::setw(14) << rs.sel_lt
                << std::setw(14) << rs.sel_h8 << "\n";
    };
    Print("Overlay",  overlay_rates[i]);
    Print("Dirt",     dirt_rates[i]);
    Print("Beam-off", ext_rates[i]);
  }

  gSystem->Exec("mkdir -p Plots/EventRates");

  DrawFigure("Overlay: event rates per POT",
             "Events / POT",
             "Plots/EventRates/Overlay",
             run_labels, overlay_rates, /*show_signal=*/true);

  DrawFigure("Dirt: event rates per POT",
             "Events / POT",
             "Plots/EventRates/Dirt",
             run_labels, dirt_rates, /*show_signal=*/false);

  DrawFigure("Beam-off: event rates per trigger",
             "Events / trigger",
             "Plots/EventRates/BeamOff",
             run_labels, ext_rates, /*show_signal=*/false);

  std::cout << "\nPlots saved to Plots/EventRates/" << std::endl;

  // --- Variable distributions ---

  // Load template histograms for variable binning (same pattern as HistogramManager::LoadTemplate)
  std::cout << "\n=== Loading variable templates ===" << std::endl;
  std::map<std::string, TH1D*> h8_templates, truth_templates;
  for (const auto& var : vars) {
    TFile* f_h8 = TFile::Open(("Analysis/" + var + "/rootfiles/BinningTemplate.root").c_str());
    if (f_h8 && !f_h8->IsZombie()) {
      TH1D* h = (TH1D*)f_h8->Get("h_template_All");
      if (h) { h->SetDirectory(0); h->SetName(("h_template_h8_" + var).c_str()); h8_templates[var] = h; }
      f_h8->Close();
    }
    TFile* f_t = TFile::Open(("Analysis/" + var + "/rootfiles/TruthBinningTemplate.root").c_str());
    if (f_t && !f_t->IsZombie()) {
      TH1D* ht = (TH1D*)f_t->Get("h_template_All");
      if (ht) { ht->SetDirectory(0); ht->SetName(("h_template_truth_" + var).c_str()); truth_templates[var] = ht; }
      f_t->Close();
    }
  }

  // Fill histograms for each run and sample type
  std::cout << "\n=== Processing variable histograms ===" << std::endl;
  std::vector<FileVarHists> ov_vhists, dirt_vhists, ext_vhists;
  for (int i = 0; i < (int)runs.size(); ++i) {
    const std::string pfx = "r" + std::to_string(i);
    std::cout << "\n--- " << runs[i].label << " overlay ---" << std::endl;
    ov_vhists.push_back(ProcessFileVars(in_dir + runs[i].overlay,   h8_templates, truth_templates, pfx + "_ov"));
    std::cout << "--- " << runs[i].label << " dirt ---" << std::endl;
    dirt_vhists.push_back(ProcessFileVars(in_dir + runs[i].dirt,    h8_templates, truth_templates, pfx + "_dirt"));
    std::cout << "--- " << runs[i].label << " beam-off ---" << std::endl;
    ext_vhists.push_back(ProcessFileVars(in_dir + runs[i].beam_off, h8_templates, truth_templates, pfx + "_ext"));
  }

  // Helper: collect one histogram per run for a given variable and source map
  auto collect = [&](const std::vector<FileVarHists>& all_runs,
                     const std::string& var, bool use_truth) {
    std::vector<TH1D*> out;
    for (const auto& fh : all_runs) {
      const auto& src = use_truth ? fh.truth : fh.h8;
      auto it = src.find(var);
      if (it != src.end()) out.push_back(it->second);
    }
    return out;
  };

  for (const auto& kv : h8_templates) {
    const std::string& var = kv.first;
    const std::string  dir = "Analysis/" + var + "/Plots/EventRates/";
    gSystem->Exec(("mkdir -p " + dir).c_str());
    DrawVarFigure(dir + "Overlay_h8_" + var, run_labels, collect(ov_vhists,   var, false));
    DrawVarFigure(dir + "Dirt_h8_"    + var, run_labels, collect(dirt_vhists, var, false));
    DrawVarFigure(dir + "BeamOff_h8_" + var, run_labels, collect(ext_vhists,  var, false));
  }
  for (const auto& kv : truth_templates) {
    const std::string& var = kv.first;
    const std::string  dir = "Analysis/" + var + "/Plots/EventRates/";
    gSystem->Exec(("mkdir -p " + dir).c_str());
    DrawVarFigure(dir + "Overlay_truth_" + var, run_labels, collect(ov_vhists, var, true));
  }

  std::cout << "\nVariable plots saved to Analysis/<var>/Plots/EventRates/" << std::endl;
}
