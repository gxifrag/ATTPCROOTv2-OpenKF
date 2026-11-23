/// @file tb_edges.C
/// @brief Measure the micromegas and window time buckets (TB) from reconstructed hits, and compare the
///        drift calibration they imply with the one in the parameter file and with a given length
///        (e.g. the 940 mm used in the (p,d) analysis).
///
/// Written 2026-11-18. See the simulation log, section "Drift: time buckets, drift velocity and ZPadPlane".
///
/// ============================================================================================
/// THE IDEA
/// ============================================================================================
/// Every hit stores the TB of its pulse (AtHit::GetTimeStampCorr(), or GetTimeStamp() if the PSA does
/// not fill the corrected one). The TB does NOT depend on the drift velocity or the TBs written in the
/// parameter file: those are only used later to turn TB into z. So the TB distribution of all hits,
/// weighted by charge, shows the detector directly:
///   - charge created right above the pad plane arrives first  -> rising edge = MICROMEGAS TB
///   - charge created by the beam at the entrance window arrives last -> falling edge = WINDOW TB
/// With the physical length L between them:
///   v_drift = L / ((TB_window - TB_micromegas) * t_TB),   t_TB = 160 ns (SamplingRate 6)
///
/// Two edge definitions are given (they should agree within a few TB):
///   (1) maximum of |d ln(Sum Q) / dTB|   (as in the "Timebucket choice" slide)
///   (2) crossing of half the plateau height (plateau = median of TB 100-400)
/// The window edge is also measured with the central pads only (r < beamRadius), where the beam is.
///
/// ============================================================================================
/// INPUT
/// ============================================================================================
/// ROOT files written by the unpack + PSA stage (e.g. unpack_a1975.C), with the branch "AtEventH".
/// No tracking or fit needed. Use several runs spread over the H2 period (runs 104-189).
/// The first argument is a pattern understood by TChain, e.g. "/path/run_01*.root".
///
/// ============================================================================================
/// HOW TO RUN
/// ============================================================================================
///   root -l 'tb_edges.C("/path/to/run_01*.root")'
///   root -l 'tb_edges.C("/path/run_0131.root", 1.30, 480, 1000., 940.)'    // all arguments
///
/// OUTPUT: tb_edges.png / tb_edges.pdf (next to this macro) and a summary printed on screen.

void tb_edges(TString files = "run_0131.root", // TChain pattern of unpack+PSA files (branch AtEventH)
              Double_t vdPar = 1.30,           // cm/us, drift velocity in the parameter file used
              Int_t tbEntrancePar = 480,       // TBEntrance in the parameter file used
              Double_t zPadPlanePar = 1000.,   // mm, ZPadPlane in the parameter file used
              Double_t lengthTest = 940.,      // mm, a length to test (e.g. the 940 of the (p,d) analysis)
              Double_t lengthPhys = 1000.,     // mm, physical distance window -> micromegas
              Double_t beamRadius = 30.,       // mm, radius of the "beam pads" (pads the beam drifts onto), used ONLY for the window edge
              Long64_t maxEvents = 200000,
              TString outName = "")        // output file name without extension; "" = tb_edges next to this macro
{
   gSystem->Load("libAtData.so");
   gStyle->SetOptStat(0);
   const Double_t tTB = 0.160;        // us per time bucket (SamplingRate 6 -> 160 ns)
   const Double_t mmPerTB = tTB * vdPar * 10.; // mm/TB with the par drift velocity

   TChain ch("cbmsim");
   if (ch.Add(files) == 0) {
      printf("\033[1;31mno file matches %s\033[0m\n", files.Data());
      return;
   }
   TClonesArray *evArr = nullptr;
   if (ch.SetBranchAddress("AtEventH", &evArr) < 0) {
      printf("\033[1;31mno branch AtEventH in %s (run the unpack + PSA stage first)\033[0m\n", files.Data());
      return;
   }

   // Sum of hit charge per TB: all pads, and central (beam) pads only. Bins of 1 TB centred on integers.
   TH1D *hAll = new TH1D("hAll", ";time bucket;#Sigma charge", 532, -10.5, 521.5);
   TH1D *hBeam = new TH1D("hBeam", ";time bucket;#Sigma charge", 532, -10.5, 521.5);
   hAll->SetDirectory(nullptr);
   hBeam->SetDirectory(nullptr);
   Long64_t nEv = std::min(maxEvents, ch.GetEntries()), nHits = 0, nCorr = 0;
   for (Long64_t i = 0; i < nEv; ++i) {
      ch.GetEntry(i);
      if (!evArr || evArr->GetEntriesFast() == 0)
         continue;
      auto *ev = (AtEvent *)evArr->At(0);
      if (!ev)
         continue;
      for (auto &h : ev->GetHits()) {
         // Prefer the corrected (centre-of-gravity) TB; MultiFit-type PSAs only fill the integer one.
         Double_t tb = h->GetTimeStampCorr();
         if (tb != 0) ++nCorr; else tb = h->GetTimeStamp();
         const Double_t q = h->GetCharge();
         hAll->Fill(tb, q);
         const auto &pos = h->GetPosition();
         if (std::hypot(pos.X(), pos.Y()) < beamRadius)
            hBeam->Fill(tb, q);
         ++nHits;
      }
   }
   printf("\n%lld events, %lld hits (%s TB used)\n", nEv, nHits,
          nCorr > nHits / 2 ? "corrected" : "integer");
   if (nHits == 0) return;

   // ---- edges ---------------------------------------------------------------------------------
   auto edges = [&](TH1D *h, Double_t &upDeriv, Double_t &dnDeriv, Double_t &upHalf, Double_t &dnHalf) {
      std::vector<double> c, y;
      for (int b = 1; b <= h->GetNbinsX(); ++b) { c.push_back(h->GetBinCenter(b)); y.push_back(h->GetBinContent(b)); }
      std::vector<double> plat;
      for (size_t k = 0; k < c.size(); ++k) if (c[k] > 100 && c[k] < 400) plat.push_back(y[k]);
      std::nth_element(plat.begin(), plat.begin() + plat.size() / 2, plat.end());
      const double p = plat[plat.size() / 2];
      // (1) derivative of ln(sum), 3-TB running average to tame fluctuations
      double best = -1e9, worst = 1e9; upDeriv = dnDeriv = -999;
      for (size_t k = 2; k + 2 < c.size(); ++k) {
         auto L = [&](size_t j) { return std::log(std::max(y[j], 1e-3 * p)); };
         double d = (L(k + 2) + L(k + 1) - L(k - 1) - L(k - 2)) / 6.0;
         if (c[k] > -5 && c[k] < 150 && d > best) { best = d; upDeriv = c[k]; }
         if (c[k] > 300 && c[k] < 520 && d < worst) { worst = d; dnDeriv = c[k]; }
      }
      // (2) half-plateau crossings: first bin above p/2 from the left, last above p/2 from the right
      upHalf = dnHalf = -999;
      for (size_t k = 0; k < c.size(); ++k) if (c[k] > -5 && y[k] > 0.5 * p) { upHalf = c[k]; break; }
      for (size_t k = c.size(); k-- > 0;) if (c[k] < 520 && y[k] > 0.5 * p) { dnHalf = c[k]; break; }
      return p;
   };
   Double_t aUpD, aDnD, aUpH, aDnH, bUpD, bDnD, bUpH, bDnH;
   edges(hAll, aUpD, aDnD, aUpH, aDnH);
   // Central pads may be empty (e.g. a beam hole, or a simulation where the beam leaves no hits):
   // then the window edge is taken from all pads, and this is said explicitly.
   const bool beamOK = hBeam->Integral() > 0.01 * hAll->Integral();
   if (beamOK)
      edges(hBeam, bUpD, bDnD, bUpH, bDnH);
   else {
      bUpD = bDnD = bUpH = bDnH = -999;
      printf("\033[1;33mWARNING: almost no charge in the central pads (r < %.0f mm): window edge taken from ALL pads.\033[0m\n",
             beamRadius);
   }
   // A sharp window edge needs the beam (or tracks starting at the window). Without a plateau the edges
   // are not meaningful: warn if the two definitions disagree by more than 5 TB.
   if (std::fabs(aUpD - aUpH) > 5 || std::fabs(aDnD - aDnH) > 5)
      printf("\033[1;33mWARNING: the two edge definitions disagree by > 5 TB -- no clean plateau; do not trust "
             "these edges (more events, or events with beam hits, are needed).\033[0m\n");

   const Double_t tbMM = aUpD, tbWin = beamOK ? bDnD : aDnD; // adopted: derivative method
   const Double_t span = tbWin - tbMM;
   const Double_t tbMMpar = tbEntrancePar - zPadPlanePar / mmPerTB; // where the par file puts the micromegas
   printf("\n================ EDGES (time buckets) ================\n");
   printf("                         max |dln/dTB|   half plateau\n");
   printf("  micromegas, all pads      %7.1f        %7.1f\n", aUpD, aUpH);
   printf("  window,     all pads      %7.1f        %7.1f\n", aDnD, aDnH);
   printf("  window,     beam pads     %7.1f        %7.1f\n", bDnD, bDnH);
   printf("  adopted: micromegas TB %.1f, window TB %.1f  ->  span %.1f TB\n", tbMM, tbWin, span);
   printf("\n================ WHAT THEY IMPLY ================\n");
   printf("  (A) physical %.0f mm      -> v_drift = %.4f cm/us\n", lengthPhys, lengthPhys / (span * tTB) / 10.);
   printf("  (B) par v = %.4f cm/us  -> L = %.1f mm\n", vdPar, span * mmPerTB);
   printf("      to get L = %.0f mm at the par v, the span would be %.1f TB\n", lengthTest, lengthTest / mmPerTB);
   printf("\n================ THE PARAMETER FILE ================\n");
   printf("  v %.4f cm/us, TBEntrance %d, ZPadPlane %.1f mm -> assumes micromegas at TB %.1f\n", vdPar, tbEntrancePar,
          zPadPlanePar, tbMMpar);
   printf("  with the measured TBs, the reconstruction puts the window at z_lab = %.1f mm and the micromegas at %.1f mm\n",
          (tbEntrancePar - tbWin) * mmPerTB, (tbEntrancePar - tbMM) * mmPerTB);
   printf("  (true: 0 and %.0f mm). Part of the chamber with reconstructed z_lab >= 0: %.1f mm\n", lengthPhys,
          lengthPhys * (tbEntrancePar - tbMM) / span);

   // ---- picture ---------------------------------------------------------------------------------
   TCanvas *c = new TCanvas("cTB", "time-bucket edges", 1300, 520);
   c->Divide(2, 1);
   c->cd(1);
   gPad->SetLogy();
   hAll->SetLineColor(kAzure + 2);
   hAll->SetTitle("#Sigma hit charge vs TB;time bucket;#Sigma charge");
   hAll->Draw("hist");
   hBeam->SetLineColor(kOrange + 7);
   hBeam->Draw("hist same");
   auto vline = [&](Double_t x, int col, int style, const char *lab, TLegend *lg) {
      TLine *l = new TLine(x, hAll->GetMinimum() > 0 ? hAll->GetMinimum() : 1, x, hAll->GetMaximum());
      l->SetLineColor(col); l->SetLineStyle(style); l->SetLineWidth(2); l->Draw();
      if (lg) lg->AddEntry(l, lab, "l");
   };
   hAll->SetMinimum(std::max(1.0, 1e-3 * hAll->GetMaximum()));
   TLegend *lg = new TLegend(0.30, 0.15, 0.88, 0.45);
   lg->SetTextSize(0.03);
   lg->AddEntry(hAll, "all pads", "l");
   lg->AddEntry(hBeam, Form("beam pads (r < %.0f mm)", beamRadius), "l");
   vline(tbMM, kGreen + 2, 1, Form("measured micromegas TB %.1f", tbMM), lg);
   vline(tbWin, kGreen + 2, 1, Form("measured window TB %.1f", tbWin), nullptr);
   vline(tbEntrancePar, kRed + 1, 2, Form("par: TBEntrance %d, micromegas %.1f", tbEntrancePar, tbMMpar), lg);
   vline(tbMMpar, kRed + 1, 2, "", nullptr);
   vline(tbWin - lengthTest / mmPerTB, kMagenta + 1, 3,
         Form("micromegas for L = %.0f mm at %.2f cm/#mus: TB %.1f", lengthTest, vdPar, tbWin - lengthTest / mmPerTB), lg);
   lg->Draw();
   c->cd(2);
   TH1D *hd = (TH1D *)hAll->Clone("hd");
   hd->Reset();
   for (int b = 3; b <= hAll->GetNbinsX() - 2; ++b) {
      auto L = [&](int j) { return std::log(std::max(hAll->GetBinContent(j), 1e-3 * hAll->GetMaximum())); };
      hd->SetBinContent(b, (L(b + 2) + L(b + 1) - L(b - 1) - L(b - 2)) / 6.0);
   }
   hd->SetTitle("d ln(#Sigma charge) / dTB;time bucket;d ln#Sigma / dTB");
   hd->SetLineColor(kRed + 1);
   hd->SetMinimum(hd->GetMinimum() * 1.1);
   hd->Draw("hist");
   TString out = outName.IsNull() ? TString(gSystem->DirName(gInterpreter->GetCurrentMacroName())) + "/tb_edges" : outName;
   c->SaveAs(out + ".png");
   c->SaveAs(out + ".pdf");
}
