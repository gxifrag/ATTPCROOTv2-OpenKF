/// @file plot_states_C16pd.C for H2 @300Torr
/// @brief The four generated 16C(p,d)15C states: kinematics, and what the beam energy loss does.
///
/// Everything here is MC TRUTH: no digitisation, no reconstruction, no resolution. Any width
/// seen is generated physics plus the analysis convention.
///
/// The effect is the beam losing energy before it reacts. The vertex is uniform along the gas,
/// the 16C enters at 184 MeV and loses about 10.5 MeV over the 940 mm active length (LISE,
/// H2 300 torr), i.e. ~10.49 MeV per 940 mm. Ex reconstructed with a constant beam energy
/// therefore drifts with the vertex z.
///
/// The bottom row shows what correcting it costs: only the vertex position. Replacing the
/// constant Ebeam with Ebeam - dEdrift*z/940 should collapse the width and put every state on
/// the energy it was generated at.
///
/// dEdrift is NOT assumed: a scan on the g.s. finds the value that minimises the Ex width, to
/// be compared with the LISE number (10.49/940 mm). A mismatch means the z origin, the
/// beam energy or maxELoss in the generator are not what we think.
///
/// The correction uses the TRUE vertex z. In data it would use the reconstructed one.
///
/// Input: all chunks <dir>/<tag>_c*_s*_sim.root (output of run_states_par.sh).
///
///   root -b -q 'plot_states_C16pd.C("/path/data_184/")'


void plot_states_C16pd(TString dir = "./data", Double_t Ebeam = 184.0,
                       Double_t dEdrift = 10.49, Double_t Lgas = 940.0,
                       TString tag = "184")
{
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   std::cout << "here directory: " << here << " " << std::endl;

   const double u = 931.49401;
   const double mb = 16.0147 * u, mt = 1.0078250322 * u, md = 2.0141017781 * u, mr = 15.0105993 * u;

   const int NS = 4;
   const char *TAGS[NS] = {"gs", "ex1", "ex2", "ex3"};
   const double EXT[NS] = {0.0, 0.740, 3.103, 4.657};
   const char *NAME[NS] = {"g.s.", "0.740", "3.103", "4.657"};
   const int COL[NS] = {kBlack, kAzure + 2, kRed + 1, kGreen + 3};

   auto exOf = [&](double Eb, double ke, double th) {
      double E1 = Eb + mb, p1 = std::sqrt(Eb * Eb + 2 * Eb * mb);
      double E3 = ke + md, p3 = std::sqrt(E3 * E3 - md * md);
      double E4 = E1 + mt - E3, p4sq = p1 * p1 + p3 * p3 - 2 * p1 * p3 * std::cos(th);
      return std::sqrt(std::max(E4 * E4 - p4sq, 0.0)) - mr;
   };

   TH2D *hkt[NS];
   TH1D *hRaw[NS], *hCor[NS];
   TH2D *hz = new TH2D("hzgs", "g.s.: reconstructed E_{x} vs vertex z;vertex z [mm];E_{x} [MeV]", 40, 0, 1000, 60,
                       -1.5, 3.0);
   for (int i = 0; i < NS; ++i) {
      hkt[i] = new TH2D(TString::Format("hkt%d", i), "", 90, 0, 45, 120, 0, 60);
      hRaw[i] = new TH1D(TString::Format("hr%d", i), "", 220, -2, 9);
      hCor[i] = new TH1D(TString::Format("hc%d", i), "", 220, -2, 9);
   }

   std::vector<double> gsKE, gsTh, gsZ; // g.s. events kept for the dEdrift scan
   std::vector<double> sKE[NS], sTh[NS], sZ[NS];   // todos los estados, para el segundo canvas
   
   printf("\n  Ebeam = %.2f MeV, dEdrift = %.2f MeV/940 mm\n", Ebeam, dEdrift);
   printf("  state   generated   reconstructed at Ebeam     corrected with vertex z      events\n");
   for (int i = 0; i < NS; ++i) {
      TChain *t = new TChain("cbmsim");
      //TString pat = dir + "/" + TAGS[i] + "_c*_s*_sim.root";
      TString pat = dir + "/" + TAGS[i] + "_s*_sim.root";
      int nf = t->Add(pat);
      if (nf == 0) {
         printf("  %-5s MISSING %s\n", TAGS[i], pat.Data());
         continue;
      }
      TClonesArray *mc = nullptr;
      t->SetBranchAddress("MCTrack", &mc);
      Long64_t nWithD = 0;
      for (Long64_t e = 0; e < t->GetEntries(); ++e) {
         t->GetEntry(e);
         bool has = false;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020)
               continue; // deuteron
            double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
            double pp = std::sqrt(px * px + py * py + pz * pz);
            if (pp <= 0)
               continue;
            has = true;
            double ke = std::sqrt(pp * pp + md * md) - md, th = std::acos(pz / pp), z = p->GetStartZ() * 10;
            hkt[i]->Fill(th * TMath::RadToDeg(), ke);
            double exRaw = exOf(Ebeam, ke, th);
            hRaw[i]->Fill(exRaw);
            hCor[i]->Fill(exOf(Ebeam - dEdrift * z / Lgas, ke, th));
            sKE[i].push_back(ke); sTh[i].push_back(th); sZ[i].push_back(z);
            if (i == 0) {
               hz->Fill(z, exRaw);
               gsKE.push_back(ke);
               gsTh.push_back(th);
               gsZ.push_back(z);
            }
         }
         if (has)
            ++nWithD;
      }
      printf("  %-5s   %7.3f     %6.2f +- %.3f              %6.2f +- %.3f      %lld / %lld\n", NAME[i], EXT[i],
             hRaw[i]->GetMean(), hRaw[i]->GetRMS(), hCor[i]->GetMean(), hCor[i]->GetRMS(), nWithD, t->GetEntries());
      delete t;
   }

   // ---- scan of dEdrift on the g.s.: the value that minimises the width
   if (!gsZ.empty()) {
      double zmin = *std::min_element(gsZ.begin(), gsZ.end()), zmax = *std::max_element(gsZ.begin(), gsZ.end());
      printf("\n  g.s. vertex z range: %.1f - %.1f mm\n", zmin, zmax);
      double bestD = 0, bestRMS = 1e9;
      for (double d = 0; d <= 20.0; d += 0.1) {
         TH1D h("hs", "", 440, -2, 9);
         for (size_t k = 0; k < gsZ.size(); ++k)
            h.Fill(exOf(Ebeam - d * gsZ[k] / Lgas, gsKE[k], gsTh[k]));
         if (h.GetRMS() < bestRMS) {
            bestRMS = h.GetRMS();
            bestD = d;
         }
      }
      printf("  best dEdrift (min g.s. width) = %.2f MeV/%.0f mm  (RMS %.3f MeV)   LISE: %.2f\n\n",
       bestD, Lgas, bestRMS, dEdrift);
   }

   TCanvas *c = new TCanvas("cpd", "16C(p,d)15C states", 1500, 950);
   c->Divide(2, 2);

   // ---- 1. the kinematic loci
   c->cd(1);
   auto *fr = new TH2D("frk", "deuteron kinematics (MC truth);#theta_{lab} [deg];KE [MeV]", 1, 0, 45, 1, 0, 60);
   fr->Draw();
   auto *lg = new TLegend(0.62, 0.62, 0.89, 0.88);
   for (int i = 0; i < NS; ++i) {
      hkt[i]->SetMarkerStyle(1);
      hkt[i]->SetMarkerColor(COL[i]);
      hkt[i]->Draw("same");
      auto *m = new TMarker(0, 0, 20);
      m->SetMarkerColor(COL[i]);
      lg->AddEntry(m, TString::Format("%s MeV", NAME[i]), "p");
   }
   lg->Draw();

   // ---- 2. Ex as the analysis reconstructs it
   c->cd(2);
   double mx = 0;
   for (int i = 0; i < NS; ++i)
      mx = std::max(mx, hRaw[i]->GetMaximum());
   auto *f2 = new TH1D("fr2", TString::Format("E_{x} at constant %.0f MeV (no vertex correction);"
                                              "E_{x} [MeV];deuterons",
                                              Ebeam),
                       1, -2, 9);
   f2->SetMaximum(mx * 1.25);
   f2->Draw();
   for (int i = 0; i < NS; ++i) {
      hRaw[i]->SetLineColor(COL[i]);
      hRaw[i]->SetLineWidth(2);
      hRaw[i]->Draw("hist same");
      auto *l = new TLine(EXT[i], 0, EXT[i], mx * 1.25);
      l->SetLineColor(COL[i]);
      l->SetLineStyle(2);
      l->Draw();
   }
   lg->Draw();

   // ---- 3. the same, with the beam loss corrected by the vertex
   c->cd(3);
   mx = 0;
   for (int i = 0; i < NS; ++i)
      mx = std::max(mx, hCor[i]->GetMaximum());
  auto *f3 = new TH1D("fr3",
                    TString::Format("E_{x} corrected with vertex z (dE = %.2f MeV over %.0f mm);E_{x} [MeV];deuterons",
                                    dEdrift, Lgas),
                    1, -2, 9);
   f3->SetMaximum(mx * 1.25);
   f3->Draw();
   for (int i = 0; i < NS; ++i) {
      hCor[i]->SetLineColor(COL[i]);
      hCor[i]->SetLineWidth(2);
      hCor[i]->Draw("hist same");
      auto *l = new TLine(EXT[i], 0, EXT[i], mx * 1.25);
      l->SetLineColor(COL[i]);
      l->SetLineStyle(2);
      l->Draw();
   }
   lg->Draw();

   // ---- 4. where the width comes from
   c->cd(4);
   gPad->SetLogz();
   hz->Draw("colz");
   auto *pf = hz->ProfileX("pfz");
   pf->SetLineColor(kRed + 1);
   pf->SetLineWidth(3);
   pf->SetMarkerColor(kRed + 1);
   pf->SetMarkerStyle(20);
   pf->Draw("same");

   TString png = here + "/plots/pd_states" + tag + ".png";
   gSystem->mkdir(here + "/plots", kTRUE);
   c->SaveAs(png);
   printf("  wrote %s\n\n", png.Data());

      // ---- segundo canvas: ¿es válida la misma corrección para los cuatro estados?
   // RMS exacto (sin binning) y solo con eventos a < 3 MeV del Ex generado, para que un evento
   // cinemáticamente imposible (Ex = -m_r) no estropee la anchura.
   auto rmsAt = [&](int i, double d) {
      double s = 0, s2 = 0; long n = 0;
      for (size_t k = 0; k < sZ[i].size(); ++k) {
         double x = exOf(Ebeam - d * sZ[i][k] / Lgas, sKE[i][k], sTh[i][k]);
         if (std::fabs(x - EXT[i]) > 3.0) continue;
         s += x; s2 += x * x; ++n;
      }
      if (n < 2) return 0.0;
      double m = s / n; return std::sqrt(std::max(s2 / n - m * m, 0.0));
   };

   TCanvas *c2 = new TCanvas("cpd_z", "Ex vs vertex z, per state", 1500, 950);
   c2->Divide(2, 2);
   printf("  state  RMS raw  RMS LISE  best d  RMS best  slope raw  slope corr   [slopes in MeV/m]\n");
   for (int i = 0; i < NS; ++i) {
      if (sZ[i].empty()) continue;
      double lo = EXT[i] - 2, hi = EXT[i] + 2;
      auto *h2  = new TH2D(Form("hzs%d", i), Form("%s MeV: E_{x} vs vertex z;vertex z [mm];E_{x} [MeV]", NAME[i]),
                           40, 0, 1000, 80, lo, hi);
      auto *h2c = new TH2D(Form("hzc%d", i), "", 40, 0, 1000, 80, lo, hi);
      for (size_t k = 0; k < sZ[i].size(); ++k) {
         h2->Fill(sZ[i][k], exOf(Ebeam, sKE[i][k], sTh[i][k]));
         h2c->Fill(sZ[i][k], exOf(Ebeam - dEdrift * sZ[i][k] / Lgas, sKE[i][k], sTh[i][k]));
      }
      double bestD = 0, bestR = 1e9;
      for (double d = 0; d <= 20.0; d += 0.1) {
         double r = rmsAt(i, d);
         if (r > 0 && r < bestR) { bestR = r; bestD = d; }
      }
      c2->cd(i + 1);
      gPad->SetLogz();
      h2->Draw("colz");
      auto *pr = h2->ProfileX(Form("prz%d", i));
      auto *pc = h2c->ProfileX(Form("pcz%d", i));
      pr->SetLineColor(kBlack);   pr->SetLineWidth(3); pr->SetMarkerStyle(20); pr->SetMarkerColor(kBlack);
      pc->SetLineColor(kRed + 1); pc->SetLineWidth(3); pc->SetMarkerStyle(21); pc->SetMarkerColor(kRed + 1);
      pr->Fit("pol1", "Q0"); pc->Fit("pol1", "Q0");
      double sRaw  = pr->GetFunction("pol1")->GetParameter(1) * 1000;
      double sCorr = pc->GetFunction("pol1")->GetParameter(1) * 1000;
      pr->Draw("same"); pc->Draw("same");
      auto *l = new TLine(0, EXT[i], 1000, EXT[i]);
      l->SetLineStyle(2); l->SetLineColor(kBlue + 1); l->Draw();
      printf("  %-5s  %7.3f  %8.3f  %6.2f  %8.3f  %9.3f  %10.3f\n", NAME[i], rmsAt(i, 0), rmsAt(i, dEdrift),
             bestD, bestR, sRaw, sCorr);
   }
   c2->SaveAs(here + "/plots/pd_states_vs_z" + tag + ".png");
}