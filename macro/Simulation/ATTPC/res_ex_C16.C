/// @file res_ex_C16.C
/// @brief Excitation-energy resolution of 16C(p,d)15C or 16C(d,t)15C at one generated Ex, from a
///        simulation and its UKF fit, with the cuts and beam-energy convention of the data analysis.
///
/// Pairs fit and truth by ENTRY INDEX (the fit must be of the ungated reco of the same sim file;
/// the macro refuses to run if the entry counts differ). In each event the fitted track closest in
/// direction to the true primary ejectile (d or t) is taken, i.e. PID is assumed perfect; the
/// data's PID gate is not applied, so this measures the width the kinematics + reconstruction +
/// fit produce, not the gate's acceptance.
///
/// Ex is computed as the analysis does it (a1975_analysis_KF pd/pd_config.h, dt/fit_dt_definitive.C):
/// FIXED beam energy (183.5 MeV for (p,d), 184.0 for (d,t)), fitted ejectile KE and lab polar angle.
/// The data-only empirical corrections (keOff, the kc angle correction, exShift) are NOT applied:
/// they calibrate the real detector. A constant offset does not change the width; the bias is
/// reported separately.
///
/// Cuts (defaults = the analysis):
///   (p,d): chi2/ndf < 5,   KE 2-30 MeV, vertex z 0-500 mm
///   (d,t): chi2/ndf < 4.5, vertex z 25-600 mm, theta_lab >= 15 deg
/// The vertex-z cut is applied to the FITTED vertex, in whatever frame the fit reports; check the
/// printed truth-vs-fit z comparison before trusting it.
///
/// Width: gaussian fitted to the core (iterated: +-2 sigma around the peak), because failed fits
/// make tails that inflate the rms and describe neither.
///
///   root -l -b -q 'res_ex_C16.C("gs_sim.root","gs_ukf.root","pd",0.0)'
/// Prints one line beginning with RESULT for scripts to grep.

void res_ex_C16(TString simFile, TString fitFile, TString reac = "pd", Double_t resEx = 0.0, TString outRoot = "",
                Double_t Ebeam = -1, Double_t chi2Max = -1, Double_t vzLo = -1e9, Double_t vzHi = 1e9,
                Double_t keLo = -1, Double_t keHi = 1e9, Double_t thLoDeg = -1, Double_t thHiDeg = 181,
                // BEAM ENERGY AT THE VERTEX (2026-11-25). catimaE0 <= 0 (default): fixed Ebeam, as the
                // analysis does. catimaE0 > 0: T_beam(z) = CATIMA energy of 16C after z mm of gas from
                // catimaE0 at z = 0 (density catimaRho), evaluated at the FITTED vertex z for the
                // reconstructed Ex (as on data) and at the TRUE z for the true-kinematics Ex.
                Double_t catimaE0 = -1, Double_t catimaRho = 3.308e-5)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gSystem->Load("libAtTools.so");
   const double u = 931.49401;
   const double mb = 16.0147 * u, mr = 15.0105993 * u;
   const bool pd = (reac == "pd");
   if (!pd && reac != "dt") {
      printf("reaction must be pd or dt\n");
      return;
   }
   const double mt = pd ? 1.0078250322 * u : 2.0141017781 * u; // target
   const double me = pd ? 2.0141017781 * u : 3.0160492779 * u; // ejectile
   const int pdgE = pd ? 1000010020 : 1000010030;
   // analysis defaults
   if (Ebeam < 0)
      Ebeam = pd ? 183.5 : 184.0;
   if (chi2Max < 0)
      chi2Max = pd ? 5.0 : 4.5;
   if (vzLo < -1e8 && vzHi > 1e8) {
      vzLo = pd ? 0.0 : 25.0;
      vzHi = pd ? 500.0 : 600.0;
   }
   if (keLo < 0 && pd) {
      keLo = 2.0;
      keHi = 30.0;
   }
   if (thLoDeg < 0)
      thLoDeg = pd ? 0.0 : 15.0;

   TFile *fs = TFile::Open(simFile), *ff = TFile::Open(fitFile);
   if (!fs || fs->IsZombie() || !ff || ff->IsZombie()) {
      printf("cannot open %s or %s\n", simFile.Data(), fitFile.Data());
      return;
   }
   auto ts = (TTree *)fs->Get("cbmsim"), tf = (TTree *)ff->Get("cbmsim");
   if (!ts || !tf)
      return;
   if (ts->GetEntries() < tf->GetEntries()) {
      printf("ENTRY MISMATCH: sim %lld < fit %lld\n", ts->GetEntries(), tf->GetEntries());
      return;
   }
   TClonesArray *mc = nullptr, *te = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   tf->SetBranchAddress("AtTrackingEvent", &te);

   AtTools::AtELossCATIMA beamCat(catimaRho);
   beamCat.SetProjectile(16, 6, 16.0147);
   beamCat.SetMaterial({std::make_tuple(1, 1, 1)});
   const bool useCat = catimaE0 > 0;
   if (useCat)
      printf("  Ebeam(z) from CATIMA: E0 %.3f MeV at z=0, density %.4g g/cm3\n", catimaE0, catimaRho);
   auto tbAt = [&](double z) {
      if (!useCat)
         return (double)Ebeam;
      return z <= 0 ? (double)catimaE0 : beamCat.GetEnergy(catimaE0, std::min(z, 1000.0));
   };
   auto exOfT = [&](double tb, double ke, double th) {
      double E1 = tb + mb, p1 = std::sqrt(tb * tb + 2 * tb * mb);
      double E3 = ke + me, p3 = std::sqrt(E3 * E3 - me * me);
      double E4 = E1 + mt - E3, p4sq = p1 * p1 + p3 * p3 - 2 * p1 * p3 * std::cos(th);
      return std::sqrt(std::max(E4 * E4 - p4sq, 0.0)) - mr;
   };

   auto hEx = new TH1D("hEx", Form("%s Ex = %.3f MeV;E_{x}^{rec} - E_{x}^{gen} [MeV];tracks", reac.Data(), resEx), 400,
                       -8, 8);
   auto hExT = new TH1D("hExT", "Ex from TRUE KE/angle at fixed Ebeam - Ex gen;[MeV];tracks", 400, -8, 8);
   auto hKE = new TH2D("hKE", "KE fit vs true;KE_{true} [MeV];KE_{fit}-KE_{true} [MeV]", 60, 0, 60, 200, -10, 10);
   auto hTh = new TH1D("hTh", "theta fit - true;#Delta#theta [deg];tracks", 200, -10, 10);
   auto hZ = new TH2D("hZ", "vertex z;z_{true} [mm];z_{fit} [mm]", 120, -200, 1100, 120, -200, 1100);
   auto hExTh = new TH2D("hExTh", ";#theta_{lab} [deg];E_{x}^{rec}-E_{x}^{gen} [MeV]", 36, 0, 180, 200, -8, 8);

   long nEvT = 0, nFit = 0, nSel = 0;
   for (Long64_t i = 0; i < tf->GetEntries(); ++i) {
      ts->GetEntry(i);
      tf->GetEntry(i);
      // true primary ejectile
      double keT = -1, thT = 0, phT = 0, zT = 0;
      for (int k = 0; mc && k < mc->GetEntriesFast(); ++k) {
         auto p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != pdgE)
            continue;
         double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
         double pp = std::sqrt(px * px + py * py + pz * pz);
         if (pp <= 0)
            continue;
         keT = std::sqrt(pp * pp + me * me) - me;
         thT = std::acos(pz / pp);
         phT = std::atan2(py, px);
         zT = p->GetStartZ() * 10; // cm -> mm
         break;
      }
      if (keT < 0)
         continue;
      ++nEvT;
      if (!te || te->GetEntriesFast() == 0)
         continue;
      auto e = (AtTrackingEvent *)te->At(0);
      if (!e)
         continue;
      // fitted track closest in polar angle to the true ejectile = perfect PID
      AtFittedTrack *best = nullptr;
      double bestD = 1e9;
      for (auto &ft : e->GetFittedTracks()) {
         if (!ft)
            continue;
         double d = std::fabs(ft->GetKinematics().theta - thT);
         if (d < bestD) {
            bestD = d;
            best = ft.get();
         }
      }
      if (!best)
         continue;
      ++nFit;
      const auto &m = best->GetTrackMetadata();
      double ndf = m ? m->GetNdf() : 0, c2 = m ? m->GetChi2() : 0;
      double ke = best->GetKinematics().kineticEnergy, th = best->GetKinematics().theta;
      double zF = best->GetVertex(0).Z(), thDeg = th * TMath::RadToDeg();
      hZ->Fill(zT, zF);
      if (!(ndf > 0 && c2 / ndf < chi2Max))
         continue;
      if (zF < vzLo || zF > vzHi || ke < keLo || ke > keHi || thDeg < thLoDeg || thDeg > thHiDeg)
         continue;
      ++nSel;
      double ex = exOfT(tbAt(zF), ke, th);
      hEx->Fill(ex - resEx);
      hExT->Fill(exOfT(tbAt(zT), keT, thT) - resEx);
      hKE->Fill(keT, ke - keT);
      hTh->Fill((th - thT) * TMath::RadToDeg());
      hExTh->Fill(thDeg, ex - resEx);
   }

   // iterated core gaussian
   auto core = [](TH1D *h, double &mu, double &sg, double &emu, double &esg) {
      mu = sg = emu = esg = -1;
      if (h->GetEntries() < 30)
         return;
      // start from the FWHM around the peak, not the rms: the tails of failed fits inflate the rms
      // several-fold and a fit seeded with it wanders off the core
      int bm = h->GetMaximumBin();
      double half = 0.5 * h->GetBinContent(bm);
      int bl = bm, br = bm;
      while (bl > 1 && h->GetBinContent(bl) > half)
         --bl;
      while (br < h->GetNbinsX() && h->GetBinContent(br) > half)
         ++br;
      double c = h->GetBinCenter(bm), w = std::max((h->GetBinCenter(br) - h->GetBinCenter(bl)) / 2.3548, 0.05);
      TF1 g("g", "gaus", c - 2 * w, c + 2 * w);
      for (int it = 0; it < 4; ++it) {
         g.SetRange(c - 2 * w, c + 2 * w);
         h->Fit(&g, "QNR");
         c = g.GetParameter(1);
         w = std::max(std::fabs(g.GetParameter(2)), 0.02);
      }
      mu = c;
      sg = w;
      emu = g.GetParError(1);
      esg = g.GetParError(2);
   };
   double mu, sg, emu, esg, muT, sgT, e1, e2;
   core(hEx, mu, sg, emu, esg);
   core(hExT, muT, sgT, e1, e2);
   double fwhm = 2.3548 * sg;
   // fraction of selected tracks within +-3 sigma of the core (how much is tails)
   double inCore = sg > 0 ? hEx->Integral(hEx->FindBin(mu - 3 * sg), hEx->FindBin(mu + 3 * sg)) : 0;

   printf("\n  %s Ex=%.3f: events with true ejectile %ld, with a fitted track %ld, selected %ld\n", reac.Data(),
          resEx, nEvT, nFit, nSel);
   printf("  KE: mean fit-true %.3f MeV (rms %.3f)   theta: mean fit-true %.3f deg (rms %.3f)\n",
          hKE->ProjectionY()->GetMean(), hKE->ProjectionY()->GetRMS(), hTh->GetMean(), hTh->GetRMS());
   printf("  vertex z: mean fit - true %.1f mm, correlation %.3f\n",
          hZ->GetMean(2) - hZ->GetMean(1), hZ->GetCorrelationFactor());
   printf("  Ex from TRUE kinematics at fixed Ebeam: bias %.3f, sigma %.3f MeV (beam-energy spread alone)\n", muT,
          sgT);
   printf("RESULT %s Ex %.3f  sel %ld  bias %.4f +- %.4f  sigma %.4f +- %.4f  FWHM %.4f  core %.3f\n", reac.Data(),
          resEx, nSel, mu, emu, sg, esg, fwhm, nSel ? inCore / nSel : 0.);

   if (outRoot.Length()) {
      TFile fo(outRoot, "RECREATE");
      hEx->Write();
      hExT->Write();
      hKE->Write();
      hTh->Write();
      hZ->Write();
      hExTh->Write();
      TParameter<double>("resEx", resEx).Write();
      TParameter<double>("sigma", sg).Write();
      TParameter<double>("sigmaErr", esg).Write();
      TParameter<double>("bias", mu).Write();
      TParameter<double>("Ebeam", Ebeam).Write();
   }
}
