// Re-check of the beam energy chain in the (p,d) simulation, entry by entry.
void beamcheck(TString file, Long64_t nMax = 20000) {
  gSystem->Load("libAtSimulationData.so");
  TFile f(file); auto t = (TTree *)f.Get("cbmsim");
  TClonesArray *mc = nullptr, *pt = nullptr;
  t->SetBranchAddress("MCTrack", &mc); t->SetBranchAddress("AtTpcPoint", &pt);
  const double u = 0.93149401, mA = 16.0147 * u, mN = mA - 6 * 0.000510999; // GeV
  auto T = [](double p, double m) { return 1000 * (std::sqrt(p * p + m * m) - m); };
  double sP0 = 0, sDp0 = 0, sW = 0, sTgA = 0, sTgN = 0, sTwA = 0, sTwN = 0, sZ0 = 0, sZw = 0; long n = 0;
  std::map<int, long> det0, detAfterWin; double maxDp0 = 0;
  for (Long64_t i = 0; i < std::min(nMax, t->GetEntries()); ++i) {
    t->GetEntry(i);
    if (!mc || !pt || pt->GetEntries() < 3) continue;
    AtMCTrack *beam = nullptr;
    for (int k = 0; k < mc->GetEntries(); ++k) { auto *m = (AtMCTrack *)mc->At(k); if (m->GetPdgCode() == 1000060160 && m->GetMotherId() == -1) { beam = m; break; } }
    if (!beam) continue;
    // keep only beam-track points (trackID of the primary beam = 0)
    std::vector<AtMCPoint *> bp;
    for (int k = 0; k < pt->GetEntries(); ++k) { auto *q = (AtMCPoint *)pt->At(k); if (q->GetTrackID() == 0) bp.push_back(q); }
    if (bp.size() < 3) continue;
    double pMC = std::sqrt(pow(beam->GetPx(), 2) + pow(beam->GetPy(), 2) + pow(beam->GetPz(), 2));
    auto pOf = [](AtMCPoint *q) { return std::sqrt(q->GetPx() * q->GetPx() + q->GetPy() * q->GetPy() + q->GetPz() * q->GetPz()); };
    double eW = 0; int firstGas = -1;
    for (size_t k = 0; k < bp.size(); ++k) { if (bp[k]->GetDetectorID() == 4) eW += bp[k]->GetEnergyLoss(); else if (firstGas < 0) firstGas = k; }
    if (firstGas < 0) continue;
    det0[bp[0]->GetDetectorID()]++; detAfterWin[bp[firstGas]->GetDetectorID()]++;
    double p0 = pOf(bp[0]), pg = pOf(bp[firstGas]);
    sP0 += pMC; sDp0 += (p0 - pMC) * 1e6; maxDp0 = std::max(maxDp0, std::fabs(p0 - pMC) * 1e6);
    sW += eW * 1000; sTwA += T(p0, mA); sTwN += T(p0, mN); sTgA += T(pg, mA); sTgN += T(pg, mN);
    sZ0 += bp[0]->GetZ(); sZw += bp[firstGas]->GetZ(); ++n;
  }
  printf("beam entries used: %ld\n", n);
  printf("MCTrack p                       : %.6f GeV/c\n", sP0 / n);
  printf("p(first point) - p(MCTrack)     : mean %.3f keV/c, max |diff| %.3f keV/c  (0 => nothing before the window)\n", sDp0 / n, maxDp0);
  printf("z of first point / first gas pt : %.4f cm / %.5f cm\n", sZ0 / n, sZw / n);
  for (auto &[d, c] : det0) printf("detector of first point: id %d  (%ld)\n", d, c);
  for (auto &[d, c] : detAfterWin) printf("detector of first non-window point: id %d  (%ld)\n", d, c);
  printf("T at window entrance : atomic mass %.3f MeV | nucleus mass %.3f MeV\n", sTwA / n, sTwN / n);
  printf("ELoss in tpc_window  : %.3f MeV\n", sW / n);
  printf("T at gas entrance    : atomic mass %.3f MeV | nucleus mass %.3f MeV\n", sTgA / n, sTgN / n);
  printf("check (nucleus): window entrance - ELoss = %.3f MeV\n", sTwN / n - sW / n);
}
