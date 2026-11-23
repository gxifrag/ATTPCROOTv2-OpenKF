// For each primary of a reaction entry: the kinetic energy Geant4 starts it with,
// T_G4 = sqrt(p^2 + m_TDB^2) - m_TDB (G4PrimaryParticle::Set4Momentum, with E from FairPrimaryGenerator),
// compared with the energy AtTPC2Body intended (EIni stored in the hits) and with T(p, macro mass).
void prodcheck(TString file, double mMacroD, double mMacroRes, double resEx, Long64_t nMax = 40000) {
  gSystem->Load("libAtSimulationData.so");
  TFile f(file); auto t = (TTree *)f.Get("cbmsim");
  TClonesArray *mc = nullptr, *pt = nullptr;
  t->SetBranchAddress("MCTrack", &mc); t->SetBranchAddress("AtTpcPoint", &pt);
  auto T = [](double p, double m) { return std::sqrt(p * p + m * m) - m; }; // MeV
  std::map<int, double> mTDB = {{1000060160, 14914.5336}, {1000060150, 13979.2186}, {1000010020, 1871.0596}, {1000010030, 2809.4140}};
  struct Acc { double dG4intended = 0, dEini = 0, n = 0, nE = 0, maxAbs = 0; };
  std::map<int, Acc> acc;
  for (Long64_t i = 0; i < std::min(nMax, t->GetEntries()); ++i) {
    t->GetEntry(i);
    if (!mc) continue;
    for (int k = 0; k < mc->GetEntries(); ++k) {
      auto *m = (AtMCTrack *)mc->At(k);
      int pdg = m->GetPdgCode();
      if (m->GetMotherId() != -1 || pdg == 1000060160 || !mTDB.count(pdg)) continue; // primaries, not the beam
      double p = 1000 * std::sqrt(pow(m->GetPx(), 2) + pow(m->GetPy(), 2) + pow(m->GetPz(), 2));
      double mMacro = (pdg == 1000060150) ? mMacroRes + resEx : mMacroD;
      double tG4 = T(p, mTDB[pdg]), tInt = T(p, mMacro);
      auto &a = acc[pdg]; a.dG4intended += (tG4 - tInt) * 1000; a.n++; a.maxAbs = std::max(a.maxAbs, std::fabs(tG4 - tInt) * 1000);
      for (int h = 0; h < (pt ? pt->GetEntries() : 0); ++h) { auto *q = (AtMCPoint *)pt->At(h); if (q->GetTrackID() == k) { a.dEini += (tInt - q->GetEIni()) * 1000; a.nE++; break; } }
    }
  }
  for (auto &[pdg, a] : acc)
    printf("pdg %d: %6.0f tracks | T_G4 - T_intended: mean %+8.2f keV (max |.| %.2f) | T_intended - EIni: mean %+7.3f keV (%0.f with hits)\n",
           pdg, a.n, a.dG4intended / a.n, a.maxAbs, a.nE ? a.dEini / a.nE : 0, a.nE);
}
