// Statistical comparison of the primary deuteron in two simulations (before / after the fix).
// (1) velocity from the first step: beta = dLength / (c dTime), compared with p/E for the real
//     (atomic) deuteron mass at the generated momentum -> tells which mass Geant4 really tracks.
// (2) energy deposited per cm in the first 10 cm of gas, in bins of generated momentum.
void fill(TString file, TH1D *hBeta, TProfile *pDedx) {
  gSystem->Load("libAtSimulationData.so");
  TFile F(file); auto *t = (TTree *)F.Get("cbmsim");
  TClonesArray *mc = nullptr, *pt = nullptr;
  t->SetBranchAddress("MCTrack", &mc); t->SetBranchAddress("AtTpcPoint", &pt);
  const double mD = 2.0141017781 * 931.49401, c = 29.9792458; // MeV, cm/ns
  for (Long64_t i = 1; i < t->GetEntries(); i += 2) {
    t->GetEntry(i);
    int id = -1; double p = 0;
    for (int k = 0; k < mc->GetEntries(); ++k) { auto *m = (AtMCTrack *)mc->At(k); if (m->GetPdgCode() == 1000010020 && m->GetMotherId() == -1) { id = k; p = 1000 * std::sqrt(pow(m->GetPx(),2)+pow(m->GetPy(),2)+pow(m->GetPz(),2)); } }
    if (id < 0) continue;
    std::vector<AtMCPoint *> q;
    for (int j = 0; j < pt->GetEntries(); ++j) { auto *x = (AtMCPoint *)pt->At(j); if (x->GetTrackID() == id) q.push_back(x); }
    if (q.size() < 3) continue;
    double dL = q[1]->GetLength() - q[0]->GetLength(), dT = q[1]->GetTime() - q[0]->GetTime();
    if (dT > 0 && dL > 0) hBeta->Fill((dL / (c * dT)) / (p / std::sqrt(p * p + mD * mD)) - 1);
    double e = 0, L = 0;
    for (auto *x : q) { if (x->GetLength() > 10) break; e += x->GetEnergyLoss() * 1e6; L = x->GetLength(); } // keV
    if (L > 9.5) pDedx->Fill(p, e / L);
  }
}
void compare2(TString fb, TString fa) {
  auto *bB = new TH1D("bB", "", 400, -0.01, 0.01), *bA = new TH1D("bA", "", 400, -0.01, 0.01);
  auto *dB = new TProfile("dB", "", 20, 100, 400), *dA = new TProfile("dA", "", 20, 100, 400);
  fill(fb, bB, dB); fill(fa, bA, dA);
  printf("(1) beta(Geant4) / beta(real deuteron mass) - 1:\n");
  printf("    BEFORE: mean %+.5f %% +- %.5f %%  (n = %.0f)\n", 100 * bB->GetMean(), 100 * bB->GetMeanError(), bB->GetEntries());
  printf("    AFTER : mean %+.5f %% +- %.5f %%  (n = %.0f)\n", 100 * bA->GetMean(), 100 * bA->GetMeanError(), bA->GetEntries());
  printf("(2) energy deposited per cm, first 10 cm, AFTER / BEFORE, by generated momentum:\n");
  double sw = 0, swr = 0;
  for (int b = 1; b <= dB->GetNbinsX(); ++b) {
    if (dB->GetBinEntries(b) < 50 || dA->GetBinEntries(b) < 50) continue;
    double r = dA->GetBinContent(b) / dB->GetBinContent(b);
    double er = r * std::hypot(dA->GetBinError(b) / dA->GetBinContent(b), dB->GetBinError(b) / dB->GetBinContent(b));
    printf("    p %3.0f-%3.0f MeV/c: before %6.2f keV/cm, after %6.2f, ratio %.4f +- %.4f\n", dB->GetBinLowEdge(b), dB->GetBinLowEdge(b) + dB->GetBinWidth(b), dB->GetBinContent(b), dA->GetBinContent(b), r, er);
    sw += 1 / (er * er); swr += r / (er * er);
  }
  printf("    weighted mean ratio: %.4f +- %.4f\n", swr / sw, 1 / std::sqrt(sw));
}
