// Beam (even entries, track 0): compare quantities that the beam charge state could change.
// Momenta are compared directly (no mass assumption); T is shown with the atomic mass for reading.
struct Acc { double s = 0, s2 = 0; long n = 0; void add(double x) { s += x; s2 += x * x; ++n; }
  double m() const { return s / n; } double e() const { double mm = m(); return std::sqrt(std::max(0., s2 / n - mm * mm) / n); } };
void one(TString file, std::map<TString, Acc> &A) {
  gSystem->Load("libAtSimulationData.so");
  TFile F(file); auto *t = (TTree *)F.Get("cbmsim"); TClonesArray *pt = nullptr; t->SetBranchAddress("AtTpcPoint", &pt);
  const double mA = 16.0147 * 931.49401;
  auto T = [&](double p) { return std::sqrt(p * p + mA * mA) - mA; };
  for (Long64_t i = 0; i < t->GetEntries(); i += 2) {
    t->GetEntry(i);
    std::vector<AtMCPoint *> b;
    for (int k = 0; k < pt->GetEntries(); ++k) { auto *q = (AtMCPoint *)pt->At(k); if (q->GetTrackID() == 0) b.push_back(q); }
    if (b.size() < 3) continue;
    auto pOf = [](AtMCPoint *q) { return 1000 * std::sqrt(q->GetPx()*q->GetPx() + q->GetPy()*q->GetPy() + q->GetPz()*q->GetPz()); };
    double eWin = 0, eGas = 0; AtMCPoint *gas0 = nullptr, *z50 = nullptr;
    for (auto *q : b) {
      if (q->GetDetectorID() == 4) eWin += q->GetEnergyLoss() * 1000;
      else { eGas += q->GetEnergyLoss() * 1000; if (!gas0) gas0 = q; if (!z50 && q->GetZ() >= 50) z50 = q; }
    }
    if (!gas0) continue;
    A["p window entrance [MeV/c]"].add(pOf(b[0]));
    A["T window entrance [MeV]"].add(T(pOf(b[0])));
    A["energy lost in window [MeV]"].add(eWin);
    A["T gas entrance [MeV]"].add(T(pOf(gas0)));
    if (z50) { A["T at z = 50 cm [MeV]"].add(T(pOf(z50))); A["r at z = 50 cm [mm]"].add(10 * std::hypot(z50->GetX(), z50->GetY())); }
    A["energy deposited in gas [MeV]"].add(eGas);
    A["beam end z (vertex) [cm]"].add(b.back()->GetZ());
    A["beam end r [mm]"].add(10 * std::hypot(b.back()->GetX(), b.back()->GetY()));
  }
}
void beamq(TString f0, TString f6) {
  std::map<TString, Acc> A0, A6; one(f0, A0); one(f6, A6);
  printf("%-32s %22s %22s %16s\n", "quantity (mean)", "q = 0", "q = 6", "diff / error");
  for (auto &[k, a] : A0) { auto &b = A6[k]; double d = b.m() - a.m(), e = std::hypot(a.e(), b.e());
    printf("%-32s %12.4f +- %7.4f %12.4f +- %7.4f %+8.4f (%4.1f sig)  n=%ld/%ld\n", k.Data(), a.m(), a.e(), b.m(), b.e(), d, e > 0 ? d / e : 0., a.n, b.n); }
}
