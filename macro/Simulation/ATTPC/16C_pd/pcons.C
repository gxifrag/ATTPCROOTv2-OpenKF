// Momentum conservation per event: |p(d) + p(15C)| at creation (reaction entry 2k+1, MCTrack)
// minus |p(beam)| at its last point (beam entry 2k). Target at rest, so it must be ~0.
void pcons(TString file) {
  gSystem->Load("libAtSimulationData.so");
  TFile F(file); auto *t = (TTree *)F.Get("cbmsim"); TClonesArray *pt = nullptr, *mc = nullptr;
  t->SetBranchAddress("AtTpcPoint", &pt); t->SetBranchAddress("MCTrack", &mc);
  double s = 0, s2 = 0, dz = 0; long n = 0;
  for (Long64_t k = 0; k + 1 < t->GetEntries(); k += 2) {
    t->GetEntry(k); AtMCPoint *last = nullptr;
    for (int i = 0; i < pt->GetEntries(); ++i) { auto *q = (AtMCPoint *)pt->At(i); if (q->GetTrackID() == 0) last = q; }
    if (!last) continue;
    TVector3 pb(last->GetPx(), last->GetPy(), last->GetPz()); double zb = last->GetZ();
    t->GetEntry(k + 1); TVector3 ps(0, 0, 0); double zv = -999; int np = 0;
    for (int i = 0; i < mc->GetEntries(); ++i) { auto *m = (AtMCTrack *)mc->At(i); if (m->GetMotherId() != -1) continue; ps += TVector3(m->GetPx(), m->GetPy(), m->GetPz()); zv = m->GetStartZ(); ++np; }
    if (np != 2) continue;
    double d = 1000 * (ps.Mag() - pb.Mag()); s += d; s2 += d * d; dz += std::fabs(zv - zb); ++n;
  }
  double m = s / n;
  printf("%s: |p_products| - |p_beam at vertex| = %+.3f MeV/c (rms %.3f), mean |z_vertex - z_beam end| = %.4f cm, %ld events\n", gSystem->BaseName(file), m, std::sqrt(s2 / n - m * m), dz / n, n);
}
