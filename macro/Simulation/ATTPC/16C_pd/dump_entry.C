/// @file dump_entry.C
/// @brief Print every variable of ONE entry of a simulation file (cbmsim tree), with units.
///
/// Written 2026-11-18 to understand the tree (see the simulation log, "What is in one entry of the
/// simulation tree"). Even entries are beam-only, odd entries are the reaction.
///
///   root -l -b -q 'dump_entry.C("data/gs_s1001_sim.root", 1)'          // reaction entry 1
///   root -l -b -q 'dump_entry.C("data/gs_s1001_sim.root", 0, 5)'       // beam entry, 5 points/track

void dump_entry(TString file = "data/gs_s1001_sim.root", Long64_t entry = 1, int pointsPerTrack = 3)
{
   gSystem->Load("libAtSimulationData.so");
   TFile f(file);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *mc = nullptr, *pt = nullptr;
   FairMCEventHeader *hdr = nullptr;
   t->SetBranchAddress("MCTrack", &mc);
   t->SetBranchAddress("AtTpcPoint", &pt);
   t->SetBranchAddress("MCEventHeader.", &hdr);
   t->GetEntry(entry);

   printf("\n=== %s, entry %lld (%s) ===\n", file.Data(), entry, entry % 2 ? "odd: reaction" : "even: beam");
   if (hdr)
      printf("MCEventHeader: RunID %u  EventID %d  vertex (%.3f, %.3f, %.3f) cm  NPrim %d\n", hdr->GetRunID(),
             hdr->GetEventID(), hdr->GetX(), hdr->GetY(), hdr->GetZ(), hdr->GetNPrim());

   printf("\n--- MCTrack: one row per particle in the stack (index = TrackID) ---\n");
   printf("%3s %11s %7s %11s %11s %11s %9s %9s %9s %9s %9s\n", "idx", "PdgCode", "Mother", "Px[GeV/c]",
          "Py[GeV/c]", "Pz[GeV/c]", "StartX", "StartY", "StartZ", "StartT", "Mass[GeV]");
   for (int k = 0; k < mc->GetEntries(); ++k) {
      auto *m = (AtMCTrack *)mc->At(k);
      printf("%3d %11d %7d %11.6f %11.6f %11.6f %9.3f %9.3f %9.3f %9.3f %9.5f\n", k, m->GetPdgCode(),
             m->GetMotherId(), m->GetPx(), m->GetPy(), m->GetPz(), m->GetStartX(), m->GetStartY(), m->GetStartZ(),
             m->GetStartT(), m->GetMass());
   }

   printf("\n--- AtTpcPoint: one row per Geant4 step inside a sensitive volume (first %d per track) ---\n",
          pointsPerTrack);
   printf("%4s %5s %13s %3s %4s %9s %9s %10s %10s %10s %10s %9s %9s %10s %9s %8s %4s %4s\n", "i", "Track",
          "VolName", "Det", "Copy", "X[cm]", "Y[cm]", "Z[cm]", "Px[GeV/c]", "Py", "Pz", "Time[ns]", "Len[cm]",
          "ELoss[MeV]", "EIni[MeV]", "AIni", "A", "Z");
   std::map<int, int> shown, total;
   for (int i = 0; i < pt->GetEntries(); ++i)
      total[((AtMCPoint *)pt->At(i))->GetTrackID()]++;
   for (int i = 0; i < pt->GetEntries(); ++i) {
      auto *q = (AtMCPoint *)pt->At(i);
      int tr = q->GetTrackID();
      if (shown[tr]++ >= pointsPerTrack)
         continue;
      printf("%4d %5d %13s %3d %4d %9.4f %9.4f %10.4f %10.6f %10.6f %10.6f %9.3f %9.4f %10.5f %9.4f %8.3f %4d %4d\n", i,
             tr, q->GetVolName().Data(), q->GetDetectorID(), q->GetDetCopyID(), q->GetX(), q->GetY(), q->GetZ(),
             q->GetPx(), q->GetPy(), q->GetPz(), q->GetTime(), q->GetLength(), q->GetEnergyLoss() * 1000,
             q->GetEIni(), q->GetAIni(), q->GetMassNum(), q->GetAtomicNum());
   }
   printf("\npoints per track:");
   for (auto &[tr, n] : total)
      printf("  track %d: %d", tr, n);
   printf("   (total %d)\n\n", pt->GetEntries());
}
