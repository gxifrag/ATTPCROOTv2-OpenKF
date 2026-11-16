// root -l 'mc_structure.C("data/gs_s1001_sim.root")'
void mc_structure(TString file, Long64_t nMax = 500000)
{
   gSystem->Load("libAtSimulationData.so");
   TChain t("cbmsim"); t.Add(file);
   t.SetBranchStatus("*", 0); t.SetBranchStatus("MCTrack*", 1);   // solo lo necesario
   TClonesArray *mc = nullptr; t.SetBranchAddress("MCTrack", &mc);
   std::map<std::string, std::pair<long, Long64_t>> comp;   // composición -> (cuenta, primera entrada)
   long withD[2] = {0, 0};
   Long64_t n = std::min<Long64_t>(t.GetEntries(), nMax);
   for (Long64_t e = 0; e < n; ++e) {
      t.GetEntry(e);
      std::vector<std::string> v; bool hasD = false;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *p = (AtMCTrack *)mc->At(k);
         if (!p) continue;
         bool prim = p->GetMotherId() == -1;
         if (prim && p->GetPdgCode() == 1000010020) hasD = true;
         v.push_back(std::to_string(p->GetPdgCode()) + (prim ? "" : "(sec)"));
      }
      std::sort(v.begin(), v.end());
      std::string key; for (auto &s : v) key += s + " ";
      auto &c = comp[key]; if (c.first++ == 0) c.second = e;
      if (hasD) ++withD[e % 2];
   }
   printf("%lld entradas; con d primario: pares %ld, impares %ld\n", n, withD[0], withD[1]);
   for (auto &c : comp) printf("%8ld  (primera en entrada %lld)  tracks: %s\n", c.second.first, c.second.second, c.first.c_str());
}