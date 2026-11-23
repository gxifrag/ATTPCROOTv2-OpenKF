// beam_energy.C    root -l 'beam_energy.C("data/gs_s1001_sim.root")'
// Energia del haz en las entradas de haz (pares). nMax=100000 -> ~50000 entradas de haz.
void beam_energy(TString file, Long64_t nMax = 100000)
{
   TFile *f = TFile::Open(file);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", file.Data()); return; }
   auto *t = (TTree *)f->Get("cbmsim");

   const TString m = "14.9176";                     // masa del 16C en GeV, la de plot_states_C16pd.C
   const TString even = "Entry$%2==0 && Length$(AtTpcPoint.fPz)>2";
   const TString beam = "AtTpcPoint.fTrackID==0";   // solo puntos del haz (no secundarios)

   // T [MeV] desde pz [GeV/c]; px y py del haz son despreciables
   auto Tof = [&](TString pz) {
      return TString::Format("1000*(sqrt((%s)*(%s)+%s*%s)-%s)", pz.Data(), pz.Data(), m.Data(), m.Data(), m.Data());
   };
   auto draw = [&](TString expr, TString name, TString bins, TString cut) {
      t->Draw(expr + ">>" + name + bins, cut, "goff", nMax);
      auto *h = (TH1D *)gDirectory->Get(name);
      printf("%-6s mean %.3f  RMS %.3f  entries %.0f\n", name.Data(), h->GetMean(), h->GetRMS(), h->GetEntries());
      return h;
   };

   TString zEnd = "MaxIf$(AtTpcPoint.fZ," + beam + ")";                       // z donde acaba el haz [cm]
   TString tEnd = Tof("MinIf$(AtTpcPoint.fPz," + beam + ")");                 // T donde acaba el haz
   TString dep  = "1000*Sum$(AtTpcPoint.fELoss*(AtTpcPoint.fDetectorID==3&&" + beam + "))";  // depositado en el gas

   auto *hGen = draw(Tof("MCTrack.fPz[0]"), "hGen", "(200,183.9,184.1)", even);   // 1) generado
   auto *hWin = draw(Tof("AtTpcPoint.fPz[0]"), "hWin", "(200,183.9,184.2)", even); // 2) entrada a la ventana
   auto *hGas = draw(Tof("AtTpcPoint.fPz[2]"), "hGas", "(300,183.0,183.7)", even); // 3) entrada al gas
   auto *hEnd = draw(tEnd, "hEnd", "(280,170,184)", even);                         // 4) donde acaba el haz
   auto *hDep = draw(dep, "hDep", "(300,0,12)", even);                             // 5) depositado en el gas

   // entradas en que el haz cruza todo el volumen (sin reaccion)
   t->Draw(zEnd + ">>hZ(1010,0,101)", even, "goff", nMax);
   auto *hZ = (TH1D *)gDirectory->Get("hZ");
   double nTot = hZ->GetEntries(), nPT = hZ->Integral(hZ->FindBin(99.95), hZ->GetNbinsX());
   printf("haz que cruza todo (z_fin>99.95 cm): %.0f de %.0f = %.3f%%   (esperado 0.51%% si son los vacios)\n",
          nPT, nTot, 100 * nPT / nTot);
   // energia perdida en todo el volumen, medida directamente con esos haces
   auto *hL = draw(dep, "hL", "(300,0,12)", even + " && " + zEnd + ">99.95");
   printf("perdida en 1000 mm (haces que cruzan): %.3f MeV   [LISE 10.49/940 mm; maxELoss=11; escaneo 11.06]\n", hL->GetMean());

   // T al final del haz frente a z: relacion exacta E(z)
   auto *p = new TProfile("pTz", "beam T at the vertex vs z;z [cm];T [MeV]", 100, 0, 100);
   t->Draw(tEnd + ":" + zEnd + ">>pTz", even, "goff", nMax);
   p->Fit("pol2", "Q");
   auto *fn = p->GetFunction("pol2");
   printf("T(0)=%.3f  perdida hasta 94 cm: %.3f MeV  hasta 100 cm: %.3f MeV\n", fn->Eval(0),
          fn->Eval(0) - fn->Eval(94), fn->Eval(0) - fn->Eval(100));

   auto *c = new TCanvas("cbe", "beam energy", 1500, 900);
   c->Divide(3, 2);
   c->cd(1); hGen->SetTitle("generated (MCTrack);T [MeV];beam entries"); hGen->Draw();
   c->cd(2); hWin->SetTitle("window entrance;T [MeV];beam entries"); hWin->Draw();
   c->cd(3); hGas->SetTitle("gas entrance;T [MeV];beam entries"); hGas->Draw();
   c->cd(4); hEnd->SetTitle("end of the beam track (vertex);T [MeV];beam entries"); hEnd->Draw();
   c->cd(5); hDep->SetTitle("energy deposited in the gas;MeV;beam entries"); hDep->Draw();
   c->cd(6); p->Draw();
   c->SaveAs("beam_energy.png");
}