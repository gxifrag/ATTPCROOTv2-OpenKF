/// @file plot_point_C16.C
/// @brief Control plots of ONE point of the Ex-resolution grid (run_res_grid_C16.sh): what was
///        generated (MCTrack), what Geant4 deposited in the gas (AtTpcPoint, split by particle),
///        and the Ex "delta" at each stage: generated -> fixed beam energy -> reconstructed.
///
///   root -l -b -q 'plot_point_C16.C("<work>/ex1_s5002/ex1_s5002_sim.root",
///                                   "<fit>/res/pd_ex1_s5002_std.root","pd",0.740,"out.pdf")'
///
/// Notes on the AtTpcPoint branches (AtDetectors/AtTpc/AtTpc.cxx, addHit):
///  - one AtTpcPoint per Geant4 STEP in the gas (StepLimiter: <= 1 mm), so every histogram filled
///    per point is weighted by track length: long tracks dominate.
///  - fELoss = gMC->Edep() of that step, in GeV.
///  - fEnergyIni / fAngleIni are NOT properties of the step: they are the generator's initial KE [MeV]
///    and lab angle [deg] of the TRACK the step belongs to (AtVertexPropagator), copied into every
///    point of that track. Track 0 (the beam in beam events, the heavy residue in reaction events)
///    gets 0/0 by construction, and so do secondaries the generator never set.

double omegaInv(double ke, double m) { return std::sqrt(ke * ke + 2 * ke * m); } // |p| from KE

void plot_point_C16(TString simFile, TString resFile = "", TString reac = "pd", Double_t exGen = 0.740,
                    TString outPdf = "point.pdf", Double_t Ebeam = -1)
{
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(1110);
   gStyle->SetPalette(kBird);
   const double u = 931.49401;
   const bool pd = (reac == "pd");
   const double mb = 16.0147 * u, mr = 15.0105993 * u;
   const double mt = pd ? 1.0078250322 * u : 2.0141017781 * u;
   const double me = pd ? 2.0141017781 * u : 3.0160492779 * u;
   const int pdgE = pd ? 1000010020 : 1000010030, pdgR = 1000060150, pdgB = 1000060160;
   const int aE = pd ? 2 : 3;
   if (Ebeam < 0)
      Ebeam = pd ? 183.5 : 184.0;
   const char *ej = pd ? "d" : "t";

   // Ex from ejectile KE and lab angle for a given beam KE (missing mass, target at rest)
   auto exOf = [&](double tb, double ke, double th) {
      double E1 = tb + mb, p1 = omegaInv(tb, mb), E3 = ke + me, p3 = omegaInv(ke, me);
      double E4 = E1 + mt - E3, p4sq = p1 * p1 + p3 * p3 - 2 * p1 * p3 * std::cos(th);
      return std::sqrt(std::max(E4 * E4 - p4sq, 0.0)) - mr;
   };

   TFile *fs = TFile::Open(simFile);
   if (!fs || fs->IsZombie())
      return;
   auto ts = (TTree *)fs->Get("cbmsim");
   TClonesArray *mc = nullptr, *pts = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   ts->SetBranchAddress("AtTpcPoint", &pts);

   // ---- generation (MCTrack primaries)
   auto hKeTh = new TH2D("hKeTh", Form("generated %s;#theta_{lab} [deg];KE [MeV]", ej), 90, 0, 90, 120, 0, 120);
   auto hRes = new TH2D("hRes", "generated ^{15}C;#theta_{lab} [deg];KE [MeV]", 100, 0, 10, 100, 100, 200);
   auto hVz = new TH1D("hVz", "generated vertex z (reaction events);z [mm];events", 110, -50, 1050);
   auto hTb = new TH1D("hTb", "beam KE at the vertex (from p_{15C}+p_{ej});T_{beam} [MeV];events", 120, 120, 190);
   auto hTbZ = new TH2D("hTbZ", "beam KE at the vertex vs z;z [mm];T_{beam} [MeV]", 110, -50, 1050, 120, 120, 190);
   // ---- the Ex deltas
   auto hExGen = new TH1D("hExGen", "E_{x} from TRUE kinematics, TRUE beam energy;E_{x}-E_{x}^{gen} [MeV];events",
                          400, -1, 1);
   auto hExFix = new TH1D("hExFix", Form("E_{x} from TRUE kinematics, FIXED E_{beam} = %.1f MeV;E_{x}-E_{x}^{gen} [MeV];events", Ebeam),
                          400, -1, 1);
   // ---- Geant4 in the gas (AtTpcPoint), by particle
   auto hElE = new TH1D("hElE", Form("energy deposited per step, %s;dE [keV];steps", ej), 150, 0, 300);
   auto hElR = new TH1D("hElR", "dE per step, ^{15}C/^{16}C, GAS (drift_volume);dE [keV];steps", 150, 0, 150);
   auto hElW = new TH1D("hElW", "dE per step, ^{15}C/^{16}C, WINDOW (tpc_window);dE [keV];steps", 150, 0, 1500);
   auto hSumE = new TH2D("hSumE", Form("%s: total dE in gas vs initial KE;KE_{ini} [MeV];#SigmadE [MeV]", ej), 120, 0,
                         120, 120, 0, 3);
   auto hStepZ = new TH1D("hStepZ", Form("%s steps: z;z [mm];steps", ej), 110, -50, 1050);
   auto hStepR = new TH1D("hStepR", Form("%s steps: radius;r [mm];steps", ej), 100, 0, 300);
   auto hSpec = new TH1D("hSpec", "AtTpcPoint steps by particle;;steps", 6, 0, 6);
   const char *lab[6] = {ej, "^{15}C", "^{16}C (beam)", "p", "e^{-}", "other (secondaries)"};
   for (int i = 0; i < 6; ++i)
      hSpec->GetXaxis()->SetBinLabel(i + 1, lab[i]);

   long nReac = 0, nBeam = 0;
   for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
      ts->GetEntry(i);
      // primaries
      TVector3 pE, pR;
      double keE = -1, zE = 0;
      bool hasR = false;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetMotherId() != -1)
            continue;
         TVector3 v(p->GetPx() * 1e3, p->GetPy() * 1e3, p->GetPz() * 1e3);
         if (p->GetPdgCode() == pdgE) {
            pE = v;
            keE = std::sqrt(v.Mag2() + me * me) - me;
            zE = p->GetStartZ() * 10;
         } else if (p->GetPdgCode() == pdgR) {
            pR = v;
            hasR = true;
         }
      }
      if (keE >= 0 && hasR) {
         ++nReac;
         double thE = pE.Theta();
         hKeTh->Fill(thE * TMath::RadToDeg(), keE);
         // residue KE with its EXCITED mass (the generator gives 15C mass m4 + Ex)
         double mR = mr + exGen, keR = std::sqrt(pR.Mag2() + mR * mR) - mR;
         hRes->Fill(pR.Theta() * TMath::RadToDeg(), keR);
         hVz->Fill(zE);
         // beam momentum at the vertex = p(ejectile) + p(residue) (target at rest)
         TVector3 pB = pE + pR;
         double tb = std::sqrt(pB.Mag2() + mb * mb) - mb;
         hTb->Fill(tb);
         hTbZ->Fill(zE, tb);
         // Ex with the ejectile's angle to the BEAM direction (beam not exactly along z)
         hExGen->Fill(exOf(tb, keE, pE.Angle(pB)) - exGen);
         hExFix->Fill(exOf(Ebeam, keE, thE) - exGen);
      } else
         ++nBeam;
      // steps
      std::map<int, double> sumE;
      std::map<int, double> eIni;
      for (int k = 0; k < pts->GetEntriesFast(); ++k) {
         auto q = (AtMCPoint *)pts->At(k);
         int A = q->GetMassNum(), Z = q->GetAtomicNum();
         double de = q->GetEnergyLoss() * 1e6; // GeV -> keV
         int s = 5;
         if (A == aE && Z == 1)
            s = 0;
         else if (A == 15 && Z == 6)
            s = 1;
         else if (A == 16 && Z == 6)
            s = 2;
         else if (A == 1 && Z == 1)
            s = 3;
         else if (Z == -1 || (A == 0 && Z == 0))
            s = 4;
         hSpec->Fill(s);
         if (s == 0) {
            hElE->Fill(de);
            hStepZ->Fill(q->GetZ() * 10);
            hStepR->Fill(std::hypot(q->GetX(), q->GetY()) * 10);
            sumE[q->GetTrackID()] += de * 1e-3;
            eIni[q->GetTrackID()] = q->GetEIni();
         } else if (s == 1 || s == 2)
            (q->GetVolName() == "tpc_window" ? hElW : hElR)->Fill(de);
      }
      for (auto &[id, e] : sumE)
         hSumE->Fill(eIni[id], e);
   }
   printf("  %s: %ld reaction events, %ld beam-only events\n", simFile.Data(), nReac, nBeam);

   // ---- reconstruction (res_ex_C16.C output), optional
   TH1D *hEx = nullptr, *hExT = nullptr, *hTh = nullptr;
   TH2D *hKE = nullptr, *hZ = nullptr;
   if (resFile.Length()) {
      TFile *fr = TFile::Open(resFile);
      if (fr && !fr->IsZombie()) {
         hEx = (TH1D *)fr->Get("hEx");
         hExT = (TH1D *)fr->Get("hExT");
         hKE = (TH2D *)fr->Get("hKE");
         hTh = (TH1D *)fr->Get("hTh");
         hZ = (TH2D *)fr->Get("hZ");
      }
   }

   TString pdf = outPdf;
   auto c1 = new TCanvas("c1", "generation", 1500, 900);
   c1->Divide(3, 2);
   c1->cd(1); hKeTh->Draw("colz");
   c1->cd(2); hRes->Draw("colz");
   c1->cd(3); hVz->Draw();
   c1->cd(4); hTb->Draw();
   c1->cd(5); hTbZ->Draw("colz");
   c1->cd(6); hSpec->Draw("hist"); gPad->SetLogy();
   c1->Print(pdf + "(");
   auto c2 = new TCanvas("c2", "Geant4 in the gas", 1500, 900);
   c2->Divide(3, 2);
   c2->cd(1); hElE->Draw(); gPad->SetLogy();
   c2->cd(2); hElR->Draw(); gPad->SetLogy();
   c2->cd(3); hSumE->Draw("colz");
   c2->cd(4); hStepZ->Draw();
   c2->cd(5); hStepR->Draw();
   c2->cd(6); hElW->Draw(); gPad->SetLogy();
   c2->Print(pdf);
   auto c3 = new TCanvas("c3", "Ex deltas", 1500, 900);
   c3->Divide(3, 2);
   c3->cd(1); hExGen->Draw();
   c3->cd(2); hExFix->Draw();
   if (hEx) {
      c3->cd(3);
      hEx->SetTitle("E_{x} RECONSTRUCTED (UKF, analysis cuts);E_{x}-E_{x}^{gen} [MeV];tracks");
      hEx->GetXaxis()->SetRangeUser(-2, 2);
      hEx->Draw();
   }
   if (hKE) { c3->cd(4); hKE->Draw("colz"); }
   if (hTh) { c3->cd(5); hTh->Draw(); }
   if (hZ) { c3->cd(6); hZ->Draw("colz"); }
   c3->Print(pdf + ")");
   // PNGs of each page next to the PDF
   TString base = pdf;
   base.ReplaceAll(".pdf", "");
   c1->Print(base + "_1_generation.png");
   c2->Print(base + "_2_geant4.png");
   c3->Print(base + "_3_ex_deltas.png");
   printf("  Ex true/true-beam: mean %.4f rms %.4f | fixed Ebeam: mean %.4f rms %.4f MeV\n", hExGen->GetMean(),
          hExGen->GetRMS(), hExFix->GetMean(), hExFix->GetRMS());
}
