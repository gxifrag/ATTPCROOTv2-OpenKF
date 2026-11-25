/// @file ebeam_z_C16.C
/// @brief Step 1 of the beam-energy correction study (2026-11-25): does CATIMA reproduce the beam
///        energy at the vertex, as a function of vertex z, in the simulation?
///
/// Truth: in every reaction event the beam momentum at the vertex is p(ejectile) + p(15C) (target at
/// rest), so T_beam(z) is known exactly from MCTrack. Model: T(z) = CATIMA energy of 16C after z mm of
/// gas, starting from E0 at z = 0 (AtELossCATIMA::GetEnergy). E0 and the density are arguments;
/// the default density is the one the (p,d) fits use (H2 300 torr, 293 K).
///
///   root -l -b -q 'ebeam_z_C16.C("<work>/ex1_s5002/ex1_s5002_sim.root","pd",0.740,"out.png")'

void ebeam_z_C16(TString simFile, TString reac = "pd", Double_t exGen = 0.740, TString outPng = "ebeam_z.png",
                 Double_t E0 = -1, Double_t density = 3.308e-5)
{
   gSystem->Load("libAtSimulationData.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   const double u = 931.49401;
   const bool pd = (reac == "pd");
   const double mb = 16.0147 * u, mr = 15.0105993 * u;
   const double me = pd ? 2.0141017781 * u : 3.0160492779 * u;
   const int pdgE = pd ? 1000010020 : 1000010030, pdgR = 1000060150;

   TFile *fs = TFile::Open(simFile);
   auto ts = (TTree *)fs->Get("cbmsim");
   TClonesArray *mc = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);
   std::vector<double> vz, vt;
   for (Long64_t i = 0; i < ts->GetEntries(); ++i) {
      ts->GetEntry(i);
      TVector3 pE, pR;
      double z = 0;
      int got = 0;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetMotherId() != -1)
            continue;
         TVector3 v(p->GetPx() * 1e3, p->GetPy() * 1e3, p->GetPz() * 1e3);
         if (p->GetPdgCode() == pdgE) {
            pE = v;
            z = p->GetStartZ() * 10;
            got |= 1;
         } else if (p->GetPdgCode() == pdgR) {
            pR = v;
            got |= 2;
         }
      }
      if (got != 3)
         continue;
      TVector3 pB = pE + pR;
      vz.push_back(z);
      vt.push_back(std::sqrt(pB.Mag2() + mb * mb) - mb);
   }
   // E0 from the truth itself if not given: mean T_beam of vertices in the first 20 mm, extrapolated
   if (E0 < 0) {
      TGraph g;
      for (size_t i = 0; i < vz.size(); ++i)
         if (vz[i] < 100)
            g.AddPoint(vz[i], vt[i]);
      TF1 l("l", "pol1", 0, 100);
      g.Fit(&l, "QN");
      E0 = l.GetParameter(0);
   }
   AtTools::AtELossCATIMA cat(density);
   cat.SetProjectile(16, 6, 16.0147);
   cat.SetMaterial({std::make_tuple(1, 1, 1)});
   auto Tcat = [&](double z) { return z <= 0 ? E0 : cat.GetEnergy(E0, z); };

   auto hT = new TH2D("hT", Form("beam KE at the vertex (%s);z_{vertex} [mm];T_{beam} [MeV]", reac.Data()), 105, -25,
                      1025, 200, 170, 186);
   auto hR = new TH2D("hR", "truth - CATIMA;z_{vertex} [mm];T_{true} - T_{CATIMA} [MeV]", 105, -25, 1025, 200, -1, 1);
   auto gC = new TGraph();
   for (double z = 0; z <= 1000; z += 10)
      gC->AddPoint(z, Tcat(z));
   for (size_t i = 0; i < vz.size(); ++i) {
      hT->Fill(vz[i], vt[i]);
      hR->Fill(vz[i], vt[i] - Tcat(vz[i]));
   }
   auto c = new TCanvas("c", "", 1400, 550);
   c->Divide(2, 1);
   c->cd(1);
   hT->Draw("colz");
   gC->SetLineColor(kRed);
   gC->SetLineWidth(2);
   gC->Draw("L same");
   auto lg = new TLegend(0.15, 0.15, 0.6, 0.27);
   lg->AddEntry(gC, Form("CATIMA, E_{0}=%.2f MeV, #rho=%.3g g/cm^{3}", E0, density), "l");
   lg->Draw();
   c->cd(2);
   hR->Draw("colz");
   auto pf = hR->ProfileX("pf");
   pf->SetLineColor(kRed);
   pf->SetMarkerColor(kRed);
   pf->SetMarkerStyle(20);
   pf->SetMarkerSize(0.5);
   pf->Draw("same");
   c->Print(outPng);
   printf("  %zu reaction events. E0 = %.3f MeV, density %.4g g/cm3\n", vz.size(), E0, density);
   printf("  CATIMA T(z):  z=0 %.2f  250 %.2f  500 %.2f  750 %.2f  1000 %.2f MeV\n", Tcat(0), Tcat(250), Tcat(500),
          Tcat(750), Tcat(1000));
   printf("  truth-CATIMA: mean %+.3f MeV, rms %.3f MeV;  in z 0-500: mean %+.3f rms %.3f\n", hR->GetMean(2),
          hR->GetRMS(2), [&] { hR->GetXaxis()->SetRangeUser(0, 500); return hR->GetMean(2); }(), hR->GetRMS(2));
}
