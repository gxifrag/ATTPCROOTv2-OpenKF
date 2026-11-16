// hit_radius.C
//   root -l 'hit_radius.C("reco.root","hole0")'
//   root -l 'hit_radius.C("reco.root","hole0","data/gs_s1001_sim.root")'   // + true vertex
// pad plane size : radius=27,5 cm
// nb = 54 matches the inner pad size (0.5cm).

void hit_radius(TString file, TString label = "x", TString simFile = "", double rmax = 280, int nb = 54,
                int maxEv = 20000)
{
   TFile *f = TFile::Open(file);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", file.Data()); return; }
   auto *t = (TTree *)f->Get("cbmsim");
   TClonesArray *arr = nullptr;
   t->SetBranchAddress("AtEventH", &arr);

   const char *ytitle = ";r [mm];hits / mm^{2} / event";
   auto *hAll = new TH1D("hr_" + label, ytitle, nb, 0, rmax);
   auto *hC16 = new TH1D("hr16_" + label, "", nb, 0, rmax);
   auto *h15  = new TH1D("hr15_" + label, "", nb, 0, rmax);
   auto *hD   = new TH1D("hrD_" + label, "", nb, 0, rmax);
   auto *hP = new TH1D("hrP_" + label, "", nb, 0, rmax);
   auto *hxyA = new TH2D("hxy_" + label, "all hits;x [mm];y [mm]", 120, -60, 60, 120, -60, 60);
   auto *hxy16 = new TH2D("hxy16_" + label, "hits labelled ^{16}C;x [mm];y [mm]", 120, -60, 60, 120, -60, 60);
   long nLab = 0, nHit = 0, n16 = 0, n15 = 0, nD = 0, nP = 0;

   Long64_t n = std::min<Long64_t>(t->GetEntries(), maxEv);
   for (Long64_t e = 0; e < n; ++e) {
      t->GetEntry(e);
      if (!arr || arr->GetEntriesFast() == 0) continue;
      auto *ev = (AtEvent *)arr->At(0);
      for (auto &h : ev->GetHitArray()) {
         auto p = h.GetPosition();
         double r = std::hypot(p.X(), p.Y());
         hAll->Fill(r); hxyA->Fill(p.X(), p.Y()); ++nHit;
         const auto &pts = h.GetMCSimPointArray();
         if (pts.empty()) continue;
         ++nLab;
         for (auto &q : pts) if (q.Z == 1 && q.A == 1) { ++nP; hP->Fill(r); break; }
         const AtHit::MCSimPoint *best = &pts[0];
         for (auto &q : pts) if (q.eloss > best->eloss) best = &q;
         if (best->Z == 6 && best->A == 16)      { hC16->Fill(r); hxy16->Fill(p.X(), p.Y()); ++n16; }
         else if (best->Z == 6 && best->A == 15) { h15->Fill(r); ++n15; }
         else if (best->Z == 1 && best->A == 2)  { hD->Fill(r);  ++nD; }
      }
   }
   printf("%s: %lld events, %ld hits, %ld with MC label | 16C %ld, 15C %ld, d %ld, p %ld, other %ld\n", label.Data(), n,
          nHit, nLab, n16, n15, nD, nP, nLab - n16 - n15 - nD);

   // true beam spot = reaction vertex from the sim file (optional)
   TH2D *hv = nullptr;
   if (simFile != "") {
      gSystem->Load("libAtSimulationData.so");
      TChain ts("cbmsim"); ts.Add(simFile);
      TClonesArray *mc = nullptr; ts.SetBranchAddress("MCTrack", &mc);
      hv = new TH2D("hv_" + label, "true vertex x-y (MC);x [mm];y [mm]", 120, -30, 30, 120, -30, 30);
      Long64_t ns = std::min<Long64_t>(ts.GetEntries(), maxEv);
      for (Long64_t e = 0; e < ns; ++e) {
         ts.GetEntry(e);
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *p = (AtMCTrack *)mc->At(k);
            if (!p || p->GetMotherId() != -1 || p->GetPdgCode() != 1000010020) continue;
            hv->Fill(p->GetStartX() * 10, p->GetStartY() * 10);   // cm -> mm
         }
      }
      printf("true vertex: mean (%.2f, %.2f) RMS (%.2f, %.2f) mm\n", hv->GetMean(1), hv->GetMean(2),
             hv->GetRMS(1), hv->GetRMS(2));
   }

   // hits / mm^2 / event: divide each bin by the ring area and by the number of events. No reference region.
   for (auto *h : {hAll, hC16, h15, hD}) {
      for (int i = 1; i <= nb; ++i) {
         double r1 = h->GetBinLowEdge(i), r2 = r1 + h->GetBinWidth(i);
         double A = TMath::Pi() * (r2 * r2 - r1 * r1);
         h->SetBinContent(i, h->GetBinContent(i) / (A * n));
         h->SetBinError(i, h->GetBinError(i) / (A * n));
      }
   }

   auto *c = new TCanvas("c_" + label, label, 1400, 900);
   c->Divide(2, 2);
      c->cd(1);
   gPad->SetLogy(); gPad->SetGridy();
   hAll->SetMarkerStyle(20); hAll->SetMarkerSize(0.8);
   hAll->SetMarkerColor(kBlack); hAll->SetLineColor(kBlack);
   hC16->SetLineColor(kRed + 1);  hC16->SetLineWidth(2);
   h15->SetLineColor(kGreen + 2); h15->SetLineWidth(2);
   hD->SetLineColor(kAzure + 1);  hD->SetLineWidth(2);
   hP->SetLineColor(kOrange + 1); hP->SetLineWidth(2);
   hAll->SetMinimum(1e-6);
   hAll->Draw("PE");
   hC16->Draw("hist same"); h15->Draw("hist same"); hD->Draw("hist same");
   if (nP > 0) hP->Draw("hist same");

   // líneas verticales: DESPUÉS de dibujar, con los límites ya fijados
   gPad->Update();
   double y1 = TMath::Power(10, gPad->GetUymin());   // en log, Uymin/Uymax son exponentes
   double y2 = TMath::Power(10, gPad->GetUymax());
   auto *lg = new TLegend(0.55, 0.62, 0.89, 0.88);
   lg->AddEntry(hAll, "all hits", "p");
   lg->AddEntry(hC16, "^{16}C (beam)", "l");
   lg->AddEntry(h15, "^{15}C (residue)", "l");
   lg->AddEntry(hD, "deuteron", "l");
   for (double R : {20., 30.}) {
      auto *l = new TLine(R, y1, R, y2);
      l->SetLineStyle(2); l->SetLineWidth(2);
      l->SetLineColor(R == 20 ? kMagenta : kBlack);
      l->Draw();
      lg->AddEntry(l, Form("r = %.0f mm", R), "l");
   }
   lg->Draw();

   c->cd(2); hxyA->Draw("colz");
   c->cd(3); hxy16->Draw("colz");
   for (double R : {20., 30.}) {
      auto *a = new TArc(0, 0, R); a->SetFillStyle(0); a->SetLineWidth(2);
      a->SetLineColor(R == 20 ? kMagenta : kBlack); a->Draw();
   }
   if (hv) { c->cd(4); hv->Draw("colz"); }
   c->SaveAs("hit_radius_" + label + ".png");
   TFile o("hr_" + label + ".root", "RECREATE");
   hAll->Write(); hC16->Write(); h15->Write(); hD->Write(); hxyA->Write(); hxy16->Write(); if (hv) hv->Write();
}