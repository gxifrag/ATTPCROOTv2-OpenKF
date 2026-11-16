// plot_padplane_hole.C
#include "../AtBeamHole.h"

void plot_padplane_hole(double holeR = 20.0, TString mapFile = "")
{
   TString dir = getenv("VMCWORKDIR");
   if (mapFile == "") mapFile = dir + "/scripts/Lookup20150611.xml";

   auto mapping = std::make_shared<AtTpcMap>();
   mapping->ParseXMLMap(mapFile.Data());
   mapping->GeneratePadPlane();
   int nOff = AtSim::InhibitBeamHole(mapping, holeR);

   // Same rule as InhibitBeamHole (centre inside holeR). Check AtMap.h for an "IsInhibited"-type
   // getter if you want to read the real flag instead of re-applying the rule.
   TGraph *gOn = new TGraph(), *gOff = new TGraph();
   double maxOff = 0, minOn = 1e9;
   for (unsigned i = 0; i < mapping->GetNumPads(); ++i) {
      auto c = mapping->CalcPadCenter(i);
      if (std::isnan(c.X()) || std::isnan(c.Y())) continue;
      double r = std::hypot(c.X(), c.Y());
      if (r < holeR) { gOff->SetPoint(gOff->GetN(), c.X(), c.Y()); maxOff = std::max(maxOff, r); }
      else           { gOn->SetPoint(gOn->GetN(), c.X(), c.Y());  minOn  = std::min(minOn, r); }
   }
   printf("inhibited %d | largest inhibited r = %.2f mm | smallest active r = %.2f mm\n", nOff, maxOff, minOn);

   TCanvas *c = new TCanvas("cpp", "pad plane", 1500, 750);
   c->Divide(2, 1);
   for (int k = 1; k <= 2; ++k) {
      c->cd(k);
      double L = (k == 1) ? 270 : 60;
      auto *fr = new TH2D(Form("fr%d", k), k == 1 ? "full pad plane (pad centres);x [mm];y [mm]"
                                                   : "zoom on the hole;x [mm];y [mm]", 1, -L, L, 1, -L, L);
      fr->Draw();
      gOn->SetMarkerStyle(20);  gOn->SetMarkerSize(k == 1 ? 0.2 : 0.8);  gOn->SetMarkerColor(kGray + 1);
      gOff->SetMarkerStyle(20); gOff->SetMarkerSize(k == 1 ? 0.3 : 0.8); gOff->SetMarkerColor(kRed + 1);
      gOn->Draw("P same"); gOff->Draw("P same");
      auto *a = new TArc(0, 0, holeR); a->SetFillStyle(0); a->SetLineColor(kAzure + 1); a->SetLineWidth(2); a->Draw();
   }
   c->SaveAs("padplane_hole.png");
}