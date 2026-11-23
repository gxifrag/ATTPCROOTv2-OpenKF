// window_loss.C   root -l 'window_loss.C("data/gs_s1001_sim.root")'
void window_loss(TString file, Long64_t nMax = 500000)
{
   TFile *f = TFile::Open(file);
   auto *t = (TTree *)f->Get("cbmsim");
   const char *cut = "Entry$%2==0 && Length$(AtTpcPoint.fPz)>2";

   // 1) energy deposited in tpc_window per beam entry [MeV]
   t->Draw("1000*Sum$(AtTpcPoint.fELoss*(AtTpcPoint.fDetectorID==4))>>hw(200,0,2)",
           "Entry$%2==0 && Length$(AtTpcPoint.fELoss)>2", "goff", nMax);
   auto *hw = (TH1D *)gDirectory->Get("hw");
   printf("window loss (sum fELoss): mean %.3f MeV, RMS %.3f, entries %.0f\n", hw->GetMean(), hw->GetRMS(),
          hw->GetEntries());

   // 2) kinetic energy at the first gas point, from the momentum [MeV]
   const char *m = "14.9176"; // 16C mass in GeV (16.0147 u), same as plot_states_C16pd.C
   t->Draw(Form("1000*(sqrt(AtTpcPoint.fPz[2]*AtTpcPoint.fPz[2]+%s*%s)-%s)>>hE(300,182,185)", m, m, m),
           "Entry$%2==0 && Length$(AtTpcPoint.fPz)>2", "goff", nMax);
   auto *hE = (TH1D *)gDirectory->Get("hE");
   printf("T at first gas point (from p): mean %.3f MeV, RMS %.3f, entries %.0f\n", hE->GetMean(), hE->GetRMS(),
          hE->GetEntries());

   auto *c = new TCanvas("cwin", "window loss", 1200, 500);
   c->Divide(2, 1);
   c->cd(1); hw->SetTitle("energy lost in tpc_window;MeV;beam entries"); hw->Draw();
   c->cd(2); hE->SetTitle("beam T at the first gas point;MeV;beam entries"); hE->Draw();
   c->SaveAs("window_loss_ex3_Ebeam184MeV.png");

    t->Draw("AtTpcPoint.fZ[2]>>hz2(200,-0.01,1.0)", cut, "goff", nMax);        // z (cm) del indice 2
   auto *hz = (TH1D *)gDirectory->Get("hz2");
   printf("z del indice 2: media %.5f cm, RMS %.5f, entradas %.0f\n", hz->GetMean(), hz->GetRMS(), hz->GetEntries());

   t->Draw("AtTpcPoint.fDetectorID[2]>>hid(10,0,10)", cut, "goff", nMax);     // 3 = gas
   auto *hi = (TH1D *)gDirectory->Get("hid");
   printf("detector del indice 2: media %.3f (3.000 = todos en gas)\n", hi->GetMean());

   t->Draw(Form("1000*(sqrt(AtTpcPoint.fPz[0]*AtTpcPoint.fPz[0]+%s*%s)-%s)>>hT0(400,183,185)", m, m, m), cut, "goff", nMax);
   auto *h0 = (TH1D *)gDirectory->Get("hT0");                                  // T a la entrada de la ventana
   printf("T en el indice 0: media %.3f MeV, RMS %.3f\n", h0->GetMean(), h0->GetRMS());
}