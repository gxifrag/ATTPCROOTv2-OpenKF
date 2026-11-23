/// @file acceptance_C16pd_allstates.C
/// @brief 16C(p,d)15C acceptance vs theta_cm for ALL states in ONE run, with every cut and every
///        hard-coded number saved inside the output file.
///
/// ============================================================================================
/// WHY THIS FILE EXISTS (created 2026-11-18, copy of acceptance_C16pd.C, which is left untouched)
/// ============================================================================================
/// The original macro makes ONE state per run and writes ONE file per state. Nothing guarantees
/// the four runs used the same cuts, and the files do not say which cuts they used, so it cannot
/// be checked afterwards. This version:
///
///   1. Loops over ALL states in a single run  ->  the same cuts for every state, by construction.
///   2. Writes ONE file with one histogram per state (hAcc_gs, hAcc_ex1, ...), like the (d,t)
///      acceptance file.
///   3. Saves every cut and hard-coded parameter in the file (folder "settings/"), so the analysis
///      (or you) can read them back and compare them with the data selection.
///   4. Reads Ex and Ebeam from each simulation file (C16_pd_sim.C stores them) instead of
///      trusting numbers typed by hand, and stops if the states were simulated differently.
///
/// Also corrected: the original header said "14C(p,p')" and "proton". This is 16C(p,d)15C and the
/// detected particle is the DEUTERON (the code already looked for deuterons; only comments were old).
///
/// ============================================================================================
/// WHAT "ACCEPTANCE" MEANS HERE
/// ============================================================================================
///   acceptance(theta_cm) = reconstructed / generated        (per theta_cm bin, per state)
///
///   generated     (denominator) = simulated reactions whose TRUE vertex z is inside [zMin, zMax].
///                                 Taken from MC truth, so reconstruction problems cannot affect it.
///   reconstructed (numerator)   = of those, the events that have a fitted track that
///                                   (a) passes the same cuts as the data (chi2/ndf, KE, vertex z), and
///                                   (b) IS the true deuteron (angle and energy close to the truth).
///
///   The data are corrected by dividing the counts by this number: N_true = N_measured / acceptance.
///   So the cuts here MUST be the cuts used on the data -- that is the whole point of saving them.
///
/// ============================================================================================
/// INPUT FILES
/// ============================================================================================
/// For each state, file names are built from a PATTERN in which "%s" is replaced by
/// "<state>_<seed>". With seed = "s1001":
///   simPattern "data/%s_sim.root"      ->  data/gs_s1001_sim.root, data/ex1_s1001_sim.root, ...
///   fitPattern "data/fit/%s_fit.root"  ->  data/fit/gs_s1001_fit.root, ...
/// The fit file must be the fit of the UNGATED reconstruction (one entry per simulated event);
/// the macro checks that the entry counts agree. Relative paths are relative to THIS macro's folder.
/// !! Check fitPattern against where your fitter actually writes its output before running. !!
///
/// ============================================================================================
/// HOW TO RUN
/// ============================================================================================
///   root -b -q 'acceptance_C16pd_allstates.C'
///   root -b -q 'acceptance_C16pd_allstates.C("s1001","data/%s_sim.root","data/fit/%s_fit.root")'
///
/// ============================================================================================
/// HOW TO READ THE SAVED SETTINGS BACK
/// ============================================================================================
///   TFile f("diagnostics/acceptance_pd_allstates.root");
///   f.Get("settings/summary")->Print();                             // everything, as readable text
///   ((TParameter<double>*)f.Get("settings/chi2Cut"))->GetVal();     // one value, as a number
///   f.ls(); f.cd("settings"); gDirectory->ls();                     // list what is inside

#include <TParameter.h>
#include <TNamed.h>
#include <TDatime.h>
#include <tuple>
#include <vector>

// ---------------------------------------------------------------------------------------------
// THE STATES. One row per 15C level: the name used in file and histogram names, and the
// excitation energy (MeV) we EXPECT. The energy actually used is read from each simulation file;
// if it differs from this table by more than 5 keV the macro stops, because then the file is not
// the state we think it is. Edit the states HERE and nowhere else.
// ---------------------------------------------------------------------------------------------
struct AccState {
   const char *name;
   double ExExpected;
};
static const AccState kStates[] = {
   {"gs", 0.000},  // 1/2+ ground state
   {"ex1", 0.740}, // 5/2+
   {"ex2", 3.100}, // (simulated at 3.100 or 3.103 depending on the script; both pass the check)
   {"ex3", 4.660}, // (4.660 or 4.657, same remark)
};
static const int kNStates = sizeof(kStates) / sizeof(kStates[0]);

// ---------------------------------------------------------------------------------------------
// HARD-CODED CONSTANTS. Every number the physics depends on is named here (no "magic numbers"
// hidden in the code below) and every one of them is saved in the output file.
// ---------------------------------------------------------------------------------------------
static const double kAmu = 931.49401;                 // MeV per atomic mass unit
static const double kMassC16 = 16.0147 * kAmu;        // beam, 16C
static const double kMassP = 1.007825 * kAmu;         // target, proton
static const double kMassD = 2.0141017781 * kAmu;     // detected particle, deuteron
static const double kMassC15 = 15.0105993 * kAmu;     // residual, 15C (ground state)
static const int kPdgDeuteron = 1000010020;           // how Geant4 labels a deuteron
static const double kExTolerance = 0.005;             // MeV: allowed |Ex(file) - Ex(table)|

static double acc_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// Two-body kinematics, identical to the original macro: from the beam energy and the detected
/// particle's (KE, theta_lab) it returns {Ex, theta_cm[deg]}.
static std::tuple<double, double> acc_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                           double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * acc_omega2(s, m1 * m1, m2 * m2) * acc_omega2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (acc_omega2(s, m1 * m1, m2 * m2) * acc_omega2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

/// Read a number that C16_pd_sim.C stored in the simulation file (TParameter<Double_t>).
/// Returns NaN if it is not there, so the caller can tell "missing" from "zero".
static double readSimParam(TFile *f, const char *name)
{
   auto *p = (TParameter<Double_t> *)f->Get(name);
   return p ? p->GetVal() : std::numeric_limits<double>::quiet_NaN();
}

// =============================================================================================
// THE MACRO. Every argument is a cut or a setting; ALL of them are saved in the output.
// The defaults for the cuts are the (p,d) DATA selection in a1975_analysis_KF/pd/macros/pd_config.h.
// If you change the data selection there, change it here too and rerun.
// =============================================================================================
void acceptance_C16pd_allstates(
   TString seed = "s1001",                                   // which simulation set (random seed tag)
   TString simPattern = "data/%s_sim.root",                  // see INPUT FILES above
   TString fitPattern = "data/fit/%s_fit.root",              // see INPUT FILES above -- CHECK IT
   TString outFile = "diagnostics/acceptance_pd_allstates.root",
   Double_t zMin = 0.0, Double_t zMax = 500.0, // mm, vertex window  (= pd_config.h vzLo, vzHi)
   Double_t chi2Cut = 5.0,                     // fitted chi2/ndf must be below (= pd_config.h chi2max)
   Double_t keLo = 2.0, Double_t keHi = 30.0,  // MeV, fitted deuteron KE window (= pd_config.h keLo, keHi)
   Bool_t useXtr = kFALSE,                     // which fitted kinematics to use (see the numerator)
   Double_t dThetaMax = 10.0,                  // deg, truth match on the angle (<=0 disables)
   Double_t keRatioMin = 0.5, Double_t keRatioMax = 2.0, // truth match on KE_fit / KE_true
   Int_t nBins = 36, Double_t cmMax = 180.0,   // theta_cm binning of the acceptance histograms
   Double_t EbeamOverride = -1.0)              // MeV. < 0: use the Ebeam stored in each sim file
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   // Relative paths are taken relative to the folder this macro lives in, so it works from any
   // directory. An absolute path (starting with / or ~) is used as it is.
   const TString macroDir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   auto inMacroDir = [&](TString p) {
      if (p.BeginsWith("/") || p.BeginsWith("~")) {
         gSystem->ExpandPathName(p);
         return p;
      }
      return macroDir + "/" + p;
   };

   // ------------------------------------------------------------------------------------------
   // Per-state results, kept until the end so we write the file only if EVERY state worked.
   // "All states or nothing": a file with some states missing would defeat the purpose.
   // ------------------------------------------------------------------------------------------
   struct StateResult {
      TString name, simFile, fitFile;
      double ExFile = 0, EbeamFile = 0, EbeamUsed = 0;
      long nEntries = 0, nGen = 0, nRec = 0;
      TH1D *hGen = nullptr, *hRec = nullptr, *hAcc = nullptr, *hLab = nullptr, *hLabR = nullptr;
   };
   std::vector<StateResult> results;

   for (int s = 0; s < kNStates; ++s) {
      StateResult r;
      r.name = kStates[s].name;
      const TString job = r.name + "_" + seed; // e.g. "gs_s1001"
      r.simFile = inMacroDir(Form(simPattern.Data(), job.Data()));
      r.fitFile = inMacroDir(Form(fitPattern.Data(), job.Data()));

      // ---- open the two files and check they belong together --------------------------------
      TFile *fs = TFile::Open(r.simFile);
      TFile *ff = TFile::Open(r.fitFile);
      if (!fs || fs->IsZombie() || !ff || ff->IsZombie()) {
         printf("\033[1;31m[%s] cannot open %s or %s -- STOPPING, no file written\033[0m\n", r.name.Data(),
                r.simFile.Data(), r.fitFile.Data());
         return;
      }
      // What was simulated: read it from the file instead of assuming it.
      r.ExFile = readSimParam(fs, "resEx");
      r.EbeamFile = readSimParam(fs, "Ebeam");
      if (std::isnan(r.ExFile) || std::isnan(r.EbeamFile)) {
         printf("\033[1;31m[%s] %s has no resEx/Ebeam parameters -- STOPPING\033[0m\n", r.name.Data(),
                r.simFile.Data());
         return;
      }
      if (std::fabs(r.ExFile - kStates[s].ExExpected) > kExTolerance) {
         printf("\033[1;31m[%s] the file was simulated at Ex = %.3f MeV, the table expects %.3f -- STOPPING\033[0m\n",
                r.name.Data(), r.ExFile, kStates[s].ExExpected);
         return;
      }
      // The beam energy for the kinematics: by default the energy the events were GENERATED with,
      // because theta_cm is computed from the TRUE deuteron, which was generated at that energy.
      r.EbeamUsed = EbeamOverride > 0 ? EbeamOverride : r.EbeamFile;

      TTree *ts = (TTree *)fs->Get("cbmsim");
      TTree *tf = (TTree *)ff->Get("cbmsim");
      if (!ts || !tf) {
         printf("\033[1;31m[%s] missing cbmsim tree -- STOPPING\033[0m\n", r.name.Data());
         return;
      }
      if (ts->GetEntries() != tf->GetEntries()) {
         printf("\033[1;31m[%s] ENTRY MISMATCH: sim %lld vs fit %lld -- the fit must be of the UNGATED reco,\n"
                "otherwise entries do not correspond and the acceptance would be nonsense. STOPPING\033[0m\n",
                r.name.Data(), ts->GetEntries(), tf->GetEntries());
         return;
      }
      r.nEntries = ts->GetEntries();

      TClonesArray *mc = nullptr, *te = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tf->SetBranchAddress("AtTrackingEvent", &te);

      // SetDirectory(nullptr): keep the histograms in memory after the input files are closed.
      r.hGen = new TH1D("hGen_" + r.name, "generated;#theta_{cm} [deg];reactions", nBins, 0, cmMax);
      r.hRec = new TH1D("hRec_" + r.name, "reconstructed;#theta_{cm} [deg];reactions", nBins, 0, cmMax);
      r.hLab = new TH1D("hLab_" + r.name, "generated;#theta_{lab} [deg];reactions", 90, 0, 180);
      r.hLabR = new TH1D("hLabR_" + r.name, "reconstructed;#theta_{lab} [deg];reactions", 90, 0, 180);
      for (TH1D *h : {r.hGen, r.hRec, r.hLab, r.hLabR}) {
         h->SetDirectory(nullptr);
         h->Sumw2();
      }

      const double massResidual = kMassC15 + r.ExFile; // 15C left in this excited state

      for (Long64_t i = 0; i < r.nEntries; ++i) {
         ts->GetEntry(i);
         if (!mc)
            continue;

         // ---- DENOMINATOR: find the true deuteron (primary particle, mother = -1) --------------
         double keT = -1, thT = -1, zT = -1e9;
         for (int k = 0; k < mc->GetEntriesFast(); ++k) {
            auto *t = (AtMCTrack *)mc->At(k);
            if (!t || t->GetPdgCode() != kPdgDeuteron || t->GetMotherId() != -1)
               continue;
            zT = t->GetStartZ() * 10.0; // FairRoot stores cm; the analysis works in mm
            double px = t->GetPx() * 1000, py = t->GetPy() * 1000, pz = t->GetPz() * 1000; // GeV -> MeV
            double p = std::sqrt(px * px + py * py + pz * pz);
            if (p <= 0)
               continue;
            keT = std::sqrt(p * p + kMassD * kMassD) - kMassD;
            thT = std::acos(pz / p);
            break;
         }
         if (keT <= 0)
            continue; // beam-only entry: no reaction happened in it

         // VERTEX WINDOW on the TRUE z: the denominator counts only reactions that really happened
         // inside the window. The luminosity used for the cross section must then be scaled to the
         // same window (in fit_pd_definitive.C: lumiFull * (vzHi - vzLo) / zFull).
         if (zT < zMin || zT > zMax)
            continue;
         auto [exT, cmT] = acc_kine(kMassC16, kMassP, kMassD, massResidual, r.EbeamUsed, thT, keT);
         if (std::isnan(cmT))
            continue;
         ++r.nGen;
         r.hGen->Fill(cmT);
         r.hLab->Fill(thT * TMath::RadToDeg());

         // ---- NUMERATOR: a fitted track in the SAME entry that passes the data cuts -----------
         tf->GetEntry(i);
         bool good = false;
         if (te && te->GetEntriesFast() > 0) {
            auto *ev = (AtTrackingEvent *)te->At(0);
            if (ev)
               for (auto &ft : ev->GetFittedTracks()) {
                  if (!ft)
                     continue;
                  const auto &md = ft->GetTrackMetadata();
                  double ndf = md ? md->GetNdf() : 0, chi2 = md ? md->GetChi2() : 0;
                  double c2n = ndf > 0 ? chi2 / ndf : 1e9;
                  // THE TWO FITTERS DEFINE THESE SLOTS OPPOSITELY:
                  //   UKF    GetKinematics() = at the vertex,        Xtr = at the first cluster
                  //   GENFIT GetKinematics() = at the first cluster, Xtr = vertex-corrected
                  // Use the SAME energy the data analysis uses: genfit needs useXtr = kTRUE.
                  const auto &kin = useXtr ? ft->GetKinematicsXtr() : ft->GetKinematics();
                  double ke = kin.kineticEnergy;
                  // CUT: fitted KE in range and chi2/ndf below chi2Cut (the data cut)
                  if (!(ke > keLo && ke < keHi && c2n < chi2Cut))
                     continue;
                  // TRUTH MATCH: the fitted track must BE the deuteron, not another track in the
                  // event (e.g. the scattered beam once it leaves the beam hole).
                  if (dThetaMax > 0) {
                     double dth = std::fabs(kin.theta * TMath::RadToDeg() - thT * TMath::RadToDeg());
                     double ratio = keT > 0 ? ke / keT : -1;
                     if (dth > dThetaMax || ratio < keRatioMin || ratio > keRatioMax)
                        continue;
                  }
                  // VERTEX WINDOW on the RECONSTRUCTED z: the data can only cut on that one.
                  double zR = ft->GetVertex(0).Z();
                  if (zR < zMin || zR > zMax)
                     continue;
                  good = true;
                  break;
               }
         }
         if (good) {
            ++r.nRec;
            r.hRec->Fill(cmT);
            r.hLabR->Fill(thT * TMath::RadToDeg());
         }
      }
      fs->Close();
      ff->Close();

      if (!r.nGen) {
         printf("\033[1;31m[%s] no true deuterons inside the vertex window -- STOPPING\033[0m\n", r.name.Data());
         return;
      }
      // acceptance = reconstructed / generated, with binomial errors ("B")
      r.hAcc = (TH1D *)r.hRec->Clone("hAcc_" + r.name);
      r.hAcc->SetDirectory(nullptr);
      r.hAcc->Divide(r.hRec, r.hGen, 1, 1, "B");
      r.hAcc->SetTitle(Form("16C(p,d) acceptance, %s, Ex = %.3f MeV;#theta_{cm} [deg];acceptance", r.name.Data(),
                            r.ExFile));
      printf("[%s] Ex %.3f MeV, Ebeam used %.3f MeV: generated %ld, reconstructed %ld, overall %.3f\n",
             r.name.Data(), r.ExFile, r.EbeamUsed, r.nGen, r.nRec, double(r.nRec) / r.nGen);
      results.push_back(r);
   }

   // ---- all states must have been simulated with the same beam energy ------------------------
   for (auto &r : results)
      if (std::fabs(r.EbeamFile - results[0].EbeamFile) > 1e-6) {
         printf("\033[1;31mstates were simulated with DIFFERENT beam energies (%s %.3f vs %s %.3f) -- STOPPING\033[0m\n",
                results[0].name.Data(), results[0].EbeamFile, r.name.Data(), r.EbeamFile);
         return;
      }

   // ============================================================================================
   // WRITE THE OUTPUT: histograms at the top level, settings in the folder "settings/".
   // ============================================================================================
   const TString outPath = inMacroDir(outFile);
   gSystem->mkdir(gSystem->DirName(outPath), kTRUE);
   TFile fo(outPath, "RECREATE");
   if (fo.IsZombie()) {
      printf("\033[1;31mcannot write %s\033[0m\n", outPath.Data());
      return;
   }
   for (auto &r : results) {
      r.hGen->Write();
      r.hRec->Write();
      r.hAcc->Write();
      r.hLab->Write();
      r.hLabR->Write();
   }

   // Each setting is saved TWICE: as a number (TParameter, for code to read) and as a line in a
   // human-readable text (TNamed "summary", for you to read with ->Print()).
   TDirectory *dSet = fo.mkdir("settings");
   dSet->cd();
   TString summary;
   auto saveNum = [&](const char *name, double value, const char *meaning) {
      TParameter<double>(name, value).Write();
      summary += Form("%-22s = %-12.10g  %s\n", name, value, meaning);
   };
   auto saveText = [&](const char *name, const TString &value, const char *meaning) {
      TNamed(name, value.Data()).Write();
      summary += Form("%-22s = %s  (%s)\n", name, value.Data(), meaning);
   };
   TDatime now;
   summary += "---- provenance ----\n";
   saveText("macro", gInterpreter->GetCurrentMacroName(), "macro that wrote this file");
   saveText("date", now.AsSQLString(), "when it was written");
   saveText("seed", seed, "simulation set");
   saveText("simPattern", simPattern, "simulation file pattern");
   saveText("fitPattern", fitPattern, "fit file pattern");
   summary += "---- cuts (must equal the data selection) ----\n";
   saveNum("zMin", zMin, "mm, vertex window low (truth on gen, reco on rec)");
   saveNum("zMax", zMax, "mm, vertex window high");
   saveNum("chi2Cut", chi2Cut, "fitted chi2/ndf must be below");
   saveNum("keLo", keLo, "MeV, fitted KE must be above (data applies keOff -0.1 first)");
   saveNum("keHi", keHi, "MeV, fitted KE must be below");
   saveNum("useXtr", useXtr, "1 = GetKinematicsXtr (genfit), 0 = GetKinematics (UKF)");
   saveNum("dThetaMax", dThetaMax, "deg, truth match on angle (<=0 = off)");
   saveNum("keRatioMin", keRatioMin, "truth match: KE_fit/KE_true above");
   saveNum("keRatioMax", keRatioMax, "truth match: KE_fit/KE_true below");
   summary += "---- binning ----\n";
   saveNum("nBins", nBins, "theta_cm bins");
   saveNum("cmMax", cmMax, "deg, theta_cm upper edge (lower edge 0)");
   summary += "---- kinematics constants ----\n";
   saveNum("amu", kAmu, "MeV/u");
   saveNum("massC16", kMassC16, "MeV, beam");
   saveNum("massP", kMassP, "MeV, target");
   saveNum("massD", kMassD, "MeV, detected deuteron");
   saveNum("massC15", kMassC15, "MeV, residual ground state");
   saveNum("pdgDeuteron", kPdgDeuteron, "PDG code searched in MC truth");
   saveNum("ExTolerance", kExTolerance, "MeV, allowed Ex(file) - Ex(table)");
   saveNum("EbeamOverride", EbeamOverride, "MeV, <0 means: Ebeam from the sim file");
   summary += "---- per state ----\n";
   for (auto &r : results) {
      const char *n = r.name.Data();
      saveNum(Form("Ex_%s", n), r.ExFile, "MeV, read from the sim file");
      saveNum(Form("EbeamSim_%s", n), r.EbeamFile, "MeV, generation energy from the sim file");
      saveNum(Form("EbeamUsed_%s", n), r.EbeamUsed, "MeV, used in the kinematics");
      saveNum(Form("nEntries_%s", n), r.nEntries, "entries in sim and fit trees");
      saveNum(Form("nGen_%s", n), r.nGen, "denominator: true reactions in the window");
      saveNum(Form("nRec_%s", n), r.nRec, "numerator: matched good fits");
      saveText(Form("simFile_%s", n), r.simFile, "simulation file");
      saveText(Form("fitFile_%s", n), r.fitFile, "fit file");
   }
   TNamed("summary", summary.Data()).Write();
   fo.Close();

   printf("\n===== SETTINGS SAVED IN %s (settings/summary) =====\n%s", outPath.Data(), summary.Data());

   // ---- picture: all states on one plot -----------------------------------------------------
   TCanvas *c = new TCanvas("cAcc", "acceptance, all states", 900, 600);
   const int colors[] = {kAzure + 2, kOrange + 7, kGreen + 2, kMagenta + 1, kRed + 1};
   TLegend *l = new TLegend(0.65, 0.15, 0.88, 0.40);
   for (size_t k = 0; k < results.size(); ++k) {
      TH1D *h = results[k].hAcc;
      h->SetMinimum(0);
      h->SetMaximum(1.05);
      h->SetMarkerStyle(20);
      h->SetLineColor(colors[k % 5]);
      h->SetMarkerColor(colors[k % 5]);
      h->SetTitle("16C(p,d) acceptance, all states;#theta_{cm} [deg];acceptance");
      h->Draw(k == 0 ? "E1" : "E1 same");
      l->AddEntry(h, Form("%s (E_{x} = %.3f)", results[k].name.Data(), results[k].ExFile), "lp");
   }
   l->Draw();
   TString png = outPath;
   png.ReplaceAll(".root", ".png");
   c->SaveAs(png);
   printf("\nwrote %s and %s\n\n", outPath.Data(), png.Data());
}
