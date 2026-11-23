/// @file acceptance_C16dt_allstates.C
/// @brief 16C(d,t)15C acceptance vs theta_cm for ALL states in ONE run, with every cut and every
///        hard-coded number saved inside the output file.
///
/// ============================================================================================
/// WHY THIS FILE EXISTS (created 2026-11-18, copy of acceptance_C16dt.C, which is left untouched)
/// ============================================================================================
/// Same strategy as 16C_pd/acceptance_C16pd_allstates.C, so both channels are built the same way:
///
///   1. ALL states in a single run, ONE output file, one histogram per state (hAcc_gs, ...).
///      The original already looped over states, but SKIPPED a missing state silently; this one
///      stops and writes nothing ("all states or nothing").
///   2. Every cut and hard-coded parameter is saved in the file (folder "settings/"), including
///      the PID gate polygon itself.
///   3. Ex and Ebeam are read from each simulation file when it stores them. C16_dt_sim.C does NOT
///      store them yet (C16_pd_sim.C does) -- see "SIMULATION PARAMETERS" below.
///   4. The default cuts are the DATA cuts of a1975_analysis_KF/dt/macros/fit_dt_definitive.C.
///      The original defaults were NOT: it used vz 10-500 (data: 25-600), no theta_lab cut (data:
///      corrected theta_lab > 15 deg), no chi2 cut (data: chi2/ndf < 4.5) and Ebeam 184.25 (data: 184.0).
///
/// ============================================================================================
/// WHAT "ACCEPTANCE" MEANS HERE
/// ============================================================================================
///   acceptance(theta_cm) = reconstructed / generated       (per theta_cm bin, per state)
///
///   generated     (denominator) = simulated reactions with the TRUE triton inside the theta_lab
///                                 window and the TRUE vertex inside [vzLo, vzHi]. MC truth only.
///   reconstructed (numerator)   = of those, the events where
///                                   (a) a PATTERN track has a valid AtSpyralPID estimate,
///                                   (b) it is inside the PID gate (sqrt(dE/dx), Brho),
///                                   (c) its RECONSTRUCTED vertex z (un-mirrored) is inside [vzLo, vzHi],
///                                   (d) more than half of its hits are triton hits (truth match), and
///                                   (e) if a fit file is given: a FITTED track in the same event with
///                                       chi2/ndf < chi2Cut, KE in [keLo, keHi], matching the truth.
///   Both are binned at the TRUE theta_cm (see the original's explanation: the data yields are
///   already binned at the reconstructed angle, so binning the numerator there too would apply the
///   angular migration twice).
///
///   About (e): the original macro deliberately did NOT require a fit ("the fitter is a moving
///   target"). But the data DO require chi2/ndf < 4.5, so without (e) the fit efficiency (about 16%,
///   see xsec_dt.C) is missing from the correction. Give fitPattern = "" to switch (e) off; the
///   file then records fitRequired = 0 so nobody can mistake it for a full correction.
///
/// ============================================================================================
/// SIMULATION PARAMETERS
/// ============================================================================================
/// The macro wants to READ Ex and Ebeam from each _sim.root, so it cannot use a state that was
/// simulated differently from what we think. C16_pd_sim.C writes them as TParameter<Double_t>
/// "resEx" and "Ebeam"; C16_dt_sim.C should do the same. Until it does, run with
/// requireSimParams = kFALSE: the macro then uses the table below and EbeamDefault, and records
/// ExSource = "table" in the file so you can see it.
///
/// ============================================================================================
/// INPUT FILES (per state; "%s" is replaced by "<state>_<seed>", e.g. "gs_s3001")
/// ============================================================================================
///   simPattern   "data/%s_sim.root"       MC truth
///   recoPattern  "data/%s_reco.root"      pattern recognition (AtPatternEvent), UNGATED
///   fitPattern   "data/fit/%s_fit.root"   fitter output (AtTrackingEvent); "" = no fit requirement
/// All three must have one entry per simulated event (checked). Relative paths are relative to
/// THIS macro's folder. !! Check the patterns against where your files really are. !!
///
/// HOW TO RUN
///   root -b -q 'acceptance_C16dt_allstates.C'
///
/// HOW TO READ THE SAVED SETTINGS BACK
///   TFile f("data/acceptance_dt_allstates.root");
///   f.Get("settings/summary")->Print();                           // everything, readable
///   ((TParameter<double>*)f.Get("settings/vzLo"))->GetVal();      // one value, as a number
///   ((TCutG*)f.Get("settings/pidGate"))->Draw("l");               // the gate that was used

#include <TParameter.h>
#include <TNamed.h>
#include <TDatime.h>
#include <TCutG.h>
#include <fstream>
#include <tuple>
#include <vector>

// ---------------------------------------------------------------------------------------------
// THE STATES. Name used in file/histogram names, and the Ex (MeV) the state was simulated at.
// ex4 is generated at the literature 6.358 MeV (run_states_dt.sh); the data put the feature at
// about 6.65-6.70 and fit_dt_definitive.C interpolates the acceptance with 6.700 for this state.
// ---------------------------------------------------------------------------------------------
struct AccState {
   const char *name;
   double ExExpected;
};
static const AccState kStates[] = {
   {"gs", 0.000}, {"ex1", 0.740}, {"ex2", 3.103}, {"ex3", 4.657}, {"ex4", 6.358},
};
static const int kNStates = sizeof(kStates) / sizeof(kStates[0]);

// ---------------------------------------------------------------------------------------------
// HARD-CODED CONSTANTS -- all saved in the output file.
// ---------------------------------------------------------------------------------------------
static const double kAmu = 931.49401;              // MeV per atomic mass unit
static const double kMassC16 = 16.0147013 * kAmu;  // beam, 16C
static const double kMassD = 2.0135532 * kAmu;     // target, deuteron
static const double kMassT = 3.01550072 * kAmu;    // detected particle, triton
static const double kMassC15 = 15.0105993 * kAmu;  // residual, 15C
static const int kPdgTriton = 1000010030;          // how Geant4 labels a triton
static const double kBField = 2.85;                // T, given to AtSpyralPID
static const double kTruthHitFraction = 0.5;       // pattern track: > this fraction of triton hits
static const double kExTolerance = 0.005;          // MeV, allowed |Ex(file) - Ex(table)|

static double dt_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics -> {Ex, theta_cm [deg]}, identical to the original macro
static std::tuple<double, double> dt_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                          double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * dt_omega2(s, m1 * m1, m2 * m2) * dt_omega2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (dt_omega2(s, m1 * m1, m2 * m2) * dt_omega2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

static double readSimParam(TFile *f, const char *name)
{
   auto *p = (TParameter<Double_t> *)f->Get(name);
   return p ? p->GetVal() : std::numeric_limits<double>::quiet_NaN();
}

// =============================================================================================
// THE MACRO. Defaults of the CUTS = the data selection of fit_dt_definitive.C (2026-11-18).
// If you change the data selection there, change it here too and rerun.
// =============================================================================================
void acceptance_C16dt_allstates(
   TString seed = "s3001",
   TString simPattern = "data/%s_sim.root",
   TString recoPattern = "data/%s_reco.root",
   TString fitPattern = "data/fit/%s_fit.root", // "" = do not require a fit (records fitRequired = 0)
   TString outFile = "data/acceptance_dt_allstates.root",
   TString gateFile = "triton_dt_sim.json",     // PID gate drawn on the SIMULATED plane
   Double_t vzLo = 25.0, Double_t vzHi = 600.0, // mm   (= fit_dt_definitive.C vzLo, vzHi)
   Double_t thLabMin = 15.0, Double_t thLabMax = 180.0, // deg, on TRUE theta_lab (= thLo, thHi)
   Double_t chi2Cut = 4.5,                      // fitted chi2/ndf below (= chi2Max)
   Double_t keLo = 0.0, Double_t keHi = 1000.0, // MeV, fitted KE window (= keLo, keHi)
   Bool_t useXtr = kFALSE,                      // UKF: kFALSE; genfit: kTRUE (see the (p,d) macro)
   Double_t dThetaMax = 10.0,                   // deg, fit truth match on angle (<= 0 disables)
   Double_t keRatioMin = 0.5, Double_t keRatioMax = 2.0, // fit truth match on KE_fit / KE_true
   // THE SIMULATION'S DRIFT z IS MIRRORED against the experiment (measured on gs_s3001:
   // z_true + z_reco = 981.8 mm). The reconstructed z is un-mirrored before the vertex cut.
   Double_t zMirror = 981.8,
   Int_t nBins = 28, Double_t cmMax = 70.0,     // theta_cm binning
   Bool_t requireSimParams = kTRUE,             // stop if a sim file does not store resEx/Ebeam
   Double_t EbeamDefault = 184.0)               // MeV, used ONLY if the sim file has no Ebeam
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);

   const TString macroDir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   auto inMacroDir = [&](TString p) {
      if (p.BeginsWith("/") || p.BeginsWith("~")) {
         gSystem->ExpandPathName(p);
         return p;
      }
      return macroDir + "/" + p;
   };
   const bool fitRequired = !fitPattern.IsNull();

   // ---- the PID gate, read from its JSON (same parser as the original) ----------------------
   const TString gatePath = inMacroDir(gateFile);
   std::ifstream in(gatePath.Data());
   if (!in) {
      printf("\033[1;31mcannot open gate %s -- STOPPING\033[0m\n", gatePath.Data());
      return;
   }
   std::string all((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
   std::vector<double> vx, vy;
   size_t vp = all.find("\"vertices\"");
   size_t p = (vp == std::string::npos) ? std::string::npos : all.find('[', vp + 10);
   while (p != std::string::npos) {
      size_t a = all.find('[', p + 1);
      if (a == std::string::npos) break;
      size_t b = all.find(']', a);
      if (b == std::string::npos) break;
      double gx, gy;
      if (sscanf(all.substr(a + 1, b - a - 1).c_str(), "%lf , %lf", &gx, &gy) == 2) {
         vx.push_back(gx);
         vy.push_back(gy);
      }
      p = b;
      size_t nxt = all.find_first_not_of(" \n\r\t,", b + 1);
      if (nxt == std::string::npos || all[nxt] == ']') break;
   }
   if (vx.size() < 3) {
      printf("\033[1;31mgate has only %zu vertices -- STOPPING\033[0m\n", vx.size());
      return;
   }
   TCutG gate("pidGate", vx.size());
   gate.SetTitle("PID gate used in the numerator;sqrt(dE/dx);B#rho");
   for (size_t i = 0; i < vx.size(); ++i)
      gate.SetPoint(i, vx[i], vy[i]);

   AtTools::AtSpyralPID spy;
   spy.SetBField(kBField);

   struct StateResult {
      TString name, simFile, recoFile, fitFile, ExSource;
      double ExUsed = 0, EbeamUsed = 0;
      long nEntries = 0, nGen = 0, nRec = 0, nRecPattern = 0;
      TH1D *hGen = nullptr, *hRec = nullptr, *hAcc = nullptr;
   };
   std::vector<StateResult> results;

   for (int s = 0; s < kNStates; ++s) {
      StateResult r;
      r.name = kStates[s].name;
      const TString job = r.name + "_" + seed;
      r.simFile = inMacroDir(Form(simPattern.Data(), job.Data()));
      r.recoFile = inMacroDir(Form(recoPattern.Data(), job.Data()));
      r.fitFile = fitRequired ? inMacroDir(Form(fitPattern.Data(), job.Data())) : TString("(none)");

      TFile *Fs = TFile::Open(r.simFile), *Fr = TFile::Open(r.recoFile);
      TFile *Ff = fitRequired ? TFile::Open(r.fitFile) : nullptr;
      if (!Fs || Fs->IsZombie() || !Fr || Fr->IsZombie() || (fitRequired && (!Ff || Ff->IsZombie()))) {
         printf("\033[1;31m[%s] cannot open one of %s, %s, %s -- STOPPING, no file written\033[0m\n",
                r.name.Data(), r.simFile.Data(), r.recoFile.Data(), r.fitFile.Data());
         return;
      }

      // ---- what was simulated: from the file if stored, otherwise the table (recorded) -------
      double exFile = readSimParam(Fs, "resEx"), ebFile = readSimParam(Fs, "Ebeam");
      if (std::isnan(exFile) || std::isnan(ebFile)) {
         if (requireSimParams) {
            printf("\033[1;31m[%s] %s does not store resEx/Ebeam. Add them in C16_dt_sim.C (as C16_pd_sim.C "
                   "does), or rerun with requireSimParams = kFALSE. STOPPING\033[0m\n",
                   r.name.Data(), r.simFile.Data());
            return;
         }
         r.ExUsed = kStates[s].ExExpected;
         r.EbeamUsed = EbeamDefault;
         r.ExSource = "table (sim file has no resEx/Ebeam)";
         printf("\033[1;33m[%s] WARNING: sim file has no resEx/Ebeam; using the table Ex %.3f and Ebeam %.3f\033[0m\n",
                r.name.Data(), r.ExUsed, r.EbeamUsed);
      } else {
         if (std::fabs(exFile - kStates[s].ExExpected) > kExTolerance) {
            printf("\033[1;31m[%s] simulated at Ex = %.3f MeV, the table expects %.3f -- STOPPING\033[0m\n",
                   r.name.Data(), exFile, kStates[s].ExExpected);
            return;
         }
         r.ExUsed = exFile;
         r.EbeamUsed = ebFile;
         r.ExSource = "sim file";
      }

      TTree *ts = (TTree *)Fs->Get("cbmsim");
      TTree *tr = (TTree *)Fr->Get("cbmsim");
      TTree *tf = fitRequired ? (TTree *)Ff->Get("cbmsim") : nullptr;
      if (!ts || !tr || (fitRequired && !tf)) {
         printf("\033[1;31m[%s] missing cbmsim tree -- STOPPING\033[0m\n", r.name.Data());
         return;
      }
      // ONE ENTRY PER GENERATED EVENT IN EVERY FILE, or the event-by-event match is a guess.
      if (ts->GetEntries() != tr->GetEntries() || (tf && ts->GetEntries() != tf->GetEntries())) {
         printf("\033[1;31m[%s] ENTRY MISMATCH: sim %lld, reco %lld, fit %lld -- STOPPING\033[0m\n", r.name.Data(),
                ts->GetEntries(), tr->GetEntries(), tf ? tf->GetEntries() : -1LL);
         return;
      }
      r.nEntries = ts->GetEntries();

      TClonesArray *mc = nullptr, *pe = nullptr, *te = nullptr;
      ts->SetBranchAddress("MCTrack", &mc);
      tr->SetBranchAddress("AtPatternEvent", &pe);
      if (tf)
         tf->SetBranchAddress("AtTrackingEvent", &te);

      r.hGen = new TH1D("hGen_" + r.name, "generated;#theta_{cm} [deg];reactions", nBins, 0, cmMax);
      r.hRec = new TH1D("hRec_" + r.name, "reconstructed;#theta_{cm} [deg];reactions", nBins, 0, cmMax);
      for (TH1D *h : {r.hGen, r.hRec}) {
         h->SetDirectory(nullptr); // keep in memory after the input files close
         h->Sumw2();
      }

      for (Long64_t i = 0; i < r.nEntries; ++i) {
         ts->GetEntry(i);
         if (!mc)
            continue;
         // ---- truth: the triton, found by SPECIES (never by index) ----------------------------
         double cmT = -1, thLabT = -1, zTrue = -1e9, keT = -1, thT = -1;
         for (int k = 0; k < mc->GetEntries(); ++k) {
            auto *t = (AtMCTrack *)mc->At(k);
            if (!t || t->GetPdgCode() != kPdgTriton)
               continue;
            zTrue = t->GetStartZ() * 10.0; // cm -> mm
            double px = t->GetPx() * 1000, py = t->GetPy() * 1000, pz = t->GetPz() * 1000; // GeV -> MeV
            double pm = std::sqrt(px * px + py * py + pz * pz);
            if (pm <= 0)
               continue;
            keT = std::sqrt(pm * pm + kMassT * kMassT) - kMassT;
            thT = std::acos(pz / pm);
            auto [ex, cm] = dt_kine(kMassC16, kMassD, kMassT, kMassC15 + r.ExUsed, r.EbeamUsed, thT, keT);
            if (!std::isfinite(cm))
               continue;
            cmT = cm;
            thLabT = thT * TMath::RadToDeg();
            break;
         }
         if (cmT < 0)
            continue;
         // CUT on TRUE theta_lab: the data's theta_lab window, applied to truth on both sides.
         if (thLabT < thLabMin || thLabT > thLabMax)
            continue;
         // CUT on TRUE vertex z: the denominator counts only reactions inside the window.
         if (zTrue < vzLo || zTrue > vzHi)
            continue;
         ++r.nGen;
         r.hGen->Fill(cmT);

         // ---- numerator, part 1: a gated, Spyral-valid, truth-matched PATTERN track ------------
         tr->GetEntry(i);
         if (!pe || pe->GetEntries() == 0)
            continue;
         auto *P = (AtPatternEvent *)pe->At(0);
         if (!P)
            continue;
         bool patternOK = false;
         for (auto &trk : P->GetTrackCand()) {
            AtTrack &T = const_cast<AtTrack &>(trk);
            auto est = spy.Estimate(T);
            if (!est.valid)
               continue;
            if (!gate.IsInside(est.sqrtdEdx, est.brho))
               continue;
            // CUT on RECONSTRUCTED vertex z, un-mirrored first (see zMirror above).
            double zr = (zMirror > 0) ? (zMirror - est.vertex.Z()) : est.vertex.Z();
            if (zr < vzLo || zr > vzHi)
               continue;
            int nT = 0, nTot = 0;
            for (const auto &hit : T.GetHitArray()) {
               if (!hit)
                  continue;
               const auto &mcs = hit->GetMCSimPointArray();
               if (mcs.empty())
                  continue;
               ++nTot;
               if (mcs[0].A == 3 && mcs[0].Z == 1)
                  ++nT;
            }
            if (nTot > 0 && nT > kTruthHitFraction * nTot) {
               patternOK = true;
               break;
            }
         }
         if (!patternOK)
            continue;
         ++r.nRecPattern;

         // ---- numerator, part 2 (if a fit file is given): the data's fit cuts ------------------
         if (fitRequired) {
            tf->GetEntry(i);
            bool fitOK = false;
            if (te && te->GetEntriesFast() > 0) {
               auto *ev = (AtTrackingEvent *)te->At(0);
               if (ev)
                  for (auto &ft : ev->GetFittedTracks()) {
                     if (!ft)
                        continue;
                     const auto &md = ft->GetTrackMetadata();
                     double ndf = md ? md->GetNdf() : 0, chi2 = md ? md->GetChi2() : 0;
                     double c2n = ndf > 0 ? chi2 / ndf : 1e9;
                     const auto &kin = useXtr ? ft->GetKinematicsXtr() : ft->GetKinematics();
                     double ke = kin.kineticEnergy;
                     // CUT: chi2/ndf below chi2Cut and fitted KE inside [keLo, keHi] (the data cuts)
                     if (!(c2n < chi2Cut && ke > keLo && ke < keHi))
                        continue;
                     // TRUTH MATCH of the fitted track to the true triton
                     if (dThetaMax > 0) {
                        double dth = std::fabs(kin.theta - thT) * TMath::RadToDeg();
                        double ratio = keT > 0 ? ke / keT : -1;
                        if (dth > dThetaMax || ratio < keRatioMin || ratio > keRatioMax)
                           continue;
                     }
                     fitOK = true;
                     break;
                  }
            }
            if (!fitOK)
               continue;
         }
         ++r.nRec;
         r.hRec->Fill(cmT); // binned at the TRUE angle
      }
      Fs->Close();
      Fr->Close();
      if (Ff)
         Ff->Close();

      if (!r.nGen) {
         printf("\033[1;31m[%s] no true tritons inside the windows -- STOPPING\033[0m\n", r.name.Data());
         return;
      }
      r.hAcc = (TH1D *)r.hRec->Clone("hAcc_" + r.name);
      r.hAcc->SetDirectory(nullptr);
      r.hAcc->Divide(r.hRec, r.hGen, 1, 1, "B"); // binomial errors
      r.hAcc->SetTitle(Form("16C(d,t) acceptance, %s, Ex = %.3f MeV;#theta_{cm} [deg];acceptance", r.name.Data(),
                            r.ExUsed));
      printf("[%s] Ex %.3f (%s), Ebeam %.3f: generated %ld, pattern+gate %ld, accepted %ld -> %.3f\n",
             r.name.Data(), r.ExUsed, r.ExSource.Data(), r.EbeamUsed, r.nGen, r.nRecPattern, r.nRec,
             double(r.nRec) / r.nGen);
      results.push_back(r);
   }

   for (auto &r : results)
      if (std::fabs(r.EbeamUsed - results[0].EbeamUsed) > 1e-6) {
         printf("\033[1;31mstates have DIFFERENT beam energies (%s %.3f vs %s %.3f) -- STOPPING\033[0m\n",
                results[0].name.Data(), results[0].EbeamUsed, r.name.Data(), r.EbeamUsed);
         return;
      }

   // ============================================================================================
   // WRITE THE OUTPUT
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
   }

   TDirectory *dSet = fo.mkdir("settings");
   dSet->cd();
   gate.Write(); // the gate polygon itself
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
   saveText("recoPattern", recoPattern, "pattern-recognition file pattern");
   saveText("fitPattern", fitRequired ? fitPattern : TString("(none)"), "fit file pattern");
   saveText("gateFile", gatePath, "PID gate JSON (polygon saved as settings/pidGate)");
   summary += "---- cuts (must equal the data selection of fit_dt_definitive.C) ----\n";
   saveNum("vzLo", vzLo, "mm, vertex window low (truth on gen, un-mirrored reco on rec)");
   saveNum("vzHi", vzHi, "mm, vertex window high");
   saveNum("thLabMin", thLabMin, "deg, TRUE theta_lab low (both sides)");
   saveNum("thLabMax", thLabMax, "deg, TRUE theta_lab high (both sides)");
   saveNum("fitRequired", fitRequired, "1 = numerator requires a fit passing chi2/KE cuts");
   saveNum("chi2Cut", chi2Cut, fitRequired ? "fitted chi2/ndf below" : "NOT APPLIED (no fit file)");
   saveNum("keLo", keLo, fitRequired ? "MeV, fitted KE above" : "NOT APPLIED (no fit file)");
   saveNum("keHi", keHi, fitRequired ? "MeV, fitted KE below" : "NOT APPLIED (no fit file)");
   saveNum("useXtr", useXtr, "1 = GetKinematicsXtr (genfit), 0 = GetKinematics (UKF)");
   saveNum("dThetaMax", dThetaMax, "deg, fit truth match on angle (<=0 = off)");
   saveNum("keRatioMin", keRatioMin, "fit truth match: KE_fit/KE_true above");
   saveNum("keRatioMax", keRatioMax, "fit truth match: KE_fit/KE_true below");
   saveNum("truthHitFraction", kTruthHitFraction, "pattern track: triton hits above this fraction");
   saveNum("zMirror", zMirror, "mm, z_reco is un-mirrored as zMirror - z (0 = off)");
   summary += "---- binning ----\n";
   saveNum("nBins", nBins, "theta_cm bins");
   saveNum("cmMax", cmMax, "deg, theta_cm upper edge (lower edge 0)");
   summary += "---- constants ----\n";
   saveNum("amu", kAmu, "MeV/u");
   saveNum("massC16", kMassC16, "MeV, beam");
   saveNum("massD", kMassD, "MeV, target");
   saveNum("massT", kMassT, "MeV, detected triton");
   saveNum("massC15", kMassC15, "MeV, residual ground state");
   saveNum("pdgTriton", kPdgTriton, "PDG code searched in MC truth");
   saveNum("BField", kBField, "T, for AtSpyralPID");
   saveNum("ExTolerance", kExTolerance, "MeV, allowed Ex(file) - Ex(table)");
   saveNum("requireSimParams", requireSimParams, "1 = stop if sim file lacks resEx/Ebeam");
   saveNum("EbeamDefault", EbeamDefault, "MeV, used only if the sim file has no Ebeam");
   summary += "---- per state ----\n";
   for (auto &r : results) {
      const char *n = r.name.Data();
      saveNum(Form("Ex_%s", n), r.ExUsed, "MeV, used in the kinematics");
      saveText(Form("ExSource_%s", n), r.ExSource, "where Ex and Ebeam came from");
      saveNum(Form("EbeamUsed_%s", n), r.EbeamUsed, "MeV, used in the kinematics");
      saveNum(Form("nEntries_%s", n), r.nEntries, "entries in each input tree");
      saveNum(Form("nGen_%s", n), r.nGen, "denominator");
      saveNum(Form("nRecPattern_%s", n), r.nRecPattern, "pattern + gate + vz + truth match");
      saveNum(Form("nRec_%s", n), r.nRec, "numerator (also passing the fit cuts if required)");
      saveText(Form("simFile_%s", n), r.simFile, "simulation file");
      saveText(Form("recoFile_%s", n), r.recoFile, "pattern-recognition file");
      saveText(Form("fitFile_%s", n), r.fitFile, "fit file");
   }
   TNamed("summary", summary.Data()).Write();
   fo.Close();

   printf("\n===== SETTINGS SAVED IN %s (settings/summary) =====\n%s", outPath.Data(), summary.Data());
   if (!fitRequired)
      printf("\033[1;31mNOTE: no fit requirement -- the data's chi2/ndf < %.2f cut is NOT in this acceptance.\033[0m\n",
             chi2Cut);

   TCanvas *c = new TCanvas("cAcc", "acceptance, all states", 900, 600);
   const int colors[] = {kBlack, kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1};
   TLegend *l = new TLegend(0.14, 0.14, 0.44, 0.40);
   l->SetBorderSize(0);
   for (size_t k = 0; k < results.size(); ++k) {
      TH1D *h = results[k].hAcc;
      h->SetMinimum(0);
      h->SetMaximum(1.05);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(0.8);
      h->SetLineColor(colors[k % 5]);
      h->SetMarkerColor(colors[k % 5]);
      h->SetTitle("16C(d,t) acceptance, all states;#theta_{cm} [deg];acceptance");
      h->Draw(k == 0 ? "E1" : "E1 same");
      l->AddEntry(h, Form("%s (E_{x} = %.3f)", results[k].name.Data(), results[k].ExUsed), "lp");
   }
   l->Draw();
   TString png = outPath;
   png.ReplaceAll(".root", ".png");
   c->SaveAs(png);
   printf("\nwrote %s and %s\n\n", outPath.Data(), png.Data());
}
