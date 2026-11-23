#include <TFile.h>
#include <TParameter.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

/// @file C16_pd_sim.C
/// @brief a1975 16C(p,d)15C simulation, matched to experimental running conditions.
///
/// Purpose: Acceptance of the deuteron channel measured per FINAL STATE. The (p,d)
/// analysis extracts four states of 15C (ground state, 0.740, 3.10, and 4.66 MeV).
/// Each state has distinct two-body kinematics, shifting theta_lab and deuteron energy.
///
/// EXPERIMENTAL CONDITIONS & BEAM ENERGETICS (a1975, proton-target runs):
///     - Magnetic field: B = 2.85 T (field convention: bFieldkG = -28.5 kG)
///     - Target gas: H2 at 300 torr, room temperature (293 K), rho = 3.308e-5 g/cm3 
///       using geometry file ATTPC_H300torr_RT.root.
///     - Active drift length: 1000 mm.
///
/// BEAM ENERGY PROFILE (16C on 1H2):
///     - Accelerator beam energy (before IC):           196.0 MeV
///     - Beam energy after IC (2 windows + IC gas):     188.34 MeV
///     - Beam energy range inside AT-TPC active volume: 186.39 MeV (entrance after tpc_window) 
///       down to 175.44 MeV (exit at z = 1000 mm).
///     - Default simulation input (Ebeam at entrance):  186.4 MeV (subject to fine-tuning 
///       and calibration against Ex spectra).
///     
/// NOTE: on Generator: The ion generator multiplies pz by A, so pz is in GeV/c per nucleon.
/// NOTE: to match the excitation energy of the dround state in the analysis, it was decided that the beam after the IC but BEFORE the tpc_window was set to 184 MeV.

///   root -l 'C16_pd_sim.C(2000)'                      // ground state
///   root -l 'C16_pd_sim.C(8000,2.,178.,"TGeant4",-28.5,"./data/pd_ex1.root",0.740,1001,184.)'
///
/// MASS CHECK (added 2026-11-18) -- see the simulation log, "Masses: generator vs FairRoot vs Geant4".
/// FairPrimaryGenerator computes every pushed particle's total energy E from its momentum and the
/// mass in TDatabasePDG; Geant4 then starts it with T = E - m (G4PrimaryParticle::Set4Momentum).
/// If the momentum was computed with a DIFFERENT mass, Geant4 starts the particle at a different
/// kinetic energy than intended. This prints, for every pushed particle, the mass used for its
/// momentum, the TDatabasePDG mass, and the resulting shift of the starting kinetic energy.
/// Returns false if the BEAM shift is above 1 keV: the beam is fixed in this macro and must stay so.
/// The PRODUCTS cannot be fixed here (their momenta come from AtTPC2Body); they are only reported.
static bool checkPushedMasses(double Ebeam, double mBeamUsedMeV, const std::vector<Int_t> &Zp,
                              const std::vector<Int_t> &Ap, const std::vector<Double_t> &massAmu,
                              const std::vector<Double_t> &ExMeV)
{
   const double u = 931.49401; // MeV/u, the same factor AtTPC2Body uses
   auto pdgOf = [](int Z, int A) { return 1000000000 + Z * 10000 + A * 10; };
   auto mTDB = [](int pdg) {
      auto *p = TDatabasePDG::Instance()->GetParticle(pdg);
      return p ? p->Mass() * 1000. : -1.; // MeV
   };
   auto tG4 = [](double T, double mUsed, double mFair) { // Geant4's starting T for intended T
      double p = std::sqrt(T * T + 2 * T * mUsed);
      return std::sqrt(p * p + mFair * mFair) - mFair;
   };
   printf("\n===== MASS CHECK: mass used for the momentum vs mass FairRoot uses for E =====\n");
   printf("  %-9s %11s %13s %13s %10s   %s\n", "particle", "pdg", "used [MeV]", "TDatabasePDG", "diff [keV]",
          "Geant4 T - intended T");
   const int pb = pdgOf(Zp[0], Ap[0]);
   const double mb = mTDB(pb), dBeam = (tG4(Ebeam, mBeamUsedMeV, mb) - Ebeam) * 1000;
   printf("  %-9s %11d %13.4f %13.4f %10.1f   %+.3f keV at T = %.2f MeV\n", "beam", pb, mBeamUsedMeV, mb,
          (mBeamUsedMeV - mb) * 1000, dBeam, Ebeam);
   for (size_t i = 2; i < Zp.size(); ++i) { // 0 = beam, 1 = target (not pushed), 2.. = products
      const int pdg = pdgOf(Zp[i], Ap[i]);
      const double mUsed = massAmu[i] * u + ExMeV[i], mf = mTDB(pdg);
      printf("  %-9s %11d %13.4f %13.4f %10.1f   %+.1f keV at T = 10 MeV, %+.1f keV at T = 30 MeV\n",
             Form("product%zu", i - 1), pdg, mUsed, mf, (mUsed - mf) * 1000, (tG4(10, mUsed, mf) - 10) * 1000,
             (tG4(30, mUsed, mf) - 30) * 1000);
   }
   const bool beamOK = std::fabs(dBeam) < 1.0;
   if (beamOK)
      printf("\033[1;32m  beam OK: Geant4 starts it within 1 keV of Ebeam.\033[0m\n");
   else
      printf("\033[1;31m  BEAM MASS MISMATCH: Geant4 would start the beam %.1f keV away from Ebeam.\033[0m\n", dBeam);
   printf("\033[1;33m  Products: any nonzero shift above is NOT corrected by this macro (AtTPC2Body pushes only the\n"
          "  momentum). See the simulation log.\033[0m\n\n");
   return beamOK;
}

void C16_pd_sim(Int_t nEvents = 20000, Double_t thetaMinCM = 2.0, Double_t thetaMaxCM = 178.0,
                TString mcEngine = "TGeant4", Double_t bFieldkG = -28.5, TString outFile = "./data/attpcsim.root",
                Double_t resEx = 0.0, UInt_t seed = 0, Double_t Ebeam = 184.0) //192.0
{
   // Parallel jobs with no seed produce byte-identical events, so "more statistics" would be the
   // same sample copied. seed = 0 keeps ROOT's time-based default; pass a distinct value per job.
   if (seed != 0)
      gRandom->SetSeed(seed);
   // Print the REQUESTED seed: gRandom->GetSeed() returns TRandom3's internal state counter, not
   // the seed, so it cannot be used to verify that a parallel run was actually seeded.
   std::cout << "RNG seed requested: " << seed << std::endl;
   TString dir = getenv("VMCWORKDIR");
   TString parFile = "./data/attpcpar.root";

   TStopwatch timer;
   timer.Start();

   FairRunSim *run = new FairRunSim();
   run->SetName(mcEngine);
   run->SetOutputFile(outFile);
   FairRuntimeDb *rtdb = run->GetRuntimeDb();

   run->SetMaterials("media.geo");

   FairModule *cave = new AtCave("CAVE");
   cave->SetGeometryFileName("cave.geo");
   run->AddModule(cave);

   FairDetector *ATTPC = new AtTpc("ATTPC", kTRUE);
   ATTPC->SetGeometryFileName("ATTPC_H300torr_RT.root"); // H2 at 300 torr, rho = 3.308e-5 g/cm3 (293 K)
   run->AddModule(ATTPC);

   // Field sign follows the a1954 lesson: generate in the DATA convention so that the sense of
   // rotation matches, otherwise anything inferring direction from curvature (AtSpyralPID) rejects
   // most of the simulated tracks and the acceptance measures the sign error instead of the
   // detector.
   AtConstField *fMagField = new AtConstField();
   fMagField->SetField(0., 0., bFieldkG);                 // kG
   fMagField->SetFieldRegion(-50, 50, -50, 50, -10, 230); // cm
   run->SetField(fMagField);

   FairPrimaryGenerator *primGen = new FairPrimaryGenerator();

   // -----   Beam : 16C   ---------------------------------------------------
   Int_t z = 6;  // atomic number
   Int_t a = 16; // mass number
   // CHARGE STATE: keep q = 0 together with ATOMIC beam masses. Tested 2026-11-18 (simulation log,
   // "Beam charge state"): q does not change the beam's energy loss, trajectory or vertex, but Geant4
   // tracks the beam with ITS OWN nucleus mass + (Z - q) electrons (it ignores Bmass), while AtTpc
   // computes the vertex energy as Etot - Bmass*u. q = 6 with atomic Bmass: reaction beam 3.07 MeV too
   // low (momentum off by 19.9 MeV/c). q = 6 with nucleus Bmass: off by 0.25 MeV/c, because the beam
   // mass in the reaction list (Mass.push_back below) is still atomic. To use q = Z, BOTH Bmass and the
   // beam entry of the mass list must be nucleus masses.
   Int_t q = 0;  // charge state
   Int_t m = 1;  // multiplicity
   // MASS FOR THE BEAM MOMENTUM (fixed 2026-11-18). FairPrimaryGenerator computes the beam's total
   // energy with the TDatabasePDG mass, which for 16C is the NUCLEUS mass (atomic - 6 electrons),
   // and Geant4 starts the beam with T = E - m. So the momentum must be computed with that same
   // nucleus mass; with the atomic mass Geant4 started the beam at 184.037 MeV instead of 184.000.
   // Bmass below stays ATOMIC: AtTpc uses it to compute the beam energy at the vertex, and Geant4
   // tracks the q = 0 beam with the atomic mass. checkPushedMasses() verifies this after Init().
   const Double_t u = 931.49401, me = 0.51099895;               // MeV
   const Double_t mBeamAtomicMeV = 16.0147 * u;                 // atomic mass (neutral 16C)
   const Double_t mBeamNucleusMeV = mBeamAtomicMeV - z * me;    // nucleus mass = what FairRoot uses
   Double_t pTot = std::sqrt(Ebeam * Ebeam + 2 * Ebeam * mBeamNucleusMeV) / 1000.0; // GeV/c
   Double_t px = 0.000 / a;
   Double_t py = 0.000 / a;
   Double_t pz = pTot / a; // GeV/c per nucleon
   Double_t BExcEner = 0.0;
   Double_t Bmass = 16.0147;   // amu (repo convention -- see the note above)
   Double_t NomEnergy = Ebeam; // MeV

   // maxELoss CONTROLS THE REACTION RATE and is not cosmetic. AtTPCIonGenerator samples a
   // threshold uniformly in [0, maxELoss] and AtTpc::reactionOccursHere() fires only once the
   // accumulated loss exceeds it. Set it far above the real traversal loss and most events never
   // react; set it far below and every reaction happens near the entrance, so the vertex is no
   // longer uniform in z and the acceptance is measured on the wrong vertex distribution.
   //
   // MEASURED, not estimated: solving the elastic truth for the beam energy that gives Ex = 0
   // event by event showed 25.2 MeV of loss across the chamber when the gas was (wrongly) set to
   // 1 bar. Scaling by the density ratio 3.308/8.27 gives about 10 MeV at 300 torr. Re-verify the
   // same way after any change: the vertex must stay flat across the drift length, and the solved
   // beam energy should fall linearly from the entrance value.
   Double_t maxELoss = 11.; // MeV 10.49

   AtTPCIonGenerator *ionGen =
      new AtTPCIonGenerator("Ion", z, a, q, m, px, py, pz, BExcEner, Bmass, NomEnergy, maxELoss);
   ionGen->SetSpotRadius(0, -100, 0); //(r, z, offx) so -100 is initial z, and it is a point-like beam (r=0)
   primGen->AddGenerator(ionGen);

   // -----   Two-body reaction 16C(p,d)15C   --------------------------------
   std::vector<Int_t> Zp, Ap, Qp;
   std::vector<Double_t> Pxp, Pyp, Pzp, Mass, ExE;
   Int_t mult = 4;         // beam, target, heavy residual, light ejectile -- must be 4
   Double_t ResEner = 0.0; // unused

   // ---- Beam : 16C ----
   Zp.push_back(z);
   Ap.push_back(a);
   Qp.push_back(q);
   Pxp.push_back(px);
   Pyp.push_back(py);
   Pzp.push_back(pz);
   Mass.push_back(16.0147); // uma
   ExE.push_back(BExcEner);

   // ---- Target : p ----
   Zp.push_back(1);
   Ap.push_back(1);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(1.0078250322); // uma
   ExE.push_back(0.0);

   // ---- Heavy residual : 15C, left in the state under study ----
   Zp.push_back(6);
   Ap.push_back(15);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(15.0105993); // uma
   // The four states the (p,d) analysis extracts: 0 (g.s., 1/2+), 0.740 (5/2+), 3.10, 4.66.
   // Note that 3.10 and 4.66 lie above the neutron separation energy of 15C (1.218 MeV), so the
   // residual is unbound and will decay. That is not simulated here: the residual is transported
   // as a stable ion. It does not affect the DEUTERON kinematics, which is what the acceptance is
   // measured on, but it does mean this simulation must not be used to say anything about the
   // recoil or about coincidences with the neutron.
   ExE.push_back(resEx);

   // ---- Light ejectile : d, the particle that is detected ----
   Zp.push_back(1);
   Ap.push_back(2);
   Qp.push_back(0);
   Pxp.push_back(0.0);
   Pyp.push_back(0.0);
   Pzp.push_back(0.0);
   Mass.push_back(2.0141017781); // uma
   ExE.push_back(0.0);

   Double_t ThetaMinCMS = thetaMinCM;
   Double_t ThetaMaxCMS = thetaMaxCM;
   std::cout << "16C(p,d)15C  Ebeam = " << Ebeam << " MeV (" << Ebeam / 16 << " MeV/u),  p = " << pTot * 1000
             << " MeV/c" << std::endl;
   std::cout << "residual 15C excitation = " << resEx << " MeV" << std::endl;
   std::cout << "CM angular range: " << ThetaMinCMS << " - " << ThetaMaxCMS << " deg" << std::endl;

   AtTPC2Body *TwoBody = new AtTPC2Body("TwoBody", &Zp, &Ap, &Qp, mult, &Pxp, &Pyp, &Pzp, &Mass, &ExE, ResEner,
                                        ThetaMinCMS, ThetaMaxCMS);
   primGen->AddGenerator(TwoBody); //unreacted tracks, punch-through, are 'empty entries' in the exit channel (odd entries)

   run->SetGenerator(primGen);
   run->SetStoreTraj(kTRUE);

   run->Init();

   // Verify the beam momentum against FairRoot's mass and report the products (see checkPushedMasses).
   if (!checkPushedMasses(Ebeam, mBeamNucleusMeV, Zp, Ap, Mass, ExE)) {
      std::cout << "STOPPING: fix the beam mass before simulating." << std::endl;
      return;
   }

   Bool_t kParameterMerged = kTRUE;
   FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
   parOut->open(parFile.Data());
   rtdb->setOutput(parOut);
   rtdb->saveOutput();
   rtdb->print();

   run->Run(nEvents);
   
   // Stamp the simulated state into the output file, so the file says by itself what it contains.
   {
      TFile f(outFile, "UPDATE");
      if (f.IsZombie()) {
         std::cout << "WARNING: could not reopen " << outFile << " to store resEx" << std::endl;
      } else {
         TParameter<Double_t> pEx("resEx", resEx);
         TParameter<Double_t> pEb("Ebeam", Ebeam);
         TParameter<Double_t> pTmin("thetaMinCM", thetaMinCM);
         TParameter<Double_t> pTmax("thetaMaxCM", thetaMaxCM);
         TParameter<Int_t> pSeed("seed", (Int_t)seed);
         pEx.Write();
         pEb.Write();
         pTmin.Write();
         pTmax.Write();
         pSeed.Write();
         f.Close();
      }
   }

   run->CreateGeometryFile("./data/geofile_C16_pd_full.root");

   timer.Stop();
   std::cout << std::endl
             << "Macro finished successfully." << std::endl
             << "Output file is " << outFile << std::endl
             << "Real time " << timer.RealTime() << " s, CPU time " << timer.CpuTime() << " s" << std::endl;
}
