/// @file run_reco_a1975_sim.C
/// @brief Digitise + reconstruct an a1975 simulation ((p,d), (p,p), (d,t)) with a configurable
///        PSA and track finder, so the chain can be set to whatever the data reconstruction used.
///
/// 2026-11-20. The new a1975 data reconstruction uses AtPSAMultiFit (the exact parameters are
/// still to be received). Until then the defaults are the MultiFit settings used for the run-106
/// PSA comparison (management/psa_cmp_2026-11-20/unpackReco_psa.C): threshold 20, peaking time
/// 0.720 us, 4 peaks, min separation 4, prominence on primary off. Every one of them is an
/// argument so the data's values can be dropped in without editing the macro.
///
///   psaType  "multifit" (default) | "max"
///   praType  "tc" (default, TC 15 / 7.5 as the run-106 reco) | "hdbscan"
///   holeR    30 mm -- the a1975 beam-hole radius, definitive.
///
/// No cleaning stage, no space-charge correction (the sim has no space charge).
///
///   root -b -q -l 'run_reco_a1975_sim.C("data/gs_s1001_sim.root","out.root","ATTPC.a1975_C16_sim.par",2000)'
#include "../AtBeamHole.h"

void run_reco_a1975_sim(TString mcFile = "./data/attpcsim.root", TString outputFile = "./data/sim_reco.root",
                        TString paramFile = "ATTPC.a1975_C16_sim.par", Int_t nEvents = 0, Int_t evStart = 0,
                        TString psaType = "multifit", Double_t thr = 20, Double_t peakingTime = 0.720,
                        Int_t maxPeaks = 4, Int_t minSep = 4, TString praType = "tc",
                        Double_t clusterRadius = 15.0, Double_t clusterDistance = 7.5, Int_t hdMcs = 20,
                        Int_t hdMs = 8, Double_t holeR = 30.0)
{
   TString scriptfile = "Lookup20150611.xml";
   TString dir = getenv("VMCWORKDIR");
   TString digiParFile = dir + "/parameters/" + paramFile;
   TString mapParFile = dir + "/scripts/" + scriptfile;

   TStopwatch timer;

   FairRunAna *fRun = new FairRunAna();
   FairFileSource *source = new FairFileSource(mcFile);
   fRun->SetSource(source);
   fRun->SetOutputFile(outputFile);

   FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
   FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
   parIo1->open(digiParFile.Data(), "in");
   rtdb->setFirstInput(parIo1);

   auto mapping = std::make_shared<AtTpcMap>();
   mapping->ParseXMLMap(mapParFile.Data());
   mapping->GeneratePadPlane();
   AtSim::InhibitBeamHole(mapping, holeR);

   AtClusterizeTask *clusterizer = new AtClusterizeTask();
   clusterizer->SetPersistence(kFALSE);

   AtPulseTask *pulse = new AtPulseTask(std::make_shared<AtPulse>(mapping));
   pulse->SetPersistence(kFALSE); // traces not needed downstream; MC info still propagated
   pulse->SetSaveMCInfo();

   AtPSA *psa = nullptr;
   if (psaType == "max") {
      auto p = new AtPSAMax();
      p->SetThreshold(thr);
      psa = p;
      std::cout << "PSA   : AtPSAMax thr " << thr << std::endl;
   } else {
      auto p = new AtPSAMultiFit();
      p->SetThreshold(thr);
      p->SetPeakingTime(peakingTime);
      p->SetMaxPeaks(maxPeaks);
      p->SetMinSeparation(minSep);
      p->SetProminenceOnPrimary(false);
      psa = p;
      std::cout << "PSA   : AtPSAMultiFit thr " << thr << " tau " << peakingTime << " maxPeaks " << maxPeaks
                << " minSep " << minSep << std::endl;
   }
   AtPSAtask *psaTask = new AtPSAtask(psa);
   psaTask->SetPersistence(kTRUE);
   psaTask->SetOutputBranch("AtEventH");

   std::unique_ptr<AtPATTERN::AtPRA> praAlgo;
   if (praType == "hdbscan") {
      auto p = std::make_unique<AtPATTERN::AtTrackFinderHDBSCAN>();
      p->SetMinClusterSize(hdMcs);
      p->SetMinSamples(hdMs);
      p->SetClusterSelectionEpsilon(10.0);
      p->SetJoinMethod("mover");
      p->SetMinClusterSizeJoin(15);
      p->SetCircleOverlapRatio(0.25);
      p->SetMotionGapTol(40);
      p->SetMotionAngleTol(35);
      praAlgo = std::move(p);
      std::cout << "PRA   : AtTrackFinderHDBSCAN mcs " << hdMcs << " ms " << hdMs << std::endl;
   } else {
      auto p = std::make_unique<AtPATTERN::AtTrackFinderTC>();
      p->SetClusterRadius(clusterRadius);
      p->SetClusterDistance(clusterDistance);
      praAlgo = std::move(p);
      std::cout << "PRA   : AtTrackFinderTC r " << clusterRadius << " d " << clusterDistance << std::endl;
   }
   std::cout << "HOLE  : " << holeR << " mm   PAR : " << paramFile << std::endl;

   AtPRAtask *praTask = new AtPRAtask(std::move(praAlgo));
   praTask->SetInputBranch("AtEventH");
   praTask->SetOutputBranch("AtPatternEvent");
   praTask->SetPersistence(kTRUE);

   fRun->AddTask(clusterizer);
   fRun->AddTask(pulse);
   fRun->AddTask(psaTask);
   fRun->AddTask(praTask);

   fRun->Init();
   timer.Start();
   if (nEvents > 0)
      fRun->Run(evStart, evStart + nEvents);
   else
      fRun->Run();
   timer.Stop();

   std::cout << "\nsim reco done -> " << outputFile << "   par = " << paramFile << std::endl;
   std::cout << "Real " << timer.RealTime() << " s, CPU " << timer.CpuTime() << " s\n" << std::endl;
}
