#!/usr/bin/env python3
"""Read a simulation log made with AtTpc debug output and, for every deuteron, compare the kinetic
energy Geant4 gives it (E - sqrt(E^2 - p^2), from Geant4's own printout) with the intended one
(generator mass at the same momentum). Written 2026-11-18, see the simulation log, "Masses".

Make the log:  root -l -b -q -e 'FairLogger::GetLogger()->SetLogScreenLevel("debug");' \
                 'C16_pd_sim.C(20,2.,178.,"TGeant4",-28.5,"/tmp/test.root",0.0,7007,184.0)' > test.log
Read it:       python3 debug_deuteron_energies.py test.log
"""
import math, sys
M_GEN = 2.0141017781 * 0.93149401  # GeV, deuteron mass in the generator (atomic, macro's list)
for fn in sys.argv[1:]:
    lines = open(fn, errors="ignore").read().splitlines()
    print(f"=== {fn} ===\n{'p [MeV/c]':>10} {'mass tracked [MeV]':>19} {'T Geant4 [MeV]':>15} {'T intended':>11} {'diff [keV]':>11}")
    for i, l in enumerate(lines):
        if "Particle : deuteron" not in l:
            continue
        blk = lines[i:i + 15]
        mom = next((x for x in blk if "Momentum :" in x), None)
        et = next((x for x in blk if "Total relativistic energy" in x), None)
        if not (mom and et):
            continue
        px, py, pz = map(float, mom.split(":")[-1].split())
        E = float(et.split()[-1]); p2 = px * px + py * py + pz * pz; m = math.sqrt(E * E - p2)
        tG4 = (E - m) * 1000; tInt = (math.sqrt(p2 + M_GEN ** 2) - M_GEN) * 1000
        print(f"{math.sqrt(p2)*1000:10.2f} {m*1000:19.2f} {tG4:15.3f} {tInt:11.3f} {(tG4-tInt)*1000:+11.1f}")
