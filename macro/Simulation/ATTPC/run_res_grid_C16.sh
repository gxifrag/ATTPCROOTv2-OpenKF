#!/usr/bin/env bash
# Generate + digitise + reconstruct (AtPSAMultiFit) 16C(p,d)15C and 16C(d,t)15C on a grid of 15C
# excitation energies, to map the experimental Ex resolution width across the Ex range.
#
# Points: the states the analysis extracts plus a grid in between, so sigma(Ex) can be
# interpolated at any energy, including unassigned features.
#   (p,d): 0, 0.740, 1, 2, 3.103, 4.657, 6, 7, 8 MeV
#   (d,t): the same plus 6.358 (literature energy of the feature the data puts at 6.648)
#
# GAIN. Generation does not depend on it, only reco does, so the stages can run separately
# (STAGE = gen | reco | all) and reco can be repeated at another gain without regenerating.
# The reco par is chosen with PAR_PD / PAR_DT and its name is part of the reco file name, so
# reconstructions at different gains coexist. Defaults are the gains matched to the OLD a1975
# data on median hits per track (new data not yet available):
#   (d,t): ATTPC.a1975_deuterium_dv1104_g35000.par  -- 75 vs 76 hits/track in run 0016 (Yassid's
#          scan, 16C_dt/run_reco_batch_dt.sh). The production par's Gain 10000 gives half that.
#   (p,d): see PAR_PD below -- the sim par's Gain 150000 gives 2.5x the hits of run 0106.
#
# Generated with the AtTPC2Body total-energy push (2026-11-18): files from before that date
# start the products with TDatabasePDG masses and must not be mixed with these.
#
#   [PAR_PD=...] [PAR_DT=...] ./run_res_grid_C16.sh <pd|dt|both> [nEvents=20000] [outDir] [nParallel=4] [stage=all]
set -eo pipefail
WHICH=${1:?usage: $0 <pd|dt|both> [nEvents] [outDir] [nParallel] [gen|reco|all]}
NEV=${2:-20000}
OUT=${3:-/media/georgina/Aurio/Georgina/resolution_C16}
NPAR=${4:-4}
STAGE=${5:-all}
PAR_PD=${PAR_PD:-ATTPC.a1975_C16_sim.par}
PAR_DT=${PAR_DT:-ATTPC.a1975_deuterium_dv1104_g35000.par}
case $STAGE in gen|reco|all) ;; *) echo "unknown stage $STAGE"; exit 1 ;; esac
export STAGE PAR_PD PAR_DT

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
export HERE NEV OUT

nent() { root -l -b -q -e "TFile f(\"$1\");auto*t=(TTree*)f.Get(\"cbmsim\");printf(\"NENT %lld\\n\",t?t->GetEntries():-1LL);" 2>/dev/null </dev/null | awk '/^NENT/{print $2}'; }
export -f nent

# One job = one reaction at one Ex: gen -> reco. Markers make reruns skip finished stages.
job() {
  local R=$1 TAG=$2 EX=$3 SEED=$4
  local D="$OUT/$R" J="${TAG}_s${SEED}"
  mkdir -p "$D/logs"
  local PAR; [ "$R" = pd ] && PAR=$PAR_PD || PAR=$PAR_DT
  # Every file of a point lives in its own work dir (layout since 2026-11-25): the generation
  # (_sim.root), the reconstruction (_reco_<par>.root), their .done markers, and the par/geometry
  # files the sim macro writes to ./data/.
  local W="$D/work/$J"
  local SIM="$W/${J}_sim.root" RECO="$W/${J}_reco_${PAR%.par}.root"
  if [ "$R" = pd ]; then
    local M=$HERE/16C_pd GEN="C16_pd_sim.C($NEV,2.,178.,\"TGeant4\",-28.5,\"$SIM\",$EX,$SEED)"
    # Same chain as the a1975 H2 data reconstructed with MultiFit (run 0106, psa_cmp_2026-11-20):
    # MultiFit thr 20 / 0.72 us / 4 / 4, prominence off, TC 15 / 7.5, hole 30 mm. The (p,d) gain is
    # matched on that chain, so the grid must use it too. Traces are not persisted by this macro.
    local RMAC="run_reco_a1975_sim.C(\"$SIM\",\"$RECO\",\"$PAR\")" RDONE="sim reco done"
  else
    local M=$HERE/16C_dt GEN="C16_dt_sim.C($NEV,2.,70.,\"TGeant4\",-28.5,\"$SIM\",$EX,$SEED)"
    local RMAC="run_reco_C16dt.C(\"$SIM\",\"$RECO\",\"$PAR\")" RDONE="(d,t) sim reco done"
  fi
  # The sim macros write ./data/attpcpar*.root and ./data/geofile_*.root relative to the cwd, so
  # parallel jobs sharing a cwd overwrite each other's par file mid-read and crash. One work dir
  # per job, as run_states_par.sh does.
  mkdir -p "$W/data"
  cd "$W"
  if [ "$STAGE" != reco ] && [ ! -e "$SIM.done" ]; then
    echo "[$(date +%H:%M:%S)] $R $J gen  Ex=$EX"
    root -b -q -l "$M/$GEN" > "$D/logs/${J}_gen.log" 2>&1 </dev/null || true
    grep -q "RNG seed requested: $SEED" "$D/logs/${J}_gen.log" &&
      grep -q "Macro finished successfully" "$D/logs/${J}_gen.log" &&
      [ "$(nent "$SIM")" = "$NEV" ] || { echo "  FAIL $R $J gen -- $D/logs/${J}_gen.log"; return 1; }
    touch "$SIM.done"
  fi
  if [ "$STAGE" != gen ] && [ ! -e "$RECO.done" ]; then
    [ -e "$SIM.done" ] || { echo "  FAIL $R $J reco -- no finished sim ($SIM.done missing)"; return 1; }
    local RLOG="$D/logs/${J}_reco_${PAR%.par}.log"
    echo "[$(date +%H:%M:%S)] $R $J reco  par=$PAR"
    root -b -q -l "$M/$RMAC" > "$RLOG" 2>&1 </dev/null || true
    grep -qF "$RDONE" "$RLOG" && grep -qF "par = $PAR" "$RLOG" &&
      [ "$(nent "$RECO")" = "$NEV" ] || { echo "  FAIL $R $J reco -- $RLOG"; return 1; }
    touch "$RECO.done"
  fi
  echo "[$(date +%H:%M:%S)] $R $J DONE ($STAGE)"
}
export -f job

# reaction tag Ex seed -- distinct seed per point so no two points share a random sequence
PTS_PD="pd gs 0.0 5001
pd ex1 0.740 5002
pd g1 1.0 5003
pd g2 2.0 5004
pd ex2 3.103 5005
pd ex3 4.657 5006
pd g6 6.0 5007
pd g7 7.0 5008
pd g8 8.0 5009"
PTS_DT="dt gs 0.0 6001
dt ex1 0.740 6002
dt g1 1.0 6003
dt g2 2.0 6004
dt ex2 3.103 6005
dt ex3 4.657 6006
dt g6 6.0 6007
dt ex4 6.358 6008
dt g7 7.0 6009
dt g8 8.0 6010"

case $WHICH in
  pd) PTS=$PTS_PD ;; dt) PTS=$PTS_DT ;; both) PTS="$PTS_PD"$'\n'"$PTS_DT" ;;
  *) echo "unknown reaction $WHICH"; exit 1 ;;
esac
# ONLY: optional extended regex on the point lines ("pd ex3 4.657 5006"), to run a subset, e.g.
# when another runner is already working on some points (two runners on one point would collide).
[ -n "${ONLY:-}" ] && PTS=$(echo "$PTS" | grep -E "$ONLY" || true)
[ -n "$PTS" ] || { echo "no points selected"; exit 1; }
mkdir -p "$OUT"
echo "[$(date +%H:%M:%S)] $WHICH stage=$STAGE: $NEV events/point, $NPAR in parallel -> $OUT"
[ "$STAGE" = gen ] || echo "  reco pars: pd=$PAR_PD dt=$PAR_DT"
echo "$PTS" | xargs -P "$NPAR" -L 1 bash -c 'job "$@"' _
echo "[$(date +%H:%M:%S)] finished; failures (if any) are listed above as FAIL"
