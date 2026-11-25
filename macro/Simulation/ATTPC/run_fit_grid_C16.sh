#!/usr/bin/env bash
# UKF fit + Ex resolution for the run_res_grid_C16.sh grid (16C(p,d) and 16C(d,t), 15C Ex grid).
#
# GenFit is not built on this machine, so the fit is the UKF, as in the a1975 UKF productions:
#   (p,d) deuteron: UKF/pipeline/fitUKF_a1975.C, settings of UKF/pd/ukf_reftrack_batch.sh
#         (measSigma 0.5, momSigmaFrac 0.1, 2 iterations, >=10 clusters, reference track, x4)
#         EXCEPT the gas density: 3.308e-5 g/cm3 (H2 300 torr, 293 K) instead of the macro's
#         9.0e-5 default, which is "H 1 bar" inherited from a1954.
#   (d,t) triton, TWO configurations, D2_UKF/fitUKF_a1975_deuterium.C:
#     _ukf_prod  the data production (D2_UKF/fit_dt_dv1104.sh): measSigma 2.0, momSigmaFrac 0.3,
#                1 iteration, >=4 clusters, density 6.61e-5, matA default (1).  On the g.s. test
#                this returns the triton KE ~+1.9 MeV high, growing with KE (2-4 MeV above 17 MeV).
#     _ukf_alt   measSigma 0.5, momSigmaFrac 0.1, >=10 clusters: bias +0.36 MeV.
#     The production one is what describes the data AS ANALYSED; the alt one shows what the fit
#     could give. Both are cheap (~9 min/point).
# Field sign -1 for both (the data's): +1 was tested on the g.s. and fails ((p,d) 0 tracks pass).
#
# Ex resolution: res_ex_C16.C with the analysis cuts and fixed beam energy; one RESULT line per
# (point, config) collected in $OUT/fit/summary.txt.
#
#   [ONLY=regex] ./run_fit_grid_C16.sh [nParallel=10] [outDir]
#   e.g. ONLY='^pd .* prod$' ./run_fit_grid_C16.sh 9
set -eo pipefail
NPAR=${1:-10}
OUT=${2:-/media/georgina/Aurio/Georgina/resolution_C16}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
A=$REPO/macro/Unpack_HDF5/a1975
PAR_PD=${PAR_PD:-ATTPC.a1975_C16_sim_g50000.par}
PAR_DT=${PAR_DT:-ATTPC.a1975_deuterium_dv1104_g35000.par}
# the fit macros read <name>_reco.root: link the gain-tagged reco files under that name
LNK=$OUT/fit/links
mkdir -p "$LNK/pd" "$LNK/dt" "$OUT/fit/pd" "$OUT/fit/dt" "$OUT/fit/logs" "$OUT/fit/res"
for f in "$OUT"/pd/work/*/*_reco_${PAR_PD%.par}.root; do b=$(basename "$f"); ln -sfn "$f" "$LNK/pd/${b%%_reco_*}_reco.root"; done
for f in "$OUT"/dt/work/*/*_reco_${PAR_DT%.par}.root; do b=$(basename "$f"); ln -sfn "$f" "$LNK/dt/${b%%_reco_*}_reco.root"; done
export HERE A OUT LNK

one() { # reaction tag Ex seed config
  local R=$1 TAG=$2 EX=$3 SEED=$4 C=$5 J="${2}_s${4}"
  local F LOG="$OUT/fit/logs/${J}_${C}.log"
  if [ "$R" = pd ]; then
    # std  = the (p,d) UKF production (ukf_reftrack_batch.sh): 0.5/0.1, 2 iter + reference track, >=10 cl
    # prod = the (d,t) production parameters (fit_dt_dv1104.sh) on (p,d), the cross-check of 2026-11-24:
    #        2.0/0.3, 1 iter, no reference track, >=4 cl. Output suffix _prod keeps std untouched.
    local P S
    if [ "$C" = prod ]; then P="2.0,0.3,1,4,\"$OUT/fit/pd/\",false,4.0"; S="_prod"
    else P="0.5,0.1,2,10,\"$OUT/fit/pd/\",true,4.0"; S=""; fi
    F="$OUT/fit/pd/${J}_ukf${S}.root"
    [ -e "$F.done" ] || { (cd "$A/UKF" && root -b -q -l "pipeline/fitUKF_a1975.C(\"$J\",-1,\"deuteron\",-1,2.85,3.308e-5,\"$S\",\"$LNK/pd/\",$P)" > "$LOG" 2>&1 </dev/null) || true; }
  else
    local P; [ "$C" = prod ] && P="2.0,0.3,1,4" || P="0.5,0.1,1,10"
    F="$OUT/fit/dt/${J}_ukf_${C}.root"
    [ -e "$F.done" ] || { (cd "$A/D2_UKF" && root -b -q -l "fitUKF_a1975_deuterium.C(\"$J\",-1,\"triton\",-1,2.85,6.61e-5,\"$LNK/dt/\",\"$OUT/fit/dt/\",$P,kTRUE,\"_reco\",\"_ukf_$C\")" > "$LOG" 2>&1 </dev/null) || true; }
  fi
  if [ ! -e "$F.done" ]; then
    grep -q 'segmentation violation' "$LOG" && { echo "  FAIL $R $J $C fit (segv) -- $LOG"; return 1; }
    [ -s "$F" ] || { echo "  FAIL $R $J $C fit (no output) -- $LOG"; return 1; }
    touch "$F.done"
  fi
  (cd "$HERE" && root -l -b -q "res_ex_C16.C(\"$OUT/$R/work/$J/${J}_sim.root\",\"$F\",\"$R\",$EX,\"$OUT/fit/res/${R}_${J}_${C}.root\")" 2>&1 | grep '^RESULT' | sed "s/^RESULT/RESULT $C/") || true
  echo "[$(date +%H:%M:%S)] $R $J $C done" >&2
}
export -f one

PTS=""
for c in std prod; do
  PTS="$PTS
pd gs 0.0 5001 $c
pd ex1 0.740 5002 $c
pd g1 1.0 5003 $c
pd g2 2.0 5004 $c
pd ex2 3.103 5005 $c
pd ex3 4.657 5006 $c
pd g6 6.0 5007 $c
pd g7 7.0 5008 $c
pd g8 8.0 5009 $c"
done
for c in prod alt; do
  PTS="$PTS
dt gs 0.0 6001 $c
dt ex1 0.740 6002 $c
dt g1 1.0 6003 $c
dt g2 2.0 6004 $c
dt ex2 3.103 6005 $c
dt ex3 4.657 6006 $c
dt g6 6.0 6007 $c
dt ex4 6.358 6008 $c
dt g7 7.0 6009 $c
dt g8 8.0 6010 $c"
done
PTS=$(echo "$PTS" | sed '/^$/d')
# ONLY: optional extended regex on the point lines ("pd gs 0.0 5001 prod"), to run a subset.
[ -n "${ONLY:-}" ] && PTS=$(echo "$PTS" | grep -E "$ONLY" || true)
[ -n "$PTS" ] || { echo "no points selected"; exit 1; }
echo "[$(date +%H:%M:%S)] fitting $(echo "$PTS" | wc -l) (point, config) pairs, $NPAR in parallel" >&2
# Results are MERGED into summary.txt: a new RESULT line replaces the old one with the same
# (config, reaction, Ex); lines of points not rerun are kept. So a subset run never wipes the rest.
SUM="$OUT/fit/summary.txt" NEW="$OUT/fit/summary.new.$$"
echo "$PTS" | xargs -P "$NPAR" -L 1 bash -c 'one "$@"' _ > "$NEW"
touch "$SUM"
awk 'NR==FNR { k[$2" "$3" "$5]=1; print; next } !(($2" "$3" "$5) in k)' "$NEW" "$SUM" |
  sort -k3,3 -k2,2 -k5,5g > "$SUM.tmp" && mv "$SUM.tmp" "$SUM" && rm -f "$NEW"
echo "[$(date +%H:%M:%S)] finished -> $SUM" >&2
cat "$SUM"
