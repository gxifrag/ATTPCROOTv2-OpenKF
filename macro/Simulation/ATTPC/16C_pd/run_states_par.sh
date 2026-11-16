#!/usr/bin/env bash
# 16C(p,d)15C generation, parallel chunks. One directory per job (attpcpar.root and the
# geofile are written relative to the cwd, so parallel jobs must not share one).
#
#   ./run_states_par.sh [nChunks] [evPerChunk] [nParallel] [outDir] [Ebeam] [seedBase]
set -eo pipefail
NCH=${1:-10}
NEV=${2:-100}
NPAR=${3:-4}
OUT=$(realpath -m "${4:-$HOME/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_pd/data_test}")
EBEAM=${5:-184.0}
SEEDBASE=${6:-1000}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
mkdir -p "$OUT"

run_job() {
  TAG=$1; EX=$2; SEED=$3; CH=$4
  J="${TAG}_c${CH}_s${SEED}"
  [ -f "$OUT/$J.marker" ] && { echo "[$(date +%H:%M:%S)] $J already COMPLETED"; return; }
  JOBDIR="$OUT/work/$J"; mkdir -p "$JOBDIR/data"
  echo "[$(date +%H:%M:%S)] $J start: Ex = $EX MeV, $NEV events"
  ( cd "$JOBDIR" && root -b -q -l \
    "$HERE/C16_pd_sim.C($NEV,2.,178.,\"TGeant4\",-28.5,\"$OUT/${J}_sim.root\",$EX,$SEED,$EBEAM)"\
    > "$OUT/${J}_gen.log" 2>&1 ) || true
  [ -s "$OUT/${J}_sim.root" ]                                  || { echo "  $J GEN_FAILED";       return; }
  grep -q "RNG seed requested: $SEED" "$OUT/${J}_gen.log"      || { echo "  $J SEED_NOT_APPLIED"; return; }
  grep -q "Macro finished successfully" "$OUT/${J}_gen.log"    || { echo "  $J DID_NOT_FINISH";   return; }
  echo COMPLETED > "$OUT/$J.marker"
  rm -rf "$JOBDIR"
  echo "[$(date +%H:%M:%S)] $J done ($(du -h "$OUT/${J}_sim.root" | cut -f1))"
}
export -f run_job
export HERE OUT NEV EBEAM
# state index * 100 + chunk keeps every (state, chunk) seed distinct
i=0
for CFG in "gs:0.0" "ex1:0.740" "ex2:3.100" "ex3:4.660"; do
  TAG=${CFG%%:*}; EX=${CFG##*:}
  for ((c=0; c<NCH; c++)); do
    echo "$TAG $EX $((SEEDBASE + i*100 + c)) $c"
  done
  i=$((i+1))
done | xargs -P "$NPAR" -L1 bash -c 'run_job $0 $1 $2 $3'

echo "[$(date +%H:%M:%S)] done. Markers: $(ls "$OUT"/*.marker 2>/dev/null | wc -l) / $((4*NCH))"