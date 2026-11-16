#!/usr/bin/env bash
# Digitise + reconstruct the generated (p,d) states, in parallel.
#   ./run_reco_states.sh [seed] [outDir] [inDir] [states...]
#   ./run_reco_states.sh                       # the four states at once
#   ./run_reco_states.sh s1001 data/reco data ex2   # a single state
set -eo pipefail
SEED=${1:-s1001}; OUT=${2:-data/reco}; IN=${3:-data}
shift 3 2>/dev/null || shift $# 
STATES=${*:-gs ex1 ex2 ex3}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$HERE"
mkdir -p "$OUT"

reco_one() {
  local TAG=$1 J="${1}_${SEED}"
  [ -f "$OUT/${J}_reco.marker" ] && { echo "$J reco already COMPLETED"; return 0; }
  # only reconstruct what the generation step finished
  [ -f "$IN/$J.marker" ] && [ -s "$IN/${J}_sim.root" ] || { echo "  $J SIM_MISSING_OR_UNFINISHED"; return 1; }
  echo "[$(date +%H:%M:%S)] $J: reconstructing"
  root -b -q -l "run_reco_C16pd.C(\"$IN/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"ATTPC.a1975_C16_sim.par\",20,20,8,0,20.0,\"mover\")" \
       > "$OUT/${J}_reco.log" 2>&1 || true
  if [ -s "$OUT/${J}_reco.root" ] && grep -q "sim reco done" "$OUT/${J}_reco.log"; then
    echo COMPLETED > "$OUT/${J}_reco.marker"; echo "[$(date +%H:%M:%S)] $J done"
  else
    echo "  $J RECO_FAILED -- see $OUT/${J}_reco.log"; return 1
  fi
}

pids=()
for TAG in $STATES; do reco_one "$TAG" & pids+=($!); done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo "[$(date +%H:%M:%S)] all reco finished (fail=$fail)"
exit $fail