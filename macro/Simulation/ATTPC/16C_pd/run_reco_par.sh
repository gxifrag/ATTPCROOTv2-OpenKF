#!/usr/bin/env bash
# Digitise + reconstruct the generated chunks, in parallel. Resumable via .marker files.
#
#   ./run_reco_par.sh [simDir] [recoDir] [nParallel] [nChunksPerState]
set -eo pipefail
SIM=$(realpath -m "${1:-data_184}")
OUT=$(realpath -m "${2:-data_184/reco_184}")
NPAR=${3:-4}
NCH=${4:-10}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
mkdir -p "$OUT"

reco_job() {
  J=$1                                   # e.g. gs_c0_s1000
  [ -f "$OUT/${J}_reco.marker" ] && { echo "[$(date +%H:%M:%S)] $J already COMPLETED"; return; }
  echo "[$(date +%H:%M:%S)] $J: reconstructing"
  ( cd "$HERE" && root -b -q -l \
    "run_reco_C16pd.C(\"$SIM/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"ATTPC.a1975_C16_sim.par\",20,20,8,0,20.0,\"mover\")" \
    > "$OUT/${J}_reco.log" 2>&1 ) || true
  if [ -s "$OUT/${J}_reco.root" ] && grep -q "sim reco done" "$OUT/${J}_reco.log"; then
    echo COMPLETED > "$OUT/${J}_reco.marker"; echo "  $J done"
  else
    echo "  $J RECO_FAILED -- see $OUT/${J}_reco.log"
  fi
}
export -f reco_job
export SIM OUT HERE

# take the job names from the sim files that actually exist and finished generating
for TAG in gs ex1 ex2 ex3; do
  for ((c=0; c<NCH; c++)); do
    M=$(ls "$SIM"/${TAG}_c${c}_s*.marker 2>/dev/null | head -1) || true
    [ -n "$M" ] && basename "$M" .marker
  done
done | xargs -P "$NPAR" -I{} bash -c 'reco_job {}'

echo "[$(date +%H:%M:%S)] reco finished: $(ls "$OUT"/*_reco.marker 2>/dev/null | wc -l) markers"