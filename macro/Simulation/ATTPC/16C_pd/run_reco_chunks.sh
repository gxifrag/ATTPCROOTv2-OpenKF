#!/usr/bin/env bash
# Digitise + reconstruct the generated (p,d) states in CHUNKS, in parallel, until a deadline.
#   ./run_reco_chunks.sh "<deadline>" [seed] [outDir] [inDir] [states...]
#   ./run_reco_chunks.sh "2026-11-14 09:30" s1001 /media/georgina/Aurio/Georgina/recoSim16Cpd data
# Env: CHUNK = events per chunk (default 10000), NTOT = events per sim file (default 500000)
# Stop early without losing anything:  touch <outDir>/STOP   (running chunks finish, no new ones start)
set -eo pipefail
DEADLINE_STR=${1:?usage: $0 "<deadline, e.g. '2026-11-14 09:30'>" [seed] [outDir] [inDir] [states...]}
DEADLINE=$(date -d "$DEADLINE_STR" +%s)
SEED=${2:-s1001}; OUT=${3:-data/reco}; IN=${4:-data}
shift 4 2>/dev/null || shift $#
STATES=${*:-gs ex1 ex2 ex3}
CHUNK=${CHUNK:-10000}; NTOT=${NTOT:-500000}

if [ "$DEADLINE" -le "$(date +%s)" ]; then echo "Deadline $DEADLINE_STR is in the past"; exit 1; fi
echo "Deadline: $(date -d @$DEADLINE)   chunk=$CHUNK events   states: $STATES"

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$HERE"
case "$OUT" in /media/*) mountpoint -q "$(echo "$OUT" | cut -d/ -f1-4)" || { echo "Drive not mounted"; exit 1; };; esac
mkdir -p "$OUT"

worker() {
  local TAG=$1 J="${1}_${SEED}" START=0 LAST=0 N C F NOW T0
  if [ ! -f "$IN/$J.marker" ] || [ ! -s "$IN/${J}_sim.root" ]; then
    echo "  $J SIM_MISSING_OR_UNFINISHED"; return 1
  fi
  while [ "$START" -lt "$NTOT" ]; do
    N=$(( NTOT - START < CHUNK ? NTOT - START : CHUNK ))
    C=$(printf "c%06d" "$START"); F="$OUT/${J}_${C}_reco"
    if [ -f "$F.marker" ]; then START=$((START+N)); continue; fi      # already done: skip
    if [ -f "$OUT/STOP" ]; then echo "[$(date +%H:%M:%S)] $J: STOP file, stopping at event $START"; return 0; fi
    NOW=$(date +%s)
    # do not start a chunk that (judging by the previous one, +10%) would not finish before the deadline
    if [ "$NOW" -ge "$DEADLINE" ] || { [ "$LAST" -gt 0 ] && [ $(( NOW + LAST + LAST/10 )) -gt "$DEADLINE" ]; }; then
      echo "[$(date +%H:%M:%S)] $J: deadline reached, stopping at event $START"; return 0
    fi
    echo "[$(date +%H:%M:%S)] $J: events $START..$((START+N-1))"
    T0=$(date +%s)
    root -b -q -l "run_reco_C16pd.C(\"$IN/${J}_sim.root\",\"$F.root\",\"ATTPC.a1975_C16_sim.par\",20,20,8,$N,20.0,\"mover\",$START)" \
         > "$F.log" 2>&1 || true
    if [ -s "$F.root" ] && grep -q "sim reco done" "$F.log"; then
      echo COMPLETED > "$F.marker"; LAST=$(( $(date +%s) - T0 ))
      echo "[$(date +%H:%M:%S)] $J: chunk $C done in ${LAST}s"
    else
      echo "  $J chunk $C RECO_FAILED -- see $F.log"; return 1
    fi
    START=$((START+N))
  done
  echo "[$(date +%H:%M:%S)] $J: ALL $NTOT events done"
}

pids=()
for TAG in $STATES; do worker "$TAG" & pids+=($!); done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
echo "[$(date +%H:%M:%S)] finished (fail=$fail)"
exit $fail