#!/bin/bash
# Usage: ./run_all.sh [HI] [WORKERS]      e.g.  ./run_all.sh 12000000 6
# Starts WORKERS background processes that certify all split primes up to HI.  Safe to re-run: finished primes are skipped.
# Set PYTHON to choose the interpreter, e.g.  PYTHON="sage -python" ./run_all.sh
HI=${1:-12000000}; N=${2:-4}; PY=${PYTHON:-python3}
KEEP=""; command -v caffeinate >/dev/null && KEEP="caffeinate -i"
[ -f bnf_B72.bin ] || $PY gen_cert.py 2 2 cert_init.txt        # builds the class group data once (about 3 minutes)
for i in $(seq 0 $((N-1))); do
  nohup $KEEP $PY gen_cert.py 2 $HI cert_w$i.txt 1e18 $N $i > log_w$i.txt 2>&1 &
done
echo "started $N workers up to $HI.  Progress:  ./status.sh     Final check:  $PY verify_cert.py cert*.txt"
