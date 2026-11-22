#!/bin/sh
# Regenerate (optional) and verify all 13 fields.  Usage: sh rig_run_all.sh [gen]
cd "$(dirname "$0")"
if [ "$1" = gen ]; then
  for i in 0 1 2 3 4 5 6 7 8; do python3 rig_gen.py $i 100000 500; done
  for i in 9 10 11 12; do python3 rig_gen.py $i 1000000 560; done
fi
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12; do python3 rig_verify.py $i > rig_log_$i.txt 2>&1; grep RESULT rig_log_$i.txt || tail -3 rig_log_$i.txt; done | tee rig_results.txt
