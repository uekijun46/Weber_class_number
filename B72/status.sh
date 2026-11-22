#!/bin/bash
echo "certified primes: $(cat cert*.txt | grep -vc '^#')   running workers: $(pgrep -f gen_cert.py | wc -l)"
for f in cert_w*.txt; do [ -f "$f" ] && echo "$f: last prime $(tail -n 1 $f | cut -d';' -f1)"; done
