#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print("Usage: parse_pfam_results.py <PfamScan_results> <Output_hits>", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

results = set()

with open(input_file, "r") as fin:
    for line in fin:
        if line.startswith("#"):
            continue

        fields = line.split()
        if len(fields) < 3:
            continue

        results.add(fields[2].split("_", 1)[0])

with open(output_file, "w") as fout:
    for item in sorted(results):
        fout.write(item + "\n")
