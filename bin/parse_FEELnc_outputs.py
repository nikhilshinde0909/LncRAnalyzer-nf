#!/usr/bin/env python3

import sys

if len(sys.argv) != 3:
    print("Usage: parse_FEELnc_outputs.py <Input_GTF> <Output_ids>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

def extract_attribute(attribute_string, attribute_name):
    for attribute in attribute_string.split(';'):
        attribute = attribute.strip()
        if attribute.startswith(attribute_name):
            return attribute.split('"')[1]
    return None

def extract_transcript_ids_from_gtf(GTF_file):
    results = set()

    with open(GTF_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue

            fields = line.rstrip().split('\t')
            if len(fields) < 9:
                continue

            attributes = fields[8]
            transcript_id = extract_attribute(attributes, "transcript_id")
            if transcript_id:
                results.add(transcript_id)

    return results

results = extract_transcript_ids_from_gtf(input_file)

with open(output_file, "w") as fout:
    for item in sorted(results):
        fout.write(item + "\n")

