#this script was adapted from:
#Wick et al., 2017. Completing bacterial genome assemblies with multiplex MinION sequencing
#https://github.com/rrwick/Bacterial-genome-assemblies-with-multiplex-MinION-sequencing/tree/master/error_rate_estimation

import sys

query_lines = []

total_length = 0
total_error_count = 0
total_mismatches = 0
total_gaps = 0

with open(sys.argv[1], 'rt') as blast_hits:
    for line in blast_hits:
        line_parts = line.split('\t')
        percent_identity = float(line_parts[2])
        alignment_length = int(line_parts[3])
        mismatches = int(line_parts[4])
        gaps = int(line_parts[5])

        if percent_identity < 95.0:
            continue
        if alignment_length < 4000:
            continue

        error_count = round(alignment_length * (100.0 - percent_identity) / 100.0)

        total_length += alignment_length
        total_error_count += error_count
        total_mismatches += mismatches
        total_gaps += gaps

overall_percent_identity = 100.0 * (1.0 - (total_error_count / total_length))
total_mismatches = str(total_mismatches)
total_gaps = str(total_gaps)
total_length = str(total_length)

print('%.4f' % overall_percent_identity + '\t' + total_mismatches + '\t' + total_gaps + '\t' + total_length, flush=True, end='\n')
