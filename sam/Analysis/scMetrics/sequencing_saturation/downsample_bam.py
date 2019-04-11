import sys
from numpy.random import rand
import time

# Usage: 
# samtools view input.bam | python downsample_bam.py aggregate_saturation_curve.tsv reads_per_cell.tsv umis_per_cell.tsv

#################
# Inputs:
summary_filename = sys.argv[1]
read_table_filename = sys.argv[2]
umi_table_filename = sys.argv[3]
sample_rates = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

#################

# For each downsampling rate, keep track of unique (cell_barcode, gene, umi) 
# combinations. Each corresponds to a count in the final counts matrix.
umi_sets = {p: {} for p in sample_rates}

# Tally number of reads sampled for each downsampling rate.
read_counts = {p: {} for p in sample_rates}

t0 = time.time()
sampled_barcodes = set([])

for iR,read in enumerate(sys.stdin):
    # Loop over alignments in bam file input (streaming from stdin)
    if (iR+1) % 5e5 == 0:
        # print status update
        t1 = time.time()
        print('{} reads processed. Time for last 500,000: {:.2f} seconds'.format((iR+1), t1-t0))
        t0 = time.time()

    rand_value = rand()
    for column in read.strip('\n').split('\t'):
        # Find tags for cell barcode, UMI, and gene alignment
        if column.startswith('XB'):
            cell_barcode = column.split(':')[-1]
        elif column.startswith('XU'):
            umi = column.split(':')[-1]
        elif column.startswith('YG'):
            gene = column.split(':')[-1]

    for p in sample_rates:
        # Sample read for each input sampling rate
        if p > rand_value:
            if cell_barcode not in umi_sets[p]:
                umi_sets[p][cell_barcode] = set([])
                read_counts[p][cell_barcode] = 0
            umi_sets[p][cell_barcode].add((gene, umi))
            read_counts[p][cell_barcode] += 1
            sampled_barcodes.add(cell_barcode)

# Reads table:
# Each row is a cell barcode
# Each column is a different read-sampling rate
# Entry (i,j) is number of reads for barcode i and sampling rate j
read_output = open(read_table_filename, 'w')

# UMI table:
# Same as reads table, but for the number of detected UMIs
umi_output = open(umi_table_filename, 'w')

read_output.write('barcode\t' + '\t'.join(['p_' + str(p) for p in sample_rates]) + '\n')
umi_output.write('barcode\t' + '\t'.join(['p_' + str(p) for p in sample_rates]) + '\n')

for bc in sampled_barcodes:
    read_out_line = [bc]
    umi_out_line = [bc]
    for p in sample_rates:
        if bc in umi_sets[p]:
            umi_out_line.append(str(len(umi_sets[p][bc])))
            read_out_line.append(str(read_counts[p][bc]))
        else:
            umi_out_line.append('0')
            read_out_line.append('0')
    read_output.write('\t'.join(read_out_line) + '\n')
    umi_output.write('\t'.join(umi_out_line) + '\n')
read_output.close()
umi_output.close()


# Save overall saturation curve (total reads and UMIs across all cell barcodes
# for each read-sampling rate)
with open(summary_filename, 'w') as o:
    o.write('sample_rate\tn_reads\tn_umis\n')
    for p in sample_rates:
        total_reads = sum(list(read_counts[p].values()))
        total_umis = sum([len(v) for v in umi_sets[p].values()])
        o.write('{}\t{}\t{}\n'.format(p, total_reads, total_umis))
