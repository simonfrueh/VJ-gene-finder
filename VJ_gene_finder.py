#!/usr/bin/python3
########################################################################
# VJ-gene-finder
# algorithm: Simon Frueh
# implementation: Martin Frueh
########################################################################

from Bio import SeqIO

import modules.init as i
import modules.functions as f
import modules.j as j
import modules.v as v
import modules.results as r

# Parse and process command line arguments
path, filename, output_dir = i.parse_arguments()

# Loop through input records
n_seq_record = 0
for seq_record in SeqIO.parse(path + filename, "fasta"):
    # Print current record
    n_seq_record += 1
    print("\nRecord\t#%i\nHeader:\t%s\nSeq:\t%s\nLength:\t%i"
          % (n_seq_record, seq_record.id, repr(seq_record.seq),
             len(seq_record)))
    print("\nProcessing data")

    # Prepare DNA strand and reverse complement
    seq_data = seq_record.seq.upper()
    seq_data_rc = seq_data.reverse_complement()

    # Prepare total number of tasks and result list
    tasks = 8
    task_i = 0
    result_list = []

    # Identify V search candidates by RSS motif CAC
    task_i += 1
    print("(%i/%i)" % (task_i, tasks) + " identify V search candidates"
          + " by RSS motif CAC")
    rss = f.ident_rss_motif_start_position(seq_data, "CAC")
    rss_rc = f.ident_rss_motif_start_position(seq_data_rc, "CAC")

    # Search V segments with single-exon leader peptide (TRAV1 and TRGV1)
    task_i += 1
    print("(%i/%i)" % (task_i, tasks) + " search for V segments with"
          + " single-exon leader peptide")
    v.task_v1(seq_data, rss, False, result_list)

    # Search V segments with single-exon leader peptide (TRAV1 and TRGV1) (RC)
    task_i += 1
    print("(%i/%i)" % (task_i, tasks) + " search for V segments with"
          + " single-exon leader peptide in reverse complement (RC)")
    v.task_v1(seq_data_rc, rss_rc, True, result_list)

    # Search V segments with two-exon leader peptide (all others)
    task_i += 1
    print("(%i/%i)" % (task_i, tasks) + " search for V segments with"
          + " two-exon leader peptide")
    v.task_v2(seq_data, rss, False, result_list)

    # Search V segments with two-exon leader peptide (all others) (RC)
    task_i += 1
    print("(%i/%i)" % (task_i, tasks) + " search for V segments with"
          + " two-exon leader peptide in reverse complement (RC)")
    v.task_v2(seq_data_rc, rss_rc, True, result_list)

    # Identify J search candidates by RSS motif GTG
    task_i += 1
    print("(%i/%i)" % (task_i, tasks) + " identify J search candidates"
          + " by RSS motif GTG")
    rss = f.ident_rss_motif_end_position(seq_data, "GTG")
    rss_rc = f.ident_rss_motif_end_position(seq_data_rc, "GTG")

    # Search J segments
    task_i += 1
    print("(%i/%i)" % (task_i, tasks) + " search for J segments")
    j.task_j(seq_data, rss, False, result_list)

    # Search J segments (RC)
    task_i += 1
    print("(%i/%i)" % (task_i, tasks) + " search for J in reverse"
          + " complement (RC)")
    j.task_j(seq_data_rc, rss_rc, True, result_list)

    # Show results and write result files
    r.show_results(result_list)
    r.write_results(output_dir, seq_record.id, result_list)

    print("...finished\n")
