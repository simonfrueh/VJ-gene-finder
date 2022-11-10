# -*- coding: UTF-8 -*-
########################################################################
# VDJ Annotation
# algorithm: 	  Simon Früh
# implementation: Martin Früh
########################################################################

from Bio.Seq import Seq
from Bio import SeqIO
import modules.Init as mI
import modules.V as mV
import modules.J as mJ
import modules.Results as mR

# parse and process command line arguments
path, filename = mI.parse_arguments()

# Loop through input records
n_seq_record=0
for seq_record in SeqIO.parse(path + filename, "fasta"):
	
	# print current record
	n_seq_record += 1
	print("\nRecord\t#%i\nHeader:\t%s\nSeq:\t%s\nLength:\t%i" 
	      % (n_seq_record, seq_record.id, repr(seq_record.seq), 
	      len(seq_record)))
	print("\nProcessing data")
	
	# prepare DNA strand and reverse complement
	seq_data = seq_record.seq.upper()
	seq_data_rc = seq_data.reverse_complement()

	# prepare total number of tasks and result list
	tasks = 8
	task_i = 0
	result_list = []
	
	# identify V search candidates by RSS motif CAC
	task_i+= 1
	print("(%i/%i)" % (task_i, tasks) + " identify V search candidates"
	    + " by RSS motif CAC")
	rss, rss_rc = mV.ident_V_rss_motif(seq_data, seq_data_rc)

	# V segments with single-exon leader peptide (TRAV1 and TRGV1)
	task_i+= 1
	print("(%i/%i)" % (task_i, tasks) + " search for V segments with"
	      + " single-exon leader peptide (TRAV1 and TRGV1)")
	mV.task_V1(seq_data, rss, "", result_list)

	# V segments with single-exon leader peptide (TRAV1 and TRGV1) (RC)
	task_i+= 1
	print("(%i/%i)" % (task_i, tasks) + " search for V segments with"
	      + " single-exon leader peptide (TRAV1 and TRGV1) in reverse"
	      + " complement (RC)")
	mV.task_V1(seq_data_rc, rss_rc, "RC", result_list)

	# V segments with two-exon leader peptide (all others)  
	task_i+= 1
	print("(%i/%i)" % (task_i, tasks) + " search for V segments with"
	      + " two-exon leader peptide")
	mV.task_V2(seq_data, rss, "", result_list)
			
	# V segments with two-exon leader peptide (all others) (RC)
	task_i+= 1
	print("(%i/%i)" % (task_i, tasks) + " search for V segments with"
	      + " two-exon leader peptide in reverse complement (RC)")
	mV.task_V2(seq_data_rc, rss_rc, "RC", result_list)
	
	# identify J search candidates by RSS motif GTG
	task_i+= 1
	print("(%i/%i)" % (task_i, tasks) + " identify J search candidates"
	    + " by RSS motif GTG")
	rss, rss_rc = mJ.ident_J_rss_motif(seq_data, seq_data_rc)

	# J segments
	task_i+= 1
	print("(%i/%i)" % (task_i, tasks) + " search for J segments")
	mJ.task_J(seq_data, rss, "", result_list)
	
	# J segments (RC)
	task_i+= 1
	print("(%i/%i)" % (task_i, tasks) + " search for J in reverse"
	      + " complement (RC)")
	mJ.task_J(seq_data_rc, rss_rc, "RC", result_list)
	
	# Results	
	mR.show_results(result_list)
	mR.write_results(path, seq_record.id, result_list)

	print("...finished\n")
