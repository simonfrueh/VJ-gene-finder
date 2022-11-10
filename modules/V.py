from Bio.Seq import Seq
from Bio import motifs


# check motif criteria for V gene part 1
# omega_aa: translated amino acid sequence of search region
def check_V1_motif_criteria(omega_aa):
	
	# omega_nn does not contain a stop codon TAA, TAG or TGA
	# or omega_aa does not contain *
	if omega_aa.count("*") > 0:
			return False
	
	# Length Ωaa > 102 and < 161
	if (len(omega_aa) <= 102 or
		len(omega_aa) >= 161):
			return False
	
	# z1 = Cys/C between 32 - 82 (Ωaa)
	count_C = omega_aa[32-1:82].count("C")
	if count_C == 0:
		return False
	
	# z2 = Trp/W between 44 - 97 (Ωaa)
	count_W = omega_aa[44-1:97].count("W")
	if count_W == 0:
		return False
		
	# z3 = Y* motif in the last 12 amino acids
	# Ωaa (Y* = YYC, YFC, YLC, YHC, YIC, TFC)
	count_Y = omega_aa[-12:len(omega_aa)].count("YYC")
	count_Y += omega_aa[-12:len(omega_aa)].count("YFC")
	count_Y += omega_aa[-12:len(omega_aa)].count("YLC")
	count_Y += omega_aa[-12:len(omega_aa)].count("YHC")
	count_Y += omega_aa[-12:len(omega_aa)].count("YIC")
	count_Y += omega_aa[-12:len(omega_aa)].count("TFC")
	if count_Y == 0:
		return False
	
	return True


# check motif criteria for V gene part 2
# omega_aa: translated amino acid sequence of search region
def check_V2_motif_criteria(omega_aa):
	
	# omega_nn does not contain a stop codon TAA, TAG or TGA
	# or omega_aa does not contain *
	if omega_aa.count("*") > 0:
			return False
	
	# Length Ωaa > 85 and < 115
	if (len(omega_aa) <= 85 or
		len(omega_aa) >= 115):
			return False
	
	# z1 = Cys/C between 19 - 33 (Ωaa )
	count_C = omega_aa[19-1:33].count("C")
	if count_C == 0:
		return False
	
	# z2 = Trp/W between 31 - 48 (Ωaa )
	count_W = omega_aa[31-1:48].count("W")
	if count_W == 0:
		return False
		
	# z3 = Y* motif in the last 12 amino acids
	# Ωaa (Y* = YYC, YFC, YLC, YHC, YIC, TFC)
	count_Y = omega_aa[-12:len(omega_aa)].count("YYC")
	count_Y += omega_aa[-12:len(omega_aa)].count("YFC")
	count_Y += omega_aa[-12:len(omega_aa)].count("YLC")
	count_Y += omega_aa[-12:len(omega_aa)].count("YHC")
	count_Y += omega_aa[-12:len(omega_aa)].count("YIC")
	count_Y += omega_aa[-12:len(omega_aa)].count("TFC")
	if count_Y == 0:
		return False
	
	return True


# assign TRgroup based on amino acid motif
# omega_aa: amino acid sequence
# n: first n amino acids in omega_aa to be searched
# tr_list: [[group name, amino acid motiv],[...]]
def assign_tr_group(omega_aa, n, tr_list):
	tr_group = ""
	
	for t in tr_list:
		if omega_aa[0:n].count(t[1]) > 0:
			tr_group += " " + t[0]
			
	return tr_group.lstrip()
	

# check whether search region defined by s and p2 overlaps with 
# previous results from V part 1
def overlaps_with_prev_result(s, p2, rc, result_list):
	overlap = False
	
	for r in result_list:
		
		# only check for overlapt non_rc <-> non_rc and rc <-> rc
		if r[4] == rc:
			
			# s: r[5], p2: r[7]
			if s >= r[5] and s < r[7]:
				overlap = True
			if p2 > r[5] and p2 <= r[7]:
				overlap = True
	
	return overlap


# search and evaluate candidates for V gene part 1
# seq: Bio.Seq DNA sequence
# p1, p2: range of seq to examine
# offset: offset after start codon
# rc: "RC" / "" is reverse complement or not
# result_list: list to append results
def search_V1_motif(seq, p1, p2, offset, rc, result_list):
	
	candidates = []
	codon = "ATG"
	# ATG only relevant within rss_i + 163 nucleotides
	for s in range(p1, p1+163):
		if seq[s:s+len(codon)] == codon:
			candidates.append(s)
	
	for s in candidates:
		
		# limit to multiple of 3
		omega_nn = seq[s+offset:p2-(p2-(s+offset))%3]
		
		# translate with NCBI standard table
		omega_aa = omega_nn.translate(table = 1)
		
		if check_V1_motif_criteria(omega_aa) == False:
			# criteria not fulfilled, continue
			continue
			
		else:
			# motiv criteria fulfilled
			# assign TRGroup
			tr_list = [["TRAV1", "QVQQ"], 
			           ["TRGV1", "QVLLQQ"]]
			tr_group = assign_tr_group(omega_aa, len(omega_aa), tr_list)
			
			if tr_group == "": tr_group = "TRV-SEL"
			
			result_list.append([omega_nn, 
								omega_aa,
								seq[p2:p2+39],
								tr_group,
								rc, 
								s, 
								p1, 
								p2,
								"V1"])
							
			# skip ATG/AG and rssi between p1 and p2 + 39 nucleotides 
			# and continue at step 2 for rssi+1
			min_next_rss = p2 + 39
			return min_next_rss
			
	return 0


# search and evaluate candidates for V gene part 2
# seq: Bio.Seq DNA sequence
# p1, p2: range of seq to examine
# offset: offset after start codon
# rc: "RC" / "" is reverse complement or not
# result_list: list to append results
def search_V2_motif(seq, p1, p2, offset, rc, result_list):
	
	candidates = []	
	codon = "AG"
	for s in range(p1, p2-(len(codon)+1)):	
		if seq[s:s+len(codon)] == codon:
			candidates.append(s)
	
	# reverse search to prioritize shorter result regions
	for s in reversed(candidates):
		
		# check if omega_nn would overlap with previous result
		s_temp = s + offset
		p2_temp = p2 - (p2 - (s + offset)) % 3
		if overlaps_with_prev_result(s_temp, p2_temp, rc, result_list):
			continue
		
		# limit to multiple of 3 (see p2_temp)
		omega_nn = seq[s_temp:p2_temp]
		
		# translate with NCBI standard table
		omega_aa = omega_nn.translate(table = 1)
		
		if check_V2_motif_criteria(omega_aa) == False:
			# criteria not fulfilled, continue
			continue
			
		else:
			# motiv criteria fulfilled
			# assign TRGroup
			tr_list = [["TRAV2", "VSQQ"], 
			           ["TRBV1", "LQQT"], 
			           ["TRBV2", "EINQ"], 
			           ["TRBV3", "ITQW"], 
			           ["TRGV2", "PIQS"], 
			           ["TRGV3", "AQA"], 
			           ["TRGV4", "WQSP"]]
			tr_group = assign_tr_group(omega_aa, 15, tr_list)
			if tr_group == "": tr_group = "TRV"
			
			result_list.append([omega_nn, 
								omega_aa,
								seq[p2:p2+39],
								tr_group,
								rc, 
								s, 
								p1, 
								p2,
								"V2"])
							
			# skip ATG/AG and rssi between p1 and p2 + 39 nucleotides 
			# and continue at step 2 for rssi+1
			min_next_rss = p2 + 39
			return min_next_rss
			
	return 0


# identify search candidates by RSS motif CAC
# seq: Bio.Seq DNA sequence
# seq_rc: Bio.Seq DNA sequence (reverse complement)
def ident_V_rss_motif(seq, seq_rc):

	# identify search candidates by RSS motif CAC
	instances = [Seq("CAC")]
	m = motifs.create(instances)

	rss = []
	for r in m.instances.search(seq):
		# pointer at 5' end of RSS motif
		rss.append(r[0])

	rss_rc = []
	for r in m.instances.search(seq_rc):
		# pointer at 5' end of RSS motif
		rss_rc.append(r[0])
		
	return rss, rss_rc


# V1: V segments with single-exon leader peptide (TRAV1 and TRGV1)
# seq: Bio.Seq DNA sequence
# rss: index list of search candidates identified by RSS motif
# rc: "RC" / "" is reverse complement or not
# result_list: list to append results
def task_V1(seq, rss, rc, result_list):
	min_next_r = 0
	
	for r in rss:
		
		# obey minimum r and perform search
		if r >= 483 and r >= min_next_r:
			
			# define search region
			# RSS motif 5'-CAC-3' is cut off at 5'-end
			p1 = r - 483
			p2 = r
			min_next_r = search_V1_motif(seq, p1, p2, 0, rc, 
			                             result_list)


# V2: V segments with two-exon leader peptide (all others)
# seq: Bio.Seq DNA sequence
# rss: index list of search candidates identified by RSS motif
# rc: "RC" / "" is reverse complement or not
# result_list: list to append results
def task_V2(seq, rss, rc, result_list):
	min_next_r = 0
	
	for r in rss:
		
		# obey minimum r and perform search
		if r >= 345 and r >= min_next_r:
			
			# define search region
			# RSS motif 5'-CAC-3' is cut off at 5'-end
			p1 = r - 345
			p2 = r
			min_next_r = search_V2_motif(seq, p1, p2, 4, rc, 
			                             result_list)
