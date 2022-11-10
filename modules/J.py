from Bio.Seq import Seq
from Bio import motifs

# return a list of all possible base sequences with wobble 
# replaced
def list_resolve_wobble_bases(seq_list):
	
	for i in range(0, len(seq_list)):
		if seq_list[i].count("R") > 0 :
			seq_list.append(seq_list[i].replace("R","A",1))
			seq_list.append(seq_list[i].replace("R","G",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("Y") > 0 :
			seq_list.append(seq_list[i].replace("Y","C",1))
			seq_list.append(seq_list[i].replace("Y","T",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("M") > 0 :
			seq_list.append(seq_list[i].replace("M","A",1))
			seq_list.append(seq_list[i].replace("M","C",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("K") > 0 :
			seq_list.append(seq_list[i].replace("K","G",1))
			seq_list.append(seq_list[i].replace("K","T",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("S") > 0 :
			seq_list.append(seq_list[i].replace("S","C",1))
			seq_list.append(seq_list[i].replace("S","G",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("W") > 0 :
			seq_list.append(seq_list[i].replace("W","A",1))
			seq_list.append(seq_list[i].replace("W","T",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("B") > 0 :
			seq_list.append(seq_list[i].replace("B","C",1))
			seq_list.append(seq_list[i].replace("B","G",1))
			seq_list.append(seq_list[i].replace("B","T",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("D") > 0 :
			seq_list.append(seq_list[i].replace("D","A",1))
			seq_list.append(seq_list[i].replace("D","G",1))
			seq_list.append(seq_list[i].replace("D","T",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("H") > 0 :
			seq_list.append(seq_list[i].replace("H","A",1))
			seq_list.append(seq_list[i].replace("H","C",1))
			seq_list.append(seq_list[i].replace("H","T",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("V") > 0 :
			seq_list.append(seq_list[i].replace("V","A",1))
			seq_list.append(seq_list[i].replace("V","C",1))
			seq_list.append(seq_list[i].replace("V","G",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)
			
		if seq_list[i].count("N") > 0 :
			seq_list.append(seq_list[i].replace("N","A",1))
			seq_list.append(seq_list[i].replace("N","C",1))
			seq_list.append(seq_list[i].replace("N","G",1))
			seq_list.append(seq_list[i].replace("N","T",1))
			seq_list.pop(i)
			list_resolve_wobble_bases(seq_list)


# check motif criteria for V gene
# omega_aa: translated amino acid sequence of search region
def check_J_motif_criteria(omega_nn, z1, z2):
	
	# omega_nn does not contain a stop codon TAA, TAG or TGA
	if (omega_nn.count("TAA") > 0 or omega_nn.count("TAG") 
	   or omega_nn.count("TGA")):
			return False
	
	# Length Ωnn > 42 and < 71
	if (len(omega_nn) <= 42 or
		len(omega_nn) >= 71):
			return False
	
	# z1 = TTNGGNNNNGG between 13 - 48 (Ωnn)
	count_z1 = 0
	for z in z1:
		count_z1 += omega_nn[13-1:48].count(z)
	if count_z1 == 0:
		return False
		
	# z2 = TNNNNNT between 31 - 63 (Ωnn)
	count_z2 = 0
	for z in z2:
		count_z2 += omega_nn[31-1:63].count(z)
	if count_z2 == 0:
		return False
	
	return True


def search_J_motif(seq, p1, p2, rc, codon_list, z1, z2, result_list):
	
	candidates = []	
	for s in range(p1+44, p2-(len(codon_list[0])+1)):
		if codon_list.count(seq[s:s+len(codon_list[0])]) > 0:
			candidates.append(s)
			
	for s in candidates:
		omega_nn = seq[p1:s]
		
		if check_J_motif_criteria(omega_nn, z1, z2) == False:
			# criteria not fulfilled, continue
			continue
				
		else:
			# motiv criteria fulfilled			
			result_list.append([omega_nn, 
								"",
								seq[p1-28:p1],
								"TRJ",
								rc, 
								s, 
								p1, 
								p2,
								"J"])
							
			# skip GTRDGD and rssi between p1 and p3 + 6 nucleotides 
			# and continue at step 2 for rssi+1
			min_next_rss = s + 6
			return min_next_rss
			
	return 0


# J
# seq: Bio.Seq DNA sequence
# rss: index list of search candidates identified by RSS motif
# rc: "RC" / "" is reverse complement or not
# result_list: list to append results
def task_J(seq, rss, rc, result_list):
	min_next_r = 0
	
	# resolve wobble bases at highest task level to save time in loops
	codon_list = ["GTRDGD"]
	list_resolve_wobble_bases(codon_list)
	
	# resolve wobble bases at highest task level to save time in loops
	z1 = ["TTNGGNNNNGG"]
	list_resolve_wobble_bases(z1)
	
	# resolve wobble bases at highest task level to save time in loops
	z2 = ["TNNNNNT"]
	list_resolve_wobble_bases(z2)
	
	for r in rss:
		
		# obey maximum r and perform search
		if r + 78 <= len(seq) and r >= min_next_r:
			
			# define search region
			p1 = r
			p2 = r + 78
			min_next_r = search_J_motif(seq, p1, p2, rc, codon_list, z1, z2, result_list)
                           

# identify search candidates by RSS motif CAC
# seq: Bio.Seq DNA sequence
# seq_rc: Bio.Seq DNA sequence (reverse complement)
def ident_J_rss_motif(seq, seq_rc):

	# identify search candidates by RSS motif GTG
	instances = [Seq("GTG")]
	m = motifs.create(instances)

	rss = []
	for r in m.instances.search(seq):
		# pointer at 3' end of RSS motif
		rss.append(r[0] + len("GTG"))

	rss_rc = []
	for r in m.instances.search(seq_rc):
		# pointer at 3' end of RSS motif
		rss_rc.append(r[0] + len("GTG"))
		
	return rss, rss_rc
