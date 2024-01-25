import modules.functions as f

# Resolve wobble bases for usage in fuction search_j_motif (here at highest
# level to save time in loops)
codon_list = ["GTRDGD"]
f.list_resolve_wobble_bases(codon_list)

# Resolve wobble bases for usage in function check_j_motif_criteria (here at
# highest task level to save time in loops
z1 = ["TTYGGNNNNGG"]
f.list_resolve_wobble_bases(z1)
z2 = ["TNNBNRT"]
f.list_resolve_wobble_bases(z2)


# Check motif criteria for J gene
# omega_aa: list of translated amino acid sequences of search region for
#           all three forwared frames
def check_j_motif_criteria(omega_nn, omega_aa):
    # Omega_aa is the frame containing "FG"
    # return False if no "FG" is found
    # return False if any omega_aa containing "FG" also contains a stop codon
    omega_aa_temp = []
    for o in omega_aa:
        if o.count("FG") > 0:
            omega_aa_temp.append(o)
    if len(omega_aa_temp) == 0:
        return False
    for o in omega_aa_temp:
        if o.count("*") > 0:
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


# Search and evaluate candidates for V gene part one
# seq: Bio.Seq DNA sequence
# start_sr, end_sr: start and end index of search region in seq
# rc: "RC" / "" means reverse complement or not
# result_list: list to append results
def search_j_motif(seq, start_sr, end_sr, rc, result_list):
    candidates = []
    for pos in range(start_sr+44, end_sr-(len(codon_list[0])+1)):
        if codon_list.count(seq[pos:pos+len(codon_list[0])]) > 0:
            candidates.append(pos)

    for s in candidates:
        omega_nn = seq[start_sr:s]

        # Translate all three forwared frames with NCBI standard table
        # (limited to multiple of three)
        omega_aa = []
        translate_end = len(omega_nn) - len(omega_nn) % 3
        omega_aa.append(omega_nn[0:translate_end].translate(table=1))
        translate_end = len(omega_nn) - (len(omega_nn) - 1) % 3
        omega_aa.append(omega_nn[1:translate_end].translate(table=1))
        translate_end = len(omega_nn) - (len(omega_nn) - 2) % 3
        omega_aa.append(omega_nn[2:translate_end].translate(table=1))

        if check_j_motif_criteria(omega_nn, omega_aa) is False:
            continue
        else:
            # start and end position (python index) of gene result
            start = start_sr
            end = s

            # start and end postition (fasta index)
            # Reduce end by one because end is not included in omega_nn
            # Add one to start/end position as python index starts at 0
            if rc == "RC":
                # convert reverse complement position to position on complement
                start_fasta = (len(seq) - 1 - start) + 1
                end_fasta = (len(seq) - 1 - (end - 1)) + 1
            else:
                start_fasta = start + 1
                end_fasta = (end - 1) + 1

            result_list.append([
                omega_nn,
                "",
                seq[start-28:start],
                "TRJ",
                rc,
                "J",
                start,
                end,
                start_fasta,
                end_fasta
                ])

            # Skip GTRDGD and rss up to end + 6 nucleotides
            min_next_rss = end + 6
            return min_next_rss

    return 0


# J gene
# seq: Bio.Seq DNA sequence
# rss: index list of search candidates identified by RSS motif
# rc: "RC" / "" is reverse complement or not
# result_list: list to append results
def task_j(seq, rss, rc, result_list):
    min_next_r = 0

    for r in rss:
        # Obey maximum r and perform search in search region [r:r+78]
        if r + 78 <= len(seq) and r >= min_next_r:
            min_next_r = search_j_motif(seq, r, r+78, rc, result_list)
