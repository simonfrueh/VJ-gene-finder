import modules.functions as f

# define list of possible tr_group identifiers for usage in functions
# search_v1_motif and search_v2_motif
tr_list = [
    ["TRAV1", "QVQQ"],
    ["TRAV2", "VSQQ"],
    ["TRAV3", "LQYP"],
    ["TRBV1", "LQQT"],
    ["TRBV2", "EINQ"],
    ["TRBV3", "ITQW"],
    ["TRGV1", "QVLLQQ"],
    ["TRGV2", "PIQS"],
    ["TRGV3", "QAVPMQ"],
    ["TRGV3", "QAAPVQ"],
    ["TRGV4", "LWQSP"],
    ["TRDV1", "ETSGGGV"],
    ["TRDV2", "LEASGGG"],
    ["TRDV3", "VEFGGDV"],
    ["TRDV4", "RIVEAG"],
    ["TRDV5", "EIHAKKSA"],
    ["TRDVH1", "QIEMVTT"]
    ]


# Check motif criteria for V gene part 1
# omega_aa: translated amino acid sequence of search region
def check_v1_motif_criteria(omega_aa):
    omega_aa_len = len(omega_aa)

    # omega_nn does not contain a stop codon TAA, TAG or TGA
    # or omega_aa does not contain *
    if omega_aa.count("*") > 0:
        return False

    # Length Ωaa > 102 and < 161
    if (omega_aa_len <= 102 or omega_aa_len >= 161):
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
    count_Y = omega_aa[-12:omega_aa_len].count("YYC")
    count_Y += omega_aa[-12:omega_aa_len].count("YLC")
    count_Y += omega_aa[-12:omega_aa_len].count("YHC")
    count_Y += omega_aa[-12:omega_aa_len].count("YFC")
    count_Y += omega_aa[-12:omega_aa_len].count("YIC")
    count_Y += omega_aa[-12:omega_aa_len].count("TFC")
    if count_Y == 0:
        return False

    return True


# Check motif criteria for V gene part 2
# omega_aa: translated amino acid sequence of search region
def check_v2_motif_criteria(omega_aa):
    omega_aa_len = len(omega_aa)

    # omega_nn does not contain a stop codon TAA, TAG or TGA
    # or omega_aa does not contain *
    if omega_aa.count("*") > 0:
        return False

    # Length Ωaa > 85 and < 115
    if (omega_aa_len <= 85 or omega_aa_len >= 115):
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
    count_Y = omega_aa[-12:omega_aa_len].count("YYC")
    count_Y += omega_aa[-12:omega_aa_len].count("YFC")
    count_Y += omega_aa[-12:omega_aa_len].count("YLC")
    count_Y += omega_aa[-12:omega_aa_len].count("YHC")
    count_Y += omega_aa[-12:omega_aa_len].count("YIC")
    count_Y += omega_aa[-12:omega_aa_len].count("TFC")
    if count_Y == 0:
        return False

    return True


# Assign TRgroup based on amino acid motif
# omega_aa: amino acid sequence
# n: first n amino acids in omega_aa to be searched
# tr_list: [[group name, amino acid motiv],[...]]
def assign_tr_group(omega_aa, n, tr_list):
    tr_group = ""

    for t in tr_list:
        if omega_aa[0:n].count(t[1]) > 0:
            tr_group += " " + t[0]

    return tr_group.lstrip()


# Check whether search region defined by start and end overlaps with
# previous results from V part 1
def overlaps_with_prev_v1_result(start, end, is_rc, result_list):
    for r in result_list:
        # Only check for overlap with V1 gene type results
        if r["gene_type"] == "V1":
            # Only check for overlap between non_rc and non_rc
            # and between rc and rc
            if r["is_reverse_complement"] == is_rc:
                if start >= r["start_pos"] and start < r["end_pos"]:
                    return True
                if end > r["start_pos"] and end <= r["end_pos"]:
                    return True
                if r["start_pos"] >= start and r["start_pos"] < end:
                    return True

    return False


# Search and evaluate candidates for V gene part one
# seq: Bio.Seq DNA sequence
# start_sr, end_sr: start and end index of search region in seq
# rc: True means seq is a reverse complement
# result_list: list to append results
def search_v1_motif(seq, start_sr, end_sr, is_rc, result_list):
    candidates = []
    codon = "ATG"

    # ATG only relevant within first 163 nucleotides
    for pos in range(start_sr, start_sr+163):
        if seq[pos:pos+len(codon)] == codon:
            candidates.append(pos)

    for s in candidates:
        omega_nn = seq[s:end_sr]

        # Translate with NCBI standard table (limited to multiple of three)
        translate_end = len(omega_nn) - len(omega_nn) % 3
        omega_aa = omega_nn[0:translate_end].translate(table=1)

        if check_v1_motif_criteria(omega_aa) is False:
            continue
        else:
            tr_group = assign_tr_group(omega_aa, len(omega_aa), tr_list)
            if tr_group == "":
                tr_group = "TRV-SEL"

            # Start and end position (python index) of gene result
            start = s
            end = end_sr

            gene_result_dict = {
                "omega_nn": omega_nn,
                "omega_aa": omega_aa,
                "seq": seq[end:end+39],
                "tr_group": tr_group,
                "is_reverse_complement": is_rc,
                "gene_type": "V1",
                "start_pos": start,
                "end_pos": end,
                "start_pos_fasta": f.start_to_fasta(len(seq), start, is_rc),
                "end_pos_fasta": f.end_to_fasta(len(seq), end, is_rc)
            }

            result_list.append(gene_result_dict)

            # Skip ATG and rss up to end + 39 nucleotides
            min_next_rss = end + 39
            return min_next_rss

    return 0


# Search and evaluate candidates for V gene part two
# seq: Bio.Seq DNA sequence
# start_sr, end_sr: start and end index of search region in seq
# is_rc: True means seq is a reverse complement
# result_list: list to append results
def search_v2_motif(seq, start_sr, end_sr, is_rc, result_list):
    candidates = []
    splice_site = "AG"
    for pos in range(start_sr, end_sr-(len(splice_site)+1)):
        if seq[pos:pos+len(splice_site)] == splice_site:
            candidates.append(pos)

    # Reverse search to prioritize shorter result regions
    for s in reversed(candidates):
        # Check if omega_nn would overlap with previous result
        if overlaps_with_prev_v1_result(s + len(splice_site), end_sr,
                                        is_rc, result_list):
            continue

        omega_nn = seq[s+len(splice_site):end_sr]

        # Translate with NCBI standard table (limited to multiple of three)
        translate_end = len(omega_nn) - (len(omega_nn) - 2) % 3
        omega_aa = omega_nn[2:translate_end].translate(table=1)

        if check_v2_motif_criteria(omega_aa) is False:
            continue
        else:
            tr_group = assign_tr_group(omega_aa, 15, tr_list)
            if tr_group == "":
                tr_group = "TRV"

            # Start and end position (python index) of gene result
            start = s + len(splice_site)
            end = end_sr

            gene_result_dict = {
                "omega_nn": omega_nn,
                "omega_aa": omega_aa,
                "seq": seq[end:end+39],
                "tr_group": tr_group,
                "is_reverse_complement": is_rc,
                "gene_type": "V2",
                "start_pos": start,
                "end_pos": end,
                "start_pos_fasta": f.start_to_fasta(len(seq), start, is_rc),
                "end_pos_fasta": f.end_to_fasta(len(seq), end, is_rc)
            }

            result_list.append(gene_result_dict)

            # Skip AG and rss up to end + 39 nucleotides
            min_next_rss = end + 39
            return min_next_rss

    return 0


# V gene V1: V segments with single-exon leader peptide (TRAV1 and TRGV1)
# seq: Bio.Seq DNA sequence
# rss: index list of search candidates identified by RSS motif
# is_rc: True means seq is a reverse complement
# result_list: list to append results
def task_v1(seq, rss, is_rc, result_list):
    min_next_r = 0

    for r in rss:
        # Obey minimum r and perform search in search region [r-483:r]
        if r >= 483 and r >= min_next_r:
            min_next_r = search_v1_motif(seq, r-483, r, is_rc, result_list)


# V gene V2: V segments with two-exon leader peptide (all others)
# seq: Bio.Seq DNA sequence
# rss: index list of search candidates identified by RSS motif
# is_rc: True means seq is a reverse complement
# result_list: list to append results
def task_v2(seq, rss, is_rc, result_list):
    min_next_r = 0

    for r in rss:
        # Obey minimum r and perform search in search region [r-345:r]
        if r >= 345 and r >= min_next_r:
            min_next_r = search_v2_motif(seq, r-345, r, is_rc, result_list)
