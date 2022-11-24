from Bio.Seq import Seq


# Check motif criteria for V gene part 1
# omega_aa: translated amino acid sequence of search region
def check_v1_motif_criteria(omega_aa):
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


# Check motif criteria for V gene part 2
# omega_aa: translated amino acid sequence of search region
def check_v2_motif_criteria(omega_aa):
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
def overlaps_with_prev_v1_result(start, end, rc_type, result_list):
    overlap = False

    for r in result_list:
        # Only check for overlap with V1 gene type results
        if r[5] == "V1":
            # Only check for overlap between non_rc and non_rc
            # and between rc and rc
            if r[4] == rc_type:
                # start = r[6], end = r[7]
                if start >= r[6] and start < r[7]:
                    overlap = True
                if end > r[6] and end <= r[7]:
                    overlap = True
                if r[6] >= start and r[6] < end:
                    overlap = True

    return overlap


# Search and evaluate candidates for V gene part one
# seq: Bio.Seq DNA sequence
# p1, p2: range of seq to examine
# rc: "RC" / "" means reverse complement or not
# result_list: list to append results
def search_v1_motif(seq, p1, p2, rc, result_list):
    candidates = []
    codon = "ATG"

    # ATG only relevant within rss_i + 163 nucleotides
    for s in range(p1, p1+163):
        if seq[s:s+len(codon)] == codon:
            candidates.append(s)

    for s in candidates:

        # Limit to multiple of three
        omega_nn = seq[s:p2]

        # Translate with NCBI standard table (limited to multiple of three)
        translate_end = len(omega_nn) - len(omega_nn) % 3
        omega_aa = omega_nn[0:translate_end].translate(table=1)

        if check_v1_motif_criteria(omega_aa) is False:
            continue
        else:
            tr_list = [
                ["TRAV1", "QVQQ"],
                ["TRGV1", "QVLLQQ"],
                ["TRAV2", "VSQQ"],
                ["TRAV3", "LQYP"],
                ["TRBV1", "LQQT"],
                ["TRBV2", "EINQ"],
                ["TRBV3", "ITQW"],
                ["TRGV2", "PIQS"],
                ["TRGV3", "AQA"],
                ["TRGV4", "WQSP"],
                ["TRDV1", "ETSGGGV"],
                ["TRDV2", "LEASGGG"],
                ["TRDV3", "EIHAKK"]
                ]
            tr_group = assign_tr_group(omega_aa, len(omega_aa), tr_list)
            if tr_group == "":
                tr_group = "TRV-SEL"

            # Start and end position (python index)
            start = s
            end = p2

            # Start and end postition (fasta index)
            # Reduce end by one because end is not included in omega_nn
            # Add one to start/end position as python index starts at 0
            if rc == "RC":
                # Convert reverse complement position to position on complement
                start_fasta = (len(seq) - 1 - s) + 1
                end_fasta = (len(seq) - 1 - (p2 - 1)) + 1
            else:
                start_fasta = s + 1
                end_fasta = (p2 - 1) + 1

            result_list.append([
                omega_nn,
                omega_aa,
                seq[p2:p2+39],
                tr_group,
                rc,
                "V1",
                start,
                end,
                start_fasta,
                end_fasta
                ])

            # Skip ATG/AG and rssi between p1 and p2 + 39 nucleotides
            # and continue at step two for rssi+1
            min_next_rss = p2 + 39
            return min_next_rss

    return 0


# Search and evaluate candidates for V gene part two
# seq: Bio.Seq DNA sequence
# p1, p2: range of seq to examine
# rc: "RC" / "" means reverse complement or not
# result_list: list to append results
def search_v2_motif(seq, p1, p2, rc, result_list):
    candidates = []
    splice_site = "AG"

    for s in range(p1, p2-(len(splice_site)+1)):
        if seq[s:s+len(splice_site)] == splice_site:
            candidates.append(s)

    # Reverse search to prioritize shorter result regions
    for s in reversed(candidates):
        # Check if omega_nn would overlap with previous result
        if overlaps_with_prev_v1_result(s + len(splice_site), p2,
                                        rc, result_list):
            continue

        # Limit to multiple of three (see p2_temp above)
        omega_nn = seq[s+len(splice_site):p2]

        # Translate with NCBI standard table (limited to multiple of three)
        translate_end = len(omega_nn) - (len(omega_nn) - 2) % 3
        omega_aa = omega_nn[2:translate_end].translate(table=1)

        if check_v2_motif_criteria(omega_aa) is False:
            continue
        else:
            tr_list = [
                ["TRAV1", "QVQQ"],
                ["TRGV1", "QVLLQQ"],
                ["TRAV2", "VSQQ"],
                ["TRAV3", "LQYP"],
                ["TRBV1", "LQQT"],
                ["TRBV2", "EINQ"],
                ["TRBV3", "ITQW"],
                ["TRGV2", "PIQS"],
                ["TRGV3", "AQA"],
                ["TRGV4", "WQSP"],
                ["TRDV1", "ETSGGGV"],
                ["TRDV2", "LEASGGG"],
                ["TRDV3", "EIHAKK"]
                ]
            tr_group = assign_tr_group(omega_aa, 15, tr_list)
            if tr_group == "":
                tr_group = "TRV"

            # Start and end position (python index)
            start = s + len(splice_site)
            end = p2

            # Start and end postition (fasta index)
            # Reduce end by one because end is not included in omega_nn
            # Add one to start/end position as python index starts at 0
            if rc == "RC":
                # Convert reverse complement position to position on complement
                start_fasta = (len(seq) - 1 - (s + len(splice_site))) + 1
                end_fasta = (len(seq) - 1 - (p2 - 1)) + 1
            else:
                start_fasta = (s + len(splice_site)) + 1
                end_fasta = (p2 - 1) + 1

            result_list.append([
                omega_nn,
                omega_aa,
                seq[p2:p2+39],
                tr_group,
                rc,
                "V2",
                start,
                end,
                start_fasta,
                end_fasta
                ])

            # Skip ATG/AG and rssi between p1 and p2 + 39 nucleotides
            # and continue at step 2 for rssi+1
            min_next_rss = p2 + 39
            return min_next_rss

    return 0


# V gene V1: V segments with single-exon leader peptide (TRAV1 and TRGV1)
# seq: Bio.Seq DNA sequence
# rss: index list of search candidates identified by RSS motif
# rc: "RC" / "" is reverse complement or not
# result_list: list to append results
def task_v1(seq, rss, rc, result_list):
    min_next_r = 0

    for r in rss:
        # Obey minimum r and perform search
        if r >= 483 and r >= min_next_r:
            # Define search region
            # RSS motif 5'-CAC-3' is cut off at 5'-end
            p1 = r - 483
            p2 = r
            min_next_r = search_v1_motif(seq, p1, p2, rc, result_list)


# V gene V2: V segments with two-exon leader peptide (all others)
# seq: Bio.Seq DNA sequence
# rss: index list of search candidates identified by RSS motif
# rc: "RC" / "" is reverse complement or not
# result_list: list to append results
def task_v2(seq, rss, rc, result_list):
    min_next_r = 0

    for r in rss:
        # Obey minimum r and perform search
        if r >= 345 and r >= min_next_r:
            # Define search region
            # RSS motif 5'-CAC-3' is cut off at 5'-end
            p1 = r - 345
            p2 = r
            min_next_r = search_v2_motif(seq, p1, p2, rc, result_list)
