from Bio.Seq import Seq
from Bio import motifs

import modules.functions as f


# Check motif criteria for J gene
# omega_aa: translated amino acid sequence of search region
def check_j_motif_criteria(omega_nn, omega_aa_1, omega_aa_2,
                           omega_aa_3, z1, z2):
    # To Do: FGXG (X = beliebige Aminosäure)
    # Simon will das noch prüfen
    omega_aa = Seq("")
    if omega_aa_1.count("FG") > 0:
        omega_aa = omega_aa_1
    if omega_aa_2.count("FG") > 0:
        omega_aa = omega_aa_2
    if omega_aa_3.count("FG") > 0:
        omega_aa = omega_aa_3

    if omega_aa.count("*") > 0:
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
# p1, p2: range of seq to examine
# rc: "RC" / "" means reverse complement or not
# codon_list: list of base sequences identifying the end of a candidate region
# z1: list of base sequences for z1 motiv criterion
# z2: list of base sequences for z2 motiv criterion
# result_list: list to append results
def search_j_motif(seq, p1, p2, rc, codon_list, z1, z2, result_list):
    candidates = []
    for s in range(p1+44, p2-(len(codon_list[0])+1)):
        if codon_list.count(seq[s:s+len(codon_list[0])]) > 0:
            candidates.append(s)

    for s in candidates:
        omega_nn = seq[p1:s]

        # Translate with NCBI standard table (limited to multiple of three)
        translate_end = len(omega_nn) - len(omega_nn) % 3
        omega_aa_1 = omega_nn[0:translate_end].translate(table=1)
        translate_end = len(omega_nn) - (len(omega_nn) - 1) % 3
        omega_aa_2 = omega_nn[1:translate_end].translate(table=1)
        translate_end = len(omega_nn) - (len(omega_nn) - 2) % 3
        omega_aa_3 = omega_nn[2:translate_end].translate(table=1)

        if check_j_motif_criteria(omega_nn, omega_aa_1, omega_aa_2,
                                  omega_aa_3, z1, z2) is False:
            continue
        else:
            # start and end position (python index)
            start = p1
            end = s

            # start and end postition (fasta index)
            # Reduce end by one because end is not included in omega_nn
            # Add one to start/end position as python index starts at 0
            if rc == "RC":
                # convert reverse complement position to position on complement
                start_fasta = (len(seq) - 1 - p1) + 1
                end_fasta = (len(seq) - 1 - (s - 1)) + 1
            else:
                start_fasta = p1 + 1
                end_fasta = (s - 1) + 1

            result_list.append([
                omega_nn,
                "",
                seq[p1-28:p1],
                "TRJ",
                rc,
                "J",
                start,
                end,
                start_fasta,
                end_fasta
                ])

            # Skip GTRDGD and rssi between p1 and p3 + 6 nucleotides
            # and continue at step 2 for rssi+1
            min_next_rss = s + 6
            return min_next_rss

    return 0


# J gene
# seq: Bio.Seq DNA sequence
# rss: index list of search candidates identified by RSS motif
# rc: "RC" / "" is reverse complement or not
# result_list: list to append results
def task_j(seq, rss, rc, result_list):
    min_next_r = 0

    # Resolve wobble bases at highest task level to save time in loops
    codon_list = ["GTRDGD"]
    f.list_resolve_wobble_bases(codon_list)

    # Resolve wobble bases at highest task level to save time in loops
    z1 = ["TTYGGNNNNGG"]
    f.list_resolve_wobble_bases(z1)

    # Resolve wobble bases at highest task level to save time in loops
    z2 = ["TNNBNRT"]
    f.list_resolve_wobble_bases(z2)

    for r in rss:
        # Obey maximum r and perform search
        if r + 78 <= len(seq) and r >= min_next_r:
            # Define search region
            p1 = r
            p2 = r + 78
            min_next_r = search_j_motif(
                seq, p1, p2, rc, codon_list, z1, z2, result_list)
