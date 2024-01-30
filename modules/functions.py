from Bio.Seq import Seq
from Bio import motifs


# Identify search candidates by RSS motif and return start position
# seq: Bio.Seq DNA sequence
# motif: RSS motif
def ident_rss_motif_start_position(seq, motif):
    # Identify search candidates by RSS motif
    m = motifs.create([Seq(motif)])

    search_candidates = []
    for p, s in seq.search(m.alignment):
        # Pointer at 5' end of RSS motif
        search_candidates.append(p)

    return search_candidates


# Identify search candidates by RSS motif and return end position
# seq: Bio.Seq DNA sequence
# motif: RSS motif
def ident_rss_motif_end_position(seq, motif):
    # Identify search candidates by RSS motif
    m = motifs.create([Seq(motif)])

    search_candidates = []
    for pos, s in seq.search(m.alignment):
        # Pointer at 3' end of RSS motif
        search_candidates.append(pos + len(motif))

    return search_candidates


# Return a list of all possible base sequences with wobble bases
# replaced
def list_resolve_wobble_bases(seq_list):
    for i in range(0, len(seq_list)):
        if seq_list[i].count("R") > 0:
            seq_list.append(seq_list[i].replace("R", "A", 1))
            seq_list.append(seq_list[i].replace("R", "G", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("Y") > 0:
            seq_list.append(seq_list[i].replace("Y", "C", 1))
            seq_list.append(seq_list[i].replace("Y", "T", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("M") > 0:
            seq_list.append(seq_list[i].replace("M", "A", 1))
            seq_list.append(seq_list[i].replace("M", "C", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("K") > 0:
            seq_list.append(seq_list[i].replace("K", "G", 1))
            seq_list.append(seq_list[i].replace("K", "T", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("S") > 0:
            seq_list.append(seq_list[i].replace("S", "C", 1))
            seq_list.append(seq_list[i].replace("S", "G", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("W") > 0:
            seq_list.append(seq_list[i].replace("W", "A", 1))
            seq_list.append(seq_list[i].replace("W", "T", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("B") > 0:
            seq_list.append(seq_list[i].replace("B", "C", 1))
            seq_list.append(seq_list[i].replace("B", "G", 1))
            seq_list.append(seq_list[i].replace("B", "T", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("D") > 0:
            seq_list.append(seq_list[i].replace("D", "A", 1))
            seq_list.append(seq_list[i].replace("D", "G", 1))
            seq_list.append(seq_list[i].replace("D", "T", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("H") > 0:
            seq_list.append(seq_list[i].replace("H", "A", 1))
            seq_list.append(seq_list[i].replace("H", "C", 1))
            seq_list.append(seq_list[i].replace("H", "T", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("V") > 0:
            seq_list.append(seq_list[i].replace("V", "A", 1))
            seq_list.append(seq_list[i].replace("V", "C", 1))
            seq_list.append(seq_list[i].replace("V", "G", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)

        if seq_list[i].count("N") > 0:
            seq_list.append(seq_list[i].replace("N", "A", 1))
            seq_list.append(seq_list[i].replace("N", "C", 1))
            seq_list.append(seq_list[i].replace("N", "G", 1))
            seq_list.append(seq_list[i].replace("N", "T", 1))
            seq_list.pop(i)
            list_resolve_wobble_bases(seq_list)


# Convert start position to fasta index
# Add +1 to start_pos as python index starts at 0
def start_to_fasta(seq_len, start_pos, is_reverse_complement):
    if is_reverse_complement:
        # convert reverse complement position to position on complement
        start_pos_fasta = (seq_len - 1 - start_pos) + 1
    else:
        start_pos_fasta = start_pos + 1

    return start_pos_fasta


# Convert end position to fasta index
# Reduce end_pos by one because end is not included in omega_nn
# Add +1 to end_pos as python index starts at 0
def end_to_fasta(seq_len, end_pos, is_reverse_complement):
    if is_reverse_complement:
        # convert reverse complement position to position on complement
        end_pos_fasta = (seq_len - 1 - (end_pos - 1)) + 1
    else:
        end_pos_fasta = (end_pos - 1) + 1

    return end_pos_fasta
