import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# Print results with legend
# result_list: list of all results to show
def show_results(result_list):
    print("\nResults\n%i result(s) found" % len(result_list))

    if len(result_list) > 0:
        print("legend: [omega_nn, omega_aa, seq_rssi, TRgroup,"
              + "RC, gene type, start, end, start_fasta, end_fasta]")
        for n in result_list:
            print(n)


# Write results to output directory, directory will be created/overwritten
# output_dir: output directory
# sec_record_id: sec_record.id for filename/header of fasta result files
# result_list: list of all results to write
def write_results(output_dir, sec_record_id, result_list):
    if len(result_list) == 0:
        return

    print("\nWriting result files to output directory")

    # Make sure output directory exists
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    records_nn = []
    records_rss = []
    records_aa = []

    for r in result_list:
        # Determine fasta header
        header = (r["tr_group"] + "-" + str(r["start_pos_fasta"]) + "*01|"
                  + sec_record_id + "|F|" + str(r["start_pos_fasta"]) + "-"
                  + str(r["end_pos_fasta"]) + "|")
        if r["is_reverse_complement"]:
            header += "RC|"

        # Prepare record for fasta file (omega_nn)
        rec = SeqRecord(
            Seq(r["omega_nn"],),
            id=header,
            description="",
        )
        records_nn.append(rec)

        # Prepare record for fasta file (omega_aa)
        # check and skip if omega_aa is empty (possible for J genes)
        if r["omega_aa"] != "":
            rec = SeqRecord(
                Seq(r["omega_aa"],),
                id=header,
                description="",
            )
            records_aa.append(rec)

        # Prepare record for fasta file (rss)
        # Use differen header for J genes
        if r["gene_type"] == "J":
            header = ("TRJ-" + str(r["start_pos_fasta"]) + "*01|"
                      + sec_record_id + "|")
        rec = SeqRecord(
            Seq(r["seq"],),
            id="RSS-"+header,
            description="",
        )
        records_rss.append(rec)

    if len(records_nn) > 0:
        filename = (sec_record_id.replace("*", "_").replace("|", "_")
                    + ".fasta")
        print(filename)
        SeqIO.write(records_nn, output_dir + filename, "fasta")

    if len(records_aa) > 0:
        filename = (sec_record_id.replace("*", "_").replace("|", "_")
                    + "_aa" + ".fasta")
        print(filename)
        SeqIO.write(records_aa, output_dir + filename, "fasta")

    if len(records_rss) > 0:
        filename = (sec_record_id.replace("*", "_").replace("|", "_")
                    + "_RSS" + ".fasta")
        print(filename)
        SeqIO.write(records_rss, output_dir + filename, "fasta")
