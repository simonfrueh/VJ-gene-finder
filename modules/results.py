import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# Print results with legend
# result_list: list of all results to show
def show_results(result_list):
    print("\nResults\n%i result(s) found" % len(result_list))

    if len(result_list) > 0:
        print("legend: [omega_nn, omega_aa, seq_rssi, TRgroup"
              + ", RC, s, p1, p2]")
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
        # Fasta header variables for V genes
        # add +1 to start/end position as python index starts at 0
        if r[8] == "V1" or r[8] == "V2":
            start = r[5] + 1
            end = r[7] + 1
        # Fasta header variables for J genes
        elif r[8] == "J":
            start = r[6] + 1
            end = r[5] + 1
        header = (r[3] + "-" + str(start) + "*01|" + sec_record_id
                  + "|F|" + str(start) + "-" + str(end) + "|")
        if r[4] == "RC":
            header += "RC|"

        # Prepare record for fasta file (omega_nn)
        rec = SeqRecord(
            Seq(r[0],),
            id=header,
            description="",
        )
        records_nn.append(rec)

        # Prepare record for fasta file (rss_i to rss_i + 39)
        rec = SeqRecord(
            Seq(r[2],),
            id="RSS-"+header,
            description="",
        )
        records_rss.append(rec)

        # Prepare record for fasta file (omega_aa)
        # check and skip if omega_aa is empty (possible for J genes)
        if r[1] != "":
            rec = SeqRecord(
                Seq(r[1],),
                id=header,
                description="",
            )
            records_aa.append(rec)

    if len(records_nn) > 0:
        filename = (sec_record_id.replace("*", "_").replace("|", "_")
                    + ".fasta")
        print(filename)
        SeqIO.write(records_nn, output_dir + filename, "fasta")

    if len(records_rss) > 0:
        filename = (sec_record_id.replace("*", "_").replace("|", "_")
                    + "_RSS" + ".fasta")
        print(filename)
        SeqIO.write(records_rss, output_dir + filename, "fasta")

    if len(records_aa) > 0:
        filename = (sec_record_id.replace("*", "_").replace("|", "_")
                    + "_aa" + ".fasta")
        print(filename)
        SeqIO.write(records_aa, output_dir + filename, "fasta")
