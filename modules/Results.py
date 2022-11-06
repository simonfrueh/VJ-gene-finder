import os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# print results
# result_list: list of all results to show
def show_results(result_list):
	print("\nResults\n%i result(s) found" % len(result_list))
	print("legend: [omega_nn, omega_aa, seq_rssi+39, TRgroup" 
		  +", RC, s, p1, p2]")
	for n in result_list:
		print(n)


# write results to output directory (will be created / overwritten)
# path: path for output directory
# sec_record_id: sec_record.id for filename/header of fasta result files
# result_list: list of all results to write
def write_results(path, sec_record_id, result_list):
	
	print("\nWriting result files to output-directory")

	# make sure output directory exists
	output_path = path + "output/"
	if not os.path.exists(output_path):
		os.mkdir(output_path)

	records_nn = []
	records_rss = []
	records_aa = []
	
	for r in result_list:
		
		# determine fasta header
		if r[8] == "V2":
			# s = r[5] = s1 + 4 --> V2: s1+2 = s -2
			s_temp = r[5] - 2
		else:
			# V1: s1 = s
			s_temp = r[5]
			
		header = (r[3] + "-" + str(s_temp) + "*01|" + sec_record_id 
		          + "|F|" + str(s_temp) + "-" + str(r[7]) + "|")
		if r[4] == "RC": header += "RC|"
		
		# prepare record for fasta file (omega_nn)
		rec = SeqRecord(Seq(r[0],), id=header, description="",)
		records_nn.append(rec)
		
		# prepare record for fasta file (rss_i to rss_i + 39)
		rec = SeqRecord(Seq(r[2],), id="RSS-" + header, description="",)
		records_rss.append(rec)
		
		# prepare record for fasta file (omega_aa)
		rec = SeqRecord(Seq(r[1],), id=header, description="",)
		records_aa.append(rec)
		
	if len(records_nn) > 0:
		filename = (sec_record_id.replace("*","_").replace("|","_") 
					+ ".fasta")
		print(filename)
		SeqIO.write(records_nn, output_path + filename, "fasta")
		
	if len(records_rss) > 0:
		filename = (sec_record_id.replace("*","_").replace("|","_") 
					+ "_RSS" + ".fasta")
		print(filename)
		SeqIO.write(records_rss, output_path + filename, "fasta")
		
	if len(records_aa) > 0:
		filename = (sec_record_id.replace("*","_").replace("|","_") 
					+ "_aa" + ".fasta")
		print(filename)
		SeqIO.write(records_aa, output_path + filename, "fasta")

	# # TRAV1 / TRV V gene results files
	# gtypes = ["TRAV1", "TRV"]
	# for g in gtypes:
		# records = []
		
		# filename = (sec_record_id.replace("*","_").replace("|","_") 
						# + "_" + g + ".fasta")
		
		# for r in result_list:
			# if r[3] == g:
				# header=(r[3] + "-" + str(r[5]) + "*01|" + sec_record_id 
						# + "|F|" + str(r[5]) + "-" + str(r[7]) + "|")
				# if r[4] == "RC": header += "RC|"
				# rec = SeqRecord(Seq(r[0],), id=header, description="",)
				# records.append(rec)
		
		# if len(records) > 0:
			# print(filename)
			# SeqIO.write(records, output_path + filename, "fasta")

	# # TRAV1 / TRV RSS results files
	# for g in gtypes:
		# records = []
		
		# filename = (sec_record_id.replace("*","_").replace("|","_") 
						# + "_" + g + "_RSS" + ".fasta")
		
		# for r in result_list:
			# if r[3] == g:
				# header=("RSS-" + r[3] + "-" + str(r[5]) + "*01|" 
						# + sec_record_id + "|")
				# if r[4] == "RC": header += "RC|"
				# rec = SeqRecord(Seq(r[2],), id=header, description="",)
				# records.append(rec)

		# if len(records) > 0:
			# print(filename)
			# SeqIO.write(records, output_path + filename, "fasta") 
