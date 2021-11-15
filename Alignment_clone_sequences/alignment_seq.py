import os
import sys
from optparse import OptionParser
import pandas
import numpy as np
import time
from Bio.Align.Applications import ClustalwCommandline
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

#####################################################################

def read_AIRR(nomFi):
	col_list = ['sequence_id','sequence','sequence_aa','rev_comp','productive','complete_vdj','vj_in_frame','stop_codon','locus','v_call','d_call','j_call','c_call','sequence_alignment','sequence_alignment_aa','germline_alignment','germline_alignment_aa','junction','junction_aa','np1','np1_aa','np2','np2_aa','cdr1','cdr1_aa','cdr2','cdr2_aa','cdr3','cdr3_aa','fwr1','fwr1_aa','fwr2','fwr2_aa','fwr3','fwr3_aa','fwr4','fwr4_aa','v_score','v_identity','v_support','v_cigar','d_score','d_identity','d_support','d_cigar','j_score','j_identity','j_support','j_cigar','c_score','c_identity','c_support','c_cigar','v_sequence_start','v_sequence_end','v_germline_start','v_germline_end','v_alignment_start','v_alignment_end','d_sequence_start','d_sequence_end','d_germline_start','d_germline_end','d_alignment_start','d_alignment_end','j_sequence_start','j_sequence_end','j_germline_start','j_germline_end','j_alignment_start','j_alignment_end','cdr1_start','cdr1_end','cdr2_start','cdr2_end','cdr3_start','cdr3_end','fwr1_start','fwr1_end','fwr2_start','fwr2_end','fwr3_start','fwr3_end','fwr4_start','fwr4_end','v_sequence_alignment','v_sequence_alignment_aa','d_sequence_alignment','d_sequence_alignment_aa','j_sequence_alignment','j_sequence_alignment_aa','c_sequence_alignment','c_sequence_alignment_aa','v_germline_alignment','v_germline_alignment_aa','d_germline_alignment','d_germline_alignment_aa','j_germline_alignment','j_germline_alignment_aa','c_germline_alignment','c_germline_alignment_aa','junction_length','junction_aa_length','np1_length','np2_length','n1_length','n2_length','p3v_length','p5d_length','p3d_length','p5j_length','consensus_count','duplicate_count','cell_id','clone_id','rearrangement_id','repertoire_id','rearrangement_set_id','sequence_analysis_category','d_number','5prime_trimmed_n_nb','3prime_trimmed_n_nb','insertions','deletions','junction_decryption']
	df = pandas.read_csv(nomFi, sep='\t',usecols=col_list)
	df.replace(np.nan,'', inplace=True)
	df["whole_seq"] = df["fwr1"].astype(str) + df["cdr1"].astype(str) +df["fwr2"].astype(str)+ df["cdr2"].astype(str)+df["fwr3"].astype(str)+df["cdr3"].astype(str)+df["fwr4"].astype(str)
	df["germline_seq"] = df["v_germline_alignment"].astype(str) + df["np1"].astype(str) + df["d_germline_alignment"].astype(str) + df["np2"].astype(str) + df["j_germline_alignment"].astype(str)
	df['germline_seq'] = df['germline_seq'].str.replace('.','')
	df.replace('', np.nan, inplace=True)
	df.dropna(axis=0, how='any', thresh=None, subset=["cdr1_start","cdr2_start","cdr3_start","cdr1_end","cdr2_end","cdr3_end","v_sequence_start","v_identity","j_identity","d_sequence_end","d_sequence_start"], inplace=True)
	#print(df.loc[22,"v_identity"])
	#print(df.loc[1,"sequence"],"     ",df.loc[1,"germline_seq"])
	return df

#=============================================================================#

def read_seq_info(airr_df):
	sequences = airr_df.sort_values(by=['v_identity','j_identity','d_identity'], ascending=False)
	return sequences.iloc[0]['sequence_id']

#=============================================================================#

def write_clonotype_align(airr_df,repertoire_name,germline):
	file_name = repertoire_name+"_selected_seq.fasta"
	#sequence_info = {}
	filetowrite=open(file_name,"w")
	G = ">germline"+ "\n" +germline+ "\n"
	filetowrite.write(G)
	for i in airr_df.index: 
		if len(airr_df["whole_seq"][i]) != 0 :
			seq = ">"+str(airr_df['sequence_id'][i]) + "\n" + airr_df["whole_seq"][i] + "\n"
			#sequence_info[] = {"cdr1":,cdr2:""}
			filetowrite.write(seq)
	filetowrite.close()
	#return sequence_info

#=============================================================================#

def alignment(repertoire_name):
	file_name = repertoire_name+"_selected_seq.fasta"
	#clustalw= "/Users/nikaabdollahi/opt/anaconda3/envs/python3.6/bin/clustalw2"
	#cline = ClustalwCommandline(clustalw, infile=file_name, outfile= "nika.aln")
	muscle_cline = MuscleCommandline(input=file_name,out=os.path.splitext(file_name)[0]+".aln")
	muscle_cline()
	align = AlignIO.read(os.path.splitext(file_name)[0]+".aln", "fasta")
	#print(align)
	count = SeqIO.write(align, os.path.splitext(file_name)[0]+"_uniq.aln.fa", "fasta")
	aligned_seq = [(seq_record.id,seq_record.seq) for seq_record in SeqIO.parse(os.path.splitext(file_name)[0]+"_uniq.aln.fa","fasta")]
	return aligned_seq

#=============================================================================#

def write_clonotype_align_seq(airr_df,repertoire_name,aligned_seq,germline_seq_id):
	sequences = {}
	germline = {}
	file_name = repertoire_name+"_region.txt"
	filetowrite=open(file_name,"w")

	region = "sequence_id	start_V end_V	start_J end_J	start_CDR1 end_CDR1	start_CDR2 end_CDR2	start_CDR3 end_CDR3	startD endD\n"
	filetowrite.write(region)
	for seq in aligned_seq:

		if seq[0] != 'germline':
			a = airr_df.loc[airr_df['sequence_id'] == seq[0]]
			if len(a["cdr1_start"].values) != 0 and len(a["d_sequence_start"].values)!= 0:
				startCDR1 = int(a["cdr1_start"].values[0]) - int(a["v_sequence_start"].values[0])
				endCDR1	= (int(a["cdr1_end"].values[0]) - int(a["v_sequence_start"].values[0])) +1 
				startCDR2 = int(a["cdr2_start"].values[0]) - int(a["v_sequence_start"].values[0])
				endCDR2 = (int(a["cdr2_end"].values[0]) - int(a["v_sequence_start"].values[0]))+1
				startCDR3 = int(float(a["cdr3_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
				startD = int(float(a["d_sequence_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
				endD = (int(float(a["d_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
				endCDR3 = (int(a["cdr3_end"].values[0]) - int(a["v_sequence_start"].values[0]))+1
				startV = 0
				endV = (int(float(a["v_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
				startJ = int(float(a["j_sequence_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
				endJ = (int(float(a["j_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
				list_loc =[startCDR1,endCDR1,startCDR2,endCDR2,startCDR3,startD,endD,endCDR3]
				sequences[seq[0]] = seq[1]
				region = seq[0] + "\t" + str(startV) + " " + str(endV) + "\t" + str(startJ) + " " + str(endJ) + "\t" + str(startCDR1) + " " + str(endCDR1) + "\t" + str(startCDR2) + " " + str(endCDR2) + "\t" + str(startCDR3) + " " + str(endCDR3) + "\t" + str(startD) + " " + str(endD) + "\n"
				filetowrite.write(region)
		else :
			germline["germline"] = seq[1] 
			a = airr_df.loc[airr_df['sequence_id'] == germline_seq_id]
			startCDR1 = int(a["cdr1_start"].values[0]) - int(a["v_sequence_start"].values[0])
			endCDR1	= (int(a["cdr1_end"].values[0]) - int(a["v_sequence_start"].values[0])) +1 
			startCDR2 = int(a["cdr2_start"].values[0]) - int(a["v_sequence_start"].values[0])
			endCDR2 = (int(a["cdr2_end"].values[0]) - int(a["v_sequence_start"].values[0]))+1
			startCDR3 = int(float(a["cdr3_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
			startD = int(float(a["d_sequence_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
			endD = (int(float(a["d_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
			endCDR3 = (int(a["cdr3_end"].values[0]) - int(a["v_sequence_start"].values[0]))+1
			startV = 0
			endV = (int(float(a["v_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
			startJ = int(float(a["j_sequence_start"].values[0])) - int(float(a["v_sequence_start"].values[0]))
			endJ = (int(float(a["j_sequence_end"].values[0])) - int(float(a["v_sequence_start"].values[0])))+1
			region = seq[0] + "\t" + str(startV) + " " + str(endV) + "\t" + str(startJ) + " " + str(endJ) + "\t" + str(startCDR1) + " " + str(endCDR1) + "\t" + str(startCDR2) + " " + str(endCDR2) + "\t" + str(startCDR3) + " " + str(endCDR3) + "\t" + str(startD) + " " + str(endD) + "\n"
			filetowrite.write(region)

	filetowrite.close()

	return sequences, germline


#=============================================================================#

def write_fasta(repertoire_name, sequences):
	file_name = repertoire_name+".fasta"
	filetowrite = open(file_name,"w")
	seq = ""
	for seqID in sequences: 
		seq += ">"+ seqID + "\n" + str(sequences[seqID]) +"\n"
		filetowrite.write(seq)

	filetowrite.close()
	#return sequence_info


#####################################################################
def main():
	start_time = time.time()
	usage = "usage: alignment_intraconal.py -a AIRR_IMGT_annotation_output -n repertoire_name "
	parser = OptionParser(usage)
	parser.add_option("-a", "--AIRR_IMGT_annotation_output", dest="IMGT_seq_info",
	      help="read data from AIRR_IMGT_annotation_output")
	parser.add_option("-n", "--repertoire_name",dest="repertoire_name",
	      help="repertoire_name")
	(options, args) = parser.parse_args()

	if len(sys.argv) != 5:
		parser.error("incorrect number of arguments")
	
	IMGT_seq_info = options.IMGT_seq_info
	repertoire = options.repertoire_name
	airr_df = read_AIRR(IMGT_seq_info)

	seq=read_seq_info(airr_df)
	germline = airr_df.loc[airr_df['sequence_id'] == seq]["germline_seq"].values[0]
	write_clonotype_align(airr_df,repertoire,germline)
	aligned_seq = alignment(repertoire)
	sequences, germlineSeq = write_clonotype_align_seq(airr_df,repertoire,aligned_seq,seq)
	write_fasta(repertoire+"_germline", germlineSeq)
	write_fasta(repertoire+"_sequences", sequences)

	print("The intraclonal alignment step execution time : %s seconds " % (time.time() - start_time))
#####################################################################
if __name__ == "__main__":
	main()
