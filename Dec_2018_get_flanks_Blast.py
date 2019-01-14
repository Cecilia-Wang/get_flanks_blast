#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 15:02:56 2018

@author: lwoo0005
"""

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from glob import glob
import os, errno
import csv

print "This program will perform a Blast search against the customized dehydrogenase AA sequence database."
print "Indicate the folder containing the protein sequences of interest."
in_file = raw_input("Be sure not to include quotation marks!\n")
db = raw_input("Where is the local database?\n")
print "/Users/lwoo0005/Documents/Laura_stuff/Ruminant_Guts/"
#Directories will be created if the path does not yet exist
#Since file names are typically uninformative, map download numbers to
#meaningful names (bacteria, strain, etc) by making a csv file (Windows 
#format) with two columns, "Number" and "Name". Input here.
print "Ensure you have produced a CSV file (Windows format) mapping download number to strain name using headings 'Number' and 'Name'."
print "These names will be applied to results in output files."
print "Current CSV:"
print "/Users/lwoo0005/Documents/Laura_stuff/Ruminant_Guts/Hungate_Non-Hungate_DL_codes.csv"
map_csv = raw_input("Input the path to the CSV file.\n")
outtie = raw_input("Where do you want the files to go?\n")

def anti_res_blast(qfile, label, outp):
    #print label
    num = qfile.split('/')[-1]
    temp_num="temp_"+num
    temp_file = qfile.replace(num, temp_num)
    seq_counter = 0
    q_read = open(qfile, 'r')
    data = []
    for line in q_read:
        if ">" in line:
            seq_counter+=1
            str_seq_counter = str(seq_counter)
            seq_pos = ">"+str_seq_counter.zfill(9)+" "
            rep_line = line.replace(">", seq_pos)
            data.append(rep_line)
        else:
            data.append(line)
    q_write =open(temp_file, 'w')
    for line in data:
        q_write.write(line)
    q_read.close()
    q_write.close()
    output_dir = str(outp) + "/final_3_db_genes_BLAST/All_rum_gut_XMLs/"
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    op = str(output_dir) + str(label) + ".xml"
    #Next three lines perform blast search
    #If already performed and you just want to parse blast XMLs,
    #hash out these lines and ensure your output directory will
    #point to the XMLs.
    #Note Blast search set to find only the top hits per query 
    #at expect values of 1e-50 (permissive IMHO).
    blastp_cline = NcbiblastpCommandline(query = temp_file,
    db = db, outfmt = "5", out= op, evalue=1e-50, max_target_seqs=1)
    stdout, stderr = blastp_cline()
    #Keep the next lines
    result_handle = open(op)
    blast_records = NCBIXML.parse(result_handle)
    suf_t = '.txt'
    flag_list = 0
    #Feel free to change the name of the output directories.
    #However, it's probably better if you have the XML folder
    #and Results folder in the same directory (in this case, Maria_hydro_BLAST)
    result_dir = str(outp) + "/final_3_db_genes_BLAST/All_rum_gut_final3_70_id_50_cov_results--no_flanks/"
    if not os.path.exists(result_dir):
        try:
            os.makedirs(result_dir)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    for blast_record in blast_records:
        if len(blast_record.alignments) >=1: 
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    coverage =0
                    sequence = str(alignment.title)
                    length = str(hsp.align_length)
                    perc_id = str(float(hsp.identities)/float(length)*100)
                    for hydro_record in SeqIO.parse(db, "fasta"):
                        if hydro_record.id in sequence:
                            hydro_len = int(len(hydro_record))
                            coverage = float(int(length)/float(hydro_len))*100
                    #This line parses the XMLs and can be set to allow HSPs with
                    #specific coverage, percent ID, score, etc values.
                    #For no furhter parsing, use:
                    #if float(perc_id) >= 0:
                    if float(perc_id) >= 70 and float(coverage) >50:
                        query_name = blast_record.query
            #The next few lines will depend on the format of your input files,
            #specifically, the header lines. It sets up the "locators" for the query,
            #upstream, and downstream sequences
                        pos2=query_name.split(" ")[0].split(">")[0]
                        up_pos2= int(pos2)-1
                        down_pos2 = int(pos2)+1
                        #genome2 = str(query_name.split("_"))[0].split(" ")[1]
                        up_flank=">"+str(up_pos2).zfill(9)
                        down_flank=">"+str(down_pos2).zfill(9)
                        #print "up flank: "+str(up_flank)
                        #print "down flank: "+str(down_flank)
                        up_flank_id = ""
                        down_flank_id = ""
                        #print query_name
            #This is where the initial input file is re-opened and the sequences
            #pulled out using a Biopython parser. For every hit, the entire
            #file is opened and the regions of interested extracted using the
            #"locators" to match the ID. Then the sequences are collected as
            #lower case letters.
                        for seq_record in SeqIO.parse(temp_file, "fasta"):
                            if seq_record.id in query_name:
                                query_sequence = str((seq_record.seq).lower())
                            elif seq_record.id in up_flank:
                                #print "seq descrip " +str(seq_record.description)
                                up_flank_id = seq_record.description
                                up_flank_sequence = str((seq_record.seq).lower())
                            elif seq_record.id in down_flank:
                                #print "seqdescripdown: "+str(seq_record.description)
                                down_flank_id = seq_record.description
                                down_flank_sequence = str((seq_record.seq).lower())
                        flag_list+=1
                        template = """
Sequence: %s
Alignment length: %s
E-value: %s
Bit score: %s
Score: %s
Identities: %s of %s
Positives: %s of %s
Gaps: %s
%s...
%s...
%s...
Query name: %s
Start position of query: %s
\n
"""
                        #Use the Biopython blast tree to set up variables for
                        #any details you want to keep.
                        sequence = str(alignment.title)
                        length = str(alignment.length)
                        exp = str(hsp.expect)
                        bit = str(hsp.bits)
                        sco = str(hsp.score)
                        ids = str(hsp.identities)
                        pos = str(hsp.positives)
                        gaps = str(hsp.gaps)
                        quer = hsp.query[0:100]
                        match = hsp.match[0:100]
                        subj = hsp.sbjct[0:100]
                        quer_start = int(hsp.query_start)
                        quer_end = quer_start+int(len(hsp.query))
                    #The next few lines set up variables for the parts of the sequences I wanted in upper or lower case
                    #The part of the query that turns up a hit is uppercase,
                    #while the upstream, downstream, and non-hit regions of the sequence are in lower case
                    #It is important to subtract gaps from the alignment length when doing the replacement
                        lower_match = str(query_sequence[quer_start-1:quer_end-1-int(hsp.query.count("-"))])
                        upper_match = str(query_sequence[quer_start-1:quer_end-1-int(hsp.query.count("-"))]).upper()
                        mixed_seq = query_sequence.replace(lower_match, upper_match)
                        if not os.path.exists(result_dir):
                            try:
                                os.makedirs(result_dir)
                            except OSError as exc:
                                if exc.errno != errno.EEXIST:
                                    raise 
                    #This is where you set up the output file for the specific strain.
                        resultz = open(str(result_dir) + label + suf_t, 'a')
                        resultz.write(str(template %(sequence, length, exp, bit, sco, ids, length, pos, length, gaps, quer, match, subj, query_name, quer_start))+"\n")
                        resultz.close()
                    #A new file is opened to collect all of the alignments+US+DS sequences as fasta sequences in a single file.
                    #This file is constantly appended ('a').
                        resultss = open(str(result_dir)+"Collated_blast_results"+str(suf_t), 'a')
                        resultss.write(">Strain: "+str(label)+"|ID: "+str(query_name)+"|Aligned to: "+sequence+" indicated in UPPERCASE|Identity: "+ids+" of "+length+" ("+perc_id+"%) "+" |Coverage (alignment/subject): "+str(coverage)+"%\n")
                        resultss.write(str(mixed_seq)+"\n")
                        #resultss.close()
                    #Note in the next few lines that it may be no US or DS region was found. This is reasonable in some cases.
                    #If this happens, it will be printed to the screen, but NOT to a file.  
                        if len(up_flank_id)>1:
                            resultss.write(">Upstream region for "+str(query_name)+" |Strain: "+str(label)+" |ID: "+str(up_flank_id)+"\n")
                            resultss.write(str(up_flank_sequence)+"\n")
                        else:
                            print "No upstream region located for "+str(temp_file)+" at "+ str(query_name)
                        if len(down_flank_id)>1:
                            resultss.write(">Downstream region for "+str(query_name)+" |Strain: "+str(label)+" |ID: "+str(down_flank_id)+"\n")
                            resultss.write(str(down_flank_sequence)+"\n")
                        else:
                            print "No downstream region located for "+str(temp_file)+" at "+ str(query_name)
                        resultss.close()
    os.remove(temp_file)
    if flag_list ==0:
        print "No hits found for "+str(label)+"."
    else:
        print "Hits found for "+str(label)+"."
    return
"""
    q_read = open(qfile, 'r')
    data = []
    for line in q_read:
        if ">" in line:
            old_header = ">"+" ".join(line.split(" ")[1:])
            data.append(old_header)
        else:
            data.append(line)
    q_write =open(qfile, 'w')
    for line in data:
        q_write.write(line)
    q_read.close()
    q_write.close()
"""
    #General presence/absence info is printed to the screen
    #This can easily be made into a chart by copying to excel,
    #replacing "for" with ":", then splitting by text to columns
    #by ":".
    

#set up dictionary for number_id as keys and bacteria name as value
#If your download names are just numbers, it's very important to make
#a number-to-strain chart and then input the path to it here.
#Be sure that your chart has two columns labeled "Number" for the ID
#and "Name" for the bacterial strain
id_dic = {}
with open(str(map_csv), 'rb') as csvfile:
    reader = csv.DictReader(csvfile, dialect='excel')
    for row in reader:
        id_dic[row['Number']]=str(row['Name'])
used_codes = []
number_ids = []
used_ids = []
#The next line will pick up anything ending with "faa" in the input file.
#This can be changed according to what else is in the file, the extension, etc.
globber = glob(str(in_file) +"/*.faa")
for globet in globber:
    #This is specific to the "non-unique" part of the download name
    number_id=(globet.split(str(in_file))[1]).split(".faa")[0]
    number_ids.append(number_id)
#This part submits each of the files through the function defined above
    for code, label in id_dic.items():    
        if number_id == code:
            anti_res_blast(globet, label, outtie)
            used_codes.append(code)
            used_ids.append(number_id)    
csvfile.close()
#This bit just checks to be sure that all files in the input folder are in 
#agreement with the code-name chart.
for code in id_dic.keys():
    if code not in used_codes:
        print "In CSV but not files: "+str(code)
for number_id in number_ids:
    if number_id not in used_ids:
        print "In files but not CSV: "+str(number_id)