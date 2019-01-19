#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 19 15:08:26 2019

@author: lwoo0005
"""

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
#The script assumes gene order in the input file is reflective of
#biological gene order.
print "Indicate the folder containing the protein sequences of interest."
in_file = raw_input("Be sure not to include quotation marks and add the final '/'!\n")
db = raw_input("Where is the local database?\n")
print "/Users/lwoo0005/Documents/Laura_stuff/Ruminant_Guts/"
#Directories will be created if the path does not yet exist
#Since file names are typically uninformative, map download numbers to
#meaningful names (bacteria, strain, etc) by making a csv file (Windows 
#format) with two columns, "Number" and "Name". Input here.
print "Ensure you have produced a CSV file (Windows format) mapping download number to strain name using headings 'Number' and 'Name'."
print "These names will be applied to results in output files."
map_csv = raw_input("Input the path to the CSV file.\n")
outtie = raw_input("Where do you want the files to go?\n")

def anti_res_blast(qfile, label, outp):
    print "Working on: "+label
    #Data from faa file written to temporary file.
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
    #Outut directories can be renamed as required
    output_dir = str(outp) + "/hydro_BLAST/XMLs/"
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
    #point to the XMLs
    #Note Blast search set to find only the top hits per query
    #Change e-value based on database and other arguments based on preference
    #https://www.ncbi.nlm.nih.gov/books/NBK279684/
    blastp_cline = NcbiblastpCommandline(query = temp_file,
    db = db, outfmt = "5", out= op, evalue=1e-50, max_target_seqs=1)
    stdout, stderr = blastp_cline()
    print "Finished blast search. Now parsing, filtering, and writing output files.\n"
    #Opening Blast XML files
    result_handle = open(op)
    blast_records = NCBIXML.parse(result_handle)
    suf_t = '.txt'
    flag_list = 0
    #Feel free to change the name of the output directories.
    #However, it's probably better if you have the XML folder
    #and Results folder in the same directory (in this case, dehydro_BLAST)
    result_dir = str(outp) + "/hydro_BLAST/Results/"
    if not os.path.exists(result_dir):
        try:
            os.makedirs(result_dir)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    total_data=[]
    for blast_record in blast_records:
        strain_data=[]
        if len(blast_record.alignments) >=1: 
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    coverage =0
                    sequence = str(alignment.title)
                    #IMPORTANT: only use alignment length when query
                    #is expected to be equal or greater in length to subject.
                    #Otherwise use hsp.align_length
                    length = str(alignment.length)
                    perc_id = str(float(hsp.identities)/float(length)*100)
                    for hydro_record in SeqIO.parse(db, "fasta"):
                        if hydro_record.id in sequence:
                            hydro_len = int(len(hydro_record))
                            coverage = float(int(length)/float(hydro_len))*100
                    #This line filters the XMLs and can be set to allow HSPs with
                    #specific coverage, percent ID, score, etc values.
                    #For no further filtering, use:
                    #if float(perc_id) >= 0:
                    if float(perc_id) >= 70 and float(coverage) >50:
                        query_name = blast_record.query
                        pos=query_name.split(" ")[0].split(">")[0]
                        up_pos= int(pos)-1
                        down_pos = int(pos)+1
                        up_flank=">"+str(up_pos).zfill(9)
                        down_flank=">"+str(down_pos).zfill(9)
                        up_flank_id = ""
                        down_flank_id = ""
                        #Opening file and searching for upstream and downstream genes
                        for seq_record in SeqIO.parse(temp_file, "fasta"):
                            if seq_record.id in query_name:
                                query_sequence = str((seq_record.seq).lower())
                            elif seq_record.id in up_flank:
                                up_flank_id = seq_record.description
                                up_flank_sequence = str((seq_record.seq).lower())
                            elif seq_record.id in down_flank:
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
                        strain_data.append(str(template %(sequence, length, exp, bit, sco, ids, length, pos, length, gaps, quer, match, subj, query_name[10:], quer_start))+"\n")
                        total_data.append(">Strain: "+str(label)+"|ID: "+str(query_name)[10:]+"|Aligned to: "+sequence+" indicated in UPPERCASE|Identity: "+ids+" of "+length+" ("+perc_id+"%) "+" |Coverage (alignment/subject): "+str(coverage)+"%\n")
                        total_data.append(str(mixed_seq+"\n"))
                    #Note in the next few lines that it may be no US or DS region was found. This is reasonable in some cases.
                    #If this happens, it will be printed to the screen, but NOT to a file.  
                        if len(up_flank_id)>1:
                            total_data.append(">Upstream region for "+str(query_name)[10:]+" |Strain: "+str(label)+" |ID: "+str(up_flank_id[10:])+"\n")
                            total_data.append(str(up_flank_sequence)+"\n")
                        else:
                            print "No upstream region located for "+str(qfile)+" at "+ str(query_name)[10:]
                        if len(down_flank_id)>1:
                            total_data.append(">Downstream region for "+str(query_name[10:])+" |Strain: "+str(label)+" |ID: "+str(down_flank_id)[10:]+"\n")
                            total_data.append(str(down_flank_sequence)+"\n")
                        else:
                            print "No downstream region located for "+str(qfile)+" at "+ str(query_name)[10:]
                        resultz = open(str(result_dir) + label + suf_t, 'w')
                        for line in strain_data:
                            resultz.write(line)
                        resultz.close()
    resultss = open(str(result_dir)+"Collated_blast_results"+str(suf_t), 'w')
    for line in total_data:
        resultss.write(line)
    resultss.close()
    #Delete temporary file
    os.remove(temp_file)
    if flag_list ==0:
        print "No hits found for "+str(label)+"."
    else:
        print "Hits found for "+str(label)+"."
    return
    #General presence/absence info is printed to the screen
    

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