#This script extracts the regions of interest from .net.axt pairwise alignments
#and prints out text files contining only the regions of interest in fasta format

import re
import string
import matplotlib
import scipy
import numpy as np
from pylab import *
import glob


filelist = glob.glob("axt_files/*.net.axt")


#function to read a file -- if header line shows coordinates overlapping
#with region of interest, prints the sequences to an outfile
def read_file(filename, genome_name, region, outfile):
    infile = open(filename)
    for line in infile:
        if region[1] in line:
            header_line = line.split(" ")
            if int(header_line[2]) >= region[3] or int(header_line[3]) <= region[2]:
                pass 
            else:
                dm3_header = ">dm3_" + "_".join(header_line[1:4])
                other_header = ">" + genome_name + "_" + "_".join(header_line[4:7])
                blast_z = header_line[8]
                header_list = [dm3_header, other_header]
                group_info = "dm3 compared with " + genome_name + " blastz score: " + blast_z
                outfile.write(group_info + "\r\n")
                for i in range(2):
                    outfile.write(header_list[i] + "\r\n")
                    outfile.write(infile.next())
                    outfile.write("\r\n")
    infile.close()


#function to identify axt.net files containing relevant chromosome data
def find_files(region):
    print "doing ", region[0]
    outfile = open(region[0] + "_aligned_fragments.txt", "w")
    for i in filelist:
        if region[1] in i:
            second_genome = i.replace(".net.axt", "")
            second_genome = second_genome.replace("axt_files/chr3R.dm3.", "")
            read_file(i, second_genome, region, outfile)
    outfile.close()

InR_region = ["InR_locus","chr3R",17345970,17495043]
foxo_region = ["foxo_locus","chr3R",9857667,9938700]
Thor_region = ["Thor_locus","chr2L",3453434,3504612]
InR_TU = ["InR_transcriptional_unit","chr3R",17395970,17445043]
InR_exon = ["InR_exon","chr3R",17399113,17401638]
InR_T1 = ["InR_promoter","chr3R",17444777,17445965]
InR_1 = ["InR_1","chr3R",17442775,17444378]
InR_2 = ["InR_2","chr3R",17441278,17442794]
InR_3 = ["InR_3","chr3R",17439554,17441296]
InR_4 =  ["InR_4", "chr3R",17438003,17439578]
InR_9 = ["InR_9","chr3R",17431892,17433222]
InR_10 = ["InR_10", "chr3R",17430326,17431915]
InR_12 = ["InR_12","chr3R",17427085,17428696]
InR_15 = ["InR_15","chr3R",17422249,17423924]
InR_18 = ["InR_18","chr3R",17417848,17419060]
InR_19 = ["InR_19","chr3R",17416054,17417645]
InR_20 = ["InR_20","chr3R",17413929,17415666]
InR_22 = ["InR_22","chr3R",17410405,17412125]

find_files(InR_region)
find_files(foxo_region)
find_files(Thor_region)
find_files(InR_TU)
find_files(InR_exon)
find_files(InR_T1)
find_files(InR_1)
find_files(InR_2)
find_files(InR_3)
find_files(InR_4)
find_files(InR_9)
find_files(InR_10)
find_files(InR_12)
find_files(InR_15)
find_files(InR_18)
find_files(InR_20)
find_files(InR_22)


