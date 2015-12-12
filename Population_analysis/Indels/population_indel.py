#This script takes adjusted vcf files and assembles a dictionary of indels
#indels that occur less than twice are dropped
#remaining indels are used to calculate
#size of indels
#frequency within population

import glob
import string
from collections import Counter
from itertools import dropwhile
from decimal import *
getcontext().prec = 2
import matplotlib
import scipy
import Bio
from Bio import Phylo
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from pylab import *


filelist = glob.glob("files/adjusted*.vcf")
filelist_length = len(filelist)

InR_region = ["InR","chr3R",17345970,17495043]
foxo_region = ["foxo","chr3R",9857667,9938700]
Thor_region = ["Thor","chr2L",3453434,3504612]
InR_TU = ["InR_TU","chr3R",17395970,17445043]
InR_exon = ["InR_exon","chr3R",17399113,17401638]
InR_T1 = ["InR_T1","chr3R",17444777,17445965]
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


list_of_regions = [InR_region, foxo_region, Thor_region, InR_TU, InR_exon, InR_T1, InR_1, InR_2, InR_3, InR_4, InR_9, InR_10, InR_12, InR_15, InR_18, InR_19, InR_20, InR_22]

chr2R_len = 21100000
chr2L_len = 23000000
chr3L_len = 24500000
chr3R_len = 27900000
chrX_len =  22400000
    

chrdict = {'3': 'chr2L', '4':'chrX', '5':'chr3L', '6': 'chr4', '7':'chr2R', '8':'chr3R'}


#adds indel lines from a file to a really big list
def search_files(filename, biglist):
    thisfile = open(filename)
    for singleline in thisfile:
        biglist.append(singleline.strip())
    thisfile.close()
    return(biglist)

#iterates through all of the files using search_files and adds every line to the big list
#turns the big list to a counter class, and drops any indel variants that occur less than twice
#sorts counter (badly)
def organize_files(filelist):
    superlist = []
    count = 0
    for i in filelist:
        superlist = search_files(i, superlist)
        count = count + 1
    print "tallying"
    indel_list = Counter(superlist)
    print len(indel_list), " total indels"
    for key, count in dropwhile(lambda key_count: key_count[1] >= 2, indel_list.most_common()):
        del indel_list[key]
    print len(indel_list), " after dropping singletons"
    indel_list2 = indel_list
    for key, count in dropwhile(lambda key_count: key_count[1] >= 30, indel_list2.most_common()):
        del indel_list2[key]    
    sorted(indel_list.items())
    return(indel_list, indel_list2)

#from the counter, outputs the entire counter to a txt file
#outputs counts (frequency of indel among all genomes/files) to a bedgraph
#outputs the length of the indel (- values are deletions, + values are insertions) to a bedgraph
def all_region(list_of_indels, chrdict, outfilename, num_of_genomes):
    mylist = []
    outfile1 = open(outfilename + ".txt", "w")
    outfile2 = open(outfilename + "_size" + ".bedgraph", "w")
    outfile3 = open(outfilename + "_freq" + ".bedgraph", "w") 
    outfile2.write("browser hide all" + "\r\n" + 'browser pack refGene encodeRegions' + "\r\n" + 'browser full altGraph' + "\r\n")
    outfile3.write("browser hide all" + "\r\n" + 'browser pack refGene encodeRegions' + "\r\n" + 'browser full altGraph' + "\r\n")
    outfile2.write("track type=bedGraph name='size " + outfilename + " description='BedGraph format' visibility=full color=60,0,0 altColor=0,60,0 priority=20" + "\r\n")
    outfile3.write("track type=bedGraph name='frequency'" + outfilename + " description='BedGraph format' visibility=full color=0,60,0 altColor=0,60,0 priority=20" + "\r\n")
    chr2L_outfile = open("chr2L_indel_freq.bedgraph", "w") 
    chr2R_outfile = open("chr2R_indel_freq.bedgraph", "w")  
    chr3L_outfile = open("chr3L_indel_freq.bedgraph", "w") 
    chr3R_outfile = open("chr3R_indel_freq.bedgraph", "w") 
    chrX_outfile = open("chrX_indel_freq.bedgraph", "w")  
    for i in list_of_indels.keys():
        line = i.split("\t")[0:4]
        if line[0] in chrdict.keys():
            outfile1.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(list_of_indels[i]) + "\r\n")
            freq = list_of_indels[i]
            this_indel = line[3].split("/")
            outfile3.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str((100 *freq)/num_of_genomes) + "\r\n")
            outfile2.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(len(this_indel[1]) - len(this_indel[0])) + "\r\n")
            mylist.append([chrdict[line[0]], line[1], line[2], freq, len(this_indel[1]) - len(this_indel[0])])
            if chrdict[line[0]] == "chr2L":
                chr2L_outfile.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(freq) + "\r\n")
            elif chrdict[line[0]] == "chr2R":
                chr2R_outfile.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(freq) + "\r\n")
            elif chrdict[line[0]] == "chr3L":
                chr3L_outfile.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(freq) + "\r\n")
            elif chrdict[line[0]] == "chr3R":
                chr3R_outfile.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(freq) + "\r\n")
            elif chrdict[line[0]] == "chrX":
                chrX_outfile.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(freq) + "\r\n")
    outfile1.close()
    outfile2.close()
    outfile3.close()
    chr2L_outfile.close()
    chr2R_outfile.close()
    chr3L_outfile.close()
    chr3R_outfile.close()
    chrX_outfile.close()
    return(mylist)

def make_plots(indel_list, region_list):
    InR_list = []
    Thor_list = []
    foxo_list = []
    chrX_list = []
    chr2L_list = []
    chr2R_list = []
    chr3L_list = []
    chr3R_list = []
    InR_size = []
    Thor_size = []
    foxo_size = []
    chrX_size = []
    chr2L_size = []
    chr2R_size = []
    chr3L_size = []
    chr3R_size = []
    InR_count = 0
    Thor_count = 0
    foxo_count = 0
    chrX_count = 0
    chr2L_count = 0
    chr2R_count = 0
    chr3L_count = 0
    chr3R_count = 0
    for j in indel_list:
        if j[0] == "chrX":
            chrX_count = chrX_count + 1
            chrX_size.append(j[4])
            chrX_list.append(j[3])
        elif j[0] == "chr2L":
            chr2L_count = chrX_count + 1
            chr2L_size.append(j[4])
            chr2L_list.append(j[3])
        elif j[0] == "chr2R":
            chr2R_count = chrX_count + 1
            chr2R_size.append(j[4])
            chr2R_list.append(j[3])
        elif j[0] == "chr3L":
            chr3L_count = chrX_count + 1
            chr3L_size.append(j[4])
            chr3L_list.append(j[3])
        elif j[0] == "chr3R":
            chr3R_count = chrX_count + 1
            chr3R_size.append(j[4])
            chr3R_list.append(j[3])
        for i in region_list:
            if j[0] == i[1]:
                if (int(j[1]) > int(i[3])) or (int(i[2]) > int(j[2])):
                    pass
                else:
                    if i[0] == "InR":
                        InR_count = InR_count + 1
                        InR_size.append(j[4])
                        InR_list.append(j[3])
                    elif i[0] == "Thor":
                        Thor_count = Thor_count + 1
                        Thor_size.append(j[4])
                        Thor_list.append(j[3])
                    elif i[0] == "foxo":
                        foxo_count = foxo_count + 1
                        foxo_size.append(j[4])
                        foxo_list.append(j[3])
    count_list = [(100000 * chrX_count/chrX_len), (100000 * chr2L_count/chr2L_len),\
                  (100000 * chr2R_count/chr2R_len), (100000 * chr3L_count/chr3L_len), (100000 * chr3R_count/chr3R_len), \
                  (100000 * foxo_count/(region_list[1][3] - region_list[1][2])),(100000 * InR_count/(region_list[0][3] - region_list[0][2]))]
    freq_list = [chrX_list, chr2L_list, chr2R_list, chr3L_list, chr3R_list, foxo_list, InR_list]
    freq_list2 = [chrX_list, chr2R_list, chr3L_list, chr3R_list, foxo_list, InR_list]
    size_list = [chrX_size, chr2L_size, chr2R_size, chr3L_size, chr3R_size, foxo_size, InR_size]
    mycount = [np.array(i) for i in count_list]
    myfreq = [np.array(i) for i in freq_list]
    myfreq2 = [np.array(i) for i in freq_list2]
    mysize = [np.array(i) for i in size_list]
    colors = ['tan','red', 'lightgreen', 'pink', 'blue', 'cyan', 'lightblue'] 
    box = plt.boxplot(mysize, notch = True, patch_artist=True, sym = "")  
    #plt.ylim(-15, 20)    
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
    plt.xticks([1, 2, 3,4,5,6,7], ['chrX', 'chr2L','chr2R', 'chr3L', 'chr3R', 'InR', 'foxo'])
    savefig('Indel_size.png', bbox_inches='tight')
    plt.close()
    box = plt.boxplot(myfreq, notch=True, patch_artist=True, sym = "")
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
    plt.xticks([1, 2, 3, 4, 5, 6, 7], ['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'InR', 'foxo'])
    savefig('Indel_freq.png', bbox_inches='tight')
    plt.close()
    colors2 = ['tan','lightgreen', 'pink', 'blue', 'cyan', 'lightblue']
    box = plt.boxplot(myfreq2, notch=True, patch_artist=True, sym = "")
    for patch, color in zip(box['boxes'], colors2):
        patch.set_facecolor(color)
    plt.xticks([1, 2, 3, 4, 5, 6], ['chrX','chr2R', 'chr3L', 'chr3R', 'InR', 'foxo'])
    savefig('Indel_freq2.png', bbox_inches='tight')
    plt.close()
    ind = np.arange(7)  # the x locations for the groups
    width = 0.5       # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar((ind + 1), mycount, width, color= colors)
    plt.xticks([1, 2, 3, 4, 5, 6, 7], ['chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R', 'InR', 'foxo'])
    savefig('Indel_count.png', bbox_inches='tight')
    plt.close()




indel_list, indel_list2 = organize_files(filelist)
mylist = all_region(indel_list, chrdict, "African_indels", filelist_length)
make_plots(mylist, list_of_regions)
mylist2 = all_region(indel_list2, chrdict, "African_indels_frequent", filelist_length)
make_plots(mylist2, list_of_regions)