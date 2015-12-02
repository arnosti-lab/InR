#This script takes adjusted vcf files and assembles a dictionary of indels
#indels that occur less than twice are dropped
#remaining indels are used to calculate
#size of indels
#frequency within population

import glob
from collections import Counter
from itertools import dropwhile
from decimal import *
getcontext().prec = 2

filelist = glob.glob("files/adjusted*.vcf")
filelist_length = len(filelist)


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
    sorted(indel_list.items())
    return(indel_list)

#from the counter, outputs the entire counter to a txt file
#outputs counts (frequency of indel among all genomes/files) to a bedgraph
#outputs the length of the indel (- values are deletions, + values are insertions) to a bedgraph
def all_region(list_of_indels, chrdict, outfilename, num_of_genomes):
    outfile1 = open(outfilename + ".txt", "w")
    outfile2 = open(outfilename + "_size" + ".bedgraph", "w")
    outfile3 = open(outfilename + "_freq" + ".bedgraph", "w")    
    for i in list_of_indels.keys():
        line = i.split("\t")[0:4]
        if line[0] in chrdict.keys():
            outfile1.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(list_of_indels[i]) + "\r\n")
            freq = 100 * Decimal(list_of_indels[i])/num_of_genomes
            this_indel = line[3].split("/")
            outfile3.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(freq) + "\r\n")
            outfile2.write(chrdict[line[0]] + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(len(this_indel[1]) - len(this_indel[0])) + "\r\n")
    outfile1.close()
    outfile2.close()
    outfile3.close()
    

indel_list = organize_files(filelist)
all_region(indel_list, chrdict, "African_indels", filelist_length)
        