#read bedfile, get regions of interest
#for each region of interest, look through summary files for regions that overlap, add the blastz score and coverage
#report region, percentage of coverage, and average blastz score
#It also plots a stacked barchart and a phylogenetic tree
#this script requires installation of matplotlib, numpy, scipy, and Biopython

import re
import string
import matplotlib
import scipy
import Bio
from Bio import Phylo
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from pylab import *

temp_file = open("temp_file.txt", "w")

#This is the worst function I've ever written :-(
def get_regions(bedfile, outfilename, additional_info):
    outfile = open(outfilename, "w")
    outfile.write("region_name" + "," + "simulans" + "," + "simulans_sd" + "," + "simulans_coverage" + "," \
                   + "yakuba" + "," + "yakuba_sd" +  "," + "yakuba_coverage" + "," + "pseudoobscura" +"," \
                    + "pseudoobscura_sd" + "," + "pseudoobscura_coverage" + "\r\n")
    file_list = []
    enhancer_list = []
    sim_list = []
    sim_cov = []
    sech_list = []
    sech_cov = []
    yak_list = []
    yak_cov = []
    ere_list = []
    ere_cov = []
    ana_list = []
    ana_cov = []
    pse_list = []
    pse_cov = []
    wil_list = []
    wil_cov = []
    vir_list = []
    vir_cov = []
    gri_list = []
    gri_cov = []
    species_list = ["sim", "sech", "yak", "ere", "ana", "pse", "wil", "vir", "gri"]
    myfile = open(bedfile)
    for line in myfile:
        enhancer_coord1 = 0
        enhancer_coord2 = 0
        if "chr" in line:
            newline = line.split("\t")
            sim_file = newline[0] + ".dm3.droSim1.summary.txt"
            sech_file = newline[0] + ".dm3.droSec1.summary.txt"
            yak_file = newline[0] + ".dm3.droYak2.summary.txt"
            ere_file = newline[0] + ".dm3.droEre2.summary.txt"
            ana_file = newline[0] + ".dm3.droAna3.summary.txt"
            pse_file = newline[0] + ".dm3.dp4.summary.txt"
            wil_file = newline[0] + ".dm3.droWil1.summary.txt"
            vir_file = newline[0] + ".dm3.droVir3.summary.txt"
            gri_file = newline[0] + ".dm3.droGri2.summary.txt"
            file_list = [sim_file, sech_file, yak_file, ere_file, ana_file, pse_file, wil_file, vir_file, gri_file]
            enhancer_coord1 = int(newline[1])
            enhancer_coord2 = int(newline[2])
            enhancer_label = newline[3]
            enhancer_list.append(enhancer_label)
            temp_file.write(enhancer_label)
            temp = [1,2,3]
            for i in file_list:
                current_name = i.strip("\.summary\.txt")
                temp_file.write(current_name)
                current_name = re.sub(".+dm3.", "", current_name)
                print current_name
                blast_z_list = []
                coverage = 0
                count = 0
                aligned = 0
                thisfile = open(i)
                #if the lines overlap with the CRM, it's determined how many base pairs from the region are in this line (fragment_length)
                #the blastz score for this line is added to the blastz list
                #blastz scores are multiplied by the fraction of the region they refer to
                for line in thisfile:
                    eachline = line.split(" ")
                    if enhancer_coord1 >= int(eachline[3]) or enhancer_coord2 <= int(eachline[2]):
                        pass
                    else:
                        count = count + 1
                        sequence_length = float(int(eachline[3]) - int(eachline[2]))
                        print "sequence length: ", sequence_length
                        temp_file.write(line)
                        fragment_length = float(int(eachline[3]) - int(eachline[2]))
                        if sequence_length == 0:
                            break
                        if int(eachline[3]) > enhancer_coord2:
                           fragment_length = fragment_length - (int(eachline[3]) - enhancer_coord2)
                        if int(eachline[2]) < enhancer_coord1:
                            fragment_length = fragment_length - (enhancer_coord1 - int(eachline[2])) 
                        aligned = aligned + fragment_length
                        blast_z_list.append(float(eachline[8].strip()) * (float(fragment_length)/float(int(eachline[3]) - int(eachline[2]))))
                coverage = 100 * float(aligned)/float(enhancer_coord2 - enhancer_coord1) 
                if aligned == 0:
                    coverage = 0
                    blast_z_list = [0 for i in blast_z_list] 
                    if current_name == "droSim1":
                        sim_list.append(0)
                        sim_cov.append(0)
                    elif current_name == "droSec1":
                        sech_list.append(0)
                        sech_cov.append(0)
                    elif current_name == "droYak2":
                        yak_list.append(0)
                        yak_cov.append(0)
                    elif current_name == "droEre2":
                        ere_list.append(0)
                        ere_cov.append(0)
                    elif current_name == "droAna3":
                        ana_list.append(0)
                        ana_cov.append(0)
                    elif current_name == "dp4":
                        pse_list.append(0)
                        pse_cov.append(0)
                    elif current_name == "droWil1":
                        wil_list.append(0)
                        wil_cov.append(0)
                    elif current_name == "droVir3":
                        vir_list.append(0)
                        vir_cov.append(0)
                    elif current_name == "droGri2":
                        gri_list.append(0)
                        gri_cov.append(0)              
                else:
                    blast_z_list = [float(i)/(float((enhancer_coord2 - enhancer_coord1)/100)) for i in blast_z_list]
                    mean_blastz_list = np.asarray(blast_z_list)
                    print np.sum(mean_blastz_list), np.mean(mean_blastz_list)
                    if current_name == "droSim1":
                        sim_list.append(np.sum(mean_blastz_list))
                        sim_cov.append(coverage)
                        print len(mean_blastz_list)
                    elif current_name == "droSec1":
                        sech_list.append(np.sum(mean_blastz_list))
                        sech_cov.append(coverage)
                    elif current_name == "droYak2":
                        yak_list.append(np.sum(mean_blastz_list))
                        yak_cov.append(coverage)
                    elif current_name == "droEre2":
                        ere_list.append(np.sum(mean_blastz_list))
                        ere_cov.append(coverage)
                    elif current_name == "droAna3":
                        ana_list.append(np.sum(mean_blastz_list))
                        ana_cov.append(coverage)
                    elif current_name == "dp4":
                        pse_list.append(np.sum(mean_blastz_list))
                        pse_cov.append(coverage)
                    elif current_name == "droWil1":
                        wil_list.append(np.sum(mean_blastz_list))
                        wil_cov.append(coverage)
                    elif current_name == "droVir3":
                        vir_list.append(np.sum(mean_blastz_list))
                        vir_cov.append(coverage)
                    elif current_name == "droGri2":
                        gri_list.append(np.sum(mean_blastz_list))
                        gri_cov.append(coverage)                
                outfile.write("," + str(np.sum(mean_blastz_list)) + "," + str(coverage))
                temp_file.write(str(aligned) + "," + str(coverage) + "\r\n" )
    outfile.close()
    myfile.close()
    return_list_blastz = [sim_list, sech_list, yak_list, ere_list, ana_list, pse_list, wil_list, vir_list, gri_list]
    return_list_cov = [sim_cov, sech_cov, yak_cov, ere_cov, ana_cov, pse_cov, wil_cov, vir_cov, gri_cov]
    return(return_list_blastz, return_list_cov, enhancer_list)

def make_plots(blastz_info, cov_info, current_regions, additional_info):
    N = len(blastz_info[0])
    ind = np.arange(N)
    width = 0.75
    rot = 65
    mydata = [np.array(log10(i)) for i in blastz_info]
    #mystd = [np.array(log10(i)) for i in std_blast]
    mycov = [np.array(i) for i in cov_info]
    fig,a = plt.subplots()
    for i in range(len(mycov)):
        if len(mycov[i]) <= 10:
            rot = 20
        else:
            rot = 65
        bottom = 100 * i
        color = plt.cm.Dark2(i * 20)
        p1 = a.bar(ind, mycov[i], color = color, bottom = bottom, width = width)
        plt.xticks(ind + width/2, current_regions, rotation = rot)
        a.yaxis.set_visible(False)
        plt.ylabel('Percent coverage by species')
    with PdfPages(additional_info + "_cov.pdf") as pdf:
        pdf.savefig()
    plt.close()
    fig2,ax = plt.subplots()
    for j in range(len(mydata)):
        if len(mycov[i]) <= 10:
            rot = 20
        else:
            rot = 65
        if j == 0:
            bottom = 0
        else:
            bottom = (j * 4)
        print "bottom = ", bottom
        color = plt.cm.Dark2(j* 20)
        p2 = plt.bar(ind, mydata[j], color = color, bottom = bottom, width = 0.75)
        plt.ylim([0,36])
        ax.tick_params(axis='y', direction = 'in', length=2, width=0.4)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        ax.set_yticks([4, 8, 12, 16, 20, 24, 28, 32, 36])
        ax.set_yticklabels(['']*9)
        #for tick in ax.yaxis.get_major_ticks():
        #    tick.label.set_fontsize(7) 
        plt.xticks(ind + width/2, current_regions, rotation = rot)
        plt.ylabel('log transformed blastZ over 100bp')
    with PdfPages(additional_info + "_newblastz.pdf") as pdf:
        pdf.savefig()
    plt.close()
    print "length is ", len(mydata)


#blastz_scores, cov_scores, regions, = get_regions("InR_sequences.bedgraph", "test_InR_seqs.tsv", "genes_exons")
#make_plots(blastz_scores, cov_scores, regions, "genes_exonsTEST")
blastz_scores2, cov_scores2, regions2, = get_regions("InR_crms.bedgraph", "test_InR_crms.tsv", "InR_regions")
make_plots(blastz_scores2, cov_scores2, regions2, "InR_regionsTEST")
#get_regions("InR_crms.bedgraph", "test_InR.tsv", "InR_crms")

tree = Phylo.read("Drosophila.dnd", "newick")
tree.rooted = True
Phylo.draw_ascii(tree)

#((grimshawi:7, virilis:7):1, (wilistoni:7, (pseudoobscura:6, (ananassae:5, (erecta:4, (yakuba:3, ((sechellia:1, simulans:1):1, melanogaster:2):1):1):1):1):1):1);




tree.root.color = "k"
tree.clade[0].color = "#a6d96a"
tree.clade[0,1].color = "#66a61e"

tree.clade[1].color = "#980043"
tree.clade[1,1].color = "#e7298a"
tree.clade[1,1,1].color = "#7570b3"
tree.clade[1,1,1,1].color = "#f768a1"
tree.clade[1,1,1,1,1].color = "#d95f02"
tree.clade[1,1,1,1,1,1,1].color = "k"
tree.clade[1,1,1,1,1,1,0].color = "#a6761d"
tree.clade[1,1,1,1,1,1,0,1].color = "#1b9e77"




Phylo.draw(tree)
with PdfPages("Drosophila_phylogeny.pdf") as pdf:
    pdf.savefig()


