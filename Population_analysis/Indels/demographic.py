#This script is intended to determine relative frequencies of indels within InR CRMs by population

import glob
import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
from pylab import *

YW_regions = glob.glob("Yiliang_regions/interesting_regions/*.txt")
Indel_files = glob.glob("files/adjusted*.vcf")
RAL_files = glob.glob("/mnt/home/sonnens2/genomes/PoolAlignments/Indel_files/DGRP/VCF_original/files/adjusted*.vcf")


all_indels = Indel_files + RAL_files




#this function determines populations associated with indels of interest
def find_pops(filelist):
    population_list = {}
    for i in filelist:
        genome = i.split("adjusted")[1].split("_")[1]
        population = genome[0:2]
        if population not in population_list.keys():
            population_list[population] = 1
        else:
            population_list[population] += 1
    return(population_list)

#for each indel, in each population, this function counts whether its allele is the wildtype or variant
def make_dict_list(new_pop_list, population, count):
    if population in new_pop_list.keys():
        if count == 0:
            new_pop_list[population][1] += 1
        else:
            new_pop_list[population][0] += 1
    else:
        if count == 0:
            new_pop_list[population] = [0,1]
        else:
            new_pop_list[population] = [1,0]
    return(new_pop_list, count)

#for a given ajdusted VCF file, this adds indels to dictionaries for each population, country and region
#using function make_dict_list
def search_vcf(indel, filelist2, pop_list):
    country_dict = {}
    country_dict["CK"] = "Congo"
    country_dict["CO"] = "Cameroon"
    country_dict["EA"] = "Ethiopia"
    country_dict["EB"] = "Ethiopia"
    country_dict["ED"] = "Ethiopia"
    country_dict["EF"] = "Ethiopia"
    country_dict["EG"] = "Egypt"
    country_dict["EM"] = "Ethiopia"
    country_dict["ER"] = "Ethiopia"
    country_dict["EZ"] = "Ethiopia"
    country_dict["FR"] = "France"
    country_dict["GA"] = "Gabon"
    country_dict["GU"] = "Guinea"
    country_dict["KM"] = "Kenya"
    country_dict["KN"] = "Kenya"
    country_dict["KO"] = "Kenya"
    country_dict["KR"] = "Kenya"
    country_dict["KT"] = "Kenya"
    country_dict["MW"] = "Malawi"
    country_dict["NG"] = "Nigeria"
    country_dict["RA"] = "United States"
    country_dict["RC"] = "Rwanda"
    country_dict["RG"] = "Rwanda"
    country_dict["SB"] = "South Africa"
    country_dict["SD"] = "South Africa"
    country_dict["SE"] = "South Africa"
    country_dict["SF"] = "South Africa"
    country_dict["SP"] = "South Africa"
    country_dict["TZ"] = "Tanzania"
    country_dict["UG"] = "Uganda"
    country_dict["UK"] = "Uganda"
    country_dict["UM"] = "Uganda"
    country_dict["ZI"] = "Zambia"
    country_dict["ZK"] = "Zimbabwe"
    country_dict["ZL"] = "Zambia"
    country_dict["ZO"] = "Zambia"
    country_dict["ZS"] = "Zimbabwe"
    region_dict = {"United States": "2_US", "France": "1_FR","Guinea": "3_W","Uganda":"4_E", "Kenya":"4_E", "Rwanda":"4_E",\
                   "Tanzania":"4_E", "Zimbabwe":"5_S", "Zambia":"5_S", "South Africa":"5_S", "Nigeria": "3_W", "Gabon" : "3_W",\
	           "Cameroon" : "3_W", "Congo":"3_W","Ethiopia":"4_E","Egypt":"N", "Malawi": "4_E"}
    vari = 0
    wildtype = 0
    ral_count = 0
    new_population_list = {}
    new_country_list = {}
    new_region_list = {}
    for i in filelist2:
        genome = i.split("adjusted")[1].split("_")[1]
        population = genome[0:2]
        indel_file = open(i)
        count1 = 0
        count2 = 0
        count3 = 0
        for line in indel_file:
            if indel[0] == line[0]:
                variant = line.split("\t")
                if (int(variant[1]) != indel[1]) or (int(variant[2]) != indel[2] or (variant[3] != indel[3])):
                    pass
                else:
                    temp2.write(genome + "," + population + "," + line)
                    count1 = 1
                    count2 = 1
                    count3 = 1
        indel_file.close()
        print population, country_dict[population], region_dict[country_dict[population]]
        new_population_list,count1 = make_dict_list(new_population_list, population, count1)
        new_country_list,count2 = make_dict_list(new_country_list, country_dict[population], count2)
        new_region_list,count2 = make_dict_list(new_region_list, region_dict[country_dict[population]], count3)
    return(new_population_list, new_country_list, new_region_list)

#this function plots output for a given indel from function each_indel 
def plot_indel_demo(pop_dict, CRM):
    N = len(pop_dict)
    print pop_dict
    pop_dict_keys = sorted(pop_dict, key=pop_dict.get)
    pop_dict_vals = [pop_dict[i] for i in pop_dict_keys]
    print N
    print pop_dict.values()
    print pop_dict.keys()
    print CRM
    variant = np.asarray([i[0] for i in pop_dict_vals])
    default = np.asarray([i[1] for i in pop_dict_vals])
    names = np.asarray(pop_dict_keys)
    ind = np.arange(N)  # the x locations for the groups
    width = 0.35       # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, variant, width, color='b')
    rects2 = ax.bar(ind + width, default, width, color='g')
    #ax.set_title(CRM)
    ax.set_xticks(ind + width)
    ax.set_xticklabels(names)
    ax.legend((rects1[0], rects2[0]), ('Variant', 'Wildtype'))
    fig_name = CRM
    print CRM
    plt.savefig(fig_name)
    plt.close()


#for a given indel of interest, this function
#determines the frequency of this indel in each 
#population, country and region
def each_indel(filelist, filelist2, pop_list):
    for thisfile in filelist:
        thisname = thisfile.replace("African_indels_all", "")
        thisname = thisname.replace(".txt", "")
        infile = open(thisfile)
        print thisfile
        count = 0
        for line in infile:
            this_indel = []
            if "chr3R" in line:
                count =count + 1
                indel_list = line.split("\t")
                this_indel = ["8", int(indel_list[1]),int(indel_list[2]),indel_list[4]]
                indel_size = str(indel_list[5])
                if "-" in indel_size:
                    indel_size = indel_size.replace("-", "deletion_")
                    indel_size = indel_size.strip()
                    indel_size = indel_size + "bp"
                else:
                    indel_size = indel_size.strip()
                    indel_size = "insertion_" + indel_size + "bp"
                CRM_name = thisname + "_" + str(count) + "_" + "freq_" + str(indel_list[3]) + "_" + indel_size
                CRM_name = CRM_name.replace("Yiliang_regions/interesting_regions/_InR", "InR")
                print CRM_name
                country_plot = CRM_name + "_country"
                region_plot = CRM_name + "_region"
                pop_info, country_info, region_info = search_vcf(this_indel, filelist2, pop_list)
                plot_indel_demo(pop_info, CRM_name)
                #plot_indel_demo(country_info, country_plot)
                plot_indel_demo(region_info, region_plot)
        infile.close()


population_list = find_pops(all_indels)
each_indel(YW_regions, all_indels, population_list)
