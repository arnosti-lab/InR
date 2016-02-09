#this script gets fasta format regions based on coordinates from John Pool's Genome Nexus seq files

import glob

def get_seq(infile, index, endex,seqname, myseq):
    fastaname = ">" +  infile.strip(".seq")+ "_" + seqname + "\r\n"
    seqlength = endex - index
    myfile = open(infile)
    sequence = ""
    coord_count = 0
    for line in myfile:
        line = line.strip()
        sequence = sequence + line.strip()
    myseq = myseq + fastaname +  sequence[index-1:index + seqlength]+ "\r\n"
    myfile.close()
    return myseq

outseq = ""

def get_file(region):
    outfilename = region[3] + ".fas"
    outfile = open(outfilename, "w")
    myfasta = ""
    for name in glob.glob("*.seq"):
        if region[0] not in name:
            pass
        else:
            myfasta = get_seq(name, region[1], region[2], "_"+ region[3], myfasta)
            print "done with ", name
    outfile.write(myfasta)
    outfile.close()
    
    	

InR_region = ["3R",17345970,17495043, "Insulin_locus"]
foxo_region = ["3R",9857667,9938700, "foxo_locus"]
Thor_region = ["2L",3453434,3504612, "Thor_locus"]
InR_exon1 = ["3R",17395970,17401638, "InR_exon1"]
InR_TU = ["3R",17395970,17445043, "InR_TU"]

chr3R_0 = ["3R",9000000,9150000, "Chr3R_0"]
chr3R_1 = ["3R",10000000,10150000, "Chr3R_1"]
chr3R_2 = ["3R",12500000,12650000, "Chr3R_2"]
chr3R_3 = ["3R",13800000,13950000, "Chr3R_3"]
chr3R_4 = ["3R",14200000,14350000, "Chr3R_4"]
chr3R_5 = ["3R",15000000,15150000, "Chr3R_5"]
chr3R_6 = ["3R",16600000,16750000, "Chr3R_6"]
chr3R_7 = ["3R",17100000,17250000, "Chr3R_7"]
chr3R_8 = ["3R",18400000,18550000, "Chr3R_8"]
chr3R_9 = ["3R",19800000,19950000, "Chr3R_9"]


#get_file(InR_region)
#get_file(foxo_region)
#get_file(Thor_region)
get_file(chr3R_0)
#get_file(InR_exon1)
#get_file(InR_TU)
get_file(chr3R_1)
get_file(chr3R_2)
get_file(chr3R_3)
get_file(chr3R_4)
get_file(chr3R_5)
get_file(chr3R_6)
get_file(chr3R_7)
get_file(chr3R_8)
get_file(chr3R_9)