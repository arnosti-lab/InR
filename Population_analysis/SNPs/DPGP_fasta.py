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
    outfilename = region[3] + ".fa"
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

get_file(InR_region)
get_file(foxo_region)
get_file(Thor_region)