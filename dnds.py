from Bio.Phylo.PAML import codeml
from Bio import SeqIO
import os,sys,pdb

#read 
pe = SeqIO.to_dict(SeqIO.parse("../output/Pe/Pe_nucl.fasta.transdecoder.cds","fasta"))
pm = SeqIO.to_dict(SeqIO.parse("../output/Pm/Pm_nucl.fasta.transdecoder.cds","fasta"))

outfile = open("../output/dnds.txt","w")
scratch="/scratch/sasha"
if not os.path.exists(scratch):
    os.makedirs(scratch)
trefile = open(scratch+"/trefile","w")
trefile.write("(a,b);\n")
trefile.close()
for line in open("../ouput/rbb.txt"):
    rec1 = pe[line.split()[1]]
    rec2 = pm[line.split()[0]]
    gene = rec2.id.split("|")[0].replace("_i1","")
    rec1.id = "a"
    rec1.description=""
    rec2.id = "b"
    rec2.description=""
    SeqIO.write([rec1,rec2],scratch+"/sasha_seqs.fasta","fasta")
    os.system("prank -quiet -nomafft -codon -d={}/sasha_seqs.fasta -o={}/sasha_seqs.paml -f=paml ".format(scratch,scratch))
    #PAML analysis
    cml = codeml.Codeml()
    cml.alignment = "{}/sasha_seqs.paml.best.phy".format(scratch)
    cml.working_dir = scratch
    cml.tree = scratch+"/trefile"
    cml.out_file = scratch+"paml.out"
    cml.set_options(seqtype=1,
            verbose=1,
            noisy=0,
            model=1,
            runmode=-2,
            Mgene=0,
            NSsites=[0],
            CodonFreq=2,
            cleandata=0)
    try:
        cml.run(verbose=True)
    except:
        continue

    # parse cml.out_file
    for line in open(cml.out_file):
        if line.find("dN/dS=") > -1:
            line = line.split()
            outfile.write("%s\t%s\n" % (gene, line[line.index("dN/dS=")+1]))
outfile.close()
