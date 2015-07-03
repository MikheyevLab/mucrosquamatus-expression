# Examining co-expression of 

'''
# Software requirements
bowtie/1.1.0 
rsem RSEM v1.2.11
'''

#sbatch "snakemake -j 999 -p --cluster-config cluster.json --cluster \"sbatch  -p {cluster.partition} -n {cluster.n}\"" 

REF="../ref/Trinity_ERCC.fa" # Trinity assembly with ERRCC92 spike-ins
RSEM_REF="../ref/Trinity_ERCC.transcripts.fa" #index produced by RSEM
SAMPLES, = glob_wildcards("../data/reads/{sample}.fastq")
CALLER=["freebayes", "GATK", "samtools"]
PYTHON="/apps/free/python/2.7.8/bin/python"
RBB = {"Pe":"Pm","Pm":"Pe"}   #dictionary to match species names for reciprocal best BLAST

rule all:
	input: "../output/dnds.txt","../data/popgen/filtered.recode.vcf" # "../ouput/rbb.txt" #expand("../output/{species}/blast.txt",species=["Pe","Pm"]) #"../out/pe2pm.xml","../out/pm2pe.xml" #"../data/popgen/consensus.vcf", expand("../data/popgen/{VCFcaller}.vcf", VCFcaller=CALLER) #  expand("../data/popgen/{sample}.bam", sample=SAMPLES)

rule rsem_prepare_reference:
	input: "../ref/Trinity_ERCC.fa"
	output: "../ref/Trinity_ERCC.idx.fa"
	version: "1.0"
	shell: """module load bowtie/1.1.0  ; \
	grep ">" ../ref/Trinity_ERCC.fa  |cut  -d" " -f1 |tr -d ">" | awk -F_ 'NF>1 {{print $1"_"$2"\t"$1"_"$2"_"$3}} NF==1 {{print $1"\t"$1}}' > ../ref/transcript-to-gene-map.txt; \
	rsem-prepare-reference --transcript-to-gene-map ../ref/transcript-to-gene-map.txt --no-polyA ../ref/Trinity_ERCC.fa ../ref/Trinity_ERCC"""

rule rsemCalculateExpression:
	input: "../data/reads/{sample}.fastq"
	output: "../data/rsem/{sample}.genes.results"
	version: "1.0"
	shell: """module load bowtie/1.1.0  ; \
	rsem-calculate-expression -p 8 ../data/reads/{wildcards.sample}.fastq ../ref/Trinity_ERCC {wildcards.sample}; \
	mv {wildcards.sample}* ../data/rsem/ """


#gather rsem results and make tables of counts and fpkm (samples in columns, genes in rows)
rule collectRsem:
	input: expand("../data/rsem/{sample}.genes.results", sample=SAMPLES)
	output: "../out/counts.csv.gz", "../out/fpkm.csv.gz"
	version: "1.0"
	shell: ". python2/bin/activate ;\
	python collect_counts.py genes ../data/rsem/ | gzip > {output[0]}; \
	python collect_fpkm.py genes ../data/rsem/ | gzip > {output[1]}"


#### Population genomics

# Make protein predictions for abundant P. elegans and P. mucrosquamatus transcripts

rule transdecoder:
	# this uses only the first Trinity isoform of both assemblies
	input: "../ref/{species}_nucl.fasta"
	output: "../output/{species}/{species}_nucl.fasta.transdecoder.pep" 
	shell: "TransDecoder -t {input} -m 70 --reuse --workdir ../output/{wildcards.species} --search_pfam /apps/unit/MikheyevU/sasha/TransDecoder_r20140704/pfam/Pfam-AB.hmm.bin --CPU 16 -v && mv {wildcards.species}_nucl.fasta.transdecoder.* ../output/{wildcards.species}"

## reciprocal BLAST of Protobothrops elegans and Protobothrops mucrosquamatus
rule blastDB:	
	input: "../output/{species}/{species}_nucl.fasta.transdecoder.pep"
	output: "../output/{species}/blast.phr"
	shell: "module load ncbi-blast/2.2.30+; makeblastdb -in {input} -out ../output/{wildcards.species}/blast -dbtype prot"

rule blast:	
	input: "../output/{species}/blast.phr"
	output: "../output/{species}/blast.txt"
	params: target=lambda wildcards: RBB[wildcards.species]
	shell: "module load ncbi-blast/2.2.30+; blastp -query ../output/{wildcards.species}/{wildcards.species}_nucl.fasta.transdecoder.pep -db ../output/{params.target}/blast -num_threads 16 -evalue 1e-8 -word_size 3 -max_target_seqs 1 -outfmt 6 -out {output}" 

rule combineBlast:
	input: expand("../output/{species}/blast.txt",species=["Pe","Pm"])
	output: "../ouput/rbb.txt"
	shell: """awk '(NR == FNR) && !($1 in seen) {{seen[$1]; query[$2]=$1; next}} ($1 in query) && ($2==query[$1]) {{print $1"\t"$2}}' {input} > {output}"""

rule dnds:
	input: "../ouput/rbb.txt"
	output: "../output/dnds.txt"
	shell: "{PYTHON} dnds.py"

#### Genetic differentiation

rule addReadGroups:
	#add read groups to rsem output
	input: "../data/rsem/{sample}.transcript.sorted.bam"
	output: temp("../data/popgen/{sample}.bam")
	version: "1.0"
	shell: "java -jar /apps/unit/MikheyevU/picard-tools-1.66/AddOrReplaceReadGroups.jar I={input} O={output} ID={wildcards.sample} LB=NEXTERA PL=ILLUMINA RGPU=X RGSM={wildcards.sample}" 

rule addIndexesToRsem:
	input: "../ref/Trinity_ERCC.transcripts.fa"
	output: "../ref/Trinity_ERCC.transcripts.dict", "../ref/Trinity_ERCC.transcripts.fa.fai"
	shell: "java -jar /apps/unit/MikheyevU/picard-tools-1.66/CreateSequenceDictionary.jar R={input} O={output[0]}; \
	samtools faidx {input}"

rule mergeBams:
	#add read groups to rsem output
	input: expand("../data/popgen/{sample}.bam",sample=SAMPLES)
	output: "../data/popgen/merged.bam"
	version: "1.0"
	shell: "novosort -t 8 --ram 22G -a -o {output} -i {input}" 

rule freeBayes:
	input:	"../data/popgen/merged.bam", "../ref/Trinity_ERCC.transcripts.fa.fai"
	output: protected("../data/popgen/freebayes.vcf")
	shell:	"freebayes --use-best-n-alleles 2 -b {input[0]} -v {output} -f {RSEM_REF}"

rule GATK:
	input:	"../data/popgen/merged.bam","../ref/Trinity_ERCC.transcripts.dict", "../ref/Trinity_ERCC.transcripts.fa.fai"
	output: protected("../data/popgen/GATK.vcf")
	shell:	"java  -Xmx30g -jar $GATK -U -nct 12 -allowPotentiallyMisencodedQuals  -T HaplotypeCaller -R  {RSEM_REF} -I {input[0]} -hets 0.002  -mbq 20 -o {output} --max_alternate_alleles 2"

rule samtools:
	input:	"../data/popgen/merged.bam", "../ref/Trinity_ERCC.transcripts.fa.fai"
	output: protected("../data/popgen/samtools.vcf")
	shell: "samtools mpileup -ugf {RSEM_REF} {input[0]} | bcftools call -vc - | vcfutils.pl varFilter -D 500 > {output}"

rule allelicPrimitives:
	input: "../data/popgen/{VCFcaller}.vcf"
	output: temp("../data/popgen/{VCFcaller}.primitives.vcf")
	shell: "java -Xmx14g -jar $GATK -T VariantsToAllelicPrimitives -R {RSEM_REF} --variant {input} -o {output}"

# generate consensus SNP calls
rule BAYSIC: 	
	input: expand("../data/popgen/{VCFcaller}.primitives.vcf", VCFcaller=CALLER)
	output: "../data/popgen/consensus.vcf"
	version: "1.0"
	run: 
		infiles = "".join([" --vcf " + i for i in input])
		shell("baysic.pl --statsOutFile ../data/popgen/combined.stats --pvalCutoff 0.8 {} --countsOutFile ../data/popgen/combined.cts --vcfOutFile {{output}}".format(infiles))

rule filterVCF:
		input: rules.BAYSIC.output
		output: "../data/popgen/filtered.recode.vcf"
		version: "1.0"
		shell: """grep -m1 ^#C ../data/popgen/consensus.vcf |cut -f10-39 | tr "\t" "\n" > ../data/popgen/names.txt; \
		vcftools --vcf ../data/popgen/consensus.vcf --max-alleles 2  --max-missing 1.0 --minGQ 20  --keep test.txt --recode --out  ../data/popgen/filtered"""

