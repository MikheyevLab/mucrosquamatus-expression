# Examining co-expression of 

'''
# Software requirements
bowtie/1.1.0 
rsem RSEM v1.2.11
'''

#sbatch "snakemake -j 999 -p --cluster-config cluster.json --cluster \"sbatch  -p {cluster.partition} -n {cluster.n}\"" 

REF="../ref/Trinity_ERCC.fa" # Trinity assembly with ERRCC92 spike-ins
SAMPLES, = glob_wildcards("../data/reads/{sample}.fastq")
CALLER=["freebayes", "GATK", "samtools"]  

rule all:
	input: "../data/popgen/consensus.vcf" #  expand("../data/popgen/{sample}.bam", sample=SAMPLES)

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

#### Genetic differentiation

# rule getAbundantTranscripts:
# 	# get sequences for the most abundant transcripts
# 	input: "../out/fpkm.csv.gz"
# 	output: "../ref/abundant.txt"
# 	shell: """zcat {input} |awk -F, 'NR>1 {{sum=0; for(i=2;i<=NF;i++) if($i<20) next;  print $1"_i1"}}' |grep -v ERCC > {output}"""

rule addReadGroups:
	#add read groups to rsem output
	input: "../data/rsem/{sample}.transcript.sorted.bam"
	output: temp("../data/popgen/{sample}.bam")
	version: "1.0"
	shell: "java -jar /apps/unit/MikheyevU/picard-tools-1.66/AddOrReplaceReadGroups.jar I={input} O={output} ID={wildcards.sample} LB=NEXTERA PL=ILLUMINA RGPU=X RGSM={wildcards.sample}" 

rule mergeBams:
	#add read groups to rsem output
	input: expand("../data/popgen/{sample}.bam",sample=SAMPLES)
	output: "../data/popgen/merged.bam"
	version: "1.0"
	shell: "novosort -t 8 --ram 22G -a -o {output} -i {input}" 

rule freeBayes:
	input:	"../data/popgen/merged.bam"
	output: protected("../data/popgen/{CALLER[0]}.vcf")
	shell:	"freebayes --use-best-n-alleles 2 -b {input} -v {output} -f {REF}"

rule GATK:
	input:	"../data/popgen/merged.bam"
	output: protected("../data/popgen/{CALLER[1]}.vcf")
	shell:	"java  -Xmx30g -jar $GATK -nct 12 -allowPotentiallyMisencodedQuals  -T HaplotypeCaller -R  {REF} -I {input} -hets 0.002  -mbq 20 -o {output} --max_alternate_alleles 2"

rule samtools:
	input:	"../data/popgen/merged.bam"
	output: protected("../data/popgen/{CALLER[2]}.vcf")
	shell: "samtools mpileup -ugf {REF} {input} | bcftools call -vc - | vcfutils.pl varFilter -D 500 > {output}"

rule allelicPrimitives:
	input: "../data/popgen/{VCFcaller}.vcf"
	output: temp("../data/popgen/{VCFcaller}.primitives.vcf")
	shell: "java -Xmx14g -jar $GATK -T VariantsToAllelicPrimitives -R {REF} --variant {input} -o {output}"

# generate consensus SNP calls
rule BAYSIC: 	
	input: expand("../data/popgen/{VCFcaller}.primitives.vcf", VCFcaller=CALLER)
	output: "../data/popgen/consensus.vcf"
	version: "1.0"
	run: 
		infiles = "".join([" --vcf " + i for i in input])
		shell("baysic.pl --statsOutFile ../data/popgen/combined.stats --pvalCutoff 0.8 {} --countsOutFile ../data/popgen/combined.cts --vcfOutFile {{output}}".format(infiles))