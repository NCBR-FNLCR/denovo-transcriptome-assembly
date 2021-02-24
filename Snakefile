
###########################################################################
# Short read denovo transcriptome assembly workflow
#
# Creator: Neelam Redekar, neelam.redekar@nih.gov
# Created: April 7, 2020
# Version 2
###########################################################################

from os.path import join
from snakemake.io import expand, glob_wildcards
from snakemake.utils import R
from os import listdir

configfile: "config.yaml"

rawdata_dir= config["rawdata_dir"]
working_dir= config["result_dir"]
taxid= config["TAXID"]
Lineage_name= config["Lineage_name"]
Lineage= config["Lineage"]

SAMPLES, = glob_wildcards(join(rawdata_dir, '{samples}_R1_001.fastq.gz'))
#SAMPLES, = glob_wildcards(join(rawdata_dir, '{samples}.R1.fastq.gz'))

RN = ['R1', 'R2']
rule All:
    input:
        # Creating data links:
        expand(join(working_dir, "raw/{samples}_{rn}_001.fastq.gz"), samples=SAMPLES, rn=RN),
        #expand(join(working_dir, "raw/{samples}_{rn}.fastq.gz"), samples=SAMPLES, rn=RN),

        # FastQC (before and after pre-processing)
        expand(join(working_dir, "rawQC/{samples}_{rn}_001_fastqc.html"), samples=SAMPLES, rn=RN),
        #expand(join(working_dir, "rawQC/{samples}_{rn}.fastqc.html"), samples=SAMPLES, rn=RN),
        expand(join(working_dir, "QC/{samples}_{rn}.trim.cor.filt_fastqc.html"), samples=SAMPLES, rn=RN),

        # Preprocessing: Adapter removal, quality trimming, k-mer based error correction
        expand(join(working_dir, "trimmed_reads/{samples}_{rn}.trim.fastq.gz"), samples=SAMPLES, rn=RN),
        expand(join(working_dir, "corrected_reads/{samples}_{rn}.trim.cor.fq"), samples=SAMPLES, rn=RN),
        expand(join(working_dir, "corrected_reads/{samples}_{rn}.trim.cor.filt.fq.gz"), samples=SAMPLES, rn=RN),

        # FastScreen (Using two sets of reference databases)
        expand(join(working_dir,"FQscreen/{samples}_{rn}.trim.cor.filt_screen.txt"), samples=SAMPLES, rn=RN),
        expand(join(working_dir,"FQscreen/{samples}_{rn}.trim.cor.filt_screen.png"), samples=SAMPLES, rn=RN),
        expand(join(working_dir,"FQscreen2/{samples}_{rn}.trim.cor.filt_screen.txt"), samples=SAMPLES, rn=RN),
        expand(join(working_dir,"FQscreen2/{samples}_{rn}.trim.cor.filt_screen.png"), samples=SAMPLES, rn=RN),

        # Kraken - Krona (removing read contamination)
        expand(join(working_dir,"kraken/{samples}.trim.fastq.kraken_bacteria.taxa.txt"), samples=SAMPLES),
        expand(join(working_dir,"kraken/{samples}.trim.fastq.kraken_bacteria.krona.html"), samples=SAMPLES),

        # Trinity + Transfuse (generate sample assemblies and merge to get assemblyA)
        expand(join(working_dir, "sample_assemblies/{samples}/Trinity.fasta"), samples=SAMPLES),
        expand(join(working_dir, "sample_assemblies/{samples}.fasta"), samples=SAMPLES),
        join(working_dir, "assemblyA/merged.fa"),

        # Trinity (generate assemblyB)
        join(working_dir, "assemblyB/trinity-out-dir/Trinity.fasta"),
        join(working_dir, "assemblyB/merged.assembly-stats.txt"),

        # Transrate, Quast, BUSCO (evaluate sample assemblies, assemblyA, assemblyB)
        ## Sample assemblies
        join(working_dir,"sample-transrate/assemblies.csv"),
        join(working_dir,"sample-quast/report.html"),
        expand(join(working_dir,"busco/{samples}/short_summary.specific.glires_odb10.{samples}.txt"), samples=SAMPLES),

        ## assemblyA
        join(working_dir,"mergeA-quast/report.html"),
        join(working_dir,"mergeA-transrate/assemblies.csv"),
        join(working_dir,"busco-merged/mergedA/short_summary.specific.glires_odb10.mergedA.txt"),

        ## assemblyB
        join(working_dir,"mergeB-quast/report.html"),
        join(working_dir,"mergeB-transrate/assemblies.csv"),
        join(working_dir,"busco-merged/mergedB/short_summary.specific.glires_odb10.mergedB.txt"),

        # BUSCO summaries
        join(working_dir,"alignment/trinity_asm_paired.sorted.bam"),
        join(working_dir,"busco-summaries/busco_figure.png"),

        # Cluster transcripts, finding orthologs
        join(working_dir, "clusters/asmB_cluster_0.99.fasta"),
        join(working_dir, "clusters/asmB_cluster_0.95.fasta"),
        join(working_dir, "clusters/asmA_cluster_0.99.fasta"),
        join(working_dir, "clusters/asmA_cluster_0.95.fasta"),

    output:
        "multiqc_report.html"
    params:
        projname="practice",
        dir=directory("Reports"),
    shell:
        """
module load multiqc/1.8
mkdir -p {params.dir}
cd {params.dir}
multiqc .
        """

rule raw_data_links:
    input:
        join(rawdata_dir, "{samples}_{rn}_001.fastq.gz")
        #join(rawdata_dir, "{samples}.{rn}.fastq.gz")
    output:
        join(working_dir, "raw/{samples}_{rn}_001.fastq.gz")
        #join(working_dir, "raw/{samples}_{rn}.fastq.gz")
    params:
        rname="raw_data_links",
        dir=directory("raw"),
    shell:
        """
mkdir -p {params.dir}
ln -s {input} {output}
        """

rule raw_fastqc:
    input:
        join(working_dir, "raw/{samples}_{rn}_001.fastq.gz")
        #join(working_dir, "raw/{samples}_{rn}.fastq.gz")
    output:
        join(working_dir, "rawQC/{samples}_{rn}_001_fastqc.html")
        #join(working_dir, "rawQC/{samples}_{rn}.fastqc.html")
    params:
        rname="raw_fastqc",
        dir=directory("rawQC"),
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
    threads: 32
    shell:
        """
module load fastqc/0.11.9
mkdir -p {params.dir}
fastqc -o {params.dir} -f fastq --threads {threads} --extract {input}
        """

rule rmv_adapter:
    input:
        F1=join(working_dir, "raw/{samples}_R1_001.fastq.gz"),
        #F1=join(working_dir, "raw/{samples}_R1.fastq.gz"),
        F2=join(working_dir, "raw/{samples}_R2_001.fastq.gz"),
        #F2=join(working_dir, "raw/{samples}_R2.fastq.gz"),
    output:
        R1=join(working_dir, "trimmed_reads/{samples}_R1.trim.fastq.gz"),
        R2=join(working_dir, "trimmed_reads/{samples}_R2.trim.fastq.gz")
    params:
        rname="rmv_adapter",
        dir=directory("trimmed_reads"),
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        fastawithadaptersetd="/data/CCBR_Pipeliner/db/PipeDB/dev/TruSeq_and_nextera_adapters.consolidated.fa"
    threads: 32
    log:
        join(working_dir, "trimmed_reads/{samples}.log")
    shell:
        """
module load cutadapt/2.9
mkdir -p {params.dir}
cutadapt --pair-filter=any --nextseq-trim=2 --trim-n -n 5 -o 5 -q 10,10 -m 35:35 -b file:{params.fastawithadaptersetd} -B file:{params.fastawithadaptersetd} -j {threads} -o {output.R1} -p {output.R2} {input.F1} {input.F2} > {log}
        """

rule rcorrector:
    input:
        F1=join(working_dir, "trimmed_reads/{samples}_R1.trim.fastq.gz"),
        F2=join(working_dir, "trimmed_reads/{samples}_R2.trim.fastq.gz")
    output:
        R1un=join(working_dir, "corrected_reads/{samples}_R1.trim.cor.fq"),
        R2un=join(working_dir, "corrected_reads/{samples}_R2.trim.cor.fq"),
    params:
        rname="rcorrector",
        dir=directory(join(working_dir, "corrected_reads")),
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        prefix='{samples}',
    threads: 72
    shell:
        """
module load rcorrector/1.0.3.1
module load python/3.7
mkdir -p {params.dir}
run_rcorrector.pl -1 {input.F1} -2 {input.F2} -od {params.dir} -t {threads}
gunzip {output.R1un}.gz
gunzip {output.R2un}.gz
        """


rule filter:
    input:
        F1=join(working_dir, "corrected_reads/{samples}_R1.trim.cor.fq"),
        F2=join(working_dir, "corrected_reads/{samples}_R2.trim.cor.fq")
    output:
        F1=join(working_dir, "corrected_reads/{samples}_R1.trim.cor.filt.fq"),
        F2=join(working_dir, "corrected_reads/{samples}_R2.trim.cor.filt.fq"),
    params:
        rname="filter",
        dir=directory(join(working_dir, "corrected_reads")),
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
        prefix='{samples}',
        script=join(working_dir, "Scripts/FilterUncorrectabledPEfastq.py"),
        F1=join(working_dir, "corrected_reads/unfixrm_{samples}_R1.trim.cor.fq"),
        F2=join(working_dir, "corrected_reads/unfixrm_{samples}_R2.trim.cor.fq")
    threads: 72
    shell:
        """
cd {params.dir}
module load python/3.7
python {params.script} -1 {input.F1} -2 {input.F2} -s {params.prefix}
cp {params.F1} {output.F1}
cp {params.F2} {output.F2}
        """


rule merge:
    input:
        F1=expand(join(working_dir, "corrected_reads/{samples}_R1.trim.cor.filt.fq"), samples=SAMPLES),
        F2=expand(join(working_dir, "corrected_reads/{samples}_R2.trim.cor.filt.fq"), samples=SAMPLES)
    output:
        M1=join(working_dir, "merged/merged.R1.trim.cor.filt.fq"),
        M2=join(working_dir, "merged/merged.R2.trim.cor.filt.fq"),
        R1=join(working_dir, "merged/merged.R1.trim.cor.filt.fq.gz"),
        R2=join(working_dir, "merged/merged.R2.trim.cor.filt.fq.gz"),
    params:
        rname="merge",
        dir=join(working_dir, "merged"),
    shell:
        """
mkdir -p {params.dir}
cat {input.F1} > {output.M1}
cat {input.F2} > {output.M2}
gzip -c {output.M1} > {output.R1}
gzip -c {output.M2} > {output.R2}
        """

rule compress:
    input:
        #M1=join(working_dir, "merged/merged.R1.trim.cor.filt.fq"),
        #M2=join(working_dir, "merged/merged.R2.trim.cor.filt.fq"),
        F1=join(working_dir, "corrected_reads/{samples}_R1.trim.cor.filt.fq"),
        F2=join(working_dir, "corrected_reads/{samples}_R2.trim.cor.filt.fq"),
        #R1un=join(working_dir, "corrected_reads/{samples}_R1.trim.cor.fq"),
        #R2un=join(working_dir, "corrected_reads/{samples}_R2.trim.cor.fq"),
    output:
        T1=join(working_dir, "corrected_reads/{samples}_R1.trim.cor.filt.fq.gz"),
        T2=join(working_dir, "corrected_reads/{samples}_R2.trim.cor.filt.fq.gz"),
        #R1=join(working_dir, "merged/merged.R1.trim.cor.filt.fq.gz"),
        #R2=join(working_dir, "merged/merged.R2.trim.cor.filt.fq.gz"),
    params:
        rname="compress",
    shell:
        """
gzip -c {input.F1} > {output.T1}
gzip -c {input.F2} > {output.T2}
        """

rule fastqc:
    input:
        join(working_dir, "corrected_reads/{samples}_{rn}.trim.cor.filt.fq.gz")
    output:
        join(working_dir, "QC/{samples}_{rn}.trim.cor.filt_fastqc.html")
    params:
        rname="fastqc",
        dir=directory("QC"),
        batch='--cpus-per-task=32 --mem=110g --time=48:00:00',
    threads: 72
    shell:
        """
module load fastqc/0.11.9
mkdir -p {params.dir}
fastqc -o {params.dir} -f fastq --threads {threads} --extract {input}
module load python/3.5
python ./Scripts/get_read_length.py {params.dir} > {params.dir}/readlength.txt  2> {params.dir}/readlength.err
        """

rule fastq_screen:
    input:
        F1=join(working_dir, "corrected_reads/{samples}_R1.trim.cor.filt.fq.gz"),
        F2=join(working_dir, "corrected_reads/{samples}_R2.trim.cor.filt.fq.gz")
    output:
        out1=join(working_dir, "FQscreen/{samples}_R1.trim.cor.filt_screen.txt"),
        out2=join(working_dir, "FQscreen/{samples}_R1.trim.cor.filt_screen.png"),
        out3=join(working_dir, "FQscreen/{samples}_R2.trim.cor.filt_screen.txt"),
        out4=join(working_dir, "FQscreen/{samples}_R2.trim.cor.filt_screen.png"),
        out5=join(working_dir, "FQscreen2/{samples}_R1.trim.cor.filt_screen.txt"),
        out6=join(working_dir, "FQscreen2/{samples}_R1.trim.cor.filt_screen.png"),
        out7=join(working_dir, "FQscreen2/{samples}_R2.trim.cor.filt_screen.txt"),
        out8=join(working_dir, "FQscreen2/{samples}_R2.trim.cor.filt_screen.png")
    params:
        rname="fastq_screen",
        batch='--cpus-per-task=24 --mem=64g --time=10:00:00',
        fastq_screen="/data/CCBR_Pipeliner/db/PipeDB/bin/fastq_screen_v0.9.3/fastq_screen",
        dir1 = directory("FQscreen"),
        dir2 = directory("FQscreen2"),
        fastq_screen_config="/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen.conf",
        fastq_screen_config2="/data/CCBR_Pipeliner/db/PipeDB/lib/fastq_screen_2.conf",
    threads: 72
    shell:
        """
module load bowtie/2-2.3.4
module load perl/5.24.3
{params.fastq_screen} --conf {params.fastq_screen_config} --outdir {params.dir1} --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.F1} {input.F2}
{params.fastq_screen} --conf {params.fastq_screen_config2} --outdir {params.dir2} --threads {threads} --subset 1000000 --aligner bowtie2 --force {input.F1} {input.F2}
        """


rule kraken_pe:
    input:
        F1=join(working_dir, "corrected_reads/{samples}_R1.trim.cor.filt.fq.gz"),
        F2=join(working_dir, "corrected_reads/{samples}_R2.trim.cor.filt.fq.gz")
    output:
        krakentaxa = join(working_dir, "kraken/{samples}.trim.fastq.kraken_bacteria.taxa.txt"),
        kronahtml = join(working_dir, "kraken/{samples}.trim.fastq.kraken_bacteria.krona.html")
    params:
        rname="kraken_pe",
        prefix = "{samples}",
        dir=directory("kraken"),
        bacdb="/fdb/kraken/20170202_bacteria"
    threads: 72
    shell:
        """
module load kraken/1.1
module load kronatools/2.7
if [ ! -d {params.dir} ];then mkdir {params.dir};fi

cd /lscratch/$SLURM_JOBID;
cp -rv {params.bacdb} /lscratch/$SLURM_JOBID/;

kraken --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` --fastq-input --gzip-compressed --threads {threads} --output /lscratch/$SLURM_JOBID/{params.prefix}.krakenout --preload --paired {input.F1} {input.F2}

kraken-translate --mpa-format --db /lscratch/$SLURM_JOBID/`echo {params.bacdb}|awk -F "/" '{{print \$NF}}'` /lscratch/$SLURM_JOBID/{params.prefix}.krakenout |cut -f2|sort|uniq -c|sort -k1,1nr > /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa

cut -f 2,3 /lscratch/$SLURM_JOBID/{params.prefix}.krakenout | ktImportTaxonomy - -o /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml

mv /lscratch/$SLURM_JOBID/{params.prefix}.krakentaxa {output.krakentaxa}
mv /lscratch/$SLURM_JOBID/{params.prefix}.kronahtml {output.kronahtml}
        """

rule sample_assembly:
    input:
        F1=join(working_dir, "corrected_reads/{samples}_R1.trim.cor.filt.fq.gz"),
        F2=join(working_dir, "corrected_reads/{samples}_R2.trim.cor.filt.fq.gz")
    output:
        FA=join(working_dir, "sample_assemblies/{samples}/Trinity.fasta"),
        FA2=join(working_dir, "sample_assemblies/{samples}.fasta"),
        ST=join(working_dir, "sample_assemblies/{samples}.assembly-stats.txt"),
    params:
        rname="sample_assembly",
        dir1=directory("sample_assemblies"),
        dir=directory("sample_assemblies/{samples}-trinity"),
        mem="100G",
        prefix='{samples}'
    threads: 72
    shell:
        """
module load trinity/2.9.0
test -d {params.dir1} || mkdir -p {params.dir1}
Trinity --seqType fq --max_memory {params.mem} --left {input.F1} --right {input.F2} --include_supertranscripts --CPU {threads} --output {params.dir}

perl /usr/local/apps/trinity/trinityrnaseq-2.9.0/util/TrinityStats.pl {output.FA} > {output.ST}
cp {output.FA} {output.FA2}
        """

rule stats_transrate:
    input:
        asm=expand(join(working_dir,"sample_assemblies/{samples}.fasta"), samples=SAMPLES),
        M1=join(working_dir, "merged/merged.R1.trim.cor.filt.fq"),
        M2=join(working_dir, "merged/merged.R2.trim.cor.filt.fq")
    output:
        ST=join(working_dir,"sample-transrate/assemblies.csv"),
    params:
        rname="stats_transrate",
        batch='--cpus-per-task=50 --mem=150g --time=10:00:00',
        dir=directory("sample-transrate"),
    threads: 50
    shell:
        """
module load transrate/1.0.3
mkdir -p {params.dir}
transrate --output {params.dir} --assembly `echo '{input.asm}'| sed 's/ /,/g'` --threads {threads} --left {input.M1} --right {input.M2}
        """

rule stats_quast:
    input:
        asm=expand(join(working_dir,"sample_assemblies/{samples}.fasta"), samples=SAMPLES),
    output:
        ST=join(working_dir,"sample-quast/report.html"),
    params:
        rname="stats_quast",
        batch='--cpus-per-task=72 --mem=80g --time=10:00:00',
        dir=directory("sample-quast")
    threads: 72
    shell:
        """
module unload python
module load quast/5.0.2
module load circos/0.69-9
quast.py -o {params.dir} -t {threads} --circos -L {input.asm}
        """

rule stats_busco:
    input:
        asm=join(working_dir, "sample_assemblies/{samples}.fasta"),
    output:
        ST=join(working_dir,"busco/{samples}/short_summary.specific.glires_odb10.{samples}.txt"),
    params:
        rname="stats_busco",
        dir=directory(join(working_dir, "busco")),
        folder="{samples}",
        lineage=Lineage,
        #lineage="~/References/Odb10/Vertebrata/eutheria/euarchontoglires/glires_odb10"
    threads: 32
    shell:
        """
module load busco/4.0.2
mkdir -p {params.dir}
cd {params.dir}
busco --offline -m transcriptome -l {params.lineage} --cpu {threads} -i {input.asm} -f -o {params.folder}
        """

rule merge_sample_assemblies:
    input:
        asm=expand(join(working_dir, "sample_assemblies/{samples}.fasta"), samples=SAMPLES),
        M1=join(working_dir, "merged/merged.R1.trim.cor.filt.fq"),
        M2=join(working_dir, "merged/merged.R2.trim.cor.filt.fq")
    output:
        asm=join(working_dir, "assemblyA/merged.fa"),
    params:
        rname="merge_sample_assemblies",
        transfuse="~/transfuse-0.5.0-linux-x86_64/transfuse",
        batch='--cpus-per-task=72 --mem=50g --time=2-00:00:00, partition=largemem,norm',
        dir = directory("assemblyA"),
    threads: 72
    shell:
        """
mkdir -p {params.dir}
cd {params.dir}
{params.transfuse} -a `echo '{input.asm}'| sed 's/ /,/g'` -t {threads} --left {input.M1} --right {input.M2} -o {output.asm}
        """


rule denovo_asmB:
    input:
        F1=join(working_dir, "merged/merged.R1.trim.cor.filt.fq"),
        F2=join(working_dir, "merged/merged.R2.trim.cor.filt.fq")
    output:
        FA=join(working_dir, "assemblyB/trinity-out-dir/Trinity.fasta"),
        ST=join(working_dir, "assemblyB/merged.assembly-stats.txt"),
    params:
        rname="denovo_asmB",
        dir=directory("assemblyB"),
        dir2=directory("assemblyB/trinity-out-dir"),
        mem="100G",
    threads: 72
    shell:
        """
module load trinity/2.9.0
mkdir -p {params.dir}
Trinity --seqType fq --max_memory {params.mem} --left {input.F1} --right {input.F2} --include_supertranscripts --CPU {threads} --output {params.dir2}

perl /usr/local/apps/trinity/trinityrnaseq-2.9.0/util/TrinityStats.pl {output.FA} > {output.ST}
        """


rule mergeA_quast:
    input:
        asm=join(working_dir,"assemblyA/merged.fa"),
    output:
        ST=join(working_dir,"mergeA-quast/report.html"),
    params:
        rname="mergeA_quast",
        batch='--cpus-per-task=72 --mem=80g --time=10:00:00',
        dir=directory("mergeA-quast")
    threads: 72
    shell:
        """
module unload python
module load quast/5.0.2
module load circos/0.69-9
quast.py -o {params.dir} -t {threads} --circos -L {input.asm}
        """


rule mergeA_transrate:
    input:
        asm=join(working_dir,"assemblyA/merged.fa"),
        M1=join(working_dir,"merged/merged.R2.trim.cor.filt.fq"),
        M2=join(working_dir,"merged/merged.R2.trim.cor.filt.fq"),
    output:
        ST=join(working_dir,"mergeA-transrate/assemblies.csv"),
    params:
        rname="mergeA_transrate",
        batch='--cpus-per-task=50 --mem=150g --time=10:00:00',
        dir=directory("mergeA-transrate")
    threads: 50
    shell:
        """
module load transrate/1.0.3
mkdir -p {params.dir}
transrate --output {params.dir} --assembly `echo '{input.asm}'| sed 's/ /,/g'` --threads {threads} --left {input.M1} --right {input.M2}
        """

rule mergeA_busco:
    input:
        asm=join(working_dir,"assemblyA/merged.fa"),
    output:
        ST=join(working_dir,"busco-merged/mergedA/short_summary.specific.glires_odb10.mergedA.txt"),
    params:
        rname="mergeA_busco",
        dir=directory(join(working_dir, "busco-merged")),
        folder="mergedA",
        lineage=Lineage,
        #lineage="~/References/Odb10/Vertebrata/eutheria/euarchontoglires/glires_odb10"
    threads: 32
    shell:
        """
module load busco/4.0.2
cd {params.dir}
busco --offline -m transcriptome -l {params.lineage} --cpu {threads} -i {input.asm} -f -o {params.folder}
        """

rule mergeB_quast:
    input:
        asm=join(working_dir,"assemblyB/trinity-out-dir/Trinity.fasta"),
    output:
        ST=join(working_dir,"mergeB-quast/report.html"),
    params:
        rname="mergeB_quast",
        batch='--cpus-per-task=72 --mem=80g --time=10:00:00',
        dir=directory("mergeB-quast")
    threads: 72
    shell:
        """
module unload python
module load quast/5.0.2
module load circos/0.69-9
quast.py -o {params.dir} -t {threads} --circos -L {input.asm}
        """

rule mergeB_transrate:
    input:
        asm=join(working_dir,"assemblyB/trinity-out-dir/Trinity.fasta"),
        M1=join(working_dir,"merged/merged.R2.trim.cor.filt.fq"),
        M2=join(working_dir,"merged/merged.R2.trim.cor.filt.fq"),
    output:
        ST=join(working_dir,"mergeB-transrate/assemblies.csv"),
    params:
        rname="mergeB_transrate",
        batch='--cpus-per-task=50 --mem=150g --time=10:00:00',
        dir=directory("mergeB-transrate")
    threads: 50
    shell:
        """
module load transrate/1.0.3
mkdir -p {params.dir}
transrate --output {params.dir} --assembly `echo '{input.asm}'| sed 's/ /,/g'` --threads {threads} --left {input.M1} --right {input.M2}
        """

rule mergeB_busco:
    input:
        asm=join(working_dir,"assemblyB/trinity-out-dir/Trinity.fasta"),
    output:
        ST=join(working_dir,"busco-merged/mergedB/short_summary.specific.glires_odb10.mergedB.txt"),
    params:
        rname="mergeB_busco",
        dir=directory(join(working_dir, "busco-merged")),
        folder="mergedB",
        lineage=Lineage,
        #lineage="~/References/Odb10/Vertebrata/eutheria/euarchontoglires/glires_odb10"
    threads: 32
    shell:
        """
module load busco/4.0.2
cd {params.dir}
busco --offline -m transcriptome -l {params.lineage} --cpu {threads} -i {input.asm} -f -o {params.folder}
        """


rule alignment_to_asm:
    input:
        asm=join(working_dir,"assemblyB/trinity-out-dir/Trinity.fasta"),
        F1=join(working_dir, "merged/merged.R2.trim.cor.filt.fq"),
        F2=join(working_dir, "merged/merged.R2.trim.cor.filt.fq")
    output:
        sam=join(working_dir,"alignment/trinity_asm_paired.sam"),
        bam=join(working_dir,"alignment/trinity_asm_paired.bam"),
        sort=join(working_dir,"alignment/trinity_asm_paired.sorted.bam"),
    params:
        rname="alignment_to_asm",
        index="mergedAsm-bowtie2",
        dir=directory("alignment")
    threads: 32
    shell:
        """
module load bowtie/2-2.3.5
mkdir -p {params.dir}
cd {params.dir}
bowtie2-build {input.asm} {params.index}
bowtie2 -q -N 1 -L 25 --phred64 --minins 1 --maxins 1000 --all --mm --threads 20 -x {params.index} -U {input.F1},{input.F2} -S {output.sam} 1>bowtie_trinity_asm_backmapping.log 2>bowtie_trinity_asm_backmapping.err

module load samtools/1.10
samtools view -bS {output.sam} > {output.bam}
samtools sort {output.bam} -o {output.sort}
samtools flagstat {output.bam} > trinity_backmap_paired_flagstat.txt

module load picard/2.1.1
java -Xmx4G -jar $PICARDJARPATH/picard.jar CollectInsertSizeMetrics I={output.bam} O=insert_size_output.txt H=insert_size_output.pdf M=0.5

java -Xmx4G -jar $PICARDJARPATH/picard.jar CollectWgsMetrics I={output.bam} O=collect_wgs_metrics.txt R={input.asm}
        """


rule cluster_merged_transcripts:
    input:
        fastaB=join(working_dir, "assemblyB/trinity-out-dir/Trinity.fasta"),
        fastaA=join(working_dir, "assemblyA/merged.fa"),
    output:
        fastaB1=join(working_dir, "clusters/asmB_cluster_0.99.fasta"),
        fastaB2=join(working_dir, "clusters/asmB_cluster_0.95.fasta"),
        fastaA1=join(working_dir, "clusters/asmA_cluster_0.99.fasta"),
        fastaA2=join(working_dir, "clusters/asmA_cluster_0.95.fasta"),
    params:
        rname="cluster_merged_transcripts",
        batch='--mem=200g --cpus-per-task=72 --partition=largemem,norm',
        mem=200000,
        dir=directory(join(working_dir, "clusters"))
    threads: 72
    shell:
        """
module load cd-hit/4.8.1
mkdir -p {params.dir}
cd-hit-est -i {input.fastaA} -o clusters/asmA_cluster_0.99 -c 0.99 -n 8 -p 1 -g 1 -M {params.mem} -T {threads} -d 40 1>asmA-cdhit-est_0.99.log 2>asmA-cdhit-est_0.99.err
mv {params.dir}/asmA_cluster_0.99 {output.fastaA1}

cd-hit-est -i {input.fastaA} -o clusters/asmA_cluster_0.95 -c 0.95 -n 8 -p 1 -g 1 -M {params.mem} -T {threads} -d 40 1>asmA-cdhit-est_0.95.log 2>asmA-cdhit-est_0.95.err
mv {params.dir}/asmA_cluster_0.95 {output.fastaA2}

cd-hit-est -i {input.fastaB} -o clusters/asmB_cluster_0.99 -c 0.99 -n 8 -p 1 -g 1 -M {params.mem} -T {threads} -d 40 1>asmB-cdhit-est_0.99.log 2>asmB-cdhit-est_0.99.err
mv {params.dir}/asmB_cluster_0.99 {output.fastaB1}

cd-hit-est -i {input.fastaB} -o clusters/asmB_cluster_0.95 -c 0.95 -n 8 -p 1 -g 1 -M {params.mem} -T {threads} -d 40 1>asmB-cdhit-est_0.95.log 2>asmB-cdhit-est_0.95.err
mv {params.dir}/asmB_cluster_0.95 {output.fastaB2}

        """

rule busco_summaries:
    input:
        sampl=expand(join(working_dir,"busco/{samples}/short_summary.specific.glires_odb10.{samples}.txt"), samples=SAMPLES),
        asmA=join(working_dir,"busco-merged/mergedA/short_summary.specific.glires_odb10.mergedA.txt"),
        asmB=join(working_dir,"busco-merged/mergedB/short_summary.specific.glires_odb10.mergedB.txt"),
    output:
        join(working_dir,"busco-summaries/busco_figure.png"),
    params:
        rname="busco_summaries",
        dir=directory(join(working_dir, "busco-summaries")),
    threads: 32
    shell:
        """
module load busco/4.0.2
mkdir -p {params.dir}
#cd {params.dir}
cp {input.sampl} {params.dir}
cp {input.asmA} {params.dir}
cp {input.asmB} {params.dir}
python3 /usr/local/apps/busco/4.0.2/generate_plot.py -rt specific –wd {params.dir}
        """


rule copy:
    input:
        fastaA = join(working_dir, "assemblyA/trinity-out-dir/Trinity.fasta"),
        fastaB = join(working_dir, "assemblyB/trinity-out-dir/Trinity.fasta"),
        sampleAsm = expand(join(working_dir, "sample_assemblies/{samples}.fasta"), samples=SAMPLES),
    output:
        join(working_dir, "FASTA-assemblies", "MergedA.fa"),
        join(working_dir, "FASTA-assemblies", "MergedB.fa"),
        expand(join(working_dir, "FASTA-assemblies/{samples}.fasta"), samples=SAMPLES),
    params:
        dir=join(working_dir, "FASTA-assemblies"),
    shell:
        """
mkdir -p {params.dir}
cp {input.fastaA} {params.dir}/MergedA.fa
cp {input.fastaB} {params.dir}/MergedB.fa
cp {input.sampleAsm} {params.dir}/
        """

rule blastn_SampleAsm:
    input:
        sampleAsm=join(working_dir, "FASTA-assemblies/{samples}.fasta"),
    output:
        sampleAsm_out=join(working_dir, "BLASTn-sample/{samples}.blastn.dcmegablast.", taxid,".txt"),
    params:
        rname="blastn_SampleAsm",
        dir=join(working_dir, "BLASTn"),
        db_nt="/fdb/blastdb/nt",
        taxid=taxid,
    threads: 56
    shell:
        """
module load blast/2.10.0+

mkdir -p {params.dir}
blastn -task megablast -db {params.db_nt} -taxids {params.taxid} -query {input.sampleAsm} -out {output.sampleAsm_out} -html -max_hsps 1 -template_type coding -outfmt '6 qseqid sseqid sacc qstart qend sstart send evalue bitscore score length pident nident mimatch positive gapopen gaps qcovss qcovhsp stitle' -num_threads {threads}
        """

rule blastn_MergedAsm:
    input:
        fastaA = join(working_dir, "FASTA-assemblies/MergedA.fa"),
        fastaB = join(working_dir, "FASTA-assemblies/MergedB.fa"),
    output:
        fastaA_out=join(working_dir, "BLASTn-sample/MergedA.blastn.dcmegablast.", taxid,".txt"),
        fastaB_out=join(working_dir, "BLASTn-sample/MergedB.blastn.dcmegablast.", taxid,".txt"),
    params:
        rname="blastn_MergedAsm",
        dir=join(working_dir, "BLASTn"),
        db_nt="/fdb/blastdb/nt",
        taxid=taxid,
    threads: 56
    shell:
        """
module load blast/2.10.0+

mkdir -p {params.dir}

blastn -task megablast -db {params.db_nt} -taxids {params.taxid} -query {input.fastaA} -out {output.fastaA_out} -html -max_hsps 1 -template_type coding -outfmt '6 qseqid sseqid sacc qstart qend sstart send evalue bitscore score length pident nident mimatch positive gapopen gaps qcovss qcovhsp stitle' -num_threads {threads}

blastn -task megablast -db {params.db_nt} -taxids {params.taxid} -query {input.fastaB} -out {output.fastaB_out} -html -max_hsps 1 -template_type coding -outfmt '6 qseqid sseqid sacc qstart qend sstart send evalue bitscore score length pident nident mimatch positive gapopen gaps qcovss qcovhsp stitle' -num_threads {threads}
        """


