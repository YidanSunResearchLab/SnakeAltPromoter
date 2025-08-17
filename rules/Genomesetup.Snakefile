import os
import re
import shlex

#############################################
# Configuration
#############################################
script_dir = os.path.dirname(workflow.snakefile)
output_dir = config["output_dir"]
organism = config["organism"]
organism_fasta = config["organism_fasta"]
genes_gtf = config["genes_gtf"]
threads = config["threads"]

# Validate config
if not all([output_dir, organism, organism_fasta, genes_gtf, threads]):
    raise ValueError("Missing required config variables: output_dir, organism, organism_fasta, genes_gtf, threads")

# Function to map organism to species
def get_species(organism):
    species_map = {
        "hg": "Homo_sapiens",
        "mm": "Mus_musculus",
        "dm": "Drosophila_melanogaster",
        "rn": "Rattus_norvegicus",
        "dr": "Danio_rerio",
        "ce": "Caenorhabditis_elegans",
        "sc": "Saccharomyces_cerevisiae"
    }
    prefix = organism[:2].lower()
    species = species_map.get(prefix, "Homo_sapiens")
    if prefix not in species_map:
        print(f"Warning: Unrecognized organism prefix '{prefix}' in '{organism}'. Defaulting to 'Homo_sapiens'.")
    return species

# Ensure output directories
os.makedirs(output_dir + "/organisms/" + organism, exist_ok=True)
os.makedirs(output_dir + "/organisms/" + organism + "/Annotation", exist_ok=True)

#############################################
# Final Target
#############################################
rule all:
    input:
        expand(
            [
                output_dir + "/organisms/" + organism + "/" + organism + ".fa",
                output_dir + "/organisms/" + organism + "/" + organism + ".chrom.sizes",
                #output_dir + "/organisms/" + organism + "/HISAT2Index/genome.6.ht2",
                output_dir + "/organisms/" + organism + "/STARIndex/SAindex",
                #output_dir + "/organisms/SalmonIndex/seq.bin", change accordingly
                #output_dir + "/organisms/{organism}/SalmonIndex", couldl not create conda environment
                output_dir + "/organisms/{organism}/SalmonIndex/seq.bin",
                output_dir + "/organisms/" + organism + "/Annotation/genes_t2g.tsv",
                output_dir + "/organisms/" + organism + "/Annotation/genes.symbol",
                output_dir + "/organisms/" + organism + "/Annotation/genes.bed",
                output_dir + "/organisms/" + organism + "/Annotation/DEXSeq_flattened_exons.gff",
                output_dir + "/organisms/" + organism + "/Annotation/proActiv_promoter_annotation.rds",
                output_dir + "/organisms/" + organism + "/Annotation/proActiv_promoter.bed",
                #output_dir + "/organisms/" + organism + "/Annotation/DEXSeq_flattened_exons.gtf",
            ],
            organism=config.get("organism", ["hg38"])
         )

#############################################
# --- Genome Index Generation Rules ---
#############################################
rule copy_gtf:
    output:
        output_dir + "/organisms/{organism}/Annotation.gtf"
    params:
        gtf=genes_gtf
    conda: "envs/basic.yaml"
    log: "logs/copy_gtf_{organism}.log"
    benchmark: "benchmarks/copy_gtf_{organism}.bmk.txt"
    shell:
        """
        if [ "{params.gtf}" == "" ]; then
            echo "Error: genes_gtf is 'NA' or not provided." >&2
            exit 1
        fi
        cp {params.gtf} {output} 2> {log}
        """

rule copy_fasta:
    output:
        fasta=output_dir + "/organisms/{organism}/{organism}.fa",
        fai=output_dir + "/organisms/{organism}/{organism}.fa.fai"
    params:
        organism_dir=output_dir + "/organisms/{organism}",
        fasta=organism_fasta
    threads: threads
    conda: "envs/basic.yaml"
    log: "logs/copy_fasta_{organism}.log"
    benchmark: "benchmarks/copy_fasta_{organism}.bmk.txt"
    shell:
        """
        if [ "{params.fasta}" == "" ]; then
            echo "Error: organism_fasta is 'NA' or not provided." >&2
            exit 1
        fi
        mkdir -p {params.organism_dir}
        cp {params.fasta} {output.fasta}
        samtools faidx {output.fasta}
        python -c "with open('{output.fai}') as f: \
                main_chroms = [line.split()[0] for line in f if line.split()[0].startswith('chr')]; \
                print(' '.join(main_chroms))" > {params.organism_dir}/main_chroms.txt 2>> {log}
        samtools faidx {output.fasta} $(cat {params.organism_dir}/main_chroms.txt) > {params.organism_dir}/{wildcards.organism}_filtered.fa
        mv {params.organism_dir}/{wildcards.organism}_filtered.fa {output.fasta}
        samtools faidx {output.fasta}
        """

rule generate_chrom_sizes:
    input:
        fai=output_dir + "/organisms/{organism}/{organism}.fa.fai"
    output:
        output_dir + "/organisms/{organism}/{organism}.chrom.sizes"
    conda: "envs/basic.yaml"
    log: "logs/generate_chrom_sizes/{organism}.log"
    benchmark: "benchmarks/generate_chrom_sizes/{organism}.txt"
    shell:
        """
        cut -f 1,2 {input.fai} > {output} 2> {log}
        """

rule hisat2Index:
    input:
        output_dir + "/organisms/{organism}/{organism}.fa"
    output:
        output_dir + "/organisms/{organism}/HISAT2Index/genome.6.ht2"
    params:
        basedir=output_dir + "/organisms/{organism}/HISAT2Index"
    threads: min(10, threads)
    conda: "envs/basic.yaml"
    log: "logs/hisat2Index_{organism}.log"
    benchmark: "benchmarks/hisat2Index_{organism}.bmk.txt"
    shell:
        """
        mkdir -p {params.basedir}
        ln -sf {input} {params.basedir}/genome.fa
        hisat2-build -q -p {threads} {params.basedir}/genome.fa {params.basedir}/genome > {log} 2>&1
        """

rule starIndex:
    input:
        output_dir + "/organisms/{organism}/{organism}.fa"
    output:
        output_dir + "/organisms/{organism}/STARIndex/SAindex"
    params:
        basedir=output_dir + "/organisms/{organism}/STARIndex"
    threads: threads
    conda: "envs/basic.yaml"
    log: "logs/starIndex_{organism}.log"
    benchmark: "benchmarks/starIndex_{organism}.bmk.txt"
    shell:
        """
        mkdir -p {params.basedir}
        STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.basedir} --genomeFastaFiles {input} > {log} 2>&1
        rm -f {params.basedir}/Log.out
        """

rule salmon_index:
    input:
        genes_fasta=output_dir + "/organisms/{organism}/Annotation/transcripts.fa",
        genome_fasta=output_dir + "/organisms/{organism}/{organism}.fa"
    #Changed because Wildcards in input files cannot be determined from output files
    output:
        #decoys=output_dir + "/organisms/SalmonIndex/decoys.txt",
        decoys=output_dir + "/organisms/{organism}/SalmonIndex/decoys.txt",
        #seq=temp(output_dir + "/organisms/SalmonIndex/seq.fa"),
        seq=temp(output_dir + "/organisms/{organism}/SalmonIndex/seq.fa"),
        #index=output_dir + "/organisms/SalmonIndex/seq.bin"
        #Salmon’s index command produces a directory (e.g., seq.bin, pos.bin, info.json).
        #index=directory(output_dir + "/organisms/{organism}/SalmonIndex")
        index=output_dir + "/organisms/{organism}/SalmonIndex/seq.bin"
    params:
        #basedir=output_dir + "/organisms/SalmonIndex"
        basedir=output_dir + "/organisms/{organism}/SalmonIndex"
    threads: min(16, threads)
    conda: "envs/basic.yaml"
    #log: "salmon/logs/salmon_index.log"
    log: "salmon/logs/salmon_index_{organism}.log"
    #benchmark: "benchmarks/salmon_index.bmk"
    benchmark: "benchmarks/salmon_index_{organism}.bmk"
    shell:
        """
        mkdir -p {params.basedir}
        grep "^>" {input.genome_fasta} | cut -d ' ' -f 1 | tr -d '>' > {output.decoys}
        cat {input.genes_fasta} {input.genome_fasta} > {output.seq}
        salmon index -p {threads} -t {output.seq} -d {output.decoys} -i {params.basedir} --kmerLen 31 --keepDuplicates > {log} 2>&1
        """

rule gtf_to_files:
    input:
        gtf=output_dir + "/organisms/{organism}/Annotation.gtf",
        genome_fasta=output_dir + "/organisms/{organism}/{organism}.fa"
    output:
        t2g=output_dir + "/organisms/{organism}/Annotation/genes_t2g.tsv",
        symbol=output_dir + "/organisms/{organism}/Annotation/genes.symbol",
        bed=output_dir + "/organisms/{organism}/Annotation/genes.bed",
        transcripts=output_dir + "/organisms/{organism}/Annotation/transcripts.fa"
    conda: "envs/basic.yaml"
    log: "logs/gtf_to_files_{organism}.log"
    benchmark: "benchmarks/gtf_to_files_{organism}.bmk.txt"
    run:
        # Initialize outputs
        with open(output.t2g, "w") as t2g, open(output.symbol, "w") as symbol, open(output.bed, "w") as bed:
            GTFdict = {}
            # Parse GTF
            for line in open(input.gtf):
                if line.startswith("#"):
                    continue
                cols = line.strip().split("\t")
                annos = re.split(''';(?=(?:[^'"]|'[^']*'|"[^"]*")*$)''', cols[8])
                if cols[2] == "gene":
                    gene_id = None
                    gene_name = None
                    for anno in annos:
                        anno = shlex.split(anno.strip())
                        if not anno:
                            continue
                        if anno[0] == "gene_id":
                            gene_id = anno[1].strip('"')
                        elif anno[0] == "gene_name":
                            gene_name = anno[1].strip('"')
                    if gene_id:
                        symbol.write(f"{gene_id}\t{gene_name or ''}\n")
                elif cols[2] in ["transcript", "mRNA"]:
                    gene_id = None
                    transcript_id = None
                    gene_name = ""
                    for anno in annos:
                        anno = shlex.split(anno.strip())
                        if not anno:
                            continue
                        if anno[0] == "gene_id":
                            gene_id = anno[1].strip('"')
                        elif anno[0] == "transcript_id":
                            transcript_id = anno[1].strip('"')
                        elif anno[0] == "gene_name":
                            gene_name = anno[1].strip('"')
                    if transcript_id:
                        t2g.write(f"{transcript_id}\t{gene_id or ''}\t{gene_name}\n")
                        GTFdict[transcript_id] = [cols[0], cols[3], cols[4], cols[6], [], []]
                elif cols[2] == "exon":
                    transcript_id = None
                    for anno in annos:
                        anno = shlex.split(anno.strip())
                        if not anno:
                            continue
                        if anno[0] == "transcript_id":
                            transcript_id = anno[1].strip('"')
                    if transcript_id in GTFdict:
                        exon_width = int(cols[4]) - int(cols[3]) + 1
                        exon_offset = int(cols[3]) - int(GTFdict[transcript_id][1])
                        GTFdict[transcript_id][4].append(str(exon_width))
                        GTFdict[transcript_id][5].append(str(exon_offset))
            # Write BED file
            for transcript_id, v in GTFdict.items():
                v[5] = [int(x) for x in v[5]]
                v[4] = [int(x) for x in v[4]]
                block_sizes = [str(x) for _, x in sorted(zip(v[5], v[4]))]
                block_starts = sorted(v[5])
                block_starts = [str(x) for x in block_starts]
                bed.write(f"{v[0]}\t{v[1]}\t{v[2]}\t{transcript_id}\t.\t{v[3]}\t{v[1]}\t{v[2]}\t255,0,0\t{len(v[4])}\t{','.join(block_sizes)}\t{','.join(block_starts)}\n")
        # Generate transcript FASTA and validate
        shell("""
        gffread -w {output.transcripts} -g {input.genome_fasta} {input.gtf} > {log} 2>&1
        if [ ! -s {output.transcripts} ]; then
            echo "Error: transcripts.fa is empty or not created." >&2
            exit 1
        fi
        """)

rule prepare_promoter_annotation:
    input:
        gtf=output_dir + "/organisms/{organism}/Annotation.gtf"
    output:
        promoter_rds=output_dir + "/organisms/{organism}/Annotation/proActiv_promoter_annotation.rds",
        promoter_bed=output_dir + "/organisms/{organism}/Annotation/proActiv_promoter.bed",
        filtered_gtf=output_dir + "/organisms/{organism}/Annotation/proActiv_protein_coding.gtf",
        txdb=output_dir + "/organisms/{organism}/Annotation/proActiv_txdb.sqlite"
    params:
        species=lambda wildcards: get_species(wildcards.organism),
        output_dir=output_dir + "/organisms/{organism}",
        script=workflow.basedir + "/../scripts/proactiv_prepare_promoter_annotation.R"
    conda: "envs/r.yaml"
    log: "logs/proactiv_prepare_{organism}.log"
    benchmark: "benchmarks/proactiv_prepare_{organism}.bmk.txt"
    shell:
        """
        Rscript {params.script} {input.gtf} {params.output_dir} {params.species} {output.promoter_rds} {output.promoter_bed} > {log} 2>&1
        """

rule dexseq_prepare_annotation:
    input:
        gtf=output_dir + "/organisms/{organism}/Annotation/proActiv_protein_coding.gtf"
    output:
        gff=output_dir + "/organisms/{organism}/Annotation/DEXSeq_flattened_exons.gff",
    #conda: "envs/basic.yaml" changed to r because r.yaml has dexseq
    conda: "envs/r.yaml"
    log: "logs/dexseq_prepare_{organism}.log"
    benchmark: "benchmarks/dexseq_prepare_{organism}.bmk.txt"
    shell:
        """
        ENV_DIR=$(dirname $(dirname $(which python)))
        SCRIPT_SRC="$ENV_DIR/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py"
        SCRIPT_DST="$ENV_DIR/bin/dexseq_prepare_annotation.py"

        # Create symlink if not exists
        if [ ! -f "$SCRIPT_DST" ]; then
            ln -sf "$SCRIPT_SRC" "$SCRIPT_DST"
        fi

        python "$SCRIPT_DST" "{input.gtf}" "{output.gff}"

        if [ ! -s "{output.gff}" ]; then
            echo "Error: DEXSeq_flattened_exons.gff is empty or not created." >&2
            exit 1
        fi
        """