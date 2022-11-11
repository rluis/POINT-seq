#!/usr/bin/env nextflow

/*
 * The following pipeline parameters specify the refence genomes
 * and read pairs and can be provided as command line options
 */

params.samples_info = "samplesInfo.txt"
params.outdir = "results"
params.cpus = 10
params.genome_fasta = "/home/rluis/Projects_RSL/Genome/hg38_HISAT2_INDEX/Ensembl/GRCh38.p10.genome.fa"
params.anot_GTF = "/home/rluis/Projects_RSL/Genome/GTF_hg38/Ensembl/Homo_sapiens.GRCh38.90_Wchr.gtf"
params.index_STAR = "/home/rluis/Projects_RSL/Genome/hg38_STAR_INDEX"
params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }

workflow {
    fastqFile = Channel
                    .fromPath(params.samples_info, checkIfExists: true)
                    .splitCsv()
                    .groupTuple()

    // Conctatenate files from same sample
    READS_CONCAT(fastqFile)
    GET_ANOT()
    // Quality control of reads
    FASTQC(READS_CONCAT.out)

    // Trimming the adaptors
    TRIM(READS_CONCAT.out)
    FLASH(TRIM.out)

    // Align the reads & keep uniquely mapped reads in proper PE
    STAR_ALIGN(params.index_STAR, TRIM.out)

    SAMTOOLS_STAT(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index)
    PRESEQ_LCEXTRAP(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index)
    FEATURE_COUNTS(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index)

    FEATURE_BIOTYPE(STAR_ALIGN.out.bams.collect())
    RSEQC_JUNCT_ANOT(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out)
    RSEQC_JUNCT_SATUR(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out)
    RSEQC_INFER_EXP(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out)
    RSEQC_READ_DISTR(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out)
    RSEQC_TIN(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out)
    RSEQC_INNER_DIST(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out)
    RSEQC_GENE_BODY_COV(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out)

    MULTIQC(FLASH.out.mix(SAMTOOLS_STAT.out)
                        .mix(PRESEQ_LCEXTRAP.out)
                        .mix(FEATURE_COUNTS.out)
                        .mix(FEATURE_BIOTYPE.out)
                        .mix(RSEQC_JUNCT_ANOT.out)
                        .mix(RSEQC_JUNCT_SATUR.out)
                        .mix(RSEQC_INFER_EXP.out)
                        .mix(RSEQC_READ_DISTR.out)
                        .mix(RSEQC_TIN.out)
                        .mix(RSEQC_INNER_DIST.out)
                        .mix(RSEQC_GENE_BODY_COV.out)
                        .collect())
}

// Beginning of pipeline
process READS_CONCAT {
    tag "$ID"
    publishDir "$params.outdir/0-fastQ", mode: 'copy'

    input:
    tuple val(ID), path(fastqFile_R1), path(fastqFile_R2)

    output:
    tuple val(ID), path("${ID}_R1.fastq.gz"), path("${ID}_R2.fastq.gz")

    script:
    """
    cat $fastqFile_R1 >> "${ID}_R1.fastq.gz"
    cat $fastqFile_R2 >> "${ID}_R2.fastq.gz"
    """
}

process STAR_INDEX {
    conda 'condaPOINTseq.yml'
    tag "STAR Index"

    output:
    path "hg38_STAR_INDEX"

    script:
    """
    STAR --runThreadN $params.cpus --runMode genomeGenerate --genomeDir hg38_STAR_INDEX --genomeFastaFiles $params.genome_fasta
        --sjdbGTFfile  $params.anot_GTF --sjdbOverhang  149 --genomeSAindexNbases 14
    """
}

process GET_ANOT {
    conda 'condaPOINTseq.yml'
    tag "GET_ANOT"

    output:
    path "hg38_protein_coding.bed12"

    script:
    """
    wget https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.chr.gtf.gz
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
    wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/genePredToBed

    chmod +x gtfToGenePred genePredToBed
    gunzip Homo_sapiens.GRCh38.108.chr.gtf.gz

    ./gtfToGenePred Homo_sapiens.GRCh38.108.chr.gtf genePred
    ./genePredToBed genePred bed12
    sort -k1,1 -k2,2n bed12 > bed12_sorted
    grep Ensembl_canonical Homo_sapiens.GRCh38.108.chr.gtf | grep protein_coding | grep -o ENST[0-9]* | sort | uniq | grep -f - bed12_sorted | awk -v OFS="\t" '{\$1="chr"\$1; print \$0}' > hg38_protein_coding.bed12
    """

}

// Pre-processing of reads
process FASTQC {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/1-fastqc", mode: 'copy'
    tag "FASTQC_$ID"

    input:
    tuple val(ID), path(fastqFile_R1), path(fastqFile_R2)

    output:
    tuple path("${ID}_R1_fastqc.html"),
          path("${ID}_R1_fastqc.zip"),
          path("${ID}_R2_fastqc.html"), 
          path("${ID}_R2_fastqc.zip")

    script:
    """
    fastqc $fastqFile_R1
    fastqc $fastqFile_R2
    """
}

process TRIM {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/2-trim", mode: 'copy'
    tag "TRIM_$ID"

    input:
    tuple val(ID),  path(fastqFile_R1), path(fastqFile_R2)

    output:
    val ID
    tuple path("${ID}_R1_val_1.fq.gz"),
          path("${fastqFile_R1}_trimming_report.txt"),
          path("${ID}_R1_val_1_fastqc.html"),
          path("${ID}_R1_val_1_fastqc.zip"),
          path("${ID}_R2_val_2.fq.gz"),
          path("${fastqFile_R2}_trimming_report.txt"),
          path("${ID}_R2_val_2_fastqc.html"),
          path("${ID}_R2_val_2_fastqc.zip")

    script:
    """
    trim_galore $fastqFile_R1 $fastqFile_R2 --fastqc -j 8 --paired --length 10 -q 20 --stringency 3
    """
}

process FLASH {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/FLASH", mode: 'copy'
    tag "FASTQC_$ID"

    input:
    val ID
    tuple path(trimFastQ_R1), path(trimReport_R1), path(trimFastQC_html_R1), path(trimFastQC_zip_R1),
            path(trimFastQ_R2), path(trimReport_R2), path(trimFastQC_html_R2), path(trimFastQC_zip_R2)

    output:
    path "*"
    
    script:
    """
    flash $trimFastQ_R1 $trimFastQ_R2 --output-prefix ${ID}_flash -t $params.cpus 2>&1 | tee ${ID}_flash_logfilename.log
    """

}

// Alignment of Reads 
process STAR_ALIGN {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/3-star", mode: 'copy'
    tag "STAR_ALIGN_${ID}"

    input:
    path genomeDir
    val ID
    tuple path(trimFastQ_R1), path(trimReport_R1), path(trimFastQC_html_R1), path(trimFastQC_zip_R1),
            path(trimFastQ_R2), path(trimReport_R2), path(trimFastQC_html_R2), path(trimFastQC_zip_R2)

    output:
    val ID, emit: sampleID
    path "${ID}.bam", emit: bams
    path "${ID}.bam.bai", emit: bam_index
    path "${ID}Log.final.out", emit: final_log

    script:
    """
    STAR --genomeDir $genomeDir --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM 100000000000 --readFilesCommand zcat  \
            --readFilesIn $trimFastQ_R1 $trimFastQ_R2  --runThreadN $params.cpus --outFileNamePrefix $ID
    rename "s/Aligned.sortedByCoord.out//g" *
    samtools index ${ID}.bam
    """
}

process CREATE_BW {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/3-star/strandBAMs", mode: 'copy'
    tag "STAR_ALIGN_${ID}"

    input:
    val ID
    path bams
    path bam_index

    output:
    path "F_${ID}.bam"
    path "F_${ID}.bam"


    script:
    """
    samtools view -f 99 -@ $params.cpus > 
    samtools view -f 147
    samtools view -f 183
    samtools view -f 83
    """
}


// Post-processing of reads
process PRESEQ_LCEXTRAP {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/PRESEQ", mode: 'copy'
    tag "PRESEQ_LCEXTRAP_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM

    output:
    path "${ID}_lc_extrap.txt"

    script:
    """
    preseq lc_extrap -B -o ${ID}_lc_extrap.txt $BAM
    """
}

process SAMTOOLS_STAT {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/SAMTOOLS_STAT", mode: 'copy'
    tag "SAMTOOLS_STAT_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM

    output:
    path "${ID}_SAMTOOLS_FLAGSTATS.txt"
    path "${ID}_SAMTOOLS_idxstats.txt"

    script:
    """
    samtools flagstats $BAM > ${ID}_SAMTOOLS_FLAGSTATS.txt
    samtools idxstats $BAM > ${ID}_SAMTOOLS_idxstats.txt
    """

}

process FEATURE_COUNTS {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/FEATURE_COUNTS", mode: 'copy'
    tag "FEATURE_COUNTS_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM

    output:
    path "*"

    script:
    """
    featureCounts -p -a $params.anot_GTF -o ${ID}_FC.counts  -t transcript -g gene_id -s 2 -T $params.cpus --fracOverlap 0.9 -R CORE $BAM
    """
}

process FEATURE_BIOTYPE {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/FEATURE_BIOTYPE", mode: 'copy'
    tag "FEATURE_BIOTYPE"

    input:
    path BAM

    output: 
    path "biotype_analysis.txt"

    script:
    """
    featureCounts -p -a $params.anot_GTF -o FC.biotype  -t transcript -g gene_biotype -s 2 -T $params.cpus --fracOverlap 0.9 -R CORE $BAM
    cut -f1,7- FC.biotype | tail -n+2 | sed 's/Geneid/Status/g' > summary_biotype.txt
    awk '{ for (i=1; i<=NF; i++) a[i]= (a[i]? a[i] FS \$i: \$i) } END{ for (i in a) print a[i] }' summary_biotype.txt > biotype_analysis.txt
    """
}

process RSEQC_JUNCT_ANOT {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/RSEQC/JUNCT_ANOT", mode: 'copy'
    tag "RSEQC_JUNCT_ANOT_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM
    path ANOT 

    output:
    path "*"

    script:
    """
    junction_annotation.py -i $BAM -o ${ID}_junction_annotation -r $ANOT 2> ${ID}_junction_annotation.log
    """
}

process RSEQC_JUNCT_SATUR {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/RSEQC/JUNCT_SATUR", mode: 'copy'
    tag "RSEQC_JUNCT_SATUR_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM
    path ANOT 

    output:
    path "*"

    script:
    """
    junction_saturation.py -i $BAM -o ${ID}_junction_saturation -r $ANOT 2> ${ID}_junction_saturation.log
    """
}

process RSEQC_INFER_EXP {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/RSEQC/INFER_EXP", mode: 'copy'
    tag "RSEQC_INFER_EXP_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM
    path ANOT 

    output:
    path "*"

    script:
    """
    infer_experiment.py -r $ANOT -i $BAM > ${ID}_infer_experiment.log
    """
}

process RSEQC_READ_DISTR {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/RSEQC/READ_DISTR", mode: 'copy'
    tag "RSEQC_READ_DISTR_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM
    path ANOT 

    output:
    path "*"

    script:
    """
    read_distribution.py -r $ANOT -i $BAM > ${ID}_read_distribution.log
    """
}

process RSEQC_TIN {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/RSEQC/TIN", mode: 'copy'
    tag "RSEQC_TIN_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM
    path ANOT 

    output:
    path "*"

    script:
    """
    tin.py -r $ANOT -i $BAM > ${ID}_tin.log
    """
}

process RSEQC_INNER_DIST {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/RSEQC/INNER_DIST", mode: 'copy'
    tag "RSEQC_INNER_DIST_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM
    path ANOT 

    output:
    path "*"

    script:
    """
    inner_distance.py -r $ANOT -i $BAM -o ${ID}_inner_distance.log
    """
}

process RSEQC_GENE_BODY_COV {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir/QControl/RSEQC/GENE_BODY_COV", mode: 'copy'
    tag "RSEQC_GENE_BODY_COV_${ID}"

    input:
    val ID 
    path BAM
    path index_BAM
    path ANOT 

    output:
    path "${ID}_geneBodyCov.geneBodyCoverage.txt"
    path "${ID}_geneBodyCov.geneBodyCoverage.curves.pdf"

    script:
    """
    geneBody_coverage.py -i $BAM -r $ANOT -o ${ID}_geneBodyCov 
    """
}

process MULTIQC {
    conda 'condaPOINTseq.yml'
    publishDir "$params.outdir", mode: 'copy'
    tag "MULTIQC"

    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc $projectDir/$params.outdir -c $projectDir/multiQC.config
    """
}