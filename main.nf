#!/usr/bin/env nextflow

/*
 * POINT-seq nextflow pipeline
 * By: Rui Sousa-Luis
*/

params.samples_info = "samplesInfo.txt"
params.outdir = "results"
params.species = "homo_sapiens"
params.version = 108
params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }

workflow {
    fastqFile = Channel
                    .fromPath(params.samples_info, checkIfExists: true)
                    .splitCsv()
                    .groupTuple()
    

    GET_ANOT(params.genomes[params.species]["gtf"],params.genomes[params.species]["fasta"])
    STAR_INDEX(GET_ANOT.out)
    
    // Conctatenate files from same sample
    READS_CONCAT(fastqFile)

    // Quality control of reads
    FASTQC(READS_CONCAT.out)

    // Trimming the adaptors
    TRIM(READS_CONCAT.out)
    FLASH(TRIM.out)

    // Align the reads & keep uniquely mapped reads in proper PE
    STAR_ALIGN(STAR_INDEX.out, TRIM.out)
    STRAND_SPLIT_BAM(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index)

    SAMTOOLS_STAT(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index)
    PRESEQ_LCEXTRAP(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index)
    FEATURE_COUNTS(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index,GET_ANOT.out.gtf)
    FEATURE_BIOTYPE(STAR_ALIGN.out.bams.collect(), GET_ANOT.out.gtf)
    RSEQC_JUNCT_ANOT(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out.bed12)
    RSEQC_JUNCT_SATUR(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out.bed12)
    RSEQC_INFER_EXP(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out.bed12)
    RSEQC_READ_DISTR(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out.bed12)
    RSEQC_TIN(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out.bed12)
    RSEQC_INNER_DIST(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out.bed12)
    RSEQC_GENE_BODY_COV(STAR_ALIGN.out.sampleID, STAR_ALIGN.out.bams, STAR_ALIGN.out.bam_index, GET_ANOT.out.bed12)

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
process GET_ANOT {
    tag "GET_ANOT"
    publishDir "$params.outdir/Anot", mode: 'copy'

    input:
    val url_gtf
    val url_fasta

    output:
    path "gtfFile.gtf", emit: gtf
    path "fastafile.fa", emit: fasta
    path "Protein_coding_genes_REF.bed12", emit: bed12

    script:
    """
    wget $url_gtf -O gtfFile.gtf.gz 
    wget $url_fasta -O fastafile.fa.gz 
    gunzip gtfFile.gtf.gz
    gunzip fastafile.fa.gz

    gtfToGenePred gtfFile.gtf genePred
    genePredToBed genePred bed12
    sort -k1,1 -k2,2n bed12 > bed12_sorted
    grep Ensembl_canonical gtfFile.gtf | grep protein_coding | grep -o ENST[0-9]* | sort | uniq | grep -f - bed12_sorted | awk -v OFS="\t" '{\$1="chr"\$1; print \$0}' > Protein_coding_genes_REF.bed12
    """
}

process STAR_INDEX {
    tag "STAR Index"
    publishDir "$params.outdir/Anot", mode: 'copy'

    input:
    path gtf 
    path fasta
    path bed12

    output:
    path "STAR_INDEX"

    script:
    """
    mkdir STAR_INDEX
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir STAR_INDEX --genomeFastaFiles $fasta \
        --sjdbGTFfile  $gtf --sjdbOverhang  149 --genomeSAindexNbases 14
    """
}

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

// Pre-processing of reads
process FASTQC {
    publishDir "$params.outdir/1-fastqc", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/2-trim", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/QControl/FLASH", mode: 'copy'
    tag "${ID}"

    input:
    val ID
    tuple path(trimFastQ_R1), path(trimReport_R1), path(trimFastQC_html_R1), path(trimFastQC_zip_R1),
            path(trimFastQ_R2), path(trimReport_R2), path(trimFastQC_html_R2), path(trimFastQC_zip_R2)

    output:
    path "*"
    
    script:
    """
    flash $trimFastQ_R1 $trimFastQ_R2 --output-prefix ${ID}_flash -t ${task.cpus} 2>&1 | tee ${ID}_flash_logfilename.log
    """

}

// Alignment of Reads 
process STAR_ALIGN {
    publishDir "$params.outdir/3-star", mode: 'copy'
    tag "${ID}"

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
    --readFilesIn $trimFastQ_R1 $trimFastQ_R2  --runThreadN ${task.cpus} --outFileNamePrefix $ID
    rename "s/Aligned.sortedByCoord.out//g" *
    samtools index ${ID}.bam
    """
}

process STRAND_SPLIT_BAM {
    publishDir "$params.outdir/3-star/strandBAMs", mode: 'copy'
    tag "${ID}"

    input:
    val ID
    path bams
    path bam_index

    output:
    path "F_${ID}.bam"
    path "R_${ID}.bam"

    script:
    """
    samtools view -f 99 -@ ${task.cpus} $bams -o "99_${ID}.bam"
    samtools view -f 147 -@ ${task.cpus} $bams -o "147_${ID}.bam"
    samtools view -f 83 -@ ${task.cpus} $bams -o "83_${ID}.bam"
    samtools view -f 163 -@ ${task.cpus} $bams -o "163_${ID}.bam"
    samtools merge -@ ${task.cpus} "F_${ID}.bam" "83_${ID}.bam" "163_${ID}.bam"
    samtools merge -@ ${task.cpus} "R_${ID}.bam" "99_${ID}.bam" "147_${ID}.bam"
    samtools index "F_${ID}.bam"
    samtools index "R_${ID}.bam"
    """
}

process CREATE_BW{
    publishDir "$params.outdir/3-star/strandBAMs", mode: 'copy'
    tag "${ID}"

    input:
    val genome
    val ID
    path "F_${ID}.bam"
    path "R_${ID}.bam"

    output:
    path "F_${ID}.bw"
    path "R_${ID}.bw"

    script:
    """
    '##########\ntrack ' + ${ID} + '\ncontainer multiWig\nshortLabel ' + ${ID} + '\nlongLabel ' + \
                        ${ID} + '\ntype bigWig 0 60000\nviewLimits -10000:10000\nvisibility full \
                        \nmaxHeightPixels 160:120:11\naggregate solidOverlay\nshowSubtrackColorOnUi on\nwindowingFunction maximum \
                        \nalwaysZero on\npriority 1.4\nconfigurable on\nautoScale on\n\ntrack ' + ${ID} + '_F\nbigDataUrl ' + \
                        "F_${ID}.bw" + \
                        '\nshortLabel ' + ${ID} + ' F\nlongLabel ' + ${ID} + ' Forward\nparent ' + ${ID} + '\ntype bigWig' + \
                        '\ncolor 0,0,255\n\ntrack ' + ${ID} + '_R\nbigDataUrl ' + \
                        'R_${ID}.bw' + \
                        '\nshortLabel ' + ${ID} + ' R' + \
                        '\nlongLabel ' + ${ID} + ' Reverse\nparent ' + ${ID} + '\ntype bigWig\ncolor 255,0,0\n\n' > $genome/trackDb.txt
    """
}

// Post-processing of reads
process PRESEQ_LCEXTRAP {
    publishDir "$params.outdir/QControl/PRESEQ", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/QControl/SAMTOOLS_STAT", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/QControl/FEATURE_COUNTS", mode: 'copy'
    tag "${ID}"

    input:
    val ID 
    path BAM
    path index_BAM
    path gtf

    output:
    path "*"

    script:
    """
    featureCounts -p -a $gtf -o ${ID}_FC.counts  -t transcript -g gene_id -s 2 -T ${task.cpus} --fracOverlap 0.9 -R CORE $BAM
    """
}

process FEATURE_BIOTYPE {
    publishDir "$params.outdir/QControl/FEATURE_BIOTYPE", mode: 'copy'
    tag "${ID}"

    input:
    path BAM
    path gtf

    output: 
    path "biotype_analysis.txt"

    script:
    """
    featureCounts -p -a $gtf -o FC.biotype  -t transcript -g gene_biotype -s 2 -T ${task.cpus} --fracOverlap 0.9 -R CORE $BAM
    cut -f1,7- FC.biotype | tail -n+2 | sed 's/Geneid/Status/g' > summary_biotype.txt
    awk '{ for (i=1; i<=NF; i++) a[i]= (a[i]? a[i] FS \$i: \$i) } END{ for (i in a) print a[i] }' summary_biotype.txt > biotype_analysis.txt
    """
}

process RSEQC_JUNCT_ANOT {
    publishDir "$params.outdir/QControl/RSEQC/JUNCT_ANOT", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/QControl/RSEQC/JUNCT_SATUR", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/QControl/RSEQC/INFER_EXP", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/QControl/RSEQC/READ_DISTR", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/QControl/RSEQC/TIN", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/QControl/RSEQC/INNER_DIST", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir/QControl/RSEQC/GENE_BODY_COV", mode: 'copy'
    tag "${ID}"

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
    publishDir "$params.outdir", mode: 'copy'
    tag "MULTIQC"

    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc $launchDir/$params.outdir -c $projectDir/multiQC.config
    """
}
