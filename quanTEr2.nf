#!/usr/bin/env nextflow

/*
================================================================================
                           quanTEfication/main
================================================================================
  Nextflow pipeline for RNA-seq analysis and quantefication of the transposons
--------------------------------------------------------------------------------
 @Repository
 https://github.com/savytskanatalia/quanTEfication
 @Author
 Savytska Natalia, 2020. Genome Biology of Neurodegenerative Diseases, DZNE TU
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run quanTEr.nf -c quanTEr.config --reads "/absolute/path/to/files/*_R{1,2}.fastq.gz" ... -with-docker savytskanatalia/quantefication2 
    Mandatory arguments:
      --reads [str]                  Path to input data (must be surrounded with quotes). Default: "input/*_R{1,2}*.fastq.gz" 
      --bam_files [str]              Path to input data (must be surrounded with quotes). Default: "input/*.bam" 
      --outdir [str]                 Output directory path. Default: output/  
      --modus [str]                  Mode to run the protgram in. Valid options: 'map' (only mapping),'mapquant' (both mapping and quantification) and 'quant' (only quantification). Default: 'mapquant'     
    References                        If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta [str]                  Path to fasta reference
      --gtf [str]                    Path to GTF annotation
      --star_index [str]             Path to STAR-Index reference
      --index                        Path to TElocal TE index file.
      --staropt [str]                Additional STAR flags for processing. If any additional STAR parameters need to be included, they can be with this flag (e.g. '--alignIntronMin 20 ...').
    Options:
      --threads [int]                Number of threads to use. Default: 4
      --read_length [int]            Length of the reads. Default: 100
      --mism [int]                   outFilterMismatchNmax in STAR  alignment. Default: 999
      --misml [int]                  outFilterMismatchNoverLmax in STAR  alignment. Default: 0.1
      --strand [str]                 Strandedness for the RNAseq data. Look up option "stranded" in TElocal. Default: "yes". Other options: "reverse", "no"
      --sorted [str]                 Option for STAR and TElocal. If .bam files are sorted. Default: 1 (=unsorted, no additional parameter for TElocal run). Other valid option: 2 (sorted)
    """.stripIndent()
}


sortt=params.sorted
stran=Channel.value(params.strand)
fastq_files=Channel.fromFilePairs(params.reads)
bams=Channel.fromPath(params.bam_files)
out=params.outdir
rfasta=params.fasta
rgtf=Channel.value(params.gtf)
tindex=params.index
read_length=Channel.value(params.read_length)
threads=Channel.value(params.threads)
mismatch_n=Channel.value(params.mism)
mismatch_nl=Channel.value(params.misml)
outdir=params.outdir
modus=params.modus
staropt=params.staropt


/*
================================================================================
                                  PREPROCESSING
================================================================================
*/


// CHECKING IF WE WANT TO MAP ONLY, QUANTIFY ONLY OR MAP AND QUANTIFY



/*
 * Build STAR index
 */
if (modus == 'map'){
    if (!params.star_index && (!params.fasta && !params.gtf)) exit 1, "Roses are red, Kaiju bloods blue, You forgot your STAR index, or give .fa for build! Either specify STAR-INDEX or Fasta and GTF!"

    if (!params.star_index){
        process build_star_index {
            echo true
            container='savytskanatalia/quantefication2'
            tag "${fasta}-${gtf}"
            label 'process_medium'
            publishDir params.outdir, mode: 'copy'
            input:
            path read from rfasta
            path gtf from rgtf
            val x from threads
            val y from read_length
            output:
            file("star-index") into star_index
            when: !(params.star_index)
            script:
            
            """
            source activate telocal
            echo TODAY WE ARE CANCELLING THE APOCALYPSE
            echo START BUILDING STAR INDEX....
            echo FASTA: $read
            echo GTF: $gtf
            mkdir star-index
            STAR \\
                --runMode genomeGenerate \\
                --runThreadN $x \\
                --sjdbGTFfile $gtf \\
                --sjdbOverhang ${y - 1} \\
                --genomeDir star-index/ \\
                --genomeFastaFiles $read \\
                --limitGenomeGenerateRAM 168632718037
                """
        }}

    ch_star_index = params.star_index ? Channel.value(file(params.star_index)).ifEmpty{exit 1, "You thought there was index, but there was none. Sorry, hurts to be wrong. STAR index not found: ${params.star_index}" } : star_index

    ch_star_index = ch_star_index.dump(tag:'ch_star_index')

    process star_mapping {
        container='savytskanatalia/quantefication2'
        echo true
        publishDir "${params.outdir}/mapped", mode: 'copy'
        input:
        set sampleId, file(reads) from fastq_files
        path x from Channel.value(params.star_index)
        val y from threads
        val a from mismatch_n
        val b from mismatch_nl
        val s from staropt
        output:
        file "*.bam"
        file "*"
        script:
        if( sortt == 1)
            """
            source activate telocal
            pwd
            echo HI THERE ${reads[0]} $sampleId
            STAR \\
                --outSAMstrandField intronMotif \\
                --outSAMunmapped Within \\
                --runThreadN $y \\
                --genomeDir $x \\
                --readFilesIn <(zcat ${reads[0]}) <(zcat ${reads[1]})  \\
                --outFilterMultimapNmax 100  \\
                --winAnchorMultimapNmax 100 \\
                --outSAMtype BAM Unsorted \\
                --outFilterMismatchNmax $a \\
                --outFilterMismatchNoverLmax $b \\
                --outFileNamePrefix  ${sampleId}_MM_100 $s


            """
        else if(sortt == 2)
            """
            source activate telocal
            pwd
            echo HI THERE ${reads[0]} $sampleId
            STAR \\
                --outSAMstrandField intronMotif \\
                --outSAMunmapped Within \\
                --runThreadN $y \\
                --genomeDir $x \\
                --readFilesIn <(zcat ${reads[0]}) <(zcat ${reads[1]}) \\
                --outFilterMultimapNmax 100  \\
                --winAnchorMultimapNmax 100 \\
                --outSAMtype BAM SortedByCoordinate \\
                --outFilterMismatchNmax $a \\
                --outFilterMismatchNoverLmax $b \\
                --outFileNamePrefix  ${sampleId}_MM_100 $s


            """
        else
            error "Invalid alignment mode: ${sortt}"
    }
}

else if(modus == 'mapquant'){
    if (!params.star_index && (!params.fasta && !params.gtf)) exit 1, "Roses are red, Kaiju bloods blue, You forgot your STAR index, or give .fa for build! Either specify STAR-INDEX or Fasta and GTF!"
    if (!params.star_index){
        process build_star_index_q {
            echo true
            container='savytskanatalia/quantefication2'
            tag "${fasta}-${gtf}"
            label 'process_medium'
            publishDir params.outdir, mode: 'copy'
            input:
            path read from rfasta
            path gtf from rgtf
            val x from threads
            val y from read_length
            output:
            file("star-index") into star_index
            when: !(params.star_index)
            script:
            
            """
            source activate telocal
            echo TODAY WE ARE CANCELLING THE APOCALYPSE
            echo START BUILDING STAR INDEX....
            echo FASTA: $read
            echo GTF: $gtf
            mkdir star-index
            STAR \\
                --runMode genomeGenerate \\
                --runThreadN $x \\
                --sjdbGTFfile $gtf \\
                --sjdbOverhang ${y - 1} \\
                --genomeDir star-index/ \\
                --genomeFastaFiles $read \\
                --limitGenomeGenerateRAM 168632718037
                """
        }}

    ch_star_index = params.star_index ? Channel.value(file(params.star_index)).ifEmpty{exit 1, "You thought there was index, but there was none. Sorry, hurts to be wrong. STAR index not found: ${params.star_index}" } : star_index

    ch_star_index = ch_star_index.dump(tag:'ch_star_index')

    process star_mapping_q {
        container='savytskanatalia/quantefication2'
        echo true
        publishDir "${params.outdir}/mapped", mode: 'copy'
        input:
        set sampleId, file(reads) from fastq_files
        path x from Channel.value(params.star_index)
        val y from threads
        val a from mismatch_n
        val b from mismatch_nl
        val s from staropt
        output:
        file "*.bam" into mapped_ch
        file "*"
        script:
        if( sortt == 1)
            """
            source activate telocal
            pwd
            echo HI THERE ${reads[0]} $sampleId
            STAR \\
                --outSAMstrandField intronMotif \\
                --outSAMunmapped Within \\
                --runThreadN $y \\
                --genomeDir $x \\
                --readFilesIn <(zcat ${reads[0]}) <(zcat ${reads[1]}) \\
                --outFilterMultimapNmax 100  \\
                --winAnchorMultimapNmax 100 \\
                --outSAMtype BAM Unsorted \\
                --outFilterMismatchNmax $a \\
                --outFilterMismatchNoverLmax $b \\
                --outFileNamePrefix  ${sampleId}_MM_100 $s


            """
        else if(sortt == 2)
            """
            source activate telocal
            pwd
            echo HI THERE ${reads[0]} $sampleId
            STAR \\
                --outSAMstrandField intronMotif \\
                --outSAMunmapped Within \\
                --runThreadN $y \\
                --genomeDir $x \\
                --readFilesIn <(zcat ${reads[0]}) <(zcat ${reads[1]})  \\
                --outFilterMultimapNmax 100  \\
                --winAnchorMultimapNmax 100 \\
                --outSAMtype BAM SortedByCoordinate \\
                --outFilterMismatchNmax $a \\
                --outFilterMismatchNoverLmax $b \\
                --outFileNamePrefix  ${sampleId}_MM_100 $s


            """
        else
            error "Invalid alignment mode: ${sortt}"
    }


    process telocal_quantification {
        container='savytskanatalia/quantefication2'
        echo true
        publishDir "${params.outdir}/quantified", mode: 'move', pattern: "*cntTable"
        publishDir "${params.outdir}/quantified/logs", mode: 'move', pattern: "*.log"
        input:
        file reads from mapped_ch
        path rgtf
        path tindex
        val x from stran
        val y from sortt
        output:
        file "*cntTable"
        file "*.log"
        when: modus == 'mapquant' && modus != 'map'
        script:
        if( sortt == 1)
            """
            source activate telocal
            echo HI THERE $reads $rgtf $tindex $x
            TElocal -b $reads --GTF $rgtf --TE $tindex --project ${reads}_TElocal_UM --mode uniq --stranded $x 2>${reads}_TElocal.log
            """
        else if(sortt == 2)
            """
            source activate telocal
            echo HI THERE $reads $rgtf $tindex $x
            TElocal --sortByPos -b $reads --GTF $rgtf --TE $tindex --project ${reads}_TElocal_UM --mode uniq --stranded $x 2>${reads}_TElocal.log
            """
        else
            error "Invalid alignment mode: ${sortt}"
    }

}

// "--sortByPos" (for sorted .bam) 
else if (params.modus == 'quant') {
    if (!params.bam_files) exit 1, "To run quantification you need to provide bam_files from which to quantify!! Alternatively, provide the .fast files and run mapping with STAR first!"
    process telocal_quant_alone {
        container='savytskanatalia/quantefication2'
        echo true
        publishDir "${params.outdir}/quantified", mode: 'move', pattern: "*cntTable"
        publishDir "${params.outdir}/quantified/logs", mode: 'move', pattern: "*.log"
        input:
        file reads from bams
        path rgtf
        path tindex
        val x from stran
        val y from sortt
        output:
        file "*cntTable"
        file "*.log"
        when: params.modus == 'quant'
        script:
        if( sortt == 1)
            """
            source activate telocal
            echo HI THERE $reads $rgtf $tindex $x
            TElocal -b $reads --GTF $rgtf --TE $tindex --project ${reads}_TElocal_UM --mode uniq --stranded  $x 2>${reads}_TElocal.log
            """
        else if(sortt == 2)
            """
            source activate telocal
            echo HI THERE $reads $rgtf $tindex $x
            TElocal --sortByPos -b $reads --GTF $rgtf --TE $tindex --project ${reads}_TElocal_UM --mode uniq --stranded  $x 2>${reads}_TElocal.log
            """
        else
            error "Invalid alignment mode: ${sortt}"
    }

}

else {
    if (!params.modus) exit 1, "What would you like to do? Map? Quantify? Map and Quantify? Choose your --modus wisely."
    else error "Invalid input for modus operandi! Please, choose an appropriate parameter for --modus: 'map', 'mapquant', or 'quant'"
}
