
params.reference_data_url = 'ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.reference_name = 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.reference_dir  = 'data'
params.index_dir = 'indexed_data_files'


/* derive names at compile time (so we can use them in output globs) */
// def refGz   = params.reference_name
// def refFa   = refGz.replaceAll(/\.gz$/, '') // .fa files
// def refBase = refFa.replaceAll(/\.fa$/, '') // basename files w no extension
// Channel
//     .value(params.reference_data_name)
//     .set { reference_file_ch }

process trim_reads {
    // 
    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)


    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz")

    script:
    """
    java -jar /apps/genomics/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 \
    $read1 $read2 \
    ${sample_id}_trimmed_R1.fastq.gz ${sample_id}_unpaired_R1.fastq.gz \
    ${sample_id}_trimmed_R2.fastq.gz ${sample_id}_unpaired_R2.fastq.gz \
    ILLUMINACLIP:data/TruSeq-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process fastqc {
    // Use FASTQC to quality control trimmed sequences
    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path("*_fastqc.*")

    script:
    """
    module load FastQC/0.11.8

    fastqc $read1 $read2
    """
}

process mapping {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    module load bwa/0.7.17

    bwa mem /home/c.c24082291/DISSERTATION/exome-workflow-proj/data/reference.fa $read1 $read2 > ${sample_id}.sam
    """
}

process make_bam{
    tag "$sample_id"

    input:
    tuple val(sample_id), path(samfile)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    module load samtools/1.17

    samtools view -bS $samfile > ${sample_id}.bam
    samtools sort -o ${sample_id}.sorted.bam ${sample_id}.bam
    """
}

process mark_duplicates {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bamfile)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markdup.bam"), path("${sample_id}.sorted.metrics")

    script:
    """
    module load java/1.8
    module load picard/2.27.5

    java -jar \$PICARD MarkDuplicates \
    I=$bamfile O=${sample_id}.sorted.markdup.bam \
    M=${sample_id}.sorted.metrics REMOVE_DUPLICATES=false
    """
}

process bamstats {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(markdup_file), path(_)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markdup.stats")

    script:
    """
    module load bamtools/170119

    bamtools stats -in $markdup_file > ${sample_id}.sorted.markdup.stats
    """
}

process index_bam {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(markdup_bamfile), path(_)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markdup.bam.bai")

    script:
    """
    module load samtools/1.17

    samtools index $markdup_bamfile
    """
}

process add_readgroups {
    tag "sample_id"

    input:
    tuple val(sample_id), path(markdup_bam), path(_)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markup.rg.bam")

    script:
    """
    module load picard/2.27.5
    module load java/1.8
    java -jar \$PICARD AddOrReplaceReadGroups I=$markdup_bam O=${sample_id}.sorted.markup.rg.bam SO=coordinate RGID=1 RGLB=libl RGPL=illumina RGPU=unit1 RGSM=${sample_id} CREATE_INDEX=True
    """
}

process index_reference {
    tag "$sample_id"

    input:
    path reference_genome

    output:
    tuple val(sample_id), path(indexed_reference)
    
    script:
    """   
    module load GATK/4.1.2.0    
    
    samtools faidx /home/c.c24082291/DISSERTATION/exome-workflow-proj/data/reference.fa
    """

}

process download_var_reference {
    tag "Download variant reference panel"
    
    output:
    path 'Mills_and_1000G_gold_standard.indels.hg38.vcf.gz', emit: var_reference_gz

    script:
    """
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    """

}

process index_variant_file {
    tag "Index variant reference panel"

    input:
    path var_reference_gz

    output:
    path '*.tbi', emit: indexed_var_reference

    """
    module load GATK/4.1.2.0
    gatk IndexFeatureFile -F $var_reference_gz
    """
}

workflow {
    // Read in pair-ended fastq files
    reads_ch = Channel.fromFilePairs("data/fastq/*_R{1,2}.slim.fastq.gz")
    // Map files into correct tuple format for process
    trimmed_reads = reads_ch.map { sample_id, reads -> 
        tuple(sample_id, reads[0], reads[1])
    }   
    // Pass data to trimming process
    trimmed_output = trim_reads(trimmed_reads)

    // Pass trimmed data to fastqc
    fastqc_output = fastqc(trimmed_output)

    // Map trimmed data to reference genome
    mapped_data = mapping(trimmed_output)

    // Make BAM file
    bam_file = make_bam(mapped_data)

    // Mark duplicates
    markdup_files = mark_duplicates(bam_file)
    
    // Create bamstats report
    bamstats_report = bamstats(markdup_files)

    // Index sample BAM file
    indexed_bam = index_bam(markdup_files)

    // Add readgroups to BAM file
    bam_with_readgroups = add_readgroups(markdup_files)

    // Download variant reference file
    var_reference_gz = download_var_reference()

    // Index variant reference file
    indexed_var_reference = index_variant_file(var_reference_gz)
}
