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
    -threads 4 \
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
    path("*.sam")

    script:
    """
    module load bwa/0.7.17

    bwa mem /home/c.c24082291/DISSERTATION/exome-workflow-proj/data/reference.fa $read1 $read2 > ${sample_id}.sam
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

    trimmed_output.view { it }
    // Pass trimmed data to fastqc
    fastqc_output = fastqc(trimmed_output)

    // Map trimmed data to reference genome
    mapped_data = mapping(trimmed_output)
}
