process trim_reads {
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




workflow {
    reads_ch = Channel.fromFilePairs("data/fastq/*_R{1,2}.slim.fastq.gz")

    trimmed_reads = reads_ch.map { sample_id, reads -> 
        tuple(sample_id, reads[0], reads[1])
    }

    trimmed_reads.view { it }

    trimmed_output = trim_reads(trimmed_reads)

    trimmed_output.view { it }
}
