
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


process download_reference {
    tag "Download reference genome"

    output:
    path 'Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa', emit: reference_genome

    script:
    """
    wget "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"

    gunzip Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
    """
}

process extract_chr22 {
    tag "Extract chr22 from reference genome"

    input:
    path reference_genome

    output:
    path "chr22.fa", emit: chr22_fasta

    script:
    """
    module load samtools/1.17    
    samtools faidx $reference_genome 22 > chr22.fa
    sed -i 's/^>22/>chr22/' chr22.fa
    """
}

process index_ref_bwa {
    tag "Index reference genome"

    input:
    path reference_genome

    output:  
    tuple path(reference_genome), path("${reference_genome}.*"), emit: ref_index_bundle

    script:
    """
    module load bwa/0.7.17
    bwa index -a bwtsw $reference_genome
    """
}


process trim_reads {
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
    tuple val(sample_id), path(read1), path(read2), path(reference_genome), path(reference_files)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    module load bwa/0.7.17
    bwa mem ${reference_genome} $read1 $read2 > ${sample_id}.sam
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
    tag "$sample_id"

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
    tag "Index Reference .fai"

    input:
    path reference_genome

    output:
    path '*.fai', emit: indexed_ref_fai
    
    script:
    """   
    module load GATK/4.1.2.0    
    module load samtools/1.17
    samtools faidx $reference_genome
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
    tuple path(var_reference_gz), path('*.tbi'), emit: indexed_var_reference

    """
    module load GATK/4.1.2.0
    gatk IndexFeatureFile -F $var_reference_gz
    """
}

process index_variant_file_2 {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf_file), path(recal_bam), path(reference_genome), path(indexed_ref_fai), path(seq_dict)

    output:
    tuple val(sample_id), path(vcf_file), path('*.tbi'), emit: indexed_vcf

    """
    module load GATK/4.1.2.0
    gatk IndexFeatureFile -F $vcf_file
    """
}

process create_seq_dir {
    tag "Create sequence dictionary"
    input:
    path reference_genome
    output:
    path 'chr22.dict', emit: reference_seq_dict

    script:
    """
    module load picard/2.27.5
    java -jar \$PICARD CreateSequenceDictionary R=$reference_genome O=chr22.dict
    """
}


process create_recalibration_model {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_with_readgroups), path(reference_genome), path(var_reference_gz), path(seq_dict), path(indexed_ref_fai), path(indexed_var_reference)

    output:
    tuple val(sample_id), path("${sample_id}.recal_data.table"), path(bam_with_readgroups), path(reference_genome), path(var_reference_gz), path(seq_dict), path(indexed_ref_fai), path(indexed_var_reference)

    script:
    """
    module load GATK/4.1.2.0

    gatk BaseRecalibrator \
        -I $bam_with_readgroups \
        -R $reference_genome \
        --known-sites $var_reference_gz \
        -O ${sample_id}.recal_data.table
    """
}

process recalibrate_bam {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(recal_file), path(bam_with_readgroups), path(reference_genome), path(indexed_var_reference), path(var_reference_gz), path(indexed_ref_fai), path(seq_dict)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markdup.rg.recal.bam"), path(reference_genome), path(indexed_ref_fai), path(seq_dict)

    script:
    """
    module load GATK/4.1.2.0

    gatk ApplyBQSR -R $reference_genome -I $bam_with_readgroups --bqsr-recal-file $recal_file -O ${sample_id}.sorted.markdup.rg.recal.bam
    """

}

process call_variants {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(recal_bam), path(reference_genome), path(indexed_ref_fai), path(seq_dict)

    output:
    tuple val(sample_id), path("*.vcf.gz"), path(recal_bam), path(reference_genome), path(indexed_ref_fai), path(seq_dict)

    script:
    """
    module load GATK/4.1.2.0
    module load java/1.8

    gatk --java-options "-Xmx4g" HaplotypeCaller -R $reference_genome -I $recal_bam -O ${sample_id}.vcf.gz    
    """
}

process remove_indels {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf_file), path(indexed_vcf)

    output:
    tuple val(sample_id), path("*snp.vcf.gz"), path(vcf_file), path(indexed_vcf)

    script:
    """
    module load GATK/4.1.2.0
    module load java/1.8
    gatk SelectVariants --variant $vcf_file --select-type SNP --output ${sample_id}.snp.vcf.gz
    """
}

process filter_variants {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(snp_vcf), path(vcf_file), path(indexed_vcf)

    output:
    tuple val(sample_id), path("*.snp.filtered.vcf.gz"), path(vcf_file), path(indexed_vcf)

    script:
    """
    module load GATK/4.1.2.0
    module load java/1.8
    gatk VariantFiltration  --variant $vcf_file --filter-expression "QD<2.0" --filter-name "QD2" --filter-expression "QUAL<30.0" --filter-name "QUAL30"  --filter-expression "SOR>3.0" --filter-name "SOR3" --filter-expression "FS>60.0" --filter-name "FS60" --filter-expression "MQ<40.0" --filter-name "MQ40" --filter-expression "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" --filter-expression "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" --output ${sample_id}.snp.filtered.vcf.gz
    """
}

process download_variant_array {
    tag "Download 1000 genome vcf"
    input:    

    output:
    tuple path("*.vcf"), path("*.vcf.idx"), emit: variant_array

    script:
    """
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx
    """
}

process refine_genotypes {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(filtered_vcf), path(filtered_indexed_vcf), path(thousand_genome_vcf), path(thousand_genome_vcf_idx)

    output:
    tuple val(sample_id), path("*.snp.filtered.ref.vcf.gz")

    script:
    """
    module load GATK/4.1.2.0
    module load java/1.8
    gatk --java-options "-Xmx4g" CalculateGenotypePosteriors -V $filtered_vcf -O ${sample_id}.snp.filtered.ref.vcf.gz -supporting $thousand_genome_vcf
    """
}

process index_vcf {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(filtered_vcf_file)

    output:
    tuple val(sample_id), path(filtered_vcf_file), path('*.tbi'), emit: filtered_indexed_vcf

    """
    module load GATK/4.1.2.0
    gatk IndexFeatureFile -F $filtered_vcf_file
    """
}

workflow {
    // Read in pair-ended fastq files
    reads_ch = Channel.fromFilePairs("data/fastq/*_R{1,2}.slim.fastq.gz")
    // Map files into correct tuple format for process
    trimmed_reads = reads_ch.map { sample_id, reads -> 
        tuple(sample_id, reads[0], reads[1])
    }

    // Download reference genome
    reference_genome = download_reference()

    // Download variant array
    thousand_genome_vcf = download_variant_array()

    chr22_fasta = extract_chr22(reference_genome)

    // Index reference genome
    indexed_reference = index_ref_bwa(chr22_fasta)

    // Pass data to trimming process
    trimmed_output = trim_reads(trimmed_reads)

    // Pass trimmed data to fastqc
    fastqc_output = fastqc(trimmed_output)

    // Combine trimmed reads with reference genome
    // Combine sample data with indexed reference
    mapping_input = trimmed_output.combine(indexed_reference)

    // Map trimmed data to reference genome
    mapped_data = mapping(mapping_input)

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

    // Index reference file (again but differently)
    reference_fai_file = index_reference(chr22_fasta)

    // Download variant reference file
    var_reference_gz = download_var_reference()

    // Index variant reference file
    indexed_var_reference = index_variant_file(var_reference_gz)

    // Generate sequence dictionary from reference genome file
    // seq_dict = create_seq_dir(reference_genome)
    seq_dict = create_seq_dir(chr22_fasta)

    recalibration_inputs = bam_with_readgroups
    .combine(chr22_fasta)        // add FASTA
    .combine(indexed_var_reference)   // add known-sites VCF
    .combine(reference_fai_file) // add .fai
    .combine(seq_dict)
        .map { sample_id, bam, fasta, vcf, fai, seq_dict, indexed_ref_fai->
        tuple(sample_id, bam, fasta, vcf, fai, seq_dict, indexed_ref_fai)
        }
    // dont touch it works
    recal_file = create_recalibration_model(recalibration_inputs)
    // Recalibrate bam
    recal_bam = recalibrate_bam(recal_file)    

    // Call variants
    vcf_file = call_variants(recal_bam)
    // Index (again)
    indexed_vcf = index_variant_file_2(vcf_file)

    // Filter variants
    // Remove indels
    snp_vcf = remove_indels(indexed_vcf)

    filtered_variants = filter_variants(snp_vcf) // *.snp.filtered.vcf.gz

    indexed_filtered_vcf = index_vcf(filtered_variants) // .snp.filtered.vcf.gz and index file

    thousand_genome_vcf.view { it }

    // Refine genotypes
    refining_inputs = indexed_filtered_vcf.combine(thousand_genome_vcf)    

    
    refined_vcf = refine_genotypes(refining_inputs)

}
