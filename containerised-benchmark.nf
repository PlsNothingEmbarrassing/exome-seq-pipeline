
//params.reference_path = '/home/c.c24082291/DISSERTATION/exome-workflow-proj/work/d7/c4cceda8eae1c209f474b983b962f7/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'
params.reference_index_fai = 'Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai'
params.fasta = 'Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'
params.data_dir  = '/scratch/c.c24082291/data'
params.variant_reference_path = 'Mills.GRCh38_ensembl.vcf'
params.variant_reference_index_path = 'Mills.GRCh38_ensembl.vcf.idx'
params.data_dict = 'Homo_sapiens.GRCh38.dna_sm.primary_assembly.dict'
params.out_dir = '/scratch/c.c24082291/exome-output'
params.fasta_dir = ''
params.truseq = '/scratch/c.c24082291/data/TruSeq-PE.fa'
params.indir   = "/scratch/c.c24082291/fastq"
params.samples = ['ERR166317','ERR166320','ERR166323', 'ERR166327', 'ERR166329']

process trim_reads {
    maxRetries params.retries
    maxErrors -1
    cpus params.trimReadsCpus
    time params.trimReadsJobLength
    memory params.trimReadsMemory

    tag "$sample_id"   

    input:
    tuple val(sample_id), path(read1), path(read2)


    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz")

    script:
    """
    trimmomatic PE -phred33 \
    $read1 $read2 \
    ${sample_id}_trimmed_R1.fastq.gz ${sample_id}_unpaired_R1.fastq.gz \
    ${sample_id}_trimmed_R2.fastq.gz ${sample_id}_unpaired_R2.fastq.gz \
    ILLUMINACLIP:${params.truseq}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

process fastqc {
    // Use FASTQC to quality control trimmed sequences
    maxRetries params.retries
    maxErrors -1
    cpus params.fastqcCpus
    time params.fastqcJobLength
    memory params.fastqcMemory


    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    path("*_fastqc.*")

    script:
    """
    fastqc $read1 $read2
    """
}

process mapping {
    maxRetries params.retries
    maxErrors -1
    cpus params.mappingCpus
    time params.mappingJobLength
    memory params.mappingMemory

    tag "$sample_id"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.sam")

    script:
    """
    bwa mem ${params.data_dir}/${params.fasta} $read1 $read2 > ${sample_id}.sam
    """
}


process make_bam{
    maxRetries params.retries
    maxErrors -1
    cpus params.makeBamCpus
    time params.makeBamJobLength
    memory params.makeBamMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(samfile)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam")

    script:
    """
    samtools view -bS $samfile > ${sample_id}.bam
    samtools sort -o ${sample_id}.sorted.bam ${sample_id}.bam
    """
}

process mark_duplicates {
    maxRetries params.retries
    maxErrors -1
    cpus params.markDuplicatesCpus
    time params.markDuplicatesJobLength
    memory params.markDuplicatesMemory

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bamfile)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markdup.bam"), path("${sample_id}.sorted.metrics")

    script:
    """
    picard MarkDuplicates \
    I=$bamfile O=${sample_id}.sorted.markdup.bam \
    M=${sample_id}.sorted.metrics REMOVE_DUPLICATES=false
    """
}

process bamstats {
    maxRetries params.retries
    maxErrors -1
    cpus params.bamstatsCpus
    time params.bamstatsJobLength
    memory params.bamstatsMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(markdup_file), path(_)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markdup.stats")

    script:
    """
    bamtools stats -in $markdup_file > ${sample_id}.sorted.markdup.stats
    """
}

process index_bam {
    maxRetries params.retries
    maxErrors -1
    cpus params.indexBamCpus
    time params.indexBamJobLength
    memory params.indexBamMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(markdup_bamfile), path(_)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markdup.bam.bai")

    script:
    """
    samtools index $markdup_bamfile
    """
}

process add_readgroups {
    maxRetries params.retries
    maxErrors -1
    cpus params.addReadgroupsCpus
    time params.addReadgroupsJobLength
    memory params.addReadgroupsMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(markdup_bam), path(_)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markup.rg.bam")

    script:
    """   
    picard AddOrReplaceReadGroups I=$markdup_bam O=${sample_id}.sorted.markup.rg.bam SO=coordinate RGID=1 RGLB=libl RGPL=illumina RGPU=unit1 RGSM=${sample_id} CREATE_INDEX=True
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
    gatk IndexFeatureFile -F $var_reference_gz
    """
}



// process create_seq_dir {
//     tag "Create sequence dictionary"
//     output:
//     path 'chr22.dict', emit: reference_seq_dict

//     publishDir params.data_dir, mode: 'copy', overwrite: true

//     script:
//     """
//     picard CreateSequenceDictionary R=${params.data_dir}/${params.fasta} O=chr22.dict
//     """
// }


process create_recalibration_model {
    maxRetries params.retries
    maxErrors -1
    cpus params.createRecalibrationModelCpus
    time params.createRecalibrationModelJobLength
    memory params.createRecalibrationModelMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam_with_readgroups)

    output:
    tuple val(sample_id), path("${sample_id}.recal_data.table"), path(bam_with_readgroups)

    script:
    """
    gatk BaseRecalibrator \
        -I $bam_with_readgroups \
        -R ${params.data_dir}/${params.fasta} \
        --known-sites ${params.data_dir}/${params.variant_reference_path} \
        -O ${sample_id}.recal_data.table
    """
}

process recalibrate_bam {
    maxRetries params.retries
    maxErrors -1
    cpus params.recalibrateBamCpus
    time params.recalibrateBamJobLength
    memory params.recalibrateBamMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(recal_file), path(bam_with_readgroups)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.markdup.rg.recal.bam")

    script:
    """
    gatk ApplyBQSR -R ${params.data_dir}/${params.fasta} -I $bam_with_readgroups --bqsr-recal-file $recal_file -O ${sample_id}.sorted.markdup.rg.recal.bam
    """

}

process index_bam_2 {
    maxRetries params.retries
    maxErrors -1
    cpus params.indexBamCpus
    time params.indexBamJobLength
    memory params.indexBamMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(recal_bam_file)

    output:
    tuple val(sample_id), path(recal_bam_file), path("${sample_id}.sorted.markdup.rg.recal.bam.bai")

    script:
    """
    samtools index $recal_bam_file
    """
}

process call_variants {
    maxRetries params.retries
    maxErrors -1
    cpus params.callVariantsCpus
    time params.callVariantsJobLength
    memory params.callVariantsMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(recal_bam), path(indexed_recal_bam)

    output:
    tuple val(sample_id), path("*.vcf.gz"), path(recal_bam)

    publishDir params.data_dir, mode: 'copy', overwrite: true

    script:
    """
    gatk --java-options "-Xmx4g" HaplotypeCaller -R ${params.data_dir}/${params.fasta} -I $recal_bam -O ${sample_id}.vcf.gz    
    """
}

process index_variant_file_2 {
    maxRetries params.retries
    maxErrors -1
    cpus params.indexVariantFile2Cpus
    time params.indexVariantFile2JobLength
    memory params.indexVariantFile2Memory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf_file), path(recal_bam)

    output:
    tuple val(sample_id), path(vcf_file), path('*.tbi'), emit: indexed_vcf

    publishDir params.data_dir, mode: 'copy', overwrite: true

    """
    gatk IndexFeatureFile -F $vcf_file
    """
}

process remove_indels {
    maxRetries params.retries
    maxErrors -1
    cpus params.removeIndelsCpus
    time params.removeIndelsJobLength
    memory params.removeIndelsMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(vcf_file), path(indexed_vcf)

    output:
    tuple val(sample_id), path("*snp.vcf.gz"), path(vcf_file), path(indexed_vcf)

    publishDir params.data_dir, mode: 'copy', overwrite: true

    script:
    """
    gatk SelectVariants --variant $vcf_file --select-type SNP --output ${sample_id}.snp.vcf.gz
    """
}

process index_snp_vcf {
    maxRetries params.retries
    maxErrors -1
    cpus params.indexSnpVcfCpus
    time params.indexSnpVcfJobLength
    memory params.indexSnpVcfMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(snp_vcf), path(vcf_file), path(indexed_vcf)

    output:
    tuple val(sample_id), path(snp_vcf), path('*.tbi'), emit: indexed_snp_vcf

    script:
    """
    gatk IndexFeatureFile -F $snp_vcf
    """
}
process filter_variants {
    maxRetries params.retries
    maxErrors -1
    cpus params.filterVariantsCpus
    time params.filterVariantsJobLength
    memory params.filterVariantsMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(snp_vcf), path(indexed_snp_vcf)

    output:
    tuple val(sample_id), path("*.snp.filtered.vcf.gz"), path(snp_vcf), path(indexed_snp_vcf)

    script:
    """
    gatk VariantFiltration  --variant $snp_vcf --filter-expression "QD<2.0" --filter-name "QD2" --filter-expression "QUAL<30.0" --filter-name "QUAL30"  --filter-expression "SOR>3.0" --filter-name "SOR3" --filter-expression "FS>60.0" --filter-name "FS60" --filter-expression "MQ<40.0" --filter-name "MQ40" --filter-expression "MQRankSum<-12.5" --filter-name "MQRankSum-12.5" --filter-expression "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" --output ${sample_id}.snp.filtered.vcf.gz
    """
}

process index_vcf {
    maxRetries params.retries
    maxErrors -1
    cpus params.indexVcfCpus
    time params.indexVcfJobLength
    memory params.indexVcfMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(filtered_vcf_file), path(vcf_file), path(indexed_vcf)

    output:
    tuple val(sample_id), path(filtered_vcf_file), path('*.tbi'), emit: filtered_indexed_vcf

    """
    gatk IndexFeatureFile -F $filtered_vcf_file
    """
}

process refine_genotypes {
    maxRetries params.retries
    maxErrors -1
    cpus params.refineGenotypesCpus
    time params.refineGenotypesJobLength
    memory params.refineGenotypesMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(snp_vcf), path(indexed_snp_vcf)

    output:
    tuple val(sample_id), path("*.snp.filtered.ref.vcf.gz"), emit: refined_vcf

    script:
    """
    gatk --java-options "-Xmx4g" CalculateGenotypePosteriors -V $snp_vcf -O ${sample_id}.snp.filtered.ref.vcf.gz -supporting ${params.data_dir}/${params.variant_reference_path}
    """
}

process annotate_variants {
    maxRetries params.retries
    maxErrors -1
    cpus params.annotateVariantsCpus
    time params.annotateVariantsJobLength
    memory params.annotateVariantsMemory
    tag "$sample_id"

    input:
    tuple val(sample_id), path(refined_vcf_file)

    output:
    path "*.hg38_multianno.txt"
    path "*.hg38_multianno.vcf"

    publishDir params.out_dir, mode: 'copy', overwrite: true

    script:
    """
    # Annotate variants using ANNOVAR
    /home/c.c24082291/DISSERTATION/exome-workflow-proj/annovar/table_annovar.pl \
        ${refined_vcf_file} \
        /home/c.c24082291/DISSERTATION/exome-workflow-proj/annovar/humandb \
        -buildver hg38 \
        -out ${refined_vcf_file.baseName} \
        -remove \
        -protocol refGene,cytoBand,avsnp150,clinvar_20210501 \
        -operation g,r,f,f \
        -nastring . \
        -vcfinput
        

    """
}

workflow {
    // Read in pair-ended fastq files
    reads_ch = Channel
    .fromList(params.samples)
    .map { id -> tuple(id,
                     file("${params.indir}/${id}_1.fastq.gz"),
                     file("${params.indir}/${id}_2.fastq.gz")) }
    // Map files into correct tuple format for process
    // trimmed_reads = reads_ch.map { sample_id, reads -> 
    //     tuple(sample_id, reads[0], reads[1])
    // }

    // Pass data to trimming process - done
    trimmed_output = trim_reads(reads_ch)

    // Pass trimmed data to fastqc - done
    fastqc_output = fastqc(trimmed_output)

    // Map trimmed data to reference genome - done
    mapped_data = mapping(trimmed_output)

    // Make BAM file - done
    bam_file = make_bam(mapped_data)

    // Mark duplicates - done
    markdup_files = mark_duplicates(bam_file)
    
    // Create bamstats report - done
    bamstats_report = bamstats(markdup_files)

    // Index sample BAM file - done
    indexed_bam = index_bam(markdup_files)

    // Add readgroups to BAM file - done
    bam_with_readgroups = add_readgroups(markdup_files)

    // Index reference file (again but differently) - done
    //reference_fai_file = index_reference(chr22_fasta)

    // Download variant reference file
    //var_reference_gz = download_var_reference()

    // Index variant reference file
    //indexed_var_reference = index_variant_file(var_reference_gz)

    // Generate sequence dictionary from reference genome file
    // seq_dict = create_seq_dir(reference_genome)
    //seq_dict = create_seq_dir()

    recalibration_inputs = bam_with_readgroups
    // dont touch it works
    recal_file = create_recalibration_model(recalibration_inputs)
    // Recalibrate bam
    recal_bam = recalibrate_bam(recal_file)

    // Index recalibrated bam (again)
    indexed_recal_bam = index_bam_2(recal_bam)  

    // Call variants
    vcf_file = call_variants(indexed_recal_bam)
    // Index (again)
    indexed_vcf = index_variant_file_2(vcf_file)

    // Filter variants
    // Remove indels
    snp_vcf = remove_indels(indexed_vcf)

    indexed_snp_vcf = index_snp_vcf(snp_vcf) // .snp.vcf.gz and index file

    filtered_variants = filter_variants(indexed_snp_vcf) // .snp.filtered.vcf.gz

    indexed_filtered_vcf = index_vcf(filtered_variants) // .snp.filtered.vcf.gz and index file
    
    refined_vcf = refine_genotypes(indexed_filtered_vcf)

    annotate_variants = annotate_variants(refined_vcf)

    

}
