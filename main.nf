#!/usr/bin/env nextflow

log.info "===================================================================="
log.info "                 Select samples from multisample VCF                 "
log.info "===================================================================="

// set threadmem equal to total memory divided by number of threads
int threads = Runtime.getRuntime().availableProcessors()
threadmem = (((Runtime.getRuntime().maxMemory() * 4) / threads) as nextflow.util.MemoryUnit)

// fasta
Channel.fromPath(params.fasta)
        .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
        .into { ch_fasta ; ch_fasta_gather}

// fai
Channel.fromPath(params.fai)
        .ifEmpty { exit 1, "fasta index file not found: ${params.fai}" }
        .into { ch_fai ; ch_fai_gather }

// dict
Channel.fromPath(params.dict)
        .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
        .into { ch_dict ; ch_dict_gather }

// TODO: Multiply channel arrays by N = number of vcf files (most likely 22, from autosomal chr count)
// Create ch with [fasta, fai, dict]
// ch_reference        = ch_fasta.combine(ch_fai)
// ch_reference_bundle = ch_reference.combine(ch_dict)
// ch_reference_bundle_to_use.into { ch_reference_bundle ; ch_reference_bundle_to_inspect}
// ch_reference_bundle_to_inspect.view()

Channel.fromPath(params.multiVCF_table)
       .ifEmpty { exit 1, "File with vcf and respective index not found or not passed to --multiVCF_table" }
       .splitCsv(sep: ',',  skip: 1 )
       .map{ vcf, vcf_index -> [file(vcf), file(vcf_index)] }
       .set { ch_multiVCF_table }

Channel.fromPath("${params.list_folder}/*.list")
       .flatten()
       .into { ch_subset_lists; ch_subset_lists_view}

// Create ch with [pop.list, vcf, vcf_index]
ch_multiVCF = ch_subset_lists_view.combine(ch_multiVCF_table)

process SubsetMultiVCF {

    tag {"${sample_list.simpleName}-${vcf.baseName}"}
    container 'broadinstitute/gatk:4.1.3.0'
    publishDir "${params.outdir}/populations/${sample_list.simpleName}/individual_chr_vcfs/", mode: 'copy'

    input:
    set file(sample_list), file(vcf), file(vcf_index) from ch_multiVCF
    each file(fasta) from ch_fasta
    each file(fai) from ch_fai
    each file(dict) from ch_dict

    output:
    set val("${sample_list.simpleName}"), file("${vcf.baseName}.${sample_list.simpleName}.vcf.gz") into (ch_pops_vcfs, ch_pops_vcfs_to_inspect)

    script:
    """
    gatk SelectVariants \
    -R ${fasta} \
    -V $vcf \
    -O with_rsID.vcf \
    --sample-name ${sample_list}  \
    --restrict-alleles-to BIALLELIC \
    --select-type-to-include SNP \
    --verbosity ERROR > stderr.txt

    # Recode ID to: chr:pos:ref:alt
    bcftools annotate -x ID iberian.vcf | bcftools annotate --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' > ${vcf.baseName}.${sample_list.simpleName}.vcf

    bgzip -c ${vcf.baseName}.${sample_list.simpleName}.vcf > ${vcf.baseName}.${sample_list.simpleName}.vcf.gz
   """
}

ch_pops_vcfs_to_inspect
                        .groupTuple(by: 0)
                        //.map {pop_name, subsetted_vcf -> subsetted_vcf}
                        .view()

ch_grouped_pop_vcfs = ch_pops_vcfs.groupTuple(by: 0)


process GatherVcfs {

    tag "${pop_name}"
    publishDir "${params.outdir}/populations/${pop_name}/subsampled_multisample_vcf/", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    set val(pop_name), file (vcf_bundle) from ch_grouped_pop_vcfs
    each file(fasta) from ch_fasta_gather
    each file(fai) from ch_fai_gather
    each file(dict) from ch_dict_gather

    output:
    file("*") into ch_complete_chr_vcf
    set val("${pop_name}"), file("${pop_name}.vcf.gz") into (ch_plink_count_freqs, ch_plink_count_freqs_to_inspect)


    script:
    """
    ls *.vcf.gz | while read vcf; do tabix -fp vcf \$vcf; done

    ## make list of input variant files
    for vcf in \$(ls *vcf.gz); do
    echo \$vcf >> ${pop_name}.vcf.list
    done

    gatk GatherVcfs \
    --INPUT  ${pop_name}.vcf.list \
    --OUTPUT ${pop_name}.vcf.gz
    """
    }

ch_plink_count_freqs_to_inspect.view()

process PlinkFilterAndFreqCount {

    tag "${pop_name}"
    publishDir "${params.outdir}/${pop_name}/plink_metrics/", mode: 'copy'
    container 'alliecreason/plink:1.90'

    input:
    set val(pop_name), file(all_chr_vcf) from ch_plink_count_freqs

    output:
    file("*") into ch_plink_results
    file("${pop_name}.frq.counts") into ch_plink_frq_counts


    script:
    """
    plink \
    --vcf $all_chr_vcf \
    --snps-only \
    --biallelic-only strict list \
    --geno 0.05 \
    --maf  0.05 \
    --freq counts \
    --out $pop_name > ${pop_name}_plink.stdout.log

    rm *.vcf.gz
    """
    }


// process FreqCounts2Ref {

//     tag "${pop_name}"
//     publishDir "${params.outdir}/FreqCounts2Ref/", mode: 'copy'
//     container ' R with data.table '

//     input:
//     set value("${pop_name}"), file("${pop_name}_chr_pos_id.csv"), ("${pop_name}.frq.counts") from ch_plink_frq_counts 

//     output:
//     file("*") into ch_freq_dataframes

//     script:
//     """
//     #!/usr/bin/env Rscript

//     library(data.table)
    
//     # chr pos rs
//     snp_metadata           <- data.table::fread("${pop_name}_chr_pos_id.csv")
//     colnames(snp_metadata) <- c("chr", "pos", "rs")

//     # ref panel:  chr,rs,gpos,pos,a1,a2 
//     # CHR,SNP,A1,A2,C1,C2,G0 
//     frq_counts           <- data.table::fread("${pop_name}.frq.counts")
//     colnames(frq_counts) <- c("chr","rs","a1","a2", "c1","c2","gpos" )

    
//     """
//     }
