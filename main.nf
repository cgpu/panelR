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

// ch_reference        = ch_fasta.combine(ch_fai)
// ch_reference_bundle = ch_reference.combine(ch_dict)
// ch_reference_bundle.view()

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

process chop_multiVCF {

    tag "${sample_list.simpleName}"
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}/subsampled_multisample_vcf/${sample_list.simpleName}", mode: 'copy'

    input:
    set file(sample_list), file(vcf), file(vcf_index) from ch_multiVCF
    each file(fasta) from ch_fasta
    each file(fai) from ch_fai
    each file(dict) from ch_dict

    output:
    file("*") into ch_nowhere
    set val("${sample_list.simpleName}"), file("${vcf.baseName}.${sample_list.simpleName}.vcf") into (ch_pops_vcfs, ch_pops_vcfs_to_inspect)

    script:
    """
    gatk SelectVariants \
    -R ${fasta} \
    -V $vcf \
    -O ${vcf.baseName}.${sample_list.simpleName}.vcf \
    --sample-name ${sample_list} \
   """
}

ch_pops_vcfs_to_inspect
                        .groupTuple(by: 0)
                        //.map {pop_name, subsetted_vcf -> subsetted_vcf}
                        .view()


// process GatherVcfs {

//     tag "${pop_name}"
//     publishDir "${params.outdir}/subsampled_multisample_vcf/${sample_list.simpleName}", mode: 'copy'
//     container 'broadinstitute/gatk:latest'

//     input:
//     set val(pop_name), file ('*vcf') from ch_pops_vcfs
//     each file(fasta) from ch_fasta_gather
//     each file(fai) from ch_fai_gather
//     each file(dict) from ch_dict_gather

//     output:
//     file("*") into ch_complete_chr_vcf


//     script:
//     """
//     ## make list of input variant files
//     for vcf in \$(ls *vcf); do
//     echo \$vcf >> input_variant_files.list
//     done
//     gatk GatherVcfs \
//     --INPUT= input_variant_files.list \
//     --OUTPUT= ${pop_name}.vcf
//     """
//     }