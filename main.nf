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

process SubsetPopVCF {

    tag {"${sample_list.simpleName}-${vcf.baseName}"}
    publishDir "${params.outdir}/populations/${sample_list.simpleName}/individual_chr_vcfs/", mode: 'copy'
    echo true 

    input:
    set file(sample_list), file(vcf), file(vcf_index) from ch_multiVCF
    each file(fasta) from ch_fasta
    each file(fai) from ch_fai
    each file(dict) from ch_dict

    output:
    set val("${sample_list.simpleName}"), file("${vcf.baseName}.${sample_list.simpleName}.vcf") into (ch_pops_vcfs_to_bcftools, ch_pops_vcfs_to_bcftools_to_inspect)

    script:
    """
    gatk SelectVariants \
    -R ${fasta} \
    -V $vcf \
    -O ${vcf.baseName}.${sample_list.simpleName}.vcf \
    --sample-name ${sample_list}  \
    --restrict-alleles-to BIALLELIC \
    --select-type-to-include SNP \
    --verbosity ERROR 2> stderr.txt
    """
}

println('Inspecting ch_pops_vcfs_to_bcftools_to_inspect: (should have [pop_name, vcf])')
ch_pops_vcfs_to_bcftools_to_inspect.view()

process RecodeID {

    tag {"${pop_name}-${vcf.baseName}"}
    echo true 

    input:
    set val(pop_name), file(vcf) from ch_pops_vcfs_to_bcftools

    output:
    set val("${pop_name}"), file("${vcf.baseName}.${pop_name}.vcf.gz") into (ch_pops_vcfs, ch_pops_vcfs_to_inspect)

    script:
    """
    # Recode ID to: chr:pos:ref:alt
    bcftools annotate -x ID $vcf | bcftools annotate --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' > ${vcf.baseName}.${pop_name}.vcf

    bgzip -c ${vcf.baseName}.${pop_name}.vcf > ${vcf.baseName}.${pop_name}.vcf.gz
   """
}

println('Inspecting ch_pops_vcfs_to_inspect: (should have [pop_name, [vcf1, vcf2, vcf3]])')
ch_pops_vcfs_to_inspect
                        .groupTuple(by: 0)
                        //.map {pop_name, subsetted_vcf -> subsetted_vcf}
                        .view()

ch_grouped_pop_vcfs = ch_pops_vcfs.groupTuple(by: 0)



process GatherVCFs {

    tag "${pop_name}"
    publishDir "${params.outdir}/populations/${pop_name}/subsampled_multisample_vcf/", mode: 'copy'

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


process GetFrqCounts {

    tag "${pop_name}"
    publishDir "${params.outdir}/${pop_name}/plink_metrics/", mode: 'copy'

    input:
    set val(pop_name), file(all_chr_vcf) from ch_plink_count_freqs

    output:
    file("*") into ch_plink_results
    set val("${pop_name}"), file("${pop_name}.frq.counts") into  ch_plink_frq_counts_pop_tables
    file("${pop_name}.frq.counts") into ch_plink_frq_counts_for_panel

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

// Take one of the files
ch_panel_base = ch_plink_frq_counts_for_panel.take( 1 )

process GetPanelBase {

    tag "panel template"
    publishDir "${params.outdir}/FreqCountsDataframes/", mode: 'copy'

    input:
    file(frq_counts) from ch_panel_base

    output:
    file("template.panel.csv") into ch_panel_base_dataframe

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    data.table::setDTthreads(${task.cpus})

    frq_file <- list.files(getwd(), full.names = TRUE, pattern = ".frq.counts")
    frq_counts           <- data.table::fread(frq_file)
    colnames(frq_counts) <- c("chr","rs","a1","a2", "c1","c2","gpos" )
    
    # Extract pos from rs column
    frq_counts\$pos <- stringr::str_split_fixed(frq_counts[["rs"]], ":", n = Inf)[,2]
    frq_counts <- frq_counts[,c("chr", "rs", "gpos", "pos", "a1", "a2")]

    # Write file in .panel.csv
    data.table::fwrite(frq_counts, file = paste0("template.panel.csv"), sep = ",", col.names = TRUE)
    """
    }

process GetPopTables {

    tag "${pop_name}"
    publishDir "${params.outdir}/PopTables/", mode: 'copy'

    input:
    set val(pop_name), file(frq_counts) from ch_plink_frq_counts_pop_tables

    output:
    file("${pop_name}.pop.csv") into ch_pop_dataframes_for_panel

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    data.table::setDTthreads(${task.cpus})

    frq_file <- list.files(getwd(), full.names = TRUE, pattern = ".frq.counts")
    frq_counts           <- data.table::fread(frq_file)
    colnames(frq_counts) <- c("chr","rs","a1","a2", "c1","c2","gpos" )

    # Extract pos from rs column
    frq_counts\$pos <- stringr::str_split_fixed(frq_counts[["rs"]], ":", n = Inf)[,2]
    frq_counts\$pop <- paste0(frq_counts\$c1, ",",  frq_counts\$c2)
    pop_table       <- frq_counts [, c("rs", "pop")]
    colnames(pop_table) <- c("rs", "${pop_name}")

    # Write file in .pop.csv
    data.table::fwrite(pop_table, file = paste0("${pop_name}", ".pop.csv"), sep = ",", col.names = TRUE)
    """
    }

process JoinPanel {

    tag "Joining panel"
    publishDir "${params.outdir}/RefPanel/", mode: 'copy'

    input:
    file(panel_csv) from ch_panel_base_dataframe
    file(pop_csv)   from ch_pop_dataframes_for_panel.collect()

    output:
    file("*") into ch_freq_dataframes

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)
    data.table::setDTthreads(${task.cpus})

    panel_csv_path   <- list.files(getwd(), full.names = TRUE, pattern = ".panel.csv")

    # Collect all the .csv files that match the pattern ".pop.csv" (engineered so in previous process mua ha)
    all_pop_csv_paths <- list.files(getwd(), full.names = TRUE, pattern = ".pop.csv")
    all_csv_paths <- c(panel_csv_path, all_pop_csv_paths)

    all_csv <- lapply(all_csv_paths, data.table::fread)

    # (Redundant but I <3 this)
    # Collect all the names of the .csv files that match the pattern ".pop.csv"
    names(all_csv) <- gsub(".csv","", basename(all_csv_paths), fixed = TRUE)

    # Now all the data.tables are in a lists already - how convenient! Just in time for plyr::join()
    panel <- Reduce(function(dtf1, dtf2) plyr::join(dtf1, dtf2, by = "rs", type = "inner"), all_csv)

    # Write file in ref.panel.txt
    data.table::fwrite(panel, file = paste0("refpanel_", length(all_pop_csv_paths), "pops.txt"), sep = " ", col.names = TRUE)

    """
    }