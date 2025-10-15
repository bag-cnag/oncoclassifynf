/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_oncoclassifynf_pipeline'

//
// MODULE: Installed directly from nf-core/modules
//
include { BCFTOOLS_NORM          } from '../modules/nf-core/bcftools/norm'
include { BCFTOOLS_PLUGINFILLTAGS} from '../modules/nf-core/bcftools/pluginfilltags'
include { VCFANNO                } from '../modules/nf-core/vcfanno'

//
// SUBWORKFLOW: Installed directly from nf-core/subworkflows
//
include { VCF_ANNOTATE_ENSEMBLVEP_SNPEFF } from '../subworkflows/nf-core/vcf_annotate_ensemblvep_snpeff/main'
include { PREPARE_REFERENCES    } from '../subworkflows/local/prepare_references'
//
// MODULE: Custom local module
//
include { CLASSIFY          } from '../modules/local/classify'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ONCOCLASSIFYNF {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Prepare channels from params
    ch_genome_fasta = Channel.fromPath(params.fasta).map { it -> [[id:it.simpleName], it] }.collect()

    ch_vcfanno_extra_unprocessed = params.vcfanno_extra_resources ? Channel.fromPath(params.vcfanno_extra_resources).map { it -> [[id:it.baseName], it] }.collect()
                                                                : Channel.empty()
    ch_vcfanno_lua              = params.vcfanno_lua                        ? Channel.fromPath(params.vcfanno_lua).collect()
                                                                            : Channel.value([])
    ch_vcfanno_toml             = params.vcfanno_toml                       ? Channel.fromPath(params.vcfanno_toml).collect()
                                                                            : Channel.value([])
    ch_vcfanno_resources        = params.vcfanno_resources                  ? Channel.fromPath(params.vcfanno_resources).splitText().map{it -> it.trim()}.collect()
                                                                            : Channel.value([])

    ch_vep_cache = params.vep_cache ? Channel.fromPath(params.vep_cache) : Channel.value([])
    ch_snpeff_cache = params.snpeff_cache ? Channel.fromPath(params.snpeff_cache) : Channel.value([])


    ch_bcftools_regions = params.bcftools_regions ? Channel.fromPath(params.bcftools_regions) : Channel.value([])
    ch_bcftools_targets = params.bcftools_targets ? Channel.fromPath(params.bcftools_targets) : Channel.value([])
    ch_bcftools_samples = params.bcftools_samples ? Channel.fromPath(params.bcftools_samples) : Channel.value([])

    // Prepare VEP extra files

    vep_extra_files = []

    if (params.dbnsfp && params.dbnsfp_tbi) {
        vep_extra_files.add(file(params.dbnsfp, checkIfExists: true))
        vep_extra_files.add(file(params.dbnsfp_tbi, checkIfExists: true))
    }

    if (params.spliceai_snv && params.spliceai_snv_tbi && params.spliceai_indel && params.spliceai_indel_tbi) {
        vep_extra_files.add(file(params.spliceai_indel, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_indel_tbi, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv, checkIfExists: true))
        vep_extra_files.add(file(params.spliceai_snv_tbi, checkIfExists: true))
    }

    //
    // SUBWORKFLOW: PREPAREREFERENCES
    //
    PREPARE_REFERENCES(
        ch_vcfanno_extra_unprocessed
    )
    ch_versions = ch_versions.mix(PREPARE_REFERENCES.out.versions)
    ch_vcfanno_extra = PREPARE_REFERENCES.out.vcfanno_extra


    //
    // MODULE: BCFTOOLS_NORM
    //

    BCFTOOLS_NORM(
        ch_samplesheet,ch_genome_fasta
    )

    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)
    ch_vcf_norm = BCFTOOLS_NORM.out.vcf
        .join(BCFTOOLS_NORM.out.tbi)
        .map { 
            meta, vcf, tbi -> [meta, vcf, tbi] 
        }

    //
    // SUBWORKFLOW: VCF_ANNOTATE_ENSEMBLVEP_SNPEFF
    //
    VCF_ANNOTATE_ENSEMBLVEP_SNPEFF(
        ch_vcf_norm,
        ch_genome_fasta,
        params.vep_genome,
        params.vep_species,
        params.vep_version,
        ch_vep_cache,
        vep_extra_files,
        params.snpeff_db,
        ch_snpeff_cache,
        params.annotation_tools,
        params.sites_per_chunk
    )

    ch_versions = ch_versions.mix(VCF_ANNOTATE_ENSEMBLVEP_SNPEFF.out.versions)
    ch_in_vcfanno = VCF_ANNOTATE_ENSEMBLVEP_SNPEFF.out.vcf_tbi
        .combine(ch_vcfanno_extra)
        .map { meta, vcf, tbi, resources -> return [meta + [prefix: meta.prefix + "_vcfanno"], vcf, tbi, resources]}    


    //
    // MODULE: VCFANNO
    //
    VCFANNO(
        ch_in_vcfanno,
        ch_vcfanno_toml,
        ch_vcfanno_lua,
        ch_vcfanno_resources
    )

    ch_versions = ch_versions.mix(VCFANNO.out.versions)
    ch_vcfanno = VCFANNO.out.vcf
        .join(VCFANNO.out.tbi)
        .map { 
            meta, vcf, tbi -> [meta, vcf, tbi] 
            }

    //
    // MODULE: BCFTOOLS_FILLTAGS
    //
    BCFTOOLS_PLUGINFILLTAGS(
        ch_vcfanno,
        ch_bcftools_regions,
        ch_bcftools_targets,
        ch_bcftools_samples
        
    )
    ch_versions = ch_versions.mix(BCFTOOLS_PLUGINFILLTAGS.out.versions)
    ch_vcf_af = BCFTOOLS_PLUGINFILLTAGS.out.vcf

    //
    // MODULE: CLASSIFY
    //
    CLASSIFY(
        ch_vcf_af,
        Channel.fromPath(params.database_config, checkIfExists: true)
    )
    ch_versions = ch_versions.mix(CLASSIFY.out.versions)
    ch_classify = CLASSIFY.out.vcf

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'oncoclassifynf_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
