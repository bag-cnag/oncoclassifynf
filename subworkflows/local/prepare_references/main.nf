include { TABIX_BGZIPTABIX as TABIX_BGZIPINDEX_VCFANNOEXTRA      } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_TABIX as TABIX_VCFANNOEXTRA    } from '../../../modules/nf-core/tabix/tabix'

workflow PREPARE_REFERENCES {

    take:
    ch_vcfanno_extra_unprocessed // channel: [mandatory] [ val(meta), path(vcf) ]

    main:

    ch_versions = Channel.empty()
    ch_vcfanno_extra = Channel.empty()
    ch_vcfanno_bgzip = Channel.empty()
    ch_vcfanno_index = Channel.empty()

    ch_vcfanno_extra_unprocessed
    .branch { it ->
        bgzipindex: !it[1].toString().endsWith(".gz")
        index: it[1].toString().endsWith(".gz")
    }
    .set { ch_vcfanno_tabix_in }

    TABIX_VCFANNOEXTRA(ch_vcfanno_tabix_in.index).tbi
        .join(ch_vcfanno_tabix_in.index)
        .map { meta, tbi, vcf -> return [[vcf,tbi]]}
        .set {ch_vcfanno_index}

    TABIX_BGZIPINDEX_VCFANNOEXTRA(ch_vcfanno_tabix_in.bgzipindex)
    Channel.empty()
        .mix(TABIX_BGZIPINDEX_VCFANNOEXTRA.out.gz_tbi, TABIX_BGZIPINDEX_VCFANNOEXTRA.out.gz_csi)
        .map { meta, vcf, index -> return [[vcf,index]] }
        .set {ch_vcfanno_bgzip}

    Channel.empty()
        .mix(ch_vcfanno_bgzip, ch_vcfanno_index)
        .collect()
        .set{ch_vcfanno_extra}


    ch_versions = ch_versions.mix(TABIX_BGZIPINDEX_VCFANNOEXTRA.out.versions)
    ch_versions = ch_versions.mix(TABIX_VCFANNOEXTRA.out.versions)


    emit:
    vcfanno_extra = ch_vcfanno_extra.ifEmpty([[]])     // channel: [ [path(vcf), path(tbi)] ]    

    versions = ch_versions                                     // channel: [ versions.yml ]
}
