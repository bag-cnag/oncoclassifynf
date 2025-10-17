
process CLASSIFY {
    tag "$meta.id"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
        'biocontainers/YOUR-TOOL-HERE' }"

    input:
    tuple val(meta), path(vcf) , path(index)
    path database_config

    output:
    tuple val(meta), path("*.classify.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
   python /OncoClassify/OncoClassify.py \\
        $args \\
        -i ${vcf} \\
        -o ${prefix}.classify \\
        -d ${database_config}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncoclassify: \$oncoclassify
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo $args

    touch ${prefix}.classify.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncoclassify: \$oncoclassify
    END_VERSIONS
    """
}
