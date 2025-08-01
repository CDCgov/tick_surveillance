process BLAST_MAKEBLASTDB {
    tag "$fasta"
    label 'process_medium'

    // if using conda                                                             
    conda "${moduleDir}/environment.yml"                                          
                                                                                
    // if using singularity                                                       
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/blast:2.16.0--hc155240_2"
    } else {                                                                      
        container "quay.io/biocontainers/blast:2.16.0--hc155240_2"                
    }                                                                             

    input:
    path fasta

    output:
    path 'blast_db'     , emit: db
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    makeblastdb \\
        -in $fasta \\
        $args
    mkdir blast_db
    mv ${fasta}* blast_db
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | sed 's/^.*blastn: //; s/ .*\$//')
    END_VERSIONS
    """
}
