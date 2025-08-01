/* 
  This workflow generates phylogenetic trees that include predefined reference
  sequences and observed sequences
 */
workflow GENERATE_TREES {

  take: 
  sequencing_report
  venv_path
  targets_file

  main:

  ch_versions = Channel.empty()                           

  CREATE_FASTA_FOR_TREES(sequencing_report, venv_path, targets_file)
  ch_versions = ch_versions.mix(CREATE_FASTA_FOR_TREES.out.versions)

  MAKE_TREE_ALIGNMENT(CREATE_FASTA_FOR_TREES.out.tree_fasta.flatten())
  ch_versions = ch_versions.mix(MAKE_TREE_ALIGNMENT.out.versions)

  MAKE_ML_TREE(MAKE_TREE_ALIGNMENT.out.tree_msa)
  ch_versions = ch_versions.mix(MAKE_ML_TREE.out.versions)

  VIEW_PHYLO_TREE(MAKE_ML_TREE.out.treefile.combine(venv_path))
  ch_versions = ch_versions.mix(VIEW_PHYLO_TREE.out.versions)
  
  emit:
  tree_files = MAKE_TREE_ALIGNMENT.out.tree_msa.collect()
  tree_pdf   = VIEW_PHYLO_TREE.out.tree_pdf.collect()
  versions   = ch_versions
}

/*
   Split up assigned observed sequeces by target
   for making trees
*/
process CREATE_FASTA_FOR_TREES {
  publishDir "${params.tree_outdir}", mode: 'link'
  tag "all"

  label 'process_low'

  // if using conda 
  conda "$baseDir/conda/python_conda_environment.yaml"                               

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/python:3.10.4" 
  } else {
      container "quay.io/biocontainers/python:3.10.4"
  }

  when:
  params.make_trees

  input:
  path (sequencing_report) 
  path (venv_path) 
  path (targets_file) 

  output:
  path("*_all.fasta")  , emit: tree_fasta
  path "versions.yml"  , emit: versions                                         

  script:
  // only need to activate the venv for singularity
  def activate_venv_command = workflow.containerEngine == 'singularity' ? "source ${venv_path}/bin/activate" : ""
  """
  $activate_venv_command
  python3 ${params.script_dir}/MPAS_create_fasta.py $sequencing_report $targets_file

  cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
        python: \$(python --version 2>&1 | sed 's/Python //g')                  
  END_VERSIONS    
  """
}

/*
   Build multiple-sequencing alignments for each group of sequences using MAFFT. 

   MAFFT documentation : https://mafft.cbrc.jp/alignment/software/manual/manual.html
*/
process MAKE_TREE_ALIGNMENT {
  publishDir "${params.tree_outdir}", mode: 'link'
  tag "all"

  label 'process_medium'

  // if using conda 
  conda "${moduleDir}/mafft_environment.yml"

  // if using singularity 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/mafft:7.505--hec16e2b_0"
  } else {
      container "https://depot.galaxyproject.org/singularity/mafft:7.505--hec16e2b_0"
  }

  input:
  path(all_fasta) 
  
  output:
  path("mafft_${all_fasta}") , emit: tree_msa
  path "versions.yml"        , emit: versions                                         

  shell:
  """
  mafft --adjustdirection --quiet --auto --nuc "$all_fasta" > "mafft_${all_fasta}"

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      mafft: \$(mafft --version 2>&1 | sed 's/^v//' | sed 's/ (.*)//')
  END_VERSIONS
  """
}

/*
   Build maximum likelihood for each group of sequences using IQ-TREE. 

   IQ-TREE documentation: www.iqtree.org/doc/
*/
process MAKE_ML_TREE {

  label 'process_medium'
  tag "all"

  // if using conda 
  conda "${moduleDir}/iqtree_environment.yml"

  // if using singularity 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/iqtree:2.2.0.3--hb97b32f_1"
  } else {
      container "quay.io/biocontainers/iqtree:2.2.0.3--hb97b32f_1"
  }
  
  input:
  path(all_alignment) 

  output:
  path("tree_${all_alignment.baseName}.treefile"), emit: treefile 
  path "versions.yml"                            , emit: versions                                         
  
  shell:
  """
  iqtree -s $all_alignment -st DNA -quiet -m MFP -pre tree_${all_alignment.baseName}   

  cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
        iqtree: \$(iqtree --version 2>&1 | grep -o "version [^ ]* " | sed 's/version //')
  END_VERSIONS    
  """ 
}

/*
   Creates pdf files of each ML tree using ToyTree. 

   ToyTree documentation: https://toytree.readthedocs.io/en/latest/
*/
process VIEW_PHYLO_TREE {
  publishDir "${params.tree_outdir}", mode: 'link'
  tag "all"

  label 'process_low'

  // if using conda 
  conda "$baseDir/conda/python_conda_environment.yaml"                               

  // if using singularity 
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/python:3.10.4" 
  } else {
      container "quay.io/biocontainers/python:3.10.4"
  }

  input:
  tuple path(iqtree), path(venv_path)

  output:
  path("${iqtree.baseName}.pdf")  , emit: tree_pdf 
  path "versions.yml"             , emit: versions                                         

  shell:
  // only need to activate the venv for singularity
  def activate_venv_command = workflow.containerEngine == 'singularity' ? "source ${venv_path}/bin/activate" : ""
  """
  $activate_venv_command
  python3 ${params.script_dir}/MPAS_view_tree.py $iqtree

  cat <<-END_VERSIONS > versions.yml                                          
    "${task.process}":                                                          
        python: \$(python --version 2>&1 | sed 's/Python //g')                  
  END_VERSIONS    
  """
}


