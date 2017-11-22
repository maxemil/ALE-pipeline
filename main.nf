#!/usr/bin/env nextflow

def startup_message() {
    log.info "=========================================================="
    log.info "                       ALE pipeline"
    log.info "Author                         : Max Emil Schön"
    log.info "email                          : max-emil.schon@icm.uu.se"
    log.info "=========================================================="
    log.info "Input gene tree samples (GT)   : $params.input_files"
    log.info "Input species trees (ST)       : $params.species_tree_files"
    log.info "Map species names in GT        : $params.genes_map"
    log.info "Map species names in ST        : $params.species_map"
    log.info "Output dir for ALE ml          : $params.output_ale"
    log.info "Output dir for ALE observe     : $params.output_samples"
    log.info "Output dir for ST and pdfs     : $params.output_trees"
    log.info "Extension of GT samples        : $params.input_extension"
    log.info "Outgroup taxa use to root ST   : $params.outgroup_taxa"
    log.info "Fraction of missing genes      : $params.fraction_missing"
    log.info "Python 3 binary used           : $params.python3"
    log.info ""
    log.info "Currently looks for separator $params.separator between species and genes"
    log.info "if no extensions is given for input files, takes all files "
    log.info "in that directory"
    log.info "Also need acces to Python_lib from ettemalab's bitbucket!"
    log.info ""
}

startup_message()

//bootstrap = Channel.from(file(params.input))
fraction_missing = create_channel_fraction_missing()
bootstrap = Channel.fromPath("$params.input_files/*$params.input_extension")
Channel.fromPath(params.species_tree_files).into { species_tree; species_tree_single_cluster }
Channel.from(file(params.genes_map)).into { genes_map; genes_map_single_clusters }
Channel.from(file(params.species_map)).into { species_map; species_map_copy }
single_cluster =  Channel.fromPath(params.single_cluster)

process cleanSpeciesTree{
  input:
  file species_tree
  file 'map_species.txt' from species_map.first()

  output:
  file "${species_tree.baseName}_clean.tree" into clean_species_tree
  file "${species_tree.baseName}_root.tree" into rooted_species_tree

  publishDir params.output_trees, mode: 'copy'
  tag {"${species_tree.simpleName}"}

  script:
  template 'cleanSpeciesTree.py'
}

process cleanNames{
  input:
  file bootstrap
  file 'map_genes.txt' from genes_map.first()

  output:
  file "${bootstrap}.clean" into bootstrap_clean

  tag {"${bootstrap.simpleName}"}

  script:
  """
  cp $bootstrap "${bootstrap}.clean"
  replace_names.py -i map_genes.txt -f "${bootstrap}.clean"
  sed -r -i 's/([^_i])_([^i_])/\\1\\2/g;s/$params.separator/_/g;s/\\.[1-9]:/:/g' "${bootstrap}.clean"
  """
}

process aleObserve{
  input:
  file bootstrap_clean

  output:
  file "${bootstrap_clean}.ale" into aleObserved

  publishDir params.output_samples, mode: 'copy'
  tag {"${bootstrap_clean.simpleName}"}

  script:
  """
  ALEobserve \$PWD/$bootstrap_clean
  """
}

species_tree_vs_ale = clean_species_tree.combine(aleObserved)

process aleMlUndated{
  input:
  set file(species_tree), file(ale) from species_tree_vs_ale
  file fraction_missing from fraction_missing.first()

  output:
  set val("${species_tree.baseName}"), file("${ale}.ucons_tree") into uconsTrees
  set val("${species_tree.baseName}"), file("${ale}.uml_rec") into umlReconsiliation
  set val("${species_tree.baseName}"), file("${ale}.uTs") into uTransfers

  publishDir "${params.output_ale}/${species_tree.baseName}", mode: 'copy'
  tag {"${ale.simpleName}"}

  script:
  if (fraction_missing ==~ /.*tmp/){
      """
      ALEml_undated \$PWD/$species_tree \$PWD/$ale
      """
  } else {
      """
      ALEml_undated \$PWD/$species_tree \$PWD/$ale fraction_missing=\$PWD/$fraction_missing
      """
  }
}

process extractDTLEvents{
  input:
  set val(species_tree), file(gene) from umlReconsiliation

  output:
  file "events.txt" into events

  tag {"${gene.simpleName}"}
  // publishDir "${params.output_ale}/${species_tree}", mode: 'copy'

  script:
  template "extractDTLevents.py"
}

process includeSingleClusters {
  input:
  file fasta from single_cluster
  file 'map_genes.txt' from genes_map_single_clusters.first()
  each species_tree from species_tree_single_cluster

  output:
  file "events.txt" into single_cluster_events

  tag {"${fasta.simpleName}"}

  script:
  """
  cp $fasta ${fasta}.clean
  replace_names.py -i map_genes.txt -f ${fasta}.clean
  grep '>' ${fasta}.clean | cut -d '>' -f 2 | awk -F '$params.separator' '{print "${species_tree.baseName}_clean\\t$fasta.simpleName\\t"\$1"\\t0\\t0\\t0\\t0\\t1"}' > events.txt
  """
}

all_events = events.concat(single_cluster_events).collectFile(name: 'events.txt', storeDir: params.output_ale)

process summmarizeDTLEvents{
  input:
  file all_events
  file species_map from species_map_copy.first()

  output:
  file "*.pdf" into tree_pdfs
  file "*.csv" into events_csv

  beforeScript = "for d in \$(ls -d ${workflow.launchDir}/${params.output_ale}/*/); do grep '^S:' `ls \$d*uml_rec | shuf -n 1 | cat -` | cut -f2 > ${workflow.launchDir}/${params.output_trees}/\$(basename \$d).tree; done"
  publishDir "${params.output_trees}", mode: 'copy'

  script:
  template 'visualize_ALE_results.py'
  }


def create_channel_fraction_missing() {
  if(params.fraction_missing){
    return(Channel.from(file(params.fraction_missing)))
  }else{
    File temp = new File("temp.tmp")
    temp.with{
      write "hello temp file!"
      deleteOnExit()
    }
    return(Channel.from(file("temp.tmp")))
  }
}
