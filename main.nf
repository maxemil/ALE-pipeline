#!/usr/bin/env nextflow

def startup_message() {
    log.info "============================================"
    log.info "               ALE pipeline"
    log.info "============================================"
    log.info "Input gene tree samples (GT)   : $params.input_files"
    log.info "Input species trees (ST)       : $params.species_tree_files"
    log.info "Map species names in GT        : $params.genes_map"
    log.info "Map species names in ST        : $params.species_map"
    log.info "Output dir for ALE ml          : $params.output_ale"
    log.info "Output dir for ALE observe     : $params.output_samples"
    log.info "Output dir for ST and pdfs     : $params.output_trees"
    log.info "Outgroup taxa use to root ST   : $params.outgroup_taxa"
    log.info "Currently looks for separator _i_ between species and genes"
    log.info ""
}

startup_message()

//bootstrap = Channel.from(file(params.input))
if(params.fraction_missing){fraction_missing = Channel.from(file(params.fraction_missing))}else{fraction_missing = Channel.from(params.fraction_missing)}
bootstrap = Channel.fromPath(params.input_files)
species_tree = Channel.fromPath(params.species_tree_files)
genes_map = Channel.from(file(params.genes_map))
Channel.from(file(params.species_map)).into { species_map; species_map_copy }

process cleanSpeciesTree{
  input:
  file species_tree
  file 'map_species.txt' from species_map.first()

  output:
  file "${species_tree.baseName}_clean.tree" into clean_species_tree
  file "${species_tree.baseName}_root.tree" into rooted_species_tree

  publishDir params.output_trees, mode: 'copy'

  script:
  template 'cleanSpeciesTree.py'
}

process cleanNames{
  input:
  file bootstrap
  file 'map_genes.txt' from genes_map.first()

  output:
  file "${bootstrap}.clean" into bootstrap_clean

  // maxForks 1

  script:
  """
  cp $bootstrap "${bootstrap}.clean"
  $workflow.projectDir/scripts/replace_names.py -i map_genes.txt -f "${bootstrap}.clean"
  sed -r -i 's/([^_i])_([^i_])/\\1\\2/g;s/_i_/_/g;s/\\.[1|2]:/:/g' "${bootstrap}.clean"
  """
}

process aleObserve{
  input:
  file bootstrap_clean

  output:
  file "${bootstrap_clean}.ale" into aleObserved

  publishDir params.output_samples, mode: 'copy'
  validExitStatus 0,1
  container true
  errorStrategy 'retry'
  maxRetries 5

  script:
  """
  maxemil/alesuite:latest ALEobserve \$PWD/$bootstrap_clean
  """
}

species_tree_vs_ale = clean_species_tree.combine(aleObserved)

process aleMlUndated{
  input:
  set file(species_tree), file(ale) from species_tree_vs_ale
  file fraction_missing from fraction_missing.first()
  // each species_tree from clean_species_tree

  output:
  set val("${species_tree.baseName}"), file("${ale}.ucons_tree") into uconsTrees
  set val("${species_tree.baseName}"), file("${ale}.uml_rec") into umlReconsiliation
  set val("${species_tree.baseName}"), file("${ale}.uTs") into uTransfers

  publishDir "${params.output_ale}/${species_tree.baseName}", mode: 'copy'
  container true
  stageInMode 'copy'
  errorStrategy 'retry'
  maxRetries 5

  script:
  if (fraction_missing){
      """
      maxemil/alesuite:latest ALEml_undated \$PWD/$species_tree \$PWD/$ale fraction_missing=\$PWD/fractionMissingGenes.txt
      """
  } else {
      """
      #/local/two/Software/ALE/ALE-build/bin/ALEml_undated \$PWD/$species_tree $ale
      maxemil/alesuite:latest ALEml_undated \$PWD/$species_tree \$PWD/$ale
      """
  }
}

process extractDTLEvents{
  input:
  set val(species_tree), file(gene) from umlReconsiliation

  output:
  file "events.txt" into events

  // publishDir "${params.output_ale}/${species_tree}", mode: 'copy'

  script:
  template "extractDTLevents.py"
}

all_events = events.collectFile(name: 'events.txt', storeDir: params.output_ale)

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
