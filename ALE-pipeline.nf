#!/usr/bin/env nextflow

params.input_files = 'ufboots/*.ufboot'
params.species_tree_files = '*.tree'
params.genes_map = 'map_genes.txt'
params.species_map = 'map_species.txt'
params.output_ale = 'ALE_results'
params.output_trees = 'species_trees'
params.output_samples = 'ufboots'
params.outgroup_taxa = '[""]'
params.python3 = '/usr/local/bin/anapy3'
params.fraction_missing = false


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
  """
  #! ${params.python3}
  print("${species_tree.baseName}")

  from ete3 import *
  from ETE3_Utils import *
  from ETE3_styles import *

  tree = Tree('$species_tree')

  new_root = get_ancestor([tree.get_leaves_by_name(s)[0] for s in $params.outgroup_taxa])
  tree.set_outgroup(new_root)
  print(tree.write(), file=open("${species_tree.baseName}_root.tree", 'w'))

  tree.convert_to_ultrametric()
  species_map = {line.split()[0]: line.split()[1] for line in open('map_species.txt')}
  for leaf in tree.get_leaves():
    leaf.name = species_map[leaf.name]
  print(tree.write(), file=open("${species_tree.baseName}_clean.tree", 'w'))
  """
}

process cleanNames{
  input:
  file bootstrap
  file 'map_genes.txt' from genes_map.first()

  output:
  file "${bootstrap}.clean" into bootstrap_clean

  maxForks 1

  script:
  """
  cp $bootstrap "${bootstrap}.preclean"
  anapy3 /local/two/Software/IdeaProjects/Python3_scripts/replace_names.py -i map_genes.txt -f "${bootstrap}.preclean"
  sed -r 's/_i_/???/g' "${bootstrap}.preclean" | sed -r 's/_//g' | sed 's/???/_/g' | sed -r 's/\\.[1|2]:/:/g' > "${bootstrap}.clean"
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
  """
  #! ${params.python3}

  with open("$gene") as f, open('events.txt', 'w') as outhandle:
    for line in f:
      if any([line.startswith(x) for x in ['S_terminal_branch', 'S_internal_branch']]):
        line = line.split()
        print("\\t".join(['$species_tree', "${gene.baseName}"] + line[1:]), file=outhandle)
  """
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
