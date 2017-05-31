#!/usr/bin/env nextflow

params.input_files = 'ufboots/*.ufboot'
params.species_tree_file = '*.tree'
params.genes_map = 'map_genes.txt'
params.species_map = 'map_species.txt'
params.output_ale = 'ALE_results'
params.python3 = "/usr/local/bin/anapy3"


//bootstrap = Channel.from(file(params.input))
bootstrap = Channel.fromPath(params.input_files)
species_tree = Channel.fromPath(params.species_tree_file)
genes_map = Channel.from(file(params.genes_map))
species_map = Channel.from(file(params.species_map))

process cleanSpeciesTree{
  input:
  file species_tree
  file 'map_species.txt' from species_map.first()

  output:
  file "${species_tree.baseName}_clean.tree" into clean_species_tree

  publishDir "species_trees", mode: 'copy'

  script:
  """
  #! ${params.python3}
  print("${species_tree.baseName}")

  from ete3 import *

  tree = Tree('$species_tree')
  new_root = tree.get_leaves_by_name('BIN125-1')[0].up
  tree.set_outgroup(new_root)

  tree.convert_to_ultrametric()

  species_map = {line.split()[0]: line.split()[1] for line in open('map_species.txt')}

  edge = len(list(tree.traverse())) - len(list(tree.get_leaves()))
  for node in tree.traverse():
    if not node.is_leaf():
      node.name = edge
      edge -= 1
    else:
      node.name = species_map[node.name]
  print(tree.write(format=3), file=open("${species_tree.baseName}_clean.tree", 'w'))
  """
}

process cleanNames{
  input:
  file bootstrap
  file 'map_genes.txt' from genes_map.first()

  output:
  file "${bootstrap}.clean" into bootstrap_clean

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

  publishDir "ufboots", mode: 'copy'
  validExitStatus 0,1
//  stageInMode 'copy'

  script:
  """
  #docker run -v \$PWD:\$PWD alesuite ALEobserve \$PWD/$bootstrap_clean
  /local/two/Software/ALE/ALE-build/bin/ALEobserve $bootstrap_clean
  """
}

process aleMlUndated{
  input:
  file ale from aleObserved
  each species_tree from clean_species_tree

  output:
  set val("${species_tree.baseName}"), file("${ale}.ucons_tree") into uconsTrees
  set val("${species_tree.baseName}"), file("${ale}.uml_rec") into umlReconsiliation
  set val("${species_tree.baseName}"), file("${ale}.uTs") into uTransfers

  publishDir "${params.output_ale}/${species_tree.baseName}", mode: 'copy'
  //errorStrategy 'ignore'

  script:
  """
  #docker run -v \$PWD:\$PWD alesuite ALEml_undated \$PWD/$species_tree \$PWD/$ale
  /local/two/Software/ALE/ALE-build/bin/ALEml_undated $species_tree $ale
  """

}

process summmarizeDTLEvents{
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
