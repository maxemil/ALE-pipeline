#!/usr/bin/env nextflow

params.input_files = 'ufboots/*.ufboot'
params.species_tree_file = '*.tree'
params.genes_map = 'map_genes.txt'
params.species_map = 'map_species.txt'
params.output_ale = 'ALE_results'
params.output_trees = 'species_trees'
params.output_samples = 'ufboots'
params.python3 = '/usr/local/bin/anapy3'


//bootstrap = Channel.from(file(params.input))
bootstrap = Channel.fromPath(params.input_files)
species_tree = Channel.fromPath(params.species_tree_file)
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

  tree = Tree('$species_tree')
  new_root = tree.get_leaves_by_name('BIN125-1')[0].up
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
  stdout x
  file "*.pdf" into tree_pdfs

  beforeScript = "for d in \$(ls -d ${workflow.launchDir}/${params.output_ale}/*/); do grep '^S:' `ls \$d*uml_rec | shuf -n 1 | cat -` | cut -f2 > ${workflow.launchDir}/${params.output_trees}/\$(basename \$d).tree; done"
  // beforeScript = "touch ${workflow.launchDir}/DONE"
  publishDir "${params.output_trees}", mode: 'copy'

  script:
  """
  #! ${params.python3}
  import sys
  sys.path.append('/local/two/Software/python_lib/')
  import pandas as pd
  import ete3
  from ETE3_Utils import *
  from numpy import log

  def layout(node):
    if node.is_leaf():
        N = AttrFace("name", fsize=14, fgcolor="black")
        faces.add_face_to_node(N, node, 0)
    if "weight" in node.features:
        # C = CircleFace(radius=log(node.weight), color="Red", style="sphere")
        C = RectFace(width=node.weight, height=3, fgcolor="Red", bgcolor="Red")
        faces.add_face_to_node(C, node, 0, position="float")

  def add_genome_size(treeO, treeC, species_map, events, filename):
    for n in treeO.traverse():
      name = ''
      if n.is_leaf():
        name = species_map[n.name]
      else:
        name = get_ancestor([treeC.get_leaves_by_name(species_map[l.name])[0] for l in n.get_leaves()]).name
      n.add_features(weight=events.loc[name]['copies'])
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    treeO.render(filename, tree_style=ts)

  header = ['species_tree', 'cluster_ID', 'node', 'duplications',
            'transfers', 'losses', 'originations', 'copies']
  df = pd.read_csv("events.txt", sep='\\t', names=header)

  species_map = {line.split()[0]: line.split()[1] for line in open("$species_map")}

  split_df = [df[df['species_tree'] == t] for t in set(df.species_tree)]
  for f in split_df:
    add_genome_size(
        ete3.Tree("$workflow.launchDir/$params.output_trees/%s.tree" % set(f.species_tree).pop().replace('_clean', '_root')),
        ete3.Tree("$workflow.launchDir/$params.output_trees/%s.tree" % set(f.species_tree).pop(), format=1),
        species_map,
        f.groupby('node').sum(),
        set(f.species_tree).pop() + ".pdf")
  """
}

x.subscribe{print it}
