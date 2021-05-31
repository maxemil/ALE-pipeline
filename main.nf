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
    log.info "Fasta files of clusters < 4    : $params.small_cluster"
    log.info "Python 3 binary used           : $params.python3"
    log.info ""
    log.info "Currently looks for separator $params.separator between species and genes"
    log.info "if no extensions is given for input files, takes all files "
    log.info "in that directory"
    log.info ""
}

startup_message()

//bootstrap = Channel.from(file(params.input))
fraction_missing = create_channel_fraction_missing()
regular_bootstrap = Channel.fromPath("$params.input_files/*$params.input_extension")
Channel.fromPath(params.species_tree_files).into { species_tree; species_tree_single_cluster }
Channel.from(file(params.genes_map)).into { genes_map; genes_map_single_clusters }
Channel.from(file(params.species_map)).into { species_map; species_map_copy }

if (params.small_cluster) {
    small_clusters = Channel.fromPath(params.small_cluster)
}else {
    small_clusters = Channel.empty()
}

process smallClusters {
  input:
  file fasta from small_clusters

  output:
  file "${fasta.simpleName}.small.faa" into small_faa optional true
  file "${fasta.simpleName}.single.faa" into single_faa optional true

  tag {"${fasta.simpleName}"}

  script:
  """
  if [ "\$(grep -c '^>' $fasta)" -lt 4 ] && [ "\$(grep -c '^>' $fasta)" -gt 1 ]
  then
    mv $fasta ${fasta.simpleName}.small.faa
  elif [ "\$(grep -c '^>' $fasta)" -eq 1 ]
  then
      mv $fasta ${fasta.simpleName}.single.faa
  fi
  """
}

process smallClustersBootstrap {
  input:
  file fasta from small_faa

  output:
  file "${fasta.simpleName}.clean" into small_bootstraps

  script:
  """
  #! ${params.python3}

  from Bio import SeqIO
  with open("${fasta.simpleName}.clean", 'w') as outhandle:
    rec_ids = [rec.id for rec in SeqIO.parse("$fasta", 'fasta')]
    if len(rec_ids) == 2:
      newick = "({},{});".format(rec_ids[0], rec_ids[1])
    else:
      newick = "(({},{}),{});".format(rec_ids[0], rec_ids[1], rec_ids[2])
    outhandle.write(newick)
  """
}

process cleanSpeciesTree {
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

bootstrap = regular_bootstrap.concat(small_bootstraps)

process cleanNames {
  input:
  file bootstrap
  file 'map_genes.txt' from genes_map.first()

  output:
  file "${bootstrap.simpleName}.clean.ufboot" into bootstrap_clean

  tag {"${bootstrap.simpleName}"}

  script:
  template 'replace_names_tree.py'
}

process aleObserve {
  input:
  file bootstrap_clean

  output:
  file "${bootstrap_clean}.ale" into aleObserved

  publishDir params.output_samples, mode: 'copy'
  tag {"${bootstrap_clean.simpleName}"}

  script:
  """
  ALEobserve $bootstrap_clean
  """
}

species_tree_vs_ale = clean_species_tree.combine(aleObserved)

process aleMlUndated {
  input:
  set file(species_tree), file(ale) from species_tree_vs_ale
  file fraction_missing from fraction_missing.first()

  output:
  set val("${species_tree.baseName}"), file("${ale}.ucons_tree") into uconsTrees optional true
  set val("${species_tree.baseName}"), file("${ale}.uml_rec") into umlReconsiliation
  set val("${species_tree.baseName}"), file("${ale}.uTs") into uTransfers optional true

  publishDir "${params.output_ale}/${species_tree.baseName}", mode: 'copy'
  tag {"${species_tree.baseName} - ${ale.simpleName}"}

  script:
  if (fraction_missing ==~ /.*tmp/){
      """
      exitcode=0
      ALEml_undated $species_tree $ale || exitcode=\$?
      if [ \$exitcode -eq 134 -a -s ${ale}.uml_rec ]
        then
          echo "error, but output file present"
          exit 0
        else
          exit \$exitcode
      fi
      """
  } else {
      """
      exitcode=0
      ALEml_undated $species_tree $ale fraction_missing=$fraction_missing || exitcode=\$?
      if [ \$exitcode -eq 134 -a -s ${ale}.uml_rec ]
        then
          echo "error, but output file present"
          exit 0
        else
          exit \$exitcode
      fi
      """
  }
}

umlReconsiliation.into{ umlReconsiliation_speciesTree; umlReconsiliation_events }

umlReconsiliation_speciesTree_unique = umlReconsiliation_speciesTree.unique { it[0] }

process extractSpeciesTree {
  input:
  set val(species_tree), file(gene) from umlReconsiliation_speciesTree_unique

  output:
  file "${species_tree}.tree" into species_tree_output

  publishDir "${workflow.launchDir}/${params.output_trees}", mode: 'copy', overwrite: 'true'
  tag {"${species_tree}"}

  script:
  """
  grep '^S:' $gene | cut -f2 > ${species_tree}.tree
  """
}

process extractDTLEvents {
  input:
  set val(species_tree), file(gene) from umlReconsiliation_events

  output:
  file "events.txt" into events

  tag {"$species_tree - ${gene.simpleName}"}
  // publishDir "${params.output_ale}/${species_tree}", mode: 'copy'

  script:
  template "extractDTLevents.py"
}

process includeSingleClusters {
  input:
  file fasta from single_faa
  file 'map_genes.txt' from genes_map_single_clusters.first()
  each species_tree from species_tree_single_cluster

  output:
  file "events.txt" into single_cluster_events

  tag {"${species_tree.baseName} - ${fasta.simpleName}"}

  script:
  """
  #! ${params.python3}
  header = open('$fasta').readline().strip().replace('>','')
  name_map = {line.split('\t')[0]:line.split('\t')[1].strip() for line in open('map_genes.txt')}
  clean_id = name_map[header.split('$params.separator')[0]]
  with open('events.txt', 'w') as outfile:
      print("\t".join(['${species_tree.baseName}_clean', '${fasta.simpleName}', clean_id, '0', '0', '0', '0', '1']), file=outfile)
  """
}

all_events = events.concat(single_cluster_events).collectFile(name: 'events.txt', storeDir: params.output_ale)

process summmarizeDTLEvents {
  input:
  file all_events
  file species_map from species_map_copy.first()
  file species_trees from species_tree_output.collect()

  output:
  file "*.pdf" into tree_pdfs
  file "*.csv" into events_csv

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
