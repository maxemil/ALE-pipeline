## How to run an ALE analysis on a set of genomes/Proteomes
### Clustering of proteins
* prepare db of all Proteomes
* blast all vs. all to get pairwise hits
* use clustering tools such as
  * mcl or
  * silix/hifix to get clusters/OGs

```
silix $prefix.$FASTA $prefix.blastfilepath --prefix $prefix --ident $i --overlap $o --net > $prefix.fnodes
hifix -t 20 -n 150 $prefix.$FASTA $prefix.net $prefix.fnodes > $prefix.HFX.fnodes
```

### Preparation of ALE input
* To run, ALE needs
  * a (correct as it can be) species tree
  * A sample of trees for a gene cluster (e.g. a bootstrap sample or mcmc sampled trees) in one file, one tree per line

### Running ALE
* Running ALE is then simplified with the nextflow pipeline (executed for all clusters at once)

```
nextflow run ALE-pipeline/main.nf --species_tree_files "prefix*.tree" \
                                  --input_files samples \
                                  --input_extension '.ufboot' \
                                  --outgroup_taxa '["tax1", "tax2"]' \
                                  --fraction_missing fractionMissingGenes.txt \
                                  --genes_map map_genes.txt \
                                  --species_map map_species.txt \
                                  --single_cluster single_cluster_faas/* \
                                  --output_ale output_dir
                                  --separator "_i_"
```
* The pipeline will output a file 'events.txt' with the inferred ODTL events
* you can also include small clusters (<4 members) to get more correct numbers for each leaf (and with that for internal nodes as well).

#### small test for the pipeline
* you can run a test with the files in `tests`, so by doing the below command you should get the directories `species_trees` with the rooted species trees and pdfs, the `ufboots` directory with the ALE objects and the `ALE_results` with the main `ALEml_undated` output.
  ```
  # run test:
  nextflow run main.nf --species_tree tests/species_tree_test.new \
                       --input_files tests\
                       --input_extension '.ufboot' \
                       --genes_map tests/map_genes.txt \
                       --species_map tests/map_species.txt \
                       --outgroup_taxa '["BIN125-1"]'

  # clean up:
  rm -r work ALE_results species_trees ufboots .nextflow*
  ```

### interpreting ALE results
* with the nextflow pipeline, there is a script to extract events for certain nodes
* include EGGnog/COG annotations to have an idea of what the clusters are doing
* you can also get numbers for gains losses when looking at pairs of nodes

```
python3 bin/match_cluster_eggnog.py --faas clusters/*.faa \
                                  --annotations Emapper_annotations.tsv \
                                  --species_tree species_tree_name \
                                  --events events.txt \
                                  --ancestor_nodes 100 99 98 97
                                  --node_pairs 100 99 100 98 99 97 \
                                  --cog_annotation cognames2003-2014.tab \
                                  --threshold 0.5
```                               
