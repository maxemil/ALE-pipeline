# ALE-pipeline

- this is supposed to be a nice pipeline for running ALE on several gene clusters and collecting the results
- it can also test several species trees at the same time
- all parameters in the nextflow.config file can be changed on the command line, e.g. the name of the outgroup taxa
- you need to add the Python_lib repo to your Pythonpath
- typical usage would be like so:
  ```
  nextflow run ALE-pipeline/main.nf --species_tree_files "prefix*.tree" \
                                    --input_files samples \
                                    --input_extension '.extension' \
                                    --outgroup_taxa '["tax1", "tax2"]' \
                                    --fraction_missing fractionMissingGenes.txt \
                                    --genes_map map_genes.txt \
                                    --species_map map_species.txt
  ```
- I use to code the names for species both in the species and in the gene tree to avoid that source of errors
- ALE sometime simply crashes, then the pipeline can be resumed by adding `-resume` to the invocation
- you can run a test with the files in `tests`, so by doing the below command you should get the directories `species_trees` with the rooted species trees and pdfs, the `ufboots` directory with the ALE objects and the `ALE_results` with the main `ALEml_undated` output.
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

# Troubleshooting
* If you get an Error 'Can't root with myself' or similar, this usually means that the outgroup you specified for the species tree is not monophyletic in that tree. Try rerooting by hand first...


# Getting an overwie of gains and losses for certain node_pairs
assuming you have annotated all proteins in your clusters (cluster_faas) with the emapper, and the `.annotation` files are in eggnog_output, you can get a summary of losses and gains as well as general description of you clusters with the script `match_cluster_eggnog.py`:
```
python3 scripts/match_cluster_eggnog.py --faas cluster_faas/* \
                                        --annotations eggnog_output/*.annotations  \
                                        --ancestor_nodes 83 82 81 \
                                        --threshold 0.5 \
                                        --species_tree Alphaproteobacteria_species_recoded_clean \
                                        --node_pairs 84 83 83 82 82 81 \
                                        --cluster_out 'cluster_OG.tab' \
                                        --annotation_extension '.faa.emapper.annotations' \
                                        --cluster_extension '.bmge.aln.ufboot.clean.ale' \
                                        --oneline
```
* this will produce a file `cluster_OG.tab` with annotation info for all clusters, a set of files for each nodes, where only those clusters are in that are present (> threshold) in that node and finally a set of files that only list gains and losses for each node pair (OBS node pairs params has to be a multiple of two).
