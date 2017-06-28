# ALE-pipeline

- this is supposed to be a nice pipeline for running ALE on several gene clusters and collecting the results
- it can also test several species trees at the same time
- all parameters in the nextflow.config file can be changed on the command line, e.g. the name of the outgroup taxa
- typical usage would be like so:
  ```
  nextflow run ALE-pipeline/main.nf --species_tree_files "prefix*.tree"
                                    --outgroup_taxa '["tax1", "tax2"]'
                                    --fraction_missing fractionMissingGenes.txt
                                    --genes_map map_genes.txt
                                    --species_map map_species.txt
  ```
- I use to code the names for species both in the species and in the gene tree to avoid that source of errors
- ALE sometime simply crashes, then the pipeline can be resumed by adding `-resume` to the invocation
- you can run a test with the files in `tests`
- you need to add the Python_lib repo to your Pythonpath
