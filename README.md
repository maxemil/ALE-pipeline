# ALE-pipeline

- this is supposed to be a nice pipeline for running ALE on several gene clusters and collecting the results
- it can also test several species trees at the same time
- all parameters in the nextflow.config file can be changed on the command line, e.g. the name of the outgroup taxa
- you need to add the Python_lib repo to your Pythonpath
- typical usage would be like so:
  ```
  nextflow run ALE-pipeline/main.nf --species_tree_files "prefix*.tree" \
                                    --input_files ufboots \
                                    --input_extension '.ufboot' \
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
