# ALE-pipeline

- this is supposed to be a nice pipeline for running ALE on several gene clusters and collecting the results
- it can also test several species trees at the same time
- all parameters in the nextflow.config file can be changed on the command line, e.g. the name of the outgroup taxa
- you need to add the Python_lib repo to your Pythonpath
- For typical usage and a small tutorial, see TUTORIAL.md
- I use to code the names for species both in the species and in the gene tree to avoid that source of errors

# Troubleshooting
* If you get an Error 'Can't root with myself' or similar, this usually means that the outgroup you specified for the species tree is not monophyletic in that tree. Try rerooting by hand first...
* ALE sometime simply crashes, then the pipeline can be resumed by adding `-resume` to the invocation
