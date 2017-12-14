## How to run an ALE analysis on a set of genomes/Proteomes
### Clustering of proteins
* prepare db of all Proteomes
* blast all vs. all to get pairwise hits
* use clustering tools such as
  * mcl or
  * silix/hifix to get clusters/OGs

### Preparation of ALE input
* To run, ALE needs
  * a (correct as it can be) species tree
  * A sample of trees for a gene cluster (e.g. a bootstrap sample or mcmc sampled trees)

### Running ALE
* Running ALE is then simplified with the nextflow pipeline (executed for all clusters at once)
* The pipeline will output a file 'events.txt' with the inferred ODTL events
* you can also include small clusters (<4 members) to get more correct numbers for each leaf (and with that for internal nodes as well).

### interpreting ALE results
* with the nextflow pipeline, there is a script to extract events for certain nodes
* include EGGnog/COG annotations to have an idea of what the clusters are doing
* you can also get numbers for gains losses when looking at pairs of nodes
