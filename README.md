# selac_implementation

Sample code for debugging implementations

Test1: Rokas yeast data, one gene.

* gene1Yeast.fasta: FASTA file
* rokasYeast.tre: newick, with an extra taxon (Calb)
* tree.nex: Nexus tree file with extra taxon removed
* dna.nex: Nexus data file (same set of taxa as tree.nex)
* test_likelihood.R: script to run selac model
* selac.csv: csv file of selac model log likelihoods

Parameter settings:

* C.Phi.q.Ne = `4*4e-7*.5*5e6`
* alpha = `1.829272`
* beta = `0.101799`
* base frequencies = all equal
* gamma alpha = 5
