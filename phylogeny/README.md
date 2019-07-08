
### Generate fasta with each tip only once

 `python beast_singleton.py revbayes_e1dog1.nogaps.fa revbayes_e1dog1.nogaps.translated.fasta > pruned_e1_aln.fa`

### These were used to generate a BEAST xml file in BEAUTI

The xml file is `strictClock_constantPop.pruned_e1_aln.xml` and was run for 10,000,000 steps.


## below are outdated - they generated a resampled xml file, and added haplotype identifiers to a tree

#### replace xml sequences with sample
`for i in {1..100}; do mkdir -p lognormalClock_muFixed_ConstantPop/rep$i; python random_sample_for_beast_make_xml.py revbayes_e1dog1.nogaps.fa revbayes_e1dog1.nogaps.translated.fasta lognormalClock_muFixed_ConstantPop.xml> lognormalClock_muFixed_ConstantPop/rep$i/rep$i.xml; done`



### Add haplotypes to tree

- `add_haplotypes_to_tree.py` writes a new newick tree to stdout
  - input: `lyrtree.txt`, a newick tree object
  - input:  `../haplotypes/dog1.may22.phy`, a tab delimited file with SRR number in column 1, and functional haplotype in second (this is a phylip file)
