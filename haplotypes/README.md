

### Assign ecotype IDs to alignment

- `dog1_phylip_to_tabout.py` switches SRR names to ecotype IDs from a phylip file
  - input: `dog1.may22.phy` (phylip file of exon 1 alignment) and `SraRunTable_1001genomes.txt` (info downloaded from NCBI SRA matching SRR numbers to ecotype IDs)
  - writes to stdout, I used this to generate `dog1.may22.uniq.txt`, and by changing the final line, to generate `dog1.may22.extended.uniq.txt`, which differ in whether they have the 

### Find additional noncoding variants in exon 1 of dog1

- `filter_exon1_extended.R` finds sub-variants beyond that of the functional haplotype of Nakabayashi et al., 2015
  - input: `dog1.may22_e1.extended.aa.txt`, a tab delimited file with column 1 as the ecotype ID, and column 2 as the (aligned) E1 amino acid sequence
  - output:  `dog1.may22.extended.uniq.txt`, a tab delimited file with all aa's at >0.5% frequency added as a new column
