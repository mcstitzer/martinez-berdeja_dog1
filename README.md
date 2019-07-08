# martinez-berdeja_dog1
Scripts to reconstruct DOG1 haplotypes and understand their evolution for Martinez-Berdeja et al., 2019

- `dog1_history.Rmd` generates Figure 3, with the gene tree of DOG1 and pairwise divergence of haplotype groups plotted on a temperature record of the Pleistocene.
  - `dog1_history.md` shows code as executed.

- `haplotypes` contains scripts to 
  - download short reads, map, and reassemble exon 1 of DOG1
  - align exon 1 sequence from all individuals
  
- `divergence_dating` contains script to calculate pi for groups. (This is no longer used, as this approach is incorporated into `dog1_history.Rmd`)

- `phylogeny` contains BEAST xml files used to generate the gene tree of DOG1

#### Warning! Paths to files are hard coded for my computer, and likely need to be updated if trying to reproduce
