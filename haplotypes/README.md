### Map reads to DOG1 mRNA

- `map_to_phrap.array.sh` 
  - downloads SRA and converts to fastq files based on SLURM array ID from `sra_1135.tosample.txt`
  - maps to `dog1_mRNA.fa` and keeps only mapped reads and their pairs
  - assembles these reads using phrap

## concatenate all DOG1 assemblies

`cat all_dog1_assemblies.e1.may21.fa assembled_out_africa/*.contigs > all_dog1_1001_and_africa.fa`

## align all with MAFFT, arrange relative to exon 1 sequence

`srun -p bigmemh --mem=16000 --time=1-:00:00 mafft --adjustdirection --localpair --maxiterate 1000 --add all_dog1_assemblies.may21.fa --reorder dog1_e1.fa > all_dog1_assemblies.e1may21.aln.fa`


## prepare for tree building

manually remove unaligned seqs (bottom ones) in aliview, those with missing data across exon 1

also remember to remove empty seqs!

## remove duplicate sequences (where are these introduced?)
`awk '/^>/{f=!d[$1];d[$1]=1}f' all_dog1_1001_and_africa.aln.trimmed.fa > all_dog1_1001_and_africa.aln.trimmed.uniq.fa`

### Assign ecotype IDs to alignment

- `dog1_phylip_to_tabout.py` switches SRR names to ecotype IDs from a phylip file
  - input: `dog1.may22.phy` (phylip file of exon 1 alignment) and `SraRunTable_1001genomes.txt` (info downloaded from NCBI SRA matching SRR numbers to ecotype IDs)
  - writes to stdout, I used this to generate `dog1.may22.uniq.txt`, and by changing the final line, to generate `dog1.may22.extended.uniq.txt`, which differ in whether they have the 

### Find additional noncoding variants in exon 1 of dog1

- `filter_exon1_extended.R` finds sub-variants beyond that of the functional haplotype of Nakabayashi et al., 2015
  - input: `dog1.may22_e1.extended.aa.txt`, a tab delimited file with column 1 as the ecotype ID, and column 2 as the (aligned) E1 amino acid sequence
  - output:  `dog1.may22.extended.uniq.txt`, a tab delimited file with all aa's at >0.5% frequency added as a new column
