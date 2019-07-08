
a=read.table('dog1.may22_e1.extended.aa.txt', header=F, stringsAsFactor=F)

a$V3=substr(a$V2,13,16)
names(a)=c('ecotype', 'extendedseq', 'haplotype')
data.frame(rename(count(a, extendedseq, haplotype), Freq=n))
a=a[a$haplotype!='----',]
a=a[a$haplotype!='E---',]
a=a[a$haplotype!='---X',]
a=a[a$haplotype!='---Y',]
a=a[a$haplotype!='-XCY',]
a=a[a$haplotype!='-XSY',]
a=a[a$haplotype!='XFSY',]
a=a[a$haplotype!='X-SY',]

## get rid of messy contigs
#a=a[!grepl(pattern='---', x=substr(a$extendedseq, 16,nchar(a$extendedseq[1]))),]
a=a[!grepl(pattern='--', x=a$extendedseq),]

## output - this could be improved to be more meaningful!
library(stringr)

ae=str_split_fixed(a$extendedseq, pattern='', nchar(a$extendedseq[1]))
## here, 750 is quite arbitrary - based on 0.5% frequency in 1001 genomes (of minor allele) but who knows??
aa=cbind(a, ae[,sapply(1:ncol(ae),function(x) !max(table(ae[,x]))>750)])

colnames(aa)[4:17]=paste0('aa_', which(sapply(1:ncol(ae),function(x) !max(table(ae[,x]))>750)))
write.table(aa, 'dog1.may22.extended.uniq.txt', col.names=T, row.names=F, sep='\t', quote=F)
