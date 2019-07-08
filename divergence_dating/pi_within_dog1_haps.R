library(pegas)
library(plyr)
library(ape)

#a=read.FASTA('~/Downloads/revbayes_e1dog1.fa')
a=read.FASTA('~/Documents/berdeja_dog1/')

haplo=read.table('~/Documents/berdeja_dog1/haplotypes/dog1.may22.phy')
haplo$hap=substr(haplo$V2, 13,16)

adm=read.csv('~/Documents/berdeja_dog1/phylogeny/data_files/1001genomes-accessions and 1001genomes-admixture.csv')
eid=read.table('~/Documents/berdeja_dog1/haplotypes/eid_to_srrcontig.txt')
eid$admixture=mapvalues(eid$V1, from=adm$id, to=as.character(adm$group))
eid$haplo=mapvalues(eid$V3, from=haplo$V1, to=as.character(haplo$hap))

nuc.div(a, pairwise.deletion=T)
GC.content(a)
base.freq(a)
seg.sites(a)
theta.s(a)
b=dist.dna(a)


## pi of everybody
sum(dist.dna(a, 'RAW', pairwise.deletion=T), na.rm=T)/(1127*(1126)/2)

sum(dist.dna(a[names(a) %in% haplo$V1[haplo$hap=='D-SY']], 'RAW', pairwise.deletion=T), na.rm=T)/*(1126)/2)

nuc.div(a[names(a) %in% haplo$V1[haplo$hap=='D-SY']], pairwise.deletion=T)
nuc.div(a[names(a) %in% haplo$V1[haplo$hap=='D-RY']], pairwise.deletion=T)
nuc.div(a[names(a) %in% haplo$V1[haplo$hap=='ECCY']], pairwise.deletion=T)
nuc.div(a[names(a) %in% haplo$V1[haplo$hap=='ECSY']], pairwise.deletion=T)
nuc.div(a[names(a) %in% haplo$V1[haplo$hap=='EFSY']], pairwise.deletion=T)

nuc.div(a[names(a) %in% haplo$V1[haplo$hap%in%c('ECCY', 'ECSY', 'EFSY')]], pairwise.deletion=T)

mu=7.1e-9
mu.ci=c(6.3e-9,7.9e-9)
variablemut=9e-9
moi.mu.ci=c(2.4e-9,3e-9)
ne=250000
nuc.div(a[names(a) %in% haplo$V1[haplo$hap=='D-SY']], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)
nuc.div(a[names(a) %in% haplo$V1[haplo$hap=='D-RY']], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)
nuc.div(a[names(a) %in% haplo$V1[haplo$hap%in%c('D-SY', 'D-RY')]], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)

nuc.div(a[names(a) %in% haplo$V1[haplo$hap%in%c('ECCY', 'ECSY', 'EFSY')]], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)

nuc.div(a[names(a) %in% haplo$V1[haplo$hap%in%c('ECCY')]], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)
nuc.div(a[names(a) %in% haplo$V1[haplo$hap%in%c('ECSY')]], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)
nuc.div(a[names(a) %in% haplo$V1[haplo$hap%in%c('EFSY')]], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)



nuc.div(a[names(a) %in% haplo$V1[haplo$hap%in%c('ECCY')] & names(a) %in% eid$V3[eid$admixture=='south_sweden']], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)
nuc.div(a[names(a) %in% haplo$V1[haplo$hap%in%c('ECCY')] & names(a) %in% eid$V3[eid$admixture=='north_sweden']], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)



nuc.div(a[names(a) %in% haplo$V1[haplo$hap%in%c('ECCY')] & names(a) %in% eid$V3[eid$admixture=='north_sweden']], pairwise.deletion=T)/2/c(mu, mu.ci,variablemut)



