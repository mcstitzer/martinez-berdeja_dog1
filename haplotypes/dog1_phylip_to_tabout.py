file='SraRunTable_1001genomes.txt'
f=open(file, 'r')
srrdict={}
for line in f:
    fields=line.strip().split('\t')
    srrdict[fields[11]]=fields[5]

ph=open('dog1.may22.phy', 'r')

for line in ph:
	fields=line.strip().split()
	psn=fields[0].find('SRR')
	if psn != -1:
		print srrdict[fields[0][psn:psn+10]]+'\t'+fields[1][0:4] # to get just the nakabayashi haplotype
### use this one for the full sequence
#		print srrdict[fields[0][psn:psn+10]]+'\t'+fields[1] # to get entire sequence


