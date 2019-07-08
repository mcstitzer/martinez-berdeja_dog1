import sys


def readFasta(filename):
	fastafile=open(filename, 'r')
	fastadict={}
	for line in fastafile:		
		if line.startswith('>') and line.strip()[1:] not in fastadict:
                        seqname=line.strip()[1:].split(':')[0] ## added the split part, maybe bedtools getfasta is using a :: separator now?
			fastadict[seqname]=[]
#			continue
		else:
			fastadict[seqname].append(line.strip())
	for entry in fastadict:
		fastadict[entry]=''.join(fastadict[entry])
	return fastadict

def printFasta(sequence, width=70):
	return '\n'.join(sequence[i:i+width] for i in range(0,len(sequence), width))

b=readFasta(sys.argv[1])
#out=open('stupid_test.txt', 'w')
#out.write('read in all the fasta entries')

t=readFasta(sys.argv[2])

for seq in t:
	t[seq]=t[seq][12:16]
#print t

uniqueFasta={}
uniqueFastaHaplo={}
for entry in b:
	if b[entry] not in uniqueFasta:
		uniqueFasta[b[entry]]=[entry]
		uniqueFastaHaplo[b[entry]]=t[entry]
	if b[entry] in uniqueFasta:
		uniqueFasta[b[entry]].append(entry)



i=1
for seq in uniqueFasta:
	uniqueFasta[seq]='.'.join(uniqueFasta[seq])
	print '>'+uniqueFastaHaplo[seq]+'_'+str(i)+ '\n'+printFasta(''.join(seq))
	i+=1

