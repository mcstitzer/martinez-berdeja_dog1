#!/bin/bash -login

module load python3
module load biopython

## paths to aspera, sra-toolkit are hardcoded
TARGET=dog1_mrna.fa
SAMPLENAME=Cvi
READ1=
READ2=
#SRRNUM=SRR3156163
#SRRNUM=SRR492239
SRRNUM=SRR1756969

## if there is a SRR number provided, download and reassign read variables
if [ ! -z "$SRRNUM" ]
then

SAMPLENAME=$SRRNUM
PATH1=${SRRNUM:0:3}
PATH2=${SRRNUM:0:6}

echo "anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$PATH1/$PATH2/$SRRNUM/${SRRNUM}.sra"
mkdir -p /scratch/mstitzer/${SRRNUM}
if [ ! -f /scratch/mstitzer/${SRRNUM}/$SRRNUM.sra ]
then
# get file from sra
## -k2 from http://betascience.blogspot.com/2010/02/using-aspera-instead-of-ftp-to-download.html because getting timeouts
~/.aspera/connect/bin/ascp -i /home/mstitzer/.aspera/connect/etc/asperaweb_id_dsa.putty -k 1 -QTr -l 200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$PATH1/$PATH2/$SRRNUM/${SRRNUM}.sra /scratch/mstitzer/${SRRNUM}/${SRRNUM}.sra

fi

# convert to fq
~/software/sratoolkit.2.3.4-2-ubuntu64/bin/fastq-dump -O /scratch/mstitzer/${SRRNUM}/ --split-3 /scratch/mstitzer/${SRRNUM}/${SRRNUM}.sra


READ1=/scratch/mstitzer/${SRRNUM}/${SRRNUM}_1.fastq
READ2=/scratch/mstitzer/${SRRNUM}/${SRRNUM}_2.fastq

rm /scratch/mstitzer/${SRRNUM}/${SRRNUM}.sra

fi

#######################
### move on and map ###
#######################
module load bamtools

## generate index if needed
if [ ! -f $TARGET.bwt ]
then
bwa index $TARGET
fi


#### start mapping
if [ -f $READ2 ]  ## we have paired end sequencing
then

bwa mem -t 4 $TARGET $READ1 $READ2 | samtools view -Su - | bamtools filter -script filter_mapped_and_pairs.bamtools.json | samtools sort - $SAMPLENAME

rm ${SRRNUM}_2.fastq
fi

if [ ! -f $READ2 ]  ## single end sequencing
then
#bwa mem -t 16 kindr_kin11.CDS.fa $READ1 | samtools view -Su -F 4 - | samtools sort - $SRRNUM
### only need to filter for mapped reads when there is only single end data
bwa mem -t 4 $TARGET $READ1 | samtools view -Su -F 4 - | samtools sort - $SAMPLENAME
fi


## get rid of these files
rm ${SRRNUM}_1.fastq



#############################
### assemble mapped reads ###
#############################

bedtools bamtofastq -i ${SAMPLENAME}.bam -fq ${SAMPLENAME}.fq
#/share/apps/qiime-1.9.1/bin/convert_fastaqual_fastq.py -f ${SAMPLENAME}.fq -c fastq_to_fastaqual
#mv ${SAMPLENAME}.qual ${SAMPLENAME}.fna.qual
python3 fq_to_fa_qual.py ${SAMPLENAME}.fq ${SAMPLENAME}.fna ${SAMPLENAME}.fna.qual
~/software/phrap_1.090518/phrap ${SAMPLENAME}.fna -vector_bound 0 -forcelevel 0 -minscore 10 -new_ace 2> ${SAMPLENAME}.phrap.stderr > ${SAMPLENAME}.phrap.stdout

mv ${SAMPLENAME}* assembled_out_tibet
