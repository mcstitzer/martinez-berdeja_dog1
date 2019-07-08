#!/bin/bash -login
#SBATCH -D /home/mstitzer/projects/alejandra/dog1_haplo
#SBATCH -o /home/mstitzer/projects/alejandra/slurm-log/dog1-stdout-%j.txt
#SBATCH -e /home/mstitzer/projects/alejandra/slurm-log/dog1-stderr-%j.txt
#SBATCH -J dog1

module load python3
module load biopython

## paths to aspera, sra-toolkit are hardcoded
TARGET=dog1_mrna.fa
#LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" sra_1135.tosample.txt)
#SPLITLINE=(${LINE//\t/ })
#SAMPLENAME=SPLITLINE
READ1=
READ2=
#SRRNUM=SRR3156163
#SRRNUM=SRR492239
#SRRNUM=ERR2026804
#SRRNUM=ERR2026767

read SRRNUM <<< $(sed -n "${SLURM_ARRAY_TASK_ID}p" sra_yangtze.tosample.txt)
SAMPLENAME=$SRRNUM
#echo $SAMPLENAME
echo $SRRNUM
## if there is a SRR number provided, download and reassign read variables
if [ ! -z "$SRRNUM" ]
then

SAMPLENAME=$SRRNUM
PATH1=${SRRNUM:0:3}
PATH2=${SRRNUM:0:6}
fi

echo "anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$PATH1/$PATH2/$SRRNUM/${SRRNUM}.sra"
echo assembled_out_africa/${SRRNUM}.fna.contigs
if [ ! -f "assembled_out_africa/${SRRNUM}.fna.contigs" ]
then

mkdir -p /scratch/mstitzer/${SRRNUM}
if [ ! -f /scratch/mstitzer/${SRRNUM}/$SRRNUM.sra ]
then
# get file from sra
## -k2 from http://betascience.blogspot.com/2010/02/using-aspera-instead-of-ftp-to-download.html because getting timeouts
~/.aspera/connect/bin/ascp -i /home/mstitzer/.aspera/connect/etc/asperaweb_id_dsa.putty -k 1 -QTr -l 200m anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$PATH1/$PATH2/$SRRNUM/${SRRNUM}.sra /scratch/mstitzer/${SRRNUM}/${SRRNUM}.sra


# convert to fq - switched to use --split-3 because it was running into issues with uneven numbers of reads.
~/software/sratoolkit.2.3.4-2-ubuntu64/bin/fastq-dump -O /scratch/mstitzer/${SRRNUM}/ --split-3 /scratch/mstitzer/${SRRNUM}/${SRRNUM}.sra


READ1=/scratch/mstitzer/${SRRNUM}/${SRRNUM}_1.fastq
READ2=/scratch/mstitzer/${SRRNUM}/${SRRNUM}_2.fastq

rm /scratch/mstitzer/${SRRNUM}/${SRRNUM}.sra

## if there is a third unpaired file, remove it here.
if [ -f /scratch/mstitzer/${SRRNUM}/${SRRNUM}.fastq ]
then
 rm /scratch/mstitzer/${SRRNUM}/${SRRNUM}.fastq
fi

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

#rm ${SRRNUM}_2.fastq
rm ${READ2}
fi

if [ ! -f $READ2 ]  ## single end sequencing
then
#bwa mem -t 16 kindr_kin11.CDS.fa $READ1 | samtools view -Su -F 4 - | samtools sort - $SRRNUM
### only need to filter for mapped reads when there is only single end data
bwa mem -t 4 $TARGET $READ1 | samtools view -Su -F 4 - | samtools sort - $SAMPLENAME
fi


## get rid of these files
#rm ${SRRNUM}_1.fastq
rm ${READ1}


#############################
### assemble mapped reads ###
#############################

bedtools bamtofastq -i ${SAMPLENAME}.bam -fq ${SAMPLENAME}.fq
#/share/apps/qiime-1.9.1/bin/convert_fastaqual_fastq.py -f ${SAMPLENAME}.fq -c fastq_to_fastaqual
#mv ${SAMPLENAME}.qual ${SAMPLENAME}.fna.qual
python3 fq_to_fa_qual.py ${SAMPLENAME}.fq ${SAMPLENAME}.fna ${SAMPLENAME}.fna.qual
~/software/phrap_1.090518/phrap ${SAMPLENAME}.fna -vector_bound 0 -forcelevel 0 -minscore 10 -new_ace 2> ${SAMPLENAME}.phrap.stderr > ${SAMPLENAME}.phrap.stdout

mv ${SAMPLENAME}* assembled_out_yangtze

fi
