#!/bin/bash

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
output="output"
samplesheet="SampleSheet.csv"
base_path=""
lane='.*'
while getopts "i:o:s:l:" opt; do
    case "$opt" in
    i)  base_path=$OPTARG
        ;;
    s)  samplesheet=$OPTARG
        ;;
    o)  output=$OPTARG
        ;;
    l)  lane=$OPTARG
    esac
done

shift $((OPTIND-1))

echo $base_path
echo $samplesheet
echo $output
workdir=$output
#exit # haha securité

# STEP 1 DEMULTIPLEXING

#base_path=/u01/genomique/Instruments/Illumina_HiSeq_Heaven/Raw_data/160620_D00487S_0119_AC8KJMANXX
#workdir=/home/mvallee/Prostate_weekend

mkdir -p $workdir
current_date=$(date +%Y-%m-%d)
#samplesheet=$workdir/SampleSheet.csv
dir_input=$base_path/Data/Intensities/BaseCalls
dir_output=$workdir/Demultiplexed_"${current_date}"

base_mask="Y*,I6n*,Y*"

/home/nanuq-admin/nanuq-programs/software/bcl2fastq-1.8.4/bin/configureBclToFastq.pl \
--input-dir ${dir_input} \
--output-dir ${dir_output} \
--sample-sheet ${samplesheet} \
--use-bases-mask ${base_mask} \
--mismatches 1 \
--tiles 's_'$lane'_.*' \
--fastq-cluster-count 0

pushd ${dir_output}

make -j 8
popd


# STEP 2 TRIMMING

mkdir $workdir/trim
FILE_ADAPTORS=/home/mvallee/Prostate_weekend/adaptors.fa

for file in $(ls $workdir/Demultiplexed_*/Project_*/*/*_R1_001.fastq.gz)
do
	sampleID=$(basename $file _R1_001.fastq.gz)

	java -XX:ParallelGCThreads=1 -jar /home/nanuq-admin/nanuq-programs/software/trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar \
	PE -phred33 -threads 8 \
	$workdir/Demultiplexed_*/Project_*/*/${sampleID}_R1_001.fastq.gz \
    $workdir/Demultiplexed_*/Project_*/*/${sampleID}_R2_001.fastq.gz \
	$workdir/trim/${sampleID}_trim_R1.fastq.gz $workdir/trim/${sampleID}_trim_R1_unpaired.fastq.gz \
	$workdir/trim/${sampleID}_trim_R2.fastq.gz $workdir/trim/${sampleID}_trim_R2_unpaired.fastq.gz \
	ILLUMINACLIP:${FILE_ADAPTORS}:2:30:10 \
	TRAILING:30 MINLEN:36 &
done
wait


# SE TEST
# for file in $(ls $workdir/Demultiplexed_*/Project_*/*/*_R1_001.fastq.gz)
# do
# 	sampleID=$(basename $file _R1_001.fastq.gz)
# 
# 	java -XX:ParallelGCThreads=1 -jar /home/nanuq-admin/nanuq-programs/software/trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar \
# 	SE -phred33 -threads 8 \
# 	$workdir/Demultiplexed_*/Project_*/*/${sampleID}_R1_001.fastq.gz \
# 	$workdir/trim/${sampleID}_trim_R1.fastq.gz \
# 	ILLUMINACLIP:${FILE_ADAPTORS}:2:30:10 \
# 	TRAILING:30 MINLEN:36 &
# done
# wait

# STEP 3 ALIGNING

mkdir $workdir/STAR
#cd $workdir/STAR
mkdir $workdir/STAR/Log
mkdir $workdir/STAR/Within
mkdir $workdir/STAR/Unsorted

	# 1st pass

for file in $(ls $workdir/trim/*_trim_R1.fastq.gz)
do
	sampleID=$(basename $file _trim_R1.fastq.gz)

	mkdir -p $workdir/STAR/alignment_1stPass/${sampleID}/

	/software/STAR/STAR-2.5.2a/bin/Linux_x86_64/STAR --runMode alignReads \
		--genomeDir /is1/commonDatasets/STAR/Homo_sapiens.GRCh38/genome/star_index/Ensembl77.sjdbOverhang99 \
		--readFilesIn \
			$workdir/trim/${sampleID}_trim_R1.fastq.gz \
			$workdir/trim/${sampleID}_trim_R2.fastq.gz \
		--runThreadN 32 \
		--readFilesCommand zcat \
		--outStd Log \
		--outSAMunmapped Within \
		--outSAMtype BAM Unsorted \
		--outFileNamePrefix $workdir/STAR/alignment_1stPass/${sampleID}/ \
		--outSAMattrRGline ID:"${sampleID}"     PL:"ILLUMINA"             SM:"${sampleID}"     CN:"CRCHUQ"  \
		--limitGenomeGenerateRAM 70000000000 \
		--limitIObufferSize 1000000000 
done

	# Create new splice junction database containing the splice junctions of all samples

cat $workdir/STAR/alignment_1stPass/*/SJ.out.tab | \
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-"} {if($5>0){print $1,$2,$3,strChar[$4]}}' | sort -k1,1h -k2,2n > $workdir/STAR/alignment_1stPass/AllSamples.SJ.out.tab && \
mkdir -p $workdir/STAR/reference.Merged && \
/software/STAR/STAR-2.5.2a/bin/Linux_x86_64/STAR --runMode genomeGenerate \
	--genomeDir $workdir/STAR/reference.Merged \
	--genomeFastaFiles /is1/commonDatasets/STAR/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
	--runThreadN 32 \
	--limitGenomeGenerateRAM 50000000000 \
	--sjdbFileChrStartEnd $workdir/STAR/alignment_1stPass/AllSamples.SJ.out.tab \
	--limitIObufferSize 1000000000 \
	--sjdbOverhang 99


	# 2nd pass

for file in $(ls $workdir/trim/*_trim_R1.fastq.gz)
do
    sampleID=$(basename $file _trim_R1.fastq.gz)
	mkdir -p $workdir/alignment/${sampleID}

	/software/STAR/STAR-2.5.2a/bin/Linux_x86_64/STAR --runMode alignReads \
		--genomeDir $workdir/STAR/reference.Merged \
		--readFilesIn \
			$workdir/trim/${sampleID}_trim_R1.fastq.gz \
			$workdir/trim/${sampleID}_trim_R2.fastq.gz \
		--runThreadN 16 \
		--readFilesCommand zcat \
		--outStd Log \
		--outSAMunmapped Within \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix $workdir/alignment/${sampleID}/ \
		--limitGenomeGenerateRAM 70000000000 \
		--limitBAMsortRAM 70000000000 \
		--limitIObufferSize 1000000000 \
		--outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr 
done

# Counting reads

