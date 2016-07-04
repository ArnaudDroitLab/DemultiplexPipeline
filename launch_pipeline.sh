#!/bin/bash

# A POSIX variable
OPTIND=1         # Reset in case getopts has been used previously in the shell.

# Initialize variables with default values
output="output"                 # Directory where output should be stored.
samplesheet="SampleSheet.csv"   # Path to the samplesheet file.
base_path=""                    # Base path for the sequencer produced files.
lane='.*'                       # Regular expression for processing only a subset of the lanes or tiles.

# Parse parameters from the command line.
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

# STEP 1 DEMULTIPLEXING #######################################################

# Name and create demultiplexing input and output directories.
mkdir -p $output
current_date=$(date +%Y-%m-%d)
dir_input=$base_path/Data/Intensities/BaseCalls
dir_output=$output/Demultiplexed_"${current_date}"

base_mask="Y*,I6n*,Y*"

# Generate the make files which will do the demultiplexing.
/home/nanuq-admin/nanuq-programs/software/bcl2fastq-1.8.4/bin/configureBclToFastq.pl \
--input-dir ${dir_input} \
--output-dir ${dir_output} \
--sample-sheet ${samplesheet} \
--use-bases-mask ${base_mask} \
--mismatches 1 \
--tiles 's_'$lane'_.*' \
--fastq-cluster-count 0

# Preform demultiplexing.
pushd ${dir_output}
make -j 8
popd# STEP 2 TRIMMING #############################################################

mkdir $output/trim
FILE_ADAPTORS=adaptors.fa

# Trim all demultiplexed files in parrallel, then wait on the jobs to complete.
for file in $(ls $output/Demultiplexed_*/Project_*/*/*_R1_001.fastq.gz)
do
	sampleID=$(basename $file _R1_001.fastq.gz)

	java -XX:ParallelGCThreads=1 -jar /home/nanuq-admin/nanuq-programs/software/trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar \
	PE -phred33 -threads 8 \
	$output/Demultiplexed_*/Project_*/*/${sampleID}_R1_001.fastq.gz \
    $output/Demultiplexed_*/Project_*/*/${sampleID}_R2_001.fastq.gz \
	$output/trim/${sampleID}_trim_R1.fastq.gz $output/trim/${sampleID}_trim_R1_unpaired.fastq.gz \
	$output/trim/${sampleID}_trim_R2.fastq.gz $output/trim/${sampleID}_trim_R2_unpaired.fastq.gz \
	ILLUMINACLIP:${FILE_ADAPTORS}:2:30:10 \
	TRAILING:30 MINLEN:36 &
done
wait

# STEP 3 ALIGNMENT ############################################################

mkdir -p $output/STAR/Log
mkdir -p $output/STAR/Within
mkdir -p $output/STAR/Unsorted

# STAR 1st pass
for file in $(ls $output/trim/*_trim_R1.fastq.gz)
do
	sampleID=$(basename $file _trim_R1.fastq.gz)

	mkdir -p $output/STAR/alignment_1stPass/${sampleID}/

	/software/STAR/STAR-2.5.2a/bin/Linux_x86_64/STAR --runMode alignReads \
		--genomeDir /is1/commonDatasets/STAR/Homo_sapiens.GRCh38/genome/star_index/Ensembl77.sjdbOverhang99 \
		--readFilesIn \
			$output/trim/${sampleID}_trim_R1.fastq.gz \
			$output/trim/${sampleID}_trim_R2.fastq.gz \
		--runThreadN 32 \
		--readFilesCommand zcat \
		--outStd Log \
		--outSAMunmapped Within \
		--outSAMtype BAM Unsorted \
		--outFileNamePrefix $output/STAR/alignment_1stPass/${sampleID}/ \
		--outSAMattrRGline ID:"${sampleID}"     PL:"ILLUMINA"             SM:"${sampleID}"     CN:"CRCHUQ"  \
		--limitGenomeGenerateRAM 70000000000 \
		--limitIObufferSize 1000000000 
done

# Create new splice junction database containing the splice junctions of all samples
cat $output/STAR/alignment_1stPass/*/SJ.out.tab | \
awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-"} {if($5>0){print $1,$2,$3,strChar[$4]}}' | sort -k1,1h -k2,2n > $output/STAR/alignment_1stPass/AllSamples.SJ.out.tab && \
mkdir -p $output/STAR/reference.Merged && \
/software/STAR/STAR-2.5.2a/bin/Linux_x86_64/STAR --runMode genomeGenerate \
	--genomeDir $output/STAR/reference.Merged \
	--genomeFastaFiles /is1/commonDatasets/STAR/Homo_sapiens.GRCh38/genome/Homo_sapiens.GRCh38.fa \
	--runThreadN 32 \
	--limitGenomeGenerateRAM 50000000000 \
	--sjdbFileChrStartEnd $output/STAR/alignment_1stPass/AllSamples.SJ.out.tab \
	--limitIObufferSize 1000000000 \
	--sjdbOverhang 99


# STAR 2nd pass.
for file in $(ls $output/trim/*_trim_R1.fastq.gz)
do
    sampleID=$(basename $file _trim_R1.fastq.gz)
	mkdir -p $output/alignment/${sampleID}

	/software/STAR/STAR-2.5.2a/bin/Linux_x86_64/STAR --runMode alignReads \
		--genomeDir $output/STAR/reference.Merged \
		--readFilesIn \
			$output/trim/${sampleID}_trim_R1.fastq.gz \
			$output/trim/${sampleID}_trim_R2.fastq.gz \
		--runThreadN 32 \
		--readFilesCommand zcat \
		--outStd Log \
		--outSAMunmapped Within \
		--outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix $output/alignment/${sampleID}/ \
		--limitGenomeGenerateRAM 70000000000 \
		--limitBAMsortRAM 70000000000 \
		--limitIObufferSize 1000000000 \
		--outWigType wiggle read1_5p --outWigStrand Stranded --outWigReferencesPrefix chr 
done

# Perform QC ##################################################################
module add fastqc/0.11.2
mkdir $output/fastqc
for file in $(ls $output/trim/*_R?.fastq.gz)
do   
    sequenceFile=`basename $file .fastq.gz`
    mkdir -p $output/fastqc/$sequenceFile
    fastqc $file -o $output/fastqc/$sequenceFile &
done
wait

# Count reads ##################################################################
module add python/2.7.8
mkdir $output/counts
for bamFile in $output/alignment/*/Aligned.sortedByCoord.out.bam
do
    sampleDir=`dirname $bamFile`
    sampleID=`basename sampleDir`

    /home/foueri01/.local/bin/htseq-count -s reverse --order pos -f bam $bamFile /is1/commonDatasets/STAR/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl77.gtf > $output/counts/$sampleID.counts &
done

Rscript summary_table.R $output