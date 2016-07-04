#workdir=/home/mvallee/Prostate_weekend
samplesheet=$workdir/SampleSheet.csv

# Perform QC using FastQC
module add fastqc/0.11.2
mkdir $workdir/fastqc
for file in $(ls $workdir/trim/*_R?.fastq.gz)
do   
    sequenceFile=`basename $file .fastq.gz`
    mkdir -p $workdir/fastqc/$sequenceFile
    fastqc $file -o $workdir/fastqc/$sequenceFile &
done
wait

# Count reads
module add python/2.7.8
mkdir $workdir/counts
for sampleDir in $workdir/alignment/166*
do
    sampleID=`basename $sampleDir`
    bamFile=$sampleDir/Aligned.sortedByCoord.out.bam
    
    echo $sampleDir
    echo $sampleID
    echo $bamFile
    /home/foueri01/.local/bin/htseq-count -s no --order pos -f bam $bamFile /is1/commonDatasets/STAR/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl77.gtf > $workdir/counts/$sampleID.counts &
done



