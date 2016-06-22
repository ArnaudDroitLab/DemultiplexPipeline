# Count reads
workdir=/home/mvallee/Prostate_weekend
samplesheet=$workdir/SampleSheet.csv
mkdir $workdir/counts
for file in $(ls $workdir/STAR/alignment/*.bam)
do
	mkdir -p alignment/${sampleID}
    sampleID=`basename $file .bam`
    
    htseq-count -f bam $bamfile Homo_sapiens.GRCh38.Ensembl77.gtf > output/counts/$sampleID.counts &
done
wait


