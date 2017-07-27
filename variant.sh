for bam in 000000F*.bam
do
ctg=${bam%.bam}
qsub -v ctg=$ctg variant.q
done

