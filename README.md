# re-quiver
When an assembly is unzipped and a consensus is found by quiver using only allele-specific reads could be less than optimal because the vast regions outside of SNP where two alleles matches including the other alternative allele could improve quiver accuracy. The program requiver.py combine two quiver results. One is allele specific quiver produced by FALCON, another is done by requiver.py using all reads mapped to the FALCON consensus. Two quiver fastq files are aligned and where they differ, requiver.py choose the path with higher quiver QV scores.

The program assumes the existence of consensus (ctgName_ref.fa, where ctgName is the name of the h- or p- contig), allele-specific quiver (cns-ctgName.fastq.gz), allele specific reads (ctgName.bam). It produce a quiver from all reads (cns-ph-ctgName.fastq.gz), and a merging of cns-ctgName.fastq.gz and cns-ph-ctgName.fastq.gz: cns-m-ctgName.fastq.gz.

Dependencies:
pbcore.io statisfied by
module load smrtanalysis/mainline

Usage:
$ ./requiver.py

usage: requiver.py [-h] [--queue QUEUE] [--debug] [dir]

positional arguments:
  dir            the name of the directory containing contig files; if not
                 supplied current dir is used

optional arguments:
  -h, --help     show this help message and exit
  --queue QUEUE  the queue to run in
  --debug        Turn on verbose logging mode

