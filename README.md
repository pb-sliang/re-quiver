# re-quiver
When an assembly is unzipped and a consensus is found by quiver using only allele-specific reads could be less than optimal because the vast regions outside of SNP where two alleles matches including the other alternative allele could improve quiver accuracy. The program requiver.py combine two quiver results. One is allele specific quiver produced by FALCON, another is done by requiver.py using all reads mapped to the FALCON consensus. Two quiver fastq files are aligned and where they differ, requiver.py choose the path with higher quiver QV scores.


Dependencies:
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

