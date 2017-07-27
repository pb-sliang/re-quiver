#!/usr/bin/env python
'''
script for 'due-quiver' so as to increase quiver coverage of a unzipped contigs
'''
import os
import sys
import glob
import subprocess
import tempfile
from functools import partial
import logging
import argparse
import collections
import multiprocessing as mp
# to previous version change here
#import merge_fq
import merge_fq_blasr



logger = logging.getLogger(__name__)


def get_hcontigs():
    pcontigs=[pctg.replace('_ref.fa','') for pctg in glob.glob('0*F_ref.fa')]
    hcontigs={}
    for pctg in pcontigs:
        hcontigs[pctg]=[hctg.replace('_ref.fa','') for hctg in glob.glob(pctg+'_???_ref.fa')]
    return hcontigs

    
def get_parser():
    """Return an argparse instance"""
    parser = argparse.ArgumentParser(description='')
    #parser.add_argument('fastafile', action='store', help="Path to fasta file to filter.")
    parser.add_argument('dir', nargs='?', default=os.getcwd(),
                        help="the name of the directory containing contig files; if not supplied current dir is used")
    parser.add_argument('-q', '--queue', dest='queue', action='store', default='default',
                        help="the queue to run in")
    parser.add_argument('--nproc', dest='nproc', action='store', default=16,
                        help="number of cpus required")
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help="Turn on verbose logging mode")
    return parser


def qsub_cmd(cmd, queue, wdir, slots=16):
    """
    Generate a qsub command for each task
    :param cmd: command to run
    :param queue: what queue to run on
    :return:
    """

    with tempfile.NamedTemporaryFile(dir=os.getcwd()) as bashscript:
        bashscript.write(cmd)
        bashscript.flush()

        cmd = "qsub -sync y -q {q} -pe smp {s} -cwd -j y -V -o {l} -e {l} " \
              "-S /bin/bash {x}".format(q=queue,
                                     s=slots,
                                     l=wdir+'/sge_log',
                                     x=bashscript.name)
        logger.info("Submitting: %s", cmd)
        subprocess.call(cmd, shell=True)

    return


def run_cmds(commands, queue, wdir, slots=16, max_processes=48):
    """
    :param commands: list of commands to run
    :param slots: number of cores needed per task
    :param max_processes: Total number of processes allowed
    :param queue: what SGE queue to run on - default is 'default'
    :return
    """
    logger.info("Running %s grid jobs", commands)

    if not os.path.exists(wdir+'/scratch'):
        os.mkdir(wdir+'/scratch')
    if not os.path.exists(wdir+'/sge_log'):
        os.mkdir(wdir+'/sge_log')
    partial_qsub = partial(qsub_cmd, wdir=wdir, queue=queue)
    pool = mp.Pool(max_processes)
    pool.map(partial_qsub, commands)
    pool.close()
    pool.join()
    return


def quiver_cmd(hcontigs, wdir):
    commands=[]
    for pctg in hcontigs:
        cmd='''
source /mnt/software/Modules/current/init/bash
module load smrtanalysis/mainline
module load samtools

set -vex

samtools view -H {wdir}/{pctg}.bam > {wdir}/scratch/{pctg}.sam
samtools view {wdir}/{pctg}.bam >> {wdir}/scratch/{pctg}.sam

for BAM in {wdir}/{pctg}_???.bam
do
samtools view $BAM >> {wdir}/scratch/{pctg}.sam
done

samtools view -hb {wdir}/scratch/{pctg}.sam > {wdir}/scratch/{pctg}.bam
#rm {wdir}/scratch/{pctg}.sam

pbalign --tmpDir={wdir}/scratch/ --nproc={nproc} --minAccuracy=0.75 --minLength=50 --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr          --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=randombest --seed 0 {wdir}/scratch/{pctg}.bam {wdir}/{pctg}_ref.fa {wdir}/scratch/aln-ph-{pctg}.bam

samtools faidx {wdir}/{pctg}_ref.fa
variantCaller --algorithm=arrow -x 5 -X 120 -q 20 -j {nproc} -r {wdir}/{pctg}_ref.fa {wdir}/scratch/aln-ph-{pctg}.bam  -o {wdir}/cns-ph-{pctg}.fasta -o {wdir}/cns-ph-{pctg}.fastq.gz
gunzip cns-{pctg}.fasta.gz
blasr --minMatch 17 --fastMaxInterval --nproc {nproc} -m 5 --advanceExactMatches 10 cns-{pctg}.fasta cns-ph-{pctg}.fasta > fa_aln_{pctg}.m5
gzip cns-{pctg}.fasta cns-ph-{pctg}.fasta

#rm /localdisk/scratch/{pctg}.bam
#rm /localdisk/scratch/aln-ph-{pctg}.bam*
'''.format(pctg=pctg,wdir=wdir, nproc=nproc)
        commands.append(cmd)

    for pctg in hcontigs:
        for hctg in hcontigs[pctg]:
            cmd='''
source /mnt/software/Modules/current/init/bash
module load smrtanalysis/mainline
module load samtools

set -vex

samtools view -H {wdir}/{hctg}.bam > {wdir}/scratch/{hctg}.sam
samtools view {wdir}/{hctg}.bam >> {wdir}/scratch/{hctg}.sam
samtools view {wdir}/{pctg}.bam >> {wdir}/scratch/{hctg}.sam
samtools view -hb {wdir}/scratch/{hctg}.sam > {wdir}/scratch/{hctg}.bam
#rm {wdir}/scratch/{hctg}.sam

pbalign --tmpDir={wdir}/scratch/ --nproc={nproc} --minAccuracy=0.75 --minLength=50 --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr          --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=randombest --seed 0 {wdir}/scratch/{hctg}.bam {wdir}/{hctg}_ref.fa {wdir}/scratch/aln-ph-{hctg}.bam

samtools faidx {wdir}/{hctg}_ref.fa
variantCaller --algorithm=arrow -x 5 -X 120 -q 20 -j {nproc} -r {hctg}_ref.fa {wdir}/scratch/aln-ph-{hctg}.bam  -o {wdir}/cns-ph-{hctg}.fasta -o {wdir}/cns-ph-{hctg}.fastq.gz
gunzip cns-{hctg}.fasta.gz
blasr --minMatch 17 --fastMaxInterval --nproc {nproc} -m 5 --advanceExactMatches 10 cns-{hctg}.fasta cns-ph-{hctg}.fasta > fa_aln_{hctg}.m5
gzip cns-{hctg}.fasta cns-ph-{hctg}.fasta


#rm /localdisk/scratch/{hctg}.bam
#rm /localdisk/scratch/aln-ph-{hctg}.bam*
'''.format(hctg=hctg, pctg=pctg, wdir=wdir, nproc=nproc)
            commands.append(cmd)
    return commands


if __name__=='__main__':
    parser = get_parser()
    args=parser.parse_args()
    wdir=os.path.abspath(args.dir)
    queue = args.queue
    nproc=args.nproc
    if args.debug:
        logging.basicConfig(filename='log.txt',level=logging.DEBUG)
    else:
        logging.basicConfig(filename='log.txt',level=logging.INFO)


    #print 'args', args, 'local', local, 'queue', queue

    hcontigs = get_hcontigs()
    commands=quiver_cmd(hcontigs, wdir)
    run_cmds(commands, queue, wdir)

    contigs=[]
    for pctg in hcontigs:
        contigs.append(pctg)
    for pctg in hcontigs:
        contigs.extend(hcontigs[pctg])
    pool=mp.Pool()
    # to previous version change here
    #pool.map( merge_fq.merge_fq, contigs)
    pool.map( merge_fq_blasr.merge_fq, contigs)
    pool.close()
    pool.join

