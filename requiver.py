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
import merge_fq


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
    parser.add_argument('--queue', dest='queue', action='store', default='default',
                        help="the queue to run in")
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

    if not os.path.exists(wdir+'/sge_log'):
        os.mkdir(wdir+'/sge_log')
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

samtools view -H {wdir}/{pctg}.bam > /localdisk/scratch/{pctg}.sam
samtools view {wdir}/{pctg}.bam >> /localdisk/scratch/{pctg}.sam

for BAM in {wdir}/{pctg}_???.bam
do
samtools view $BAM >> /localdisk/scratch/{pctg}.sam
done

samtools view -hb /localdisk/scratch/{pctg}.sam > /localdisk/scratch/{pctg}.bam
rm /localdisk/scratch/{pctg}.sam

pbalign --tmpDir=/localdisk/scratch/ --nproc=24 --minAccuracy=0.75 --minLength=50 --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr          --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=random --seed=1 /localdisk/scratch/{pctg}.bam {wdir}/{pctg}_ref.fa /localdisk/scratch/aln-ph-{pctg}.bam

samtools faidx {wdir}/{pctg}_ref.fa
variantCaller --algorithm=quiver -x 5 -X 120 -q 20 -j 24 -r {wdir}/{pctg}_ref.fa /localdisk/scratch/aln-ph-{pctg}.bam  -o {wdir}/cns-ph-{pctg}.fasta.gz -o {wdir}/cns-ph-{pctg}.fastq.gz

rm /localdisk/scratch/{pctg}.bam
rm /localdisk/scratch/aln-ph-{pctg}.bam*
'''.format(pctg=pctg,wdir=wdir)
        commands.append(cmd)

    for pctg in hcontigs:
        for hctg in hcontigs[pctg]:
            cmd='''
source /mnt/software/Modules/current/init/bash
module load smrtanalysis/mainline
module load samtools

set -vex

samtools view -H {wdir}/{hctg}.bam > /localdisk/scratch/{hctg}.sam
samtools view {wdir}/{hctg}.bam >> /localdisk/scratch/{hctg}.sam
samtools view {wdir}/{pctg}.bam >> /localdisk/scratch/{hctg}.sam
samtools view -hb /localdisk/scratch/{hctg}.sam > /localdisk/scratch/{hctg}.bam
rm /localdisk/scratch/{hctg}.sam

pbalign --tmpDir=/localdisk/scratch/ --nproc=24 --minAccuracy=0.75 --minLength=50 --minAnchorSize=12 --maxDivergence=30 --concordant --algorithm=blasr          --algorithmOptions=--useQuality --maxHits=1 --hitPolicy=random --seed=1 /localdisk/scratch/{hctg}.bam {wdir}/{hctg}_ref.fa /localdisk/scratch/aln-ph-{hctg}.bam

samtools faidx {wdir}/{hctg}_ref.fa
variantCaller --algorithm=quiver -x 5 -X 120 -q 20 -j 24 -r {hctg}_ref.fa /localdisk/scratch/aln-ph-{hctg}.bam  -o {wdir}/cns-ph-{hctg}.fasta.gz -o {wdir}/cns-ph-{hctg}.fastq.gz

rm /localdisk/scratch/{hctg}.bam
rm /localdisk/scratch/aln-ph-{hctg}.bam*
'''.format(hctg=hctg, pctg=pctg, wdir=wdir)
            commands.append(cmd)
    return commands


if __name__=='__main__':
    parser = get_parser()
    args=parser.parse_args()
    wdir=os.path.abspath(args.dir)
    queue = args.queue
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
    pool.map( merge_fq.merge_fq, contigs)
    pool.close()
    pool.join

