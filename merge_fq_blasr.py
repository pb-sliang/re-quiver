#!/usr/bin/env python

from pbcore.io import FastqReader
from pbcore.io import FastqWriter
from pbcore.io import FastaWriter
from sys import argv
import logging


min_match_len=20  # minimum length of matches between two mismatched regions
max_indel_len=min_match_len+500 # maximun length of indel allowed for a mismatch


logger = logging.getLogger(__name__)


def next_mismatch(seq1, seq2, match_str, mseq1, mseq2):
    ''' 
    based on blasr -m5 alignment string match_str
    return the begining and the ending where seq1 and seq2 mismatch in two tuples (b1, e1) (b2, e2)
    such that seq1[b1-1] == seq2[b2-1] and seq1[b1-1] != seq2[b2-1]
    b1 and b2 are 'moved up' to include the last homopolymer. 
    This is because low QV for homopolymer can be at the either ends of the homopolymer
    seq1 and seq2 are from the fastq file, match the QV
    mseq1, match_str, mseq2 are from the -m5 output of blasr
    '''
    len_post_mismatch_homo=2   # this is needed b/c min QV near a homopolymer is sometimes one nucleotide outside homopolymer
    if len(seq1) <= len_post_mismatch_homo or len(seq2) <= len_post_mismatch_homo:
        return
    homo_pos1, homo_pos2=0, 0
    mpos=0
    pos1, pos2=0, 0
    prev_nucl=seq1[pos1]
    while mpos < len(match_str):
        if match_str[mpos] == '*':
            not_homo=False
            post_homo=0
            while mpos < len(match_str):
                if seq1[pos1] != prev_nucl and seq2[pos2] != prev_nucl:  # homopolymer breaks
                    not_homo=True
                if not_homo and match_str[mpos] == '|':
                    post_homo+=1
                    if post_homo > len_post_mismatch_homo:
                        yield (homo_pos1, pos1), (homo_pos2, pos2)
                        prev_nucl=seq1[pos1]
                        homo_pos1=pos1
                        homo_pos2=pos2
                        break
                if mseq1[mpos] != '-': pos1+=1
                if mseq2[mpos] != '-': pos2+=1
                mpos+=1
        else:
            if seq1[pos1] != prev_nucl:
                prev_nucl=seq1[pos1]
                homo_pos1=pos1   # position for the first nucleotide of a homopolymer 
                homo_pos2=pos2
            pos1+=1
            pos2+=1
            mpos+=1


def pick_path(seq1, seq2, match_str, mseq1, mseq2):
    ''' pick a path by merging sequences from seq1 and seq2 based on quality score 
        segment that is the edge of the path are from next_mismatch(seq1, seq2), 
        containing the intervals for mismatched sequences sandwitched between matched sequences
        quality score are compared and the sequence is taken from the one with higher quality score
        when quality score ties; seq1 is favored ; seq1 should be from phased reads
    '''

    seq_out=''
    quality_out=[]
    prev_seq_intv1=(0,0)
    ncorrected=0
    nuncorrected=0
    for seq_intv1, seq_intv2 in next_mismatch(seq1.sequence.upper(), seq2.sequence.upper(), match_str, mseq1, mseq2):
        print seq_intv1, seq1.sequence[seq_intv1[0]-1:seq_intv1[1]+1], seq1.quality[seq_intv1[0]-1:seq_intv1[1]+1], \
              seq_intv2, seq2.sequence[seq_intv2[0]-1:seq_intv2[1]+1], seq2.quality[seq_intv2[0]-1:seq_intv2[1]+1]
        seq_out += seq1.sequence[prev_seq_intv1[1]:seq_intv1[0]]    # take seq1 before first mismatch
        quality_out.extend(seq1.quality[prev_seq_intv1[1]:seq_intv1[0]])
        minQV1=min(seq1.quality[seq_intv1[0]:seq_intv1[1]])
        minQV2=min(seq2.quality[seq_intv2[0]:seq_intv2[1]])
        l=int(seq_intv1[1])-int(seq_intv1[0]) - (int(seq_intv2[1])-int(seq_intv2[0]))
        islowercase=any([c.islower() for c in seq1.sequence[seq_intv1[0]:seq_intv1[1]]])
        if islowercase or abs(l) > 2 or abs(int(minQV1)-int(minQV2)) < 20 or minQV1 >= minQV2:
            seq_out += seq1.sequence[seq_intv1[0]:seq_intv1[1]]
            quality_out.extend(seq1.quality[seq_intv1[0]:seq_intv1[1]])
            nuncorrected+=1
        else:
            seq_out += seq2.sequence[seq_intv2[0]:seq_intv2[1]]
            quality_out.extend(seq2.quality[seq_intv2[0]:seq_intv2[1]])
            ncorrected+=1
        prev_seq_intv1=seq_intv1

    seq_out += seq1.sequence[prev_seq_intv1[1]:]
    quality_out.extend(seq1.quality[prev_seq_intv1[1]:])

    logger.info("%s has %s differences. %s are corrected and %s are not", seq1.id, str(ncorrected+nuncorrected), str(ncorrected), str(nuncorrected))
    return seq_out, quality_out


def read_fq(fn):
    with FastqReader(fn) as fq:
        fastq=[f for f in fq]
    assert len(fastq)==1, '%r has more than one sequence'%fn
    return fastq[0]

def read_m5(fn):
    ''' read blasr -m 5 alignment output
    '''
    with open (fn) as fhm5:
        f=fhm5.readline().strip().split()
        mseq1, match_str, mseq2 = [x for x in f[16:]]
    return mseq1, match_str, mseq2


def merge_fq(ctg_name):
    ''' given a contig name (for example 000000F_001) merge cns-000000F_001.fastq.gz and cns-ph-000000F_001_fastq.gz
    '''
    seq_phased=read_fq( 'cns-'+ctg_name+'.fastq.gz' )
    seq_unphased=read_fq( 'cns-ph-'+ctg_name+'.fastq.gz' )
    mseq1, match_str, mseq2 = read_m5( 'fa_aln_'+ctg_name+'.m5')
    print ctg_name
    seq_merged, quality_merged = pick_path(seq_phased, seq_unphased, match_str, mseq1, mseq2)
    with FastqWriter('cns-m-'+ctg_name+'.fastq.gz') as fqwriter:
        fqwriter.writeRecord(ctg_name+'|quiver_merged', seq_merged, quality_merged)
    with FastaWriter('cns-m-'+ctg_name+'.fasta.gz') as fawriter:
        fawriter.writeRecord(ctg_name+'|quiver_merged', seq_merged)



if __name__ == '__main__':
    seq1=read_fq(argv[1])
    #print seq1.sequence
    #print seq1.quality
    seq2=read_fq(argv[2])
    pick_path(seq1, seq2)
