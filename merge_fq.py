#!/usr/bin/env python

from pbcore.io import FastqReader
from pbcore.io import FastqWriter
from sys import argv
import logging


min_match_len=20  # minimum length of matches between two mismatched regions
max_mismatch_len=40 # maximun length of mismatch
max_indel_len=min_match_len+500 # maximun length of indel allowed for a mismatch


logger = logging.getLogger(__name__)


'''
for f in fastq:
    print 'header', f.header
    print 'id', f.id 
    #print 'length', f.length
    print 'name', f.name
    print 'quality', len(f.quality), f.quality
    print 'reverseComplement', len(f.reverseComplement().sequence)
    print 'sequence', len(f.sequence)
'''


def next_mismatch(seq1, seq2):
    ''' return the begining and the ending where seq1 and seq2 mismatch in two tuples (b1, e1) (b2, e2)
        such that seq1[b1-1] == seq2[b2-1] and seq1[b1-1] != seq2[b2-1]
        e1, e2 are the first base where seq1[end1:end1+min_match_len] == seq2[end2:end2+min_match_len] 
    '''
    if len(seq1) < min_match_len:
        return

    for pos1 in range(len(seq1)-min_match_len+1):
        sstart=0 if pos1 < max_indel_len else pos1-max_indel_len
        pos2=seq2.find(seq1[pos1:pos1+min_match_len], sstart, sstart+min_match_len+max_indel_len)
        print 'init match', pos1, pos2
        if pos2 >= 0:
            break
    if pos2 != 0:
        yield (0, pos1), (0, pos2)
    prev_nucl=seq1[pos1]
    homo_pos1=pos1
    homo_pos2=pos2
    while pos1 < len(seq1) and pos2 < len(seq2):
        if seq1[pos1] != seq2[pos2]:
            for m1 in range(pos1+1, len(seq1)-min_match_len+1):
                #m2=kmerpos2.get(seq1[m1:m1+kmerlen], -1)
                d=m1-pos1
                #sstart=pos2 if d < max_indel_len else pos2+d-max_indel_len
                sstart=pos2
                m2=seq2.find(seq1[m1:m1+min_match_len], sstart, sstart+min_match_len+max_indel_len)
                #print 'm1', m1, 'm2', m2, 'seq1', seq1[m1:m1+min_match_len], 'seq2', seq2[sstart:sstart+min_match_len+d+max_indel_len]
                if m2 >= 0:
                    # m1+1, m2+1 because sometime lowest quility is at the edge
                    yield (homo_pos1, m1+1, pos1), (homo_pos2, m2+1, pos2)
                    pos1=m1+1
                    pos2=m2+1
                    break
            else:
                yield (homo_pos1, len(seq1)), (homo_pos2, len(seq2))
                break
        else: # skip match
            if seq1[pos1] != prev_nucl:
                prev_nucl=seq1[pos1]
                homo_pos1=pos1   # position for the first nucleotide of a homopolymer 
                homo_pos2=pos2
            pos1+=1
            pos2+=1


def pick_path(seq1, seq2):
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
    for seq_intv1, seq_intv2 in next_mismatch(seq1.sequence.upper(), seq2.sequence.upper()):
        #print seq_intv1, seq1.sequence[seq_intv1[0]-1:seq_intv1[1]+1], seq1.quality[seq_intv1[0]-1:seq_intv1[1]+1], seq1.quality[seq_intv1[2]], \
              #seq_intv2, seq2.sequence[seq_intv2[0]-1:seq_intv2[1]+1], seq2.quality[seq_intv2[0]-1:seq_intv2[1]+1], seq2.quality[seq_intv2[2]]
        seq_out += seq1.sequence[prev_seq_intv1[1]:seq_intv1[0]]    # take seq1 before first mismatch
        quality_out.extend(seq1.quality[prev_seq_intv1[1]:seq_intv1[0]])
        minQV1=min(seq1.quality[seq_intv1[0]:seq_intv1[1]])
        minQV2=min(seq2.quality[seq_intv2[0]:seq_intv2[1]])
        l=int(seq_intv1[1])-int(seq_intv1[0]) - (int(seq_intv2[1])-int(seq_intv2[0]))
        if abs(l) > 5 or abs(int(minQV1)-int(minQV2)) < 8 or minQV1 >= minQV2:
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


def merge_fq(ctg_name):
    ''' given a contig name (for example 000000F_001) merge cns-000000F_001.fastq.gz and cns-ph-000000F_001_fastq.gz
    '''
    seq_phased=read_fq( 'cns-'+ctg_name+'.fastq.gz' )
    seq_unphased=read_fq( 'cns-ph-'+ctg_name+'.fastq.gz' )
    seq_merged, quality_merged = pick_path(seq_phased, seq_unphased)
    with FastqWriter('cns-m-'+ctg_name+'.fastq.gz') as fqwriter:
        fqwriter.writeRecord(ctg_name+'|quiver_merged', seq_merged, quality_merged)



if __name__ == '__main__':
    seq1=read_fq(argv[1])
    #print seq1.sequence
    #print seq1.quality
    seq2=read_fq(argv[2])
    pick_path(seq1, seq2)
