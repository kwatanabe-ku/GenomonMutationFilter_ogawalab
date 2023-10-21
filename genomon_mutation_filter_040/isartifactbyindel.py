# -*- coding: utf-8 -*-
import os, copy
import re
import pysam
import subprocess
from auto_vivification import auto_vivification
from operator import itemgetter

target = re.compile( '([\+\-])([0-9]+)([ACGTNacgtn]+)' )
remove_chr = re.compile( '\^.' )

# at least 2 reads support indel
min_indel = 2

############################################################
def ReadRefGenome(chr, start, end, reference_genome):
    #1-indexed
    grange = "{}:{}-{}".format(chr, start, end)
    items = pysam.faidx(reference_genome, grange).strip('\n').split('\n')
    bases = [x for x in items if x[0] != '>'][0]
    
    return bases
        
############################################################
def pickupindel(chr, start, end, bamfile, search_length, reference_genome, mapping_quality_thres, exclude_sam_flags):

    indel           = auto_vivification()
    indels_multiple = auto_vivification()

    #exclude sam flags and MAPQ threshold
    reads = list(bamfile.fetch(chr, start-search_length, end+search_length))
    reads = [read for read in reads
                if read.flag & exclude_sam_flags == 0 and \
                   read.mapping_quality >= mapping_quality_thres
            ]
    reads_indel = [read for read in reads if len(read.get_blocks()) > 1]
    
    if len(reads_indel) > 1:
        for read in reads_indel:
            blocks = read.get_blocks()
            qname  = read.qname
            
            for k, block in enumerate(blocks[:-1]):
                indel_start = block[1] # the base before deletion (=0-indexed)
                
                if indel_start < start-search_length or end+search_length < indel_start:
                    continue
                
                #deletion
                if block[1] < blocks[k+1][0]:
                    indel_bases = 'N'*( blocks[k+1][0] - block[1] ) # you not need know exact bases.
                    key = '\t'.join( [chr, str(indel_start), "", indel_bases] )
                    indel_type  = '-'
                
                #insertion
                else:
                    aligned_pairs = read.get_aligned_pairs()
                    query_pos = [x[0] for x in aligned_pairs]
                    ref_pos   = [x[1] for x in aligned_pairs]
                    insertion_query_pos = ( query_pos[ref_pos.index(indel_start-1) + 1], 
                                            query_pos[ref_pos.index(indel_start)   - 1] )
                    indel_bases = read.query_sequence[
                                            insertion_query_pos[0]:insertion_query_pos[1]+1 ].upper()
                    key = '\t'.join( [chr, str(indel_start), "", indel_bases] )
                    indel_type  = '+'
                    
                if indel_type in indel.keys() and \
                   key in indel[indel_type].keys():
                    indel[indel_type][key].append(qname)
                    
                else:
                    indel[indel_type][key] = [qname]                    

        for indeltype in ( '+', '-' ):
            if indeltype in indel:
                for key in indel[ indeltype ].keys():
                    n_indel = len(set(indel[ indeltype ][ key ]))
                    if n_indel >= min_indel:
                        # Reading fasta only in deletion with reads more than min_indel, to minimize reading.
                        if indeltype == '-':
                            start = int(key.split('\t')[1]) + 1
                            end   = start + len(key.split('\t')[3]) - 1
                            bases = ReadRefGenome(chr, start, end, reference_genome)
                            key   = '\t'.join( key.split('\t')[:3] + [bases] )
                    
                        indels_multiple[ indeltype ][ key ] = n_indel
                          
    return indels_multiple

    
############################################################
def makeseq(refseq, starts, refs, alts, start_original, window):

    #sort 1.start 2.insertion is before SNV and deletion
    #in python, sort is stable.
    starts, refs, alts = zip(*sorted(
                                sorted(zip(starts, refs, alts), key=itemgetter(1), reverse=True),
                                key=itemgetter(0)
                            ))

    muttypes = []
    for ref, alt in zip(refs, alts):
        if ref == '-':   muttype = 'insertion'
        elif alt == '-': muttype = 'deletion'
        else:            muttype = 'SNV'
        
        muttypes.append(muttype)

    thresbypattern = {('SNV', 'insertion'):0,
                      ('SNV', 'deletion'):1,
                      ('SNV', 'SNV'):1,
                      ('deletion', 'insertion'):1,
                      ('deletion', 'deletion'):2,
                      ('deletion', 'SNV'):1,
                      ('insertion', 'insertion'):1,
                      ('insertion', 'deletion'):2,
                      ('insertion', 'SNV'):1}
        
    seq = copy.copy(refseq)
    shift = window-start_original
    #shift: 2つ以上変異を組み込むと、先にindelが入った場合に場所がずれていくので、そこを修正するため
    
    for n, (start, ref, alt, muttype) in enumerate(zip(starts, refs, alts, muttypes)):
        #同じサンプルで併存がみられるものだけをループのたびにpath_bamに残していく
        #併存がどれもみられなくなったら結合をやめてNoneを返す
        if n != 0:
            diff = start - ( starts[n-1] + len(refs[n-1]) - 1 )
            muttypepair = (muttypes[n-1], muttype)
            if diff < thresbypattern[muttypepair]:
                return None
            
        #start is 0-based
        end = start + len(ref)
            
        if muttype == 'SNV':
            seq = seq[:start+shift] + alt + seq[end+shift:] 
            shift += len(alt) - len(ref)
            
        elif muttype == 'deletion':
            seq = seq[:start+shift] + seq[end+shift:]
            shift -= len(ref)
            
        elif muttype == 'insertion':
            seq = seq[:start+shift+1] + alt + seq[end+shift:]
            shift += len(alt)
            
    return seq

############################################################
def makeEachReference(chr, start, end, ref, alt, indelpos, indelref, indelalt, output1, output2, window, reference_genome):
    
    ref_seq = ""
    label = ','.join([chr, str(start), str(end), ref, alt])
    range = chr + ":" + str(int(start) - window + 1) +"-"+ str(int(end) + window)
    for item in pysam.faidx(reference_genome, range):
        if item[0] == ">": continue
        ref_seq = ref_seq + item.rstrip('\n').upper()
        ref_seq = ref_seq.replace('>', '')
        ref_seq = ref_seq.replace(range, '')
        
    indelpos -= 1 #0-based
        
    seq_indelwithsnv = makeseq(ref_seq, [start, indelpos], [ref, indelref], [alt, indelalt], start, window)
    seq_indel = makeseq(ref_seq, [indelpos], [indelref], [indelalt], start, window)
    seq_snv = makeseq(ref_seq, [start], [ref], [alt], start, window)
    # seq = ref_seq
    # shift = (indelpos + 1) - (start + 1)
    
    # # for insertion
    # if indelref == "-":
        # seq = seq[0:(window + shift + 1)] + indelalt + seq[window + shift + 1:]
    # # for deletion
    # elif indelalt == "-":
        # seq = seq[0:(window + shift)] + seq[window + shift + len(indelref):]

    hOUT = open(output1, 'w')    
    
    print('>' + label + "_alt", file=hOUT)
    print(seq_snv, file=hOUT)
    print('>' + label + "_ref", file=hOUT)
    print(seq_indel, file=hOUT)
    
    hOUT.close()
    
    if seq_indelwithsnv is not None:
        
        hOUT = open(output2, 'w')    
        
        print('>' + label + "_alt", file=hOUT)
        print(seq_indelwithsnv, file=hOUT)
        print('>' + label + "_ref", file=hOUT)
        print(seq_indel, file=hOUT)
        
        hOUT.close()
    
    else:
        if os.path.isfile(output2):
            os.remove(output2)
    
    
    # # for insertion
    # if ref == "-":
        # seq     = seq[0:(window + 1)] + alt + seq[-window:]
    # # for deletion
    # elif alt == "-":
        # seq     = seq[0:window] + seq[-window:]
     # # for SNV
    # else:
        # seq     = seq[0:window] + alt + seq[-window:]

    
    
    
    

                            