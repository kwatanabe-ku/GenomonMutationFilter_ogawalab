# -*- coding: utf-8 -*-
import sys
import os
import re, copy
import pysam
import logging
import subprocess
import scipy.special
from scipy.stats import fisher_exact as fisher
import math
import isartifactbyindel
import cruciformfilter_core

md_divide = re.compile('([0-9]+|[\^ACTGNR]+)')

#
# Class definitions
#
class realignment_filter:

    def __init__(self,
                referenceGenome,
                tumor_min_mismatch,
                normal_max_mismatch,
                search_length,
                score_difference,
                blat,
                header_flag,
                max_depth,
                exclude_sam_flags,
                mapping_quality,
                base_quality,
                mutation_cluster,
                removedup_byumi,
                removedup_duplex,
                removedup_thres,
                cruciform_filter,
                cruciform_search_window,
                cruciform_min_matchrate,
                search_length_indel,
                debug_mode
                ):
                
        self.reference_genome      =  referenceGenome
        self.tumor_min_mismatch    =  tumor_min_mismatch
        self.normal_max_mismatch   =  normal_max_mismatch
        self.window                =  search_length
        self.score_difference      =  score_difference
        # Running blat takes time the most, but inevitable so far.
        self.blat_cmds             =  [blat, '-fine']
        self.header_flag           =  header_flag
        self.max_depth             =  max_depth
        self.exclude_sam_flags     =  exclude_sam_flags
        self.mapping_quality       =  mapping_quality
        self.base_quality          =  base_quality
        self.mutation_cluster      =  mutation_cluster
        self.removedup_byumi       =  removedup_byumi
        self.removedup_duplex      =  removedup_duplex
        self.removedup_thres       =  removedup_thres
        self.cruciform_filter      =  cruciform_filter
        self.cruciform_search_window =  cruciform_search_window
        self.cruciform_min_matchrate =  cruciform_min_matchrate
        self.search_length_indel   =  search_length_indel
        self.debug_mode            =  debug_mode
        
        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        self.target = re.compile( '([\+\-])([0-9]+)([ACGTNacgtn]+)' )
        self.remove_chr = re.compile( '\^.' )
        
        self.cruciformcore = cruciformfilter_core.cruciform_filter(self.reference_genome,
                                                                  self.cruciform_search_window,
                                                                  self.cruciform_min_matchrate,
                                                                  self.max_depth,
                                                                  self.exclude_sam_flags,
                                                                  self.mapping_quality,
                                                                  self.base_quality,
                                                                  self.removedup_byumi,
                                                                  self.removedup_thres,
                                                                  self.debug_mode
                                                                )
    
    ############################################################
    def math_log_fisher_pvalue(self,fisher_pvalue):

        val = float(0.0)
        if fisher_pvalue < 10**(-60):
            val = float(60.0)
        elif fisher_pvalue  > 1.0 - 10**(-10) :
            val = float(0.0)
        else:
            val = -math.log( fisher_pvalue, 10 )
                                                                
        return val

    #####Fix:v0.3.0 20221114: makeTwoReferenceが2回定義されていたのを削除

    ############################################################
    def extractRead(self, bamfile, chr, start, end):
        reads = [read for read in list(bamfile.fetch(chr,start-2,end+2))
                               if read.flag & self.exclude_sam_flags == 0]
        
        return reads

    def filterreads_mapqual(self, reads, mapping_quality_thres):
        n_reads = len(reads)
        n_readsmapq0 = sum([read.mapping_quality == 0 for read in reads])
        rate_mapqual0 = '{0:.3f}'.format(float(n_readsmapq0) / n_reads) if n_reads != 0 else '0'
        
        reads = [read for read in reads if read.mapping_quality >= mapping_quality_thres]
        
        return reads, str(n_readsmapq0), rate_mapqual0

            
    def removedupbyumi(self, reads):
        if len(reads) == 0 or not reads[0].has_tag('ZA') or not reads[0].has_tag('RX'): return reads
        
        #markdup
        if reads[0].has_tag('ZA'):
            #####Fix:v0.3.0 20221114: UMIがあってremovedup_byumiがTrueの場合、①UMIが片方一致+そちら側の端が一致②反対側の端がTHRES_SHIFT以内　の場合、PCRerrorを伴ったduplicateとみなし、重複を削除する
            isplusstrand = lambda x: x & 16 == 0
            isreverse = lambda x: x & 144 not in (0, 144)
            readsinfo = [{'start_insert': read.pos if isplusstrand(read.flag) else read.aend + read.tlen,
                          'end_insert'  : read.pos + read.tlen if isplusstrand(read.flag) else read.aend,
                          'qname'       : read.qname,
                          'meanQ'       : sum(read.query_qualities) / float(len(read.query_qualities)),
                          'umi1'        : read.get_tag('ZA') if not isreverse(read.flag) else read.get_tag('ZB'),
                          'umi2'        : read.get_tag('ZB') if not isreverse(read.flag) else read.get_tag('ZA')}
                        for read in reads]
                        
            readsinfo = sorted(readsinfo, key=lambda x: (x['umi1'], x['start_insert'], x['end_insert']))
            
            THRES_SHIFT = self.removedup_thres
            
            #startが一致してendが少しずれているもの
            tmp_umi1  = ''
            tmp_start = -100
            tmp_end   = -100
            tmp_meanQ = 0
            tmp_k     = -1
            isremain = [True] * len(readsinfo)
            for k, read in enumerate(readsinfo):
                if read['umi1']                 == tmp_umi1  and \
                   read['start_insert']         == tmp_start and \
                   read['end_insert'] - tmp_end <= THRES_SHIFT:
                   
                    if read['meanQ'] > tmp_meanQ:
                        isremain[tmp_k] = False
                        tmp_meanQ       = read['meanQ']
                        tmp_k           = k
                    else:
                        isremain[k]     = False
                   
                    tmp_end   = read['end_insert']
                else:
                    tmp_umi1  = read['umi1']
                    tmp_start = read['start_insert']
                    tmp_end   = read['end_insert']
                    tmp_meanQ = read['meanQ']
                    tmp_k     = k
                    
            readsinfo_rmdup = [x for x,y in zip(readsinfo, isremain) if y]
            readsinfo_rmdup = sorted(readsinfo_rmdup, key=lambda x: (x['umi2'], x['end_insert'], x['start_insert']))
                    
            #endが一致してstartが少しずれているもの
            tmp_umi2  = ''
            tmp_start = -100
            tmp_end   = -100
            tmp_meanQ = 0
            tmp_k     = -1
            isremain = [True] * len(readsinfo_rmdup)
            for k, read in enumerate(readsinfo_rmdup):
                if read['umi2']                     == tmp_umi2    and \
                   read['start_insert'] - tmp_start <= THRES_SHIFT and \
                   read['end_insert']               == tmp_end:
                   
                    if read['meanQ'] > tmp_meanQ:
                        isremain[tmp_k] = False
                        tmp_meanQ       = read['meanQ']
                        tmp_k           = k
                    else:
                        isremain[k]     = False
                   
                    tmp_start   = read['start_insert']
                else:
                    tmp_umi2  = read['umi2']
                    tmp_start = read['start_insert']
                    tmp_end   = read['end_insert']
                    tmp_meanQ = read['meanQ']
                    tmp_k     = k
                    
            readsinfo_rmdup = [x for x,y in zip(readsinfo_rmdup, isremain) if y]
            qnames_remain = [x['qname'] for x in readsinfo_rmdup]
            
            reads = [read for read in reads if read.qname in qnames_remain]

            return reads
            
        #consensus (note that markdup bam also has RX tag)
        elif reads[0].has_tag('RX'):
               
            # removedup_byumiがTrueの場合、
            # UMI一致+端がTHRES_SHIFT以内
            #　の場合、PCRerrorを伴ったduplicateとみなし、重複を削除する
            isplusstrand = lambda x: x & 16 == 0
            isreverse = lambda x: x & 144 not in (0, 144)
            readsinfo = [{'start_insert': read.pos if isplusstrand(read.flag) else read.aend + read.tlen,
                          'end_insert'  : read.pos + read.tlen if isplusstrand(read.flag) else read.aend,
                          'qname'       : read.qname,
                          'meanQ'       : sum(read.query_qualities) / float(len(read.query_qualities)),
                          'umipair'        : read.get_tag('RX')}
                        for read in reads]
                        
            readsinfo = sorted(readsinfo, key=lambda x: (x['umipair'], x['start_insert'], x['end_insert']))
            
            THRES_SHIFT = self.removedup_thres
            
            tmp_umi   = ''
            tmp_start = -100
            tmp_end   = -100
            tmp_meanQ = 0
            tmp_k     = -1
            isremain = [True] * len(readsinfo)
            for k, read in enumerate(readsinfo):
                if read['umipair']                  == tmp_umi  and \
                   read['start_insert'] - tmp_start <= THRES_SHIFT and \
                   read['end_insert']   - tmp_end   <= THRES_SHIFT:
                   
                    if read['meanQ'] > tmp_meanQ:
                        isremain[tmp_k] = False
                        tmp_meanQ       = read['meanQ']
                        tmp_k           = k
                    else:
                        isremain[k]     = False
                   
                    tmp_end   = read['end_insert']
                else:
                    tmp_umi   = read['umipair']
                    tmp_start = read['start_insert']
                    tmp_end   = read['end_insert']
                    tmp_meanQ = read['meanQ']
                    tmp_k     = k
                    
            readsinfo_rmdup = [x for x,y in zip(readsinfo, isremain) if y]
            
            qnames_remain = [x['qname'] for x in readsinfo_rmdup]
            
            reads = [read for read in reads if read.qname in qnames_remain]

            return reads
            

    # def removedupbyumi_duplex(self, reads):
        # if len(reads) == 0 or not reads[0].has_tag('RX'): return reads
           
        # # removedup_byumiがTrueの場合、
        # # UMI一致+端がTHRES_SHIFT以内
        # #　の場合、PCRerrorを伴ったduplicateとみなし、重複を削除する
        # isplusstrand = lambda x: x & 16 == 0
        # isreverse = lambda x: x & 144 not in (0, 144)
        # readsinfo = [{'start_insert': read.pos if isplusstrand(read.flag) else read.aend + read.tlen,
                      # 'end_insert'  : read.pos + read.tlen if isplusstrand(read.flag) else read.aend,
                      # 'qname'       : read.qname,
                      # 'meanQ'       : sum(read.query_qualities) / float(len(read.query_qualities)),
                      # 'umipair'        : read.get_tag('RX')}
                    # for read in reads]
                    
        # readsinfo = sorted(readsinfo, key=lambda x: (x['umipair'], x['start_insert'], x['end_insert']))
        
        # THRES_SHIFT = self.removedup_thres
        
        # tmp_umi   = ''
        # tmp_start = -100
        # tmp_end   = -100
        # tmp_meanQ = 0
        # tmp_k     = -1
        # isremain = [True] * len(readsinfo)
        # for k, read in enumerate(readsinfo):
            # if read['umipair']                  == tmp_umi  and \
               # read['start_insert'] - tmp_start <= THRES_SHIFT and \
               # read['end_insert']   - tmp_end   <= THRES_SHIFT:
               
                # if read['meanQ'] > tmp_meanQ:
                    # isremain[tmp_k] = False
                    # tmp_meanQ       = read['meanQ']
                    # tmp_k           = k
                # else:
                    # isremain[k]     = False
               
                # tmp_end   = read['end_insert']
            # else:
                # tmp_umi   = read['umipair']
                # tmp_start = read['start_insert']
                # tmp_end   = read['end_insert']
                # tmp_meanQ = read['meanQ']
                # tmp_k     = k
                
        # readsinfo_rmdup = [x for x,y in zip(readsinfo, isremain) if y]
        
        # qnames_remain = [x['qname'] for x in readsinfo_rmdup]
        
        # reads = [read for read in reads if read.qname in qnames_remain]

        # return reads
        
    def add_readNo_to_qname(self, read):
        if read.flag & 64 != 0:
            return read.qname + '/1'
        else:
            return read.qname + '/2'

    def writeRead(self, reads, ref, alt, output):
        hOUT = open(output, 'w')    
        
        for k, read in enumerate(reads):
            flags = format(int(read.flag), "#014b")[:1:-1]

            tempSeq = read.seq
            
            #####Fix:v0.3.0 SNVの場合はbase quality が低いものはNに(indelだと予期しない結果になるため注意)
            if ref != '-' and alt != '-':
                tempSeq = ''.join([x if q >= self.base_quality else 'N' for x, q in zip(tempSeq, read.query_qualities)])
            
            # reverse read
            if read.flag & 16 != 0:
                tempSeq = "".join(self.complement.get(base) for base in reversed(str(tempSeq)))

            print('>' + self.add_readNo_to_qname(read), file=hOUT)
            print(tempSeq, file=hOUT)
                
        hOUT.close()

    ############################################################
    def makeTwoReference(self, chr,start,end,ref,alt, output):

        hOUT = open(output, 'w')
        
        seq = ""
        label = ','.join([chr, str(start), str(end), ref, alt])
        range = chr + ":" + str(max(0,int(start) - self.window + 1)) +"-"+ str(int(end) + self.window)
        for item in pysam.faidx(self.reference_genome, range):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()
            seq = seq.replace('>', '')
            seq = seq.replace(range, '')
            
        ref_seq = seq

        print('>' + label + "_ref", file=hOUT)
        print(ref_seq, file=hOUT)

        # for insertion
        if ref == "-":
            seq     = seq[0:(self.window + 1)] + alt + seq[-self.window:]
            ref_seq = ref_seq[self.window + 1 :]
        # for deletion
        elif alt == "-":
            seq     = seq[0:self.window] + seq[-self.window:]
            ref_seq = ref_seq[self.window:]
         # for SNV
        else:
            seq     = seq[0:self.window] + alt + seq[-self.window:]
            ref_seq = ref_seq[self.window:]

        print('>' + label + "_alt", file=hOUT)
        print(seq, file=hOUT)

        hOUT.close()
        
        return ref_seq
        
    ############################################################
    def IsMicroSatellite(self, ref, alt, ref_seq):
    #特に指定なければ最初の4塩基が連続していればmicrosatelliteとする

        def GetRepeatUnit(bases):
            if len(bases) == 1:
                return bases
            
            for k in range(1, int(len(bases)/2)+1):
                if len(bases) % k != 0: continue
                
                if bases[:-k] == bases[k:]:
                    return bases[:k]
                
            return bases
            
        if ref == "-": #ins
            unit = GetRepeatUnit(alt)
        elif alt == "-": #del
            unit = GetRepeatUnit(ref)
        else: #SNV
            return '---'
            
        match = re.match('('+ unit +')+', ref_seq)
        if match is None: return '---'
        
        assert match.end() % len(unit) == 0
        
        n_repeat = unit + ':' + str( int( match.end() / len(unit) ) )
        
        if len(ref_seq) < match.end() + len(unit):
            n_repeat = '>' + n_repeat
        
        return n_repeat
            
    ############################################################
    def countmutationofread(self, start, end, ref, alt, read):
        ### read.get_aligned_pairs(True, True) and get_reference_sequence()
        ### can cause segmentation fault in old version of pysam
        # # this does not include del, ins and soft-clip
        # pairs = read.get_aligned_pairs(True, True)
        
        # # remove bases around current variant (for DBS or indel+SNV)
        # if ref == '-': #ins
            # pairs = [ pair for pair in pairs if pair[1] < start - 1 or pair[1] > end + 1 ]
        # else: #SNV, del
            # pairs = [ pair for pair in pairs if pair[1] < start - 2 or pair[1] > end + 1 ]
        
        # # mismatch base is lower-case
        # return sum([x[2].islower() for x in pairs])
        
        # pairs = [ pair for pair in pairs if pair[1] is not None] # del->remain, ins,soft-clip->remove
        
        pairs = read.get_aligned_pairs(True)               
        MD_split = md_divide.findall(read.get_tag('MD'))
        
        mut_pos = []
        pos_now = 0
        for MD in MD_split:
            if MD[0] == '^': continue #deletion
            
            if MD.isdigit():
                pos_now += int(MD)
                
            else: #variant bases
                mut_pos += [pair[1] for pair in pairs[pos_now : pos_now + len(MD)]]
                pos_now += len(MD)
                
        if ref == '-': #ins
            mut_pos = [x for x in mut_pos if x < start - 1 or x > end + 1 ]
        else: #SNV, del   
            mut_pos = [x for x in mut_pos if x < start - 2 or x > end + 1 ]
        
        # For debugging. Even if this condition is not met, there is no major problem
        assert len(pairs) == pos_now
                
        return len(mut_pos)
        
    ############################################################
    def summarizeRefAlt(self, inputPsl, start, end, ref, alt, reads):
        qname_readswithn = [self.add_readNo_to_qname(read) for read in reads]
        
        # Split for each qname first.
        qnames         = []
        lines_forqname = {}
        hIN = open(inputPsl, 'r')
        for line in hIN:
            F = line.rstrip('\n').split('\t')
            if F[0].isdigit() == False: continue
            
            qname = F[9][:-2]
            if qname not in qnames:
                qnames.append(qname)
                lines_forqname[qname] = [line]
            else:
                lines_forqname[qname].append(line)
                
        hIN.close()
                
        # In _mc, calculate reads excluding unreliable reads with many mutations.
        # If one of the pair reads is "ref" and the other is "alt", the read is counted as "other".
        readsOther = []
        readsAlt   = []
        readsRef   = []
        readsOther_mc = []
        readsAlt_mc   = []
        readsRef_mc   = []
        strand_alt = [0, 0] #( +, - )
                
        for qname in qnames:
        
            lines_qname = lines_forqname[qname]
            qnames_with_nth = [x.split('\t')[9] for x in lines_qname]
            
            n_other = 0
            n_alt   = 0
            n_ref   = 0
            n_other_mc = 0
            n_alt_mc   = 0
            n_ref_mc   = 0
            
            for qname_with_nth in set(qnames_with_nth):
                read = reads[ qname_readswithn.index( qname_with_nth ) ]
                flags = format(int(read.flag), "#014b")[:1:-1]
                    
                if self.countmutationofread(start, end, ref, alt, read) >= self.mutation_cluster:
                    flag_mc = True
                else:
                    flag_mc = False
                
                flag_reverse = True if flags[4] == "1" else False
                
                lines = [line for line, q in zip(lines_qname, qnames_with_nth) if q == qname_with_nth]
                
                if len(lines) < 2: continue
                
                scores_ref  = [-1000]
                scores_alt  = [-1000]
                
                for k, line in enumerate(lines):
                    F = line.rstrip('\n').split('\t')
                    
                    x = F[18].split(',')
                    y = F[19].split(',') # not used
                    z = F[20].split(',')
                    
                    number_of_mismatch = int(F[3])
                
                    for i in range(1, int(F[17])):
                    
                        lz = int(z[i]) - int(z[i - 1]) - int(x[i - 1]) 
                        if (lz > 0): number_of_mismatch += lz

                    score  = int(F[0]) - number_of_mismatch
                    if F[13][-3:] == 'alt':
                        scores_alt.append(score)
                    elif F[13][-3:] == 'ref':
                        scores_ref.append(score)
                
                sort_score = sorted(scores_alt + scores_ref, reverse=True)
                if len(lines) >= 3 and sort_score[1] - sort_score[2] <= self.score_difference:
                    continue
                    
                if max(scores_alt) == max(scores_ref) and max(scores_ref) != -1000:
                    n_other += 1
                    if not flag_mc:
                        n_other_mc += 1
                        
                elif max(scores_alt) >  max(scores_ref):
                    n_alt   += 1
                    if not flag_mc:
                        n_alt_mc += 1
                        
                    if flag_reverse:
                        strand_alt[1] += 1
                    else:
                        strand_alt[0] += 1
                        
                elif max(scores_alt) <  max(scores_ref):
                    n_ref   += 1
                    if not flag_mc:
                        n_ref_mc += 1
                        
            if   n_ref == 0 and n_alt == 0 and n_other > 0:
                readsOther.append(qname)
            elif n_ref >  0 and n_alt >  0:
                readsOther.append(qname)
            elif n_ref >  0 and n_alt == 0:
                readsRef.append(qname)
            elif n_ref == 0 and n_alt >  0:
                readsAlt.append(qname)
                
            if   n_ref_mc == 0 and n_alt_mc == 0 and n_other_mc > 0:
                readsOther_mc.append(qname)
            elif n_ref_mc >  0 and n_alt_mc >  0:
                readsOther_mc.append(qname)
            elif n_ref_mc >  0 and n_alt_mc == 0:
                readsRef_mc.append(qname)
            elif n_ref_mc == 0 and n_alt_mc >  0:
                readsAlt_mc.append(qname)
                
        if sum(strand_alt) != 0:
            strand_ratio = '{0:.3f}'.format(strand_alt[0] / float( strand_alt[0] + strand_alt[1] ))
        else:
            strand_ratio = '---'
            
        return( readsRef,    readsAlt,    readsOther,
                readsRef_mc, readsAlt_mc, readsOther_mc, strand_ratio )
    
    ############################################################
    def applycruciformfilter(self, reads, tumor_samfile, genome_fa, chr, start, end, ref, alt):
        mutreads = self.cruciformcore.extractmutread(tumor_samfile, chr, start+1, ref, alt)
        
        if len(mutreads) == 0:
            return "0", "0"
        
        if self.removedup_duplex:
            mutreads = self.removedupbyumi(mutreads)
        elif self.removedup_byumi:
            mutreads = self.removedupbyumi(mutreads)
            
        
        invreads, noinvreads = self.cruciformcore.checkinv(mutreads, genome_fa, chr, start, end, ref, alt, self.cruciform_search_window)
        
        n_noinv = str(len(set([x.qname for x in noinvreads])))
        n_inv   = str(len(set([x.qname for x in invreads]) - 
                            set([x.qname for x in noinvreads])))
                            
        return n_inv, n_noinv
    
    ############################################################
    def isindel(self, chr, start, end, ref, alt, reads, tumor_alt, qn_tumor_alt, tumor_samfile, output):
        if (ref == '-' or alt == '-') or len(qn_tumor_alt) == 0:
            return '---', '---'
        
        # indels_nearby = isartifactbyindel.pickupindel(
                            # chr, start, end, ref, alt, in_tumor_bam,
                            # self.search_length_indel, self.samtools_path, self.samtools_params, self.reference_genome
                       # )        
        indels_nearby = isartifactbyindel.pickupindel(
                            chr, start, end, tumor_samfile,
                            self.search_length_indel, self.reference_genome,
                            self.mapping_quality,     self.exclude_sam_flags)
                       
                       
        if len(indels_nearby) == 0:
            return '---', '---'
                       
        # only use variant read
        reads_alt = [read for read in reads if read.qname in qn_tumor_alt]
                 
        ###      
        hOUT = open(output + ".tmp.indel.fa", 'w')
        
        if len(reads_alt) == 0:
            hOUT.close()
            return '---', '---'
        
        for k, read in enumerate(reads_alt):
            flags = format(int(read.flag), "#014b")[:1:-1]

            tempSeq = read.seq # not filtered by base quality (indelはおかしくなることあり)
            
            if flags[4] == "1":
                tempSeq = "".join(self.complement.get(base) for base in reversed(str(tempSeq)))

            # the first read
            if flags[6] == "1":
                qname_forblat = '>' + read.qname + '/1'
            else:
                qname_forblat = '>' + read.qname + '/2'
 
            print(qname_forblat, file=hOUT)
            print(tempSeq, file=hOUT)
                
        hOUT.close()
        ###
        
                       
        min_n_alt = 10000
        min_key = '---'
        
        for indeltype in ( '+', '-' ):
            if indeltype in indels_nearby:
                for key in indels_nearby[ indeltype ].keys():
                
                    #indelpos is same as genomon output format
                    indelpos = int(key.split( '\t' )[ 1 ])
                    bases = key.split( '\t' )[ 3 ]
                    if indeltype == '+': #insertion
                        indelref = '-'
                        indelalt = bases
                    else: #deletion
                        indelref = bases
                        indelalt = '-'
                        indelpos += 1
                    
                    isartifactbyindel.makeEachReference(
                            chr, start, end, ref, alt, indelpos, indelref, indelalt,
                            output + ".tmp.refalt.indel1.fa", output + ".tmp.refalt.indel2.fa",
                            self.window, self.reference_genome
                        )
                        
                    FNULL = open(os.devnull, 'w')
                    
                    # indel vs snv
                    subprocess.check_call(
                        self.blat_cmds + [output + ".tmp.refalt.indel1.fa", output + ".tmp.indel.fa", output + ".tmp.indel1.psl"], 
                        stdout = FNULL, stderr = FNULL)
                        
                    _, qnames_tumor_alt_indel1, qnames_tumor_other_indel1, _, _, _, _  = \
                        self.summarizeRefAlt(output + ".tmp.indel1.psl", start, end, ref, alt, reads_alt)
                    
                    # indel+snv vs indel
                    if os.path.isfile(output + ".tmp.refalt.indel2.fa"):
                        subprocess.check_call(
                            self.blat_cmds + [output + ".tmp.refalt.indel2.fa", output + ".tmp.indel.fa", output + ".tmp.indel2.psl"], 
                            stdout = FNULL, stderr = FNULL)
                            
                        _, qnames_tumor_alt_indel2, _, _, _, _, _  = \
                            self.summarizeRefAlt(output + ".tmp.indel2.psl", start, end, ref, alt, reads_alt)
                    else:
                        qnames_tumor_alt_indel2 = []
                        
                    FNULL.close()
                    
                    tumor_alt_indel = len(set(qnames_tumor_alt_indel1 + qnames_tumor_alt_indel2))
                        
                    if tumor_alt_indel < min_n_alt and tumor_alt_indel < tumor_alt:
                        min_n_alt = tumor_alt_indel
                        min_key = '{0}:{1}{2},{3}'.format(
                                        str(indelpos), indeltype, bases, indels_nearby[ indeltype ][ key ])
                        
        if min_key == '---':
            min_n_alt = '---'
        else:
            min_n_alt = str(min_n_alt)
            
        return min_n_alt, min_key
                        
    
    ############################################################
    def filter(self, in_tumor_bam, in_normal_bam, output, in_mutation_file):

        if in_tumor_bam and in_normal_bam:
            mode = 'pair'
        elif in_tumor_bam:
            mode = 'single'
        else:
            return
            
        srcfile = open(in_mutation_file,'r')
        hResult = open(output,'w')
        
        tumor_samfile = pysam.Samfile(in_tumor_bam, "rb")
        if mode == 'pair':
            normal_samfile = pysam.Samfile(in_normal_bam, "rb")
        genome_fa     = pysam.FastaFile(self.reference_genome)

        if self.header_flag:
            if mode == 'pair':
                newheader = (
                    "readPairNum_tumor\tvariantPairNum_tumor\totherPairNum_tumor\t"+"readPairNum_normal\tvariantPairNum_normal\totherPairNum_normal\t"+
                    "Num_tumor(mutcluster)\tstrandRatio_tumor(realignment)\t"+
                    "Num_mapq0_tumor\tRate_mapq0_tumor\tvariantPairNum_normal_allmapq\t"+
                    "P-value(fisher_realignment)\tNum_Repeat\tvarNum_tumor_withindel\tindel_info")
                         
                         
            else:
                newheader = (
                    "readPairNum\tvariantPairNum\totherPairNum\t"+
                    "Num(mutcluster)\tstrandRatio(realignment)\tNum_mapq0\tRate_mapq0\t"+
                    "10%_posterior_quantile(realignment)\tposterior_mean(realignment)\t90%_posterior_quantile(realignment)\t"+
                    "Num_Repeat\tvarNum_tumor_withindel\tindel_info")
                
            header = srcfile.readline().rstrip('\n')  
            print(header +"\t"+ newheader, hOUT=hResult)
    
        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
            # annovar input file (not zero-based number)
            chr, start, end, ref, alt  = (itemlist[0], (int(itemlist[1]) - 1), int(itemlist[2]), itemlist[3], itemlist[4])
            
            if self.debug_mode:
                 print('{0}:{1}-{2}:{3}>{4}'.format(chr, start, end, ref, alt))
            
            tumor_ref = tumor_alt = tumor_other = normal_ref = normal_alt = normal_other = '---'
            tumor_ref_mc = tumor_alt_mc = tumor_other_mc = strand_ratio = '---'
            beta_01 = beta_mid = beta_09 = log10_fisher_pvalue = '---'
            n_repeat = min_n_alt_withindel = min_key_indel = normal_alt_all = '---'
            n_mapqual0 = rate_mapqual0 = '---'
            n_inv = n_noinv = '---'
            
            ref_seq = self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")
            n_repeat = self.IsMicroSatellite(ref, alt, ref_seq)

            if tumor_samfile.count(chr,start,end) < self.max_depth:

                reads = self.extractRead(tumor_samfile, chr, start, end)
                
                if len(reads) > 0:
                    reads, n_mapqual0, rate_mapqual0 = self.filterreads_mapqual(reads, self.mapping_quality)
                            
                    if self.removedup_duplex and len(reads) > 0:
                        reads = self.removedupbyumi(reads)
                    elif self.removedup_byumi and len(reads) > 0:
                        reads = self.removedupbyumi(reads)
                        
                    self.writeRead(reads, ref, alt, output + ".tmp.fa")
                
                    with open(os.devnull, 'w') as FNULL:
                        subprocess.check_call(
                            self.blat_cmds + [output + ".tmp.refalt.fa", output + ".tmp.fa", output + ".tmp.psl"], 
                            stdout = FNULL, stderr = subprocess.STDOUT)
                    
                    qn_tumor_ref,    qn_tumor_alt,    qn_tumor_other, \
                    qn_tumor_ref_mc, qn_tumor_alt_mc, qn_tumor_other_mc, strand_ratio = \
                        self.summarizeRefAlt(output + ".tmp.psl", start, end, ref, alt, reads)
                        
                    if self.debug_mode:
                        print('tumor-ref:', ','.join(qn_tumor_ref))
                        print('tumor-alt:', ','.join(qn_tumor_alt))
                        
                    if self.cruciform_filter:
                        n_inv, n_noinv = self.applycruciformfilter(reads, tumor_samfile, genome_fa, chr, start, end, ref, alt)
                        
                    tumor_ref, tumor_alt, tumor_other, \
                    tumor_ref_mc, tumor_alt_mc, tumor_other_mc = \
                        [len(set(qn)) for qn in (qn_tumor_ref,    qn_tumor_alt,    qn_tumor_other, 
                                                 qn_tumor_ref_mc, qn_tumor_alt_mc, qn_tumor_other_mc)]
                        
                    min_n_alt_withindel, min_key_indel = self.isindel(
                            chr, start, end, ref, alt, reads, tumor_alt, qn_tumor_alt, tumor_samfile, output)
                
                    if mode == 'single':
                        beta_01  = '{0:.3f}'.format(float(scipy.special.btdtri( int(tumor_alt) + 1, int(tumor_ref) + 1, 0.1 )))
                        beta_mid = '{0:.3f}'.format(float( int(tumor_alt) + 1 ) / float( int(tumor_ref) + int(tumor_alt) + 2 ))
                        beta_09  = '{0:.3f}'.format(float(scipy.special.btdtri( int(tumor_alt) + 1, int(tumor_ref) + 1, 0.9 )))
            
            if mode == 'pair' and normal_samfile.count(chr,start,end) < self.max_depth:

                reads = self.extractRead(normal_samfile, chr, start, end)
                
                if len(reads) > 0:
                    reads_all = copy.deepcopy(reads)
                    
                    reads, n_mapqual0_normal, _ = self.filterreads_mapqual(reads, self.mapping_quality)
                
                    if self.removedup_duplex and len(reads) > 0:
                        reads = self.removedupbyumi(reads)
                    if self.removedup_byumi and len(reads) > 0:
                        reads = self.removedupbyumi(reads)
                        
                    self.writeRead(reads, ref, alt, output + ".tmp.fa")
                
                    with open(os.devnull, 'w') as FNULL:
                        subprocess.check_call(
                            self.blat_cmds + [output + ".tmp.refalt.fa", output + ".tmp.fa", output + ".tmp.psl"], 
                            stdout = FNULL, stderr = subprocess.STDOUT)
                    
                    qn_normal_ref, qn_normal_alt, qn_normal_other, _, _, _, _ = \
                        self.summarizeRefAlt(output + ".tmp.psl", start, end, ref, alt, reads)

                    if self.debug_mode:
                        print('normal-ref:', ','.join(qn_normal_ref))
                        print('normal-alt:', ','.join(qn_normal_alt))
                        
                    normal_ref, normal_alt, normal_other = \
                        [len(set(qn)) for qn in (qn_normal_ref, qn_normal_alt, qn_normal_other)]
                        
                    if n_mapqual0_normal != '0':
                        if self.removedup_duplex and len(reads) > 0:
                            reads_all = self.removedupbyumi(reads_all)
                        if self.removedup_byumi and len(reads) > 0:
                            reads_all = self.removedupbyumi(reads_all)
                            
                        self.writeRead(reads_all, ref, alt, output + ".tmp.fa")
                    
                        with open(os.devnull, 'w') as FNULL:
                            subprocess.check_call(
                                self.blat_cmds + [output + ".tmp.refalt.fa", output + ".tmp.fa", output + ".tmp.psl"], 
                                stdout = FNULL, stderr = subprocess.STDOUT)
                        
                        _, qn_normal_alt_all, _, _, _, _, _ = \
                            self.summarizeRefAlt(output + ".tmp.psl", start, end, ref, alt, reads_all)

                        normal_alt_all = len(set(qn_normal_alt_all))
                        
                        
                                             
            if mode == 'pair' and all([x != '---' for x in (tumor_ref, tumor_alt, normal_ref, normal_alt)]):
                _, fisher_pvalue = fisher((
                                           (int(tumor_ref),int(normal_ref)),
                                           (int(tumor_alt),int(normal_alt))
                                          ), alternative='two-sided')
                log10_fisher_pvalue = '{0:.3f}'.format(float(self.math_log_fisher_pvalue(fisher_pvalue)))

            if mode == 'pair':
                print((line +"\t"+ str(tumor_ref)  +"\t"+ str(tumor_alt)  +"\t"+ str(tumor_other)
                    +"\t"+ str(normal_ref) +"\t"+ str(normal_alt) +"\t"+ str(normal_other)
                    +"\t"+ ','.join((str(tumor_ref_mc), str(tumor_alt_mc), str(tumor_other_mc)))
                    +"\t"+ strand_ratio
                    +"\t"+ n_mapqual0 +"\t"+ rate_mapqual0 + '\t' + str(normal_alt_all)
                    +"\t"+ log10_fisher_pvalue
                    +"\t"+ n_repeat +"\t"+ min_n_alt_withindel +"\t"+ min_key_indel
                    +"\t"+ n_inv +"\t"+ n_noinv),
                    file=hResult) 
            else:
                print((line +"\t"+ str(tumor_ref)  +"\t"+ str(tumor_alt)  +"\t"+ str(tumor_other)
                    +"\t"+ ','.join((str(tumor_ref_mc), str(tumor_alt_mc), str(tumor_other_mc)))
                    +"\t"+ strand_ratio +"\t"+ n_mapqual0 +"\t"+ rate_mapqual0
                    +"\t"+ str(beta_01) +"\t"+ str(beta_mid) +"\t"+ str(beta_09)
                    +"\t"+ n_repeat +"\t"+ min_n_alt_withindel +"\t"+ min_key_indel
                    +"\t"+ n_inv +"\t"+ n_noinv), 
                    file=hResult)
            
        ####
        tumor_samfile.close()
        if mode == 'pair':
            normal_samfile.close()
        genome_fa.close()

        hResult.close()
        srcfile.close()

        ####
        if not self.debug_mode:
            if os.path.exists(output + ".tmp.refalt.fa"): os.unlink(output + ".tmp.refalt.fa")
            if os.path.exists(output + ".tmp.refalt.indel1.fa"): os.unlink(output + ".tmp.refalt.indel1.fa")
            if os.path.exists(output + ".tmp.refalt.indel2.fa"): os.unlink(output + ".tmp.refalt.indel2.fa")
            if os.path.exists(output + ".tmp.indel.fa"): os.unlink(output + ".tmp.indel.fa")
            if os.path.exists(output + ".tmp.fa"): os.unlink(output + ".tmp.fa")
            if os.path.exists(output + ".tmp.psl"): os.unlink(output + ".tmp.psl")
            if os.path.exists(output + ".tmp.indel1.psl"): os.unlink(output + ".tmp.indel1.psl")
            if os.path.exists(output + ".tmp.indel2.psl"): os.unlink(output + ".tmp.indel2.psl")

