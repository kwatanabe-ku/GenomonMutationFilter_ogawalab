# -*- coding: utf-8 -*-

import pysam, re, subprocess, os, Levenshtein, math


class cruciform_filter:

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    def __init__(self,
                reference_genome,
                search_window,
                min_matchrate,
                max_depth,
                exclude_sam_flags,
                mapping_quality,
                base_quality,
                removedup_byumi,
                removedup_thres,
                debug_mode
                ):
                
        self.reference_genome      =  reference_genome
        self.search_window         =  search_window
        self.min_matchrate         =  min_matchrate
        self.max_depth             =  max_depth
        self.exclude_sam_flags     =  exclude_sam_flags
        self.mapping_quality       =  mapping_quality
        self.base_quality          =  base_quality
        self.removedup_byumi       =  removedup_byumi
        self.removedup_thres       =  removedup_thres
        self.debug_mode            =  debug_mode
        
    def filter(self, in_tumor_bam, output, in_mutation_file):


        srcfile = open(in_mutation_file,'r')
        hResult = open(output,'w')
        
        tumor_samfile = pysam.Samfile(in_tumor_bam, "rb")
        genome_fa     = pysam.FastaFile(self.reference_genome)
        
        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
            # annovar input file (not zero-based number)
            chr, start, end, ref, alt  = (itemlist[0], (int(itemlist[1]) - 1), int(itemlist[2]), itemlist[3], itemlist[4])
            
            if self.debug_mode:
                 print('{0}:{1}-{2}:{3}>{4}'.format(chr, start, end, ref, alt))
            
            tumor_ref = tumor_alt = '---'

            ### don't check depth, because it is slow
            # if tumor_samfile.count(chr,start,end) < self.max_depth:

            mutreads = self.extractmutread(tumor_samfile, chr, start+1, ref, alt)
            
            if self.removedup_byumi:
                mutreads = self.removedupbyumi_duplex(mutreads)
            
            # if len(softclipreads) > 0:
            if True:
                invreads, noinvreads = self.checkinv(mutreads, genome_fa, chr, start, end, ref, alt, self.search_window)
                
                tumor_alt = str(len(set([x.qname for x in noinvreads])))
                tumor_ref = str(len(set([x.qname for x in invreads]) - 
                                    set([x.qname for x in noinvreads])))
                    
                
            print(line +"\t"+ str(tumor_ref)  +"\t"+ str(tumor_alt), file=hResult)
        
        tumor_samfile.close()
        genome_fa.close()
        srcfile.close()
        hResult.close()
                
    
    def ins_bases(self, query):
        pos = query.query_position
        return query.alignment.query_sequence[ pos+1 : pos+query.indel+1]
                    
    #pos: 1-indexed = genomon output format
    def extractmutread(self, tumor_samfile, chr, pos, ref, alt):
            
        if alt == '-': #deletion
            pos_0ind = pos - 2
        else:
            pos_0ind = pos - 1
            
        try:
            pileups = [x.pileups for x in tumor_samfile.pileup(
                                            chr, pos_0ind, pos_0ind+1,
                                            truncate=True,
                                            ignore_overlaps=False,
                                            max_depth=self.max_depth,
                                            flag_filter=self.exclude_sam_flags,
                                            flag_require=2,
                                            min_base_quality=self.base_quality,
                                            min_mapping_quality=self.mapping_quality)][0]
        except IndexError:
            return []
            
        #SNV
        if ref != '-' and alt != '-':
            reads  = [x.alignment for x in pileups
                            if x.query_position is not None and 
                               (x.alignment.flag & 2) > 0 and 
                               (x.alignment.flag & self.exclude_sam_flags) == 0 and 
                               x.alignment.query_qualities[x.query_position] >= self.base_quality and
                               x.alignment.query_sequence[x.query_position] == alt]
        #insertion
        elif ref == '-':
            reads  = [x.alignment for x in pileups
                            if x.query_position is not None and 
                               (x.alignment.flag & 2) > 0 and
                               (x.alignment.flag & self.exclude_sam_flags) == 0 and 
                               x.indel == len(alt) and
                               self.ins_bases(x) == alt]
        #deletion
        else:
            reads  = [x.alignment for x in pileups
                            if x.query_position is not None and 
                               (x.alignment.flag & 2) > 0 and
                               (x.alignment.flag & self.exclude_sam_flags) == 0 and 
                               x.indel == -len(ref)]

        return reads
        
    #pos: 1-indexed = genomon output format
    def softclipread(self, mutreads, pos, thres_distance_softclip):
    
        softclipreads = []
        
        for read in mutreads:
            softclip = re.search(r'\A([0-9]+)S', read.cigarstring)
            if softclip is not None: #left soft-clip
                if pos-read.reference_start-1 <= thres_distance_softclip:
                    softclipreads.append(read)
                    continue
            
            softclip = re.search(r'([0-9]+)S\Z', read.cigarstring)
            if softclip is not None: #right soft-clip
                if read.reference_end-pos <= thres_distance_softclip:
                    softclipreads.append(read)
                
        return softclipreads
        
    def removedupbyumi_duplex(self, reads):
        if len(reads) == 0 or not reads[0].has_tag('RX'): return reads
           
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
        
    def getterminalseq(self, read, pos_0ind, ref, alt):
        pos_seq = read.get_reference_positions(full_length=True).index(pos_0ind)
        pos_seq_align = read.get_reference_positions().index(pos_0ind)
        
        #left clip
        if pos_seq_align < read.reference_length /2:
            if alt == '-': #deletion
                q_seq = read.query_sequence[:pos_seq+1]
            elif ref == '-': #insertion
                q_seq = read.query_sequence[:pos_seq+1+len(alt)]
            else: #SNV
                q_seq = read.query_sequence[:pos_seq+len(alt)]
        #right clip
        else:
            if alt == '-': #deletion
                q_seq = read.query_sequence[pos_seq+1:]
            elif ref == '-': #insertion
                q_seq = read.query_sequence[pos_seq+1:]
            else: #SNV
                q_seq = read.query_sequence[pos_seq:]
    
        return q_seq
                    
    def fuzzy_search(self, query, text):
        query_len = len(query)
        matches = []
        for i in range(len(text) - query_len + 1):
            substring = text[i:i + query_len]
            # indel is treated as ins + del
            # weights muts be integer
            distance = Levenshtein.distance(query, substring, weights=(1,1,2))/2\
                        - query.count('N') - substring.count('N')
            # if distance <= self.max_error:
            if distance <= math.ceil((1-self.min_matchrate) * query_len):
                matches.append((i, substring, distance))
        return matches
                    
    #SNPを考慮していない
    def checkinv(self, mutreads, genome_fa, chr, start, end, ref, alt, search_window):
        
        if alt == '-': #deletion
            pos_0ind = start - 1
        else:
            pos_0ind = start
        
        ref_seq = genome_fa.fetch(chr, max(0, start-search_window), end+search_window)
        ref_seq = "".join(self.complement.get(base) for base in reversed(str(ref_seq)))
    
        if self.debug_mode: print('inverse sequence: ' + ref_seq)
            
        invreads   = []
        noinvreads = []
        for read in mutreads:
        
            q_seq = self.getterminalseq(read, pos_0ind, ref, alt)
            
            # Using regex: N base -> "."
            reobj = re.search(q_seq.replace('N','.'), ref_seq)
        
            if reobj is not None:
                invreads.append(read)
                if self.debug_mode:
                    print('inverse. readname: {}. query_sequence: {}'.format(read.qname, q_seq))
            else:
                noinvreads.append(read)
                if self.debug_mode:
                    print('not inverse. readname: {}. query_sequence: {}'.format(read.qname, q_seq))
                    
        if len(noinvreads) > 0:
        
            tempreads  = noinvreads
            noinvreads = [] 
            
            for k, read in enumerate(tempreads):
                tempSeq = self.getterminalseq(read, pos_0ind, ref, alt)
                
                matches = self.fuzzy_search(tempSeq, ref_seq)
                
                if len(matches) > 0:
                    invreads.append(read)
                else:
                    noinvreads.append(read)
        
        return invreads, noinvreads
        

def run_cruciform_filter(arg):

    crucif = cruciform_filter(
                    arg.reference_genome,
                    arg.search_window,
                    arg.distance_to_softclip,
                    # arg.max_error,
                    arg.min_matchrate,
                    arg.max_depth,
                    arg.exclude_sam_flags,
                    arg.mapping_quality,
                    arg.base_quality,
                    arg.removedup_byumi,
                    arg.removedup_thres,
                    arg.debug_mode
                )
    crucif.filter(arg.bam, arg.output, arg.target_mutation_file)
