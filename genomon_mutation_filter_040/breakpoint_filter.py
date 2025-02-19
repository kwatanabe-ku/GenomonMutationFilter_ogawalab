import pysam
import sys
import os
import re
import logging

#
# Class definitions
#
class Breakpoint_filter:

    def __init__(self, max_depth, min_clip_size, junc_num_thres, mapq_thres, header_flag, exclude_sam_flags):
        self.max_depth = max_depth
        self.min_clip_size = min_clip_size
        self.junc_num_thres = junc_num_thres
        self.mapq_thres = mapq_thres
        self.header_flag = header_flag
        self.exclude_sam_flags = exclude_sam_flags


    def write_result_file(self, line, file_handle, dist, max_junc_cnt):
        print(line +"\t"+ str(max_junc_cnt) +"\t"+str(dist), file=file_handle) 
        
        
    def filter_main(self, chr, start, end, samfile):

        if samfile.count(chr, start, (start+1)) >= self.max_depth:
            return ('---', '---')

        bp_dict = {}
        max_junc_pos = int(0)
        max_junc_cnt_p = int(0)
        max_junc_cnt_m = int(0)
            
        ####
        for read in samfile.fetch(chr, max(0, int(start-(self.min_clip_size))), int(end+(self.min_clip_size))):

            # get the flag information
            read_flag = int(read.flag)

            if 0 != int(bin(self.exclude_sam_flags & read_flag),2): continue

            flags = format(read_flag, "#014b")[:1:-1]

            # no clipping
            if len(read.cigar) == 1: continue
                
            # skip low mapping quality
            if read.mapq < self.mapq_thres: continue

            #####Fix:v0.2.1 Only soft-clip is considered.
            #In prism pipeline, hard-clip simply means trimmed adaptor.
            if read.cigar[0][0] == 4:
                left_clipping = read.cigar[0][1]
            elif read.cigar[0][0] == 5 and read.cigar[1][0] == 4:
                left_clipping = read.cigar[1][1]
            else:
                left_clipping = 0
                
            if read.cigar[-1][0] == 4:
                right_clipping = read.cigar[-1][1]
            elif read.cigar[-1][0] == 5 and read.cigar[-2][0] == 4:
                right_clipping = read.cigar[-2][1]
            else:
                right_clipping = 0

            # get strand info
            strand = "-" if flags[4] == "1" else "+"

            # when the right side is clipped...
            if right_clipping >= self.min_clip_size:
                    
                juncPos_current = str(int(read.pos + 1) + read.alen - 1)
                key = juncPos_current +"\tR"
        
                if key not in bp_dict:
                    bp_dict[key] = {"+":0, "-":0 }
                bp_dict[key][strand] += 1

            # when the left side is clipped...
            if left_clipping >= self.min_clip_size:

                juncPos_current = str(int(read.pos + 1))
                juncPos_current = str(read.pos)
                key = juncPos_current +"\tL"

                if key not in bp_dict:
                    bp_dict[key] = {"+":0, "-":0 }
                bp_dict[key][strand] += 1

        dist = 0
        # if start == 4680207: 
        #    print(bp_dict)

        for key in bp_dict:
            juncPos_current = (key.split("\t")[0])
            # if start == 4680207: 
            #     print(juncPos_current)
            #     print(bp_dict[key]["+"] + bp_dict[key]["-"])

            if ((int(juncPos_current) - self.min_clip_size) <= start <= int(juncPos_current)
            or (start <= int(juncPos_current) <= end)
            or (int(juncPos_current) <= end <= int(juncPos_current) + self.min_clip_size)):

                tmp_sdist = abs(start - int(juncPos_current))
                tmp_edist = abs(start - int(juncPos_current))
                tmp_dist = tmp_sdist if tmp_sdist < tmp_edist else tmp_edist

                if (bp_dict[key]["+"] + bp_dict[key]["-"]) > (max_junc_cnt_p + max_junc_cnt_m) or ((bp_dict[key]["+"] + bp_dict[key]["-"]) == (max_junc_cnt_p + max_junc_cnt_m) and (tmp_dist < dist) ):
                    max_junc_cnt_p = bp_dict[key]["+"]
                    max_junc_cnt_m = bp_dict[key]["-"]
                    max_junc_pos = juncPos_current

                    sdist = abs(start - int(max_junc_pos))
                    edist = abs(end   - int(max_junc_pos))
                    dist = sdist if sdist < edist else edist
                    
                    if (max_junc_cnt_p + max_junc_cnt_m) == 0:
                        dist = 0

        return (dist, (max_junc_cnt_p + max_junc_cnt_m))

                
    def filter(self, in_mutation_file, in_bam, output):
   
        samfile = pysam.AlignmentFile(in_bam, "rb")

        srcfile = open(in_mutation_file,'r')
        hResult = open(output,'w')
        if self.header_flag:
            header = srcfile.readline().rstrip('\n')  
            newheader = "bp_mismatch_count\tdistance_from_breakpoint"
            print(header +"\t"+ newheader, file=hResult)
        
        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
    
            # input file is annovar format (not zero-based number)
            chr = itemlist[0]
            start = (int(itemlist[1]) - 1)
            end = int(itemlist[2])

            dist = "---"
            junction_num = "---"
            if samfile.count(chr, start, (start+1)) < self.max_depth:
                dist, junction_num = self.filter_main(chr, start, end, samfile)

            ####
            if junction_num =="---" or  junction_num >= self.junc_num_thres:
                self.write_result_file(line, hResult, dist, junction_num)

        ####
        hResult.close()
        srcfile.close()
