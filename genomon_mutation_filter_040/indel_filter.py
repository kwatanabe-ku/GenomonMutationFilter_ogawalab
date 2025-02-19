import sys
import os
import re
import pysam
import argparse
import logging
import subprocess
import auto_vivification as autov

#
# Class definitions
#
class Indel_filter:

    ############################################################
    def __init__(self, search_length, min_depth, min_mismatch, af_thres, neighbor, header_flag, samtools_path, samtools_params):
        self.search_length = search_length
        self.min_depth = min_depth
        self.min_mismatch = min_mismatch
        self.af_thres = af_thres
        self.neighbor = neighbor
        self.target = re.compile( '([\+\-])([0-9]+)([ACGTNRMacgtnrm]+)' )
        self.remove_chr = re.compile( '\^.' )
        self.header_flag = header_flag
        self.samtools_path = samtools_path
        self.samtools_params = samtools_params
   
 
    ############################################################
    def filter_main(self, chr, start, end, ref, alt, in_bam):

        max_mismatch_count = 0
        max_mismatch_rate = 0
        region = chr +":"+str(max(0, (start - self.search_length + 1))) +"-"+ str(end + self.search_length)

        cmd_list = [self.samtools_path,'mpileup']
        cmd_list.extend(self.samtools_params.split(" "))
        cmd_list.extend(['-r', region, in_bam])
        
        FNULL = open(os.devnull, 'w')

        ####
        # print region
        pileup = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr = FNULL)
        end_of_pipe = pileup.stdout
        for mpileup in end_of_pipe:
            # print mpileup.rstrip()

            #
            # Prepare mpileup data
            #
            # mp_list = str( mpileup.translate( None, '\n' ) ).split( '\t' )
            if sys.version_info.major == 3:
                mp_list = mpileup.decode().strip('\n').split( '\t' )
            else:
                mp_list = mpileup.strip('\n').split( '\t' )
            mp_list_len = len( mp_list )
            coordinate = mp_list[ 0:3 ]

            #
            # skip if depth is 0
            #
            if mp_list[ 3 ] == '0' or ( mp_list_len > 6 and mp_list[ 6 ] == '0' ):
                continue
                
            #
            # skip if depth < min_depth
            #
            if int(mp_list[ 3 ]) < self.min_depth:
                continue

            #
            depth = mp_list[ 3 ]
            read_bases = mp_list[ 4 ]
            qual_list = mp_list[ 5 ]

            #
            # position id,
            # mpileup output 4th row(number of read covering the site),
            # 5th row(read bases),
            # 6th row(base quality)
            #
            indel = autov.auto_vivification()

            #
            # Look for deletion/insertion and save info in 'indel' dictionary
            #
            #   ([\+\-])[0-9]+[ACGTNacgtn]+
            #
            # m.group(1): + or - (deletion/insertion)
            # m.group(2): number of deletion/insertion
            # m.group(3): nucleotides
            #
            deleted = 0
            iter = self.target.finditer( read_bases )
            for m in iter:
                site = m.start()
                type = m.group( 1 )
                num = m.group( 2 )
                bases = m.group( 3 )[ 0:int( num ) ]
                if bases.islower():
                    strand = ( '-', '+' )
                else:
                    strand = ( '+', '-' )

                key = '\t'.join( coordinate + [ bases.upper() ] )
                if type in indel and key in indel[ type ]:
                    indel[ type ][ key ][ strand[ 0 ] ] += 1
                else:
                    indel[ type ][ key ][ strand[ 0 ] ] = 1
                    indel[ type ][ key ][ strand[ 1 ] ] = 0

                read_bases = read_bases[ 0:site - deleted ] + read_bases[ site + int(num) + len( num ) + 1 - deleted: ]
                deleted += 1 + len( num ) + int( num )

            #
            # Remove '^.' and '$'
            #
            read_bases = self.remove_chr.sub( '', read_bases )
            read_bases = read_bases.replace( '$', '' )

            #
            # Error check
            #
            if len( read_bases ) != len( qual_list ):
                logging.error( "mpileup data is not good: {0}, {1}".format( mpileup, read_bases ) )
                exit(1)

            #
            for type in ( '+', '-' ):
                if type in indel:
                    for key in indel[ type ].keys():
                        start_pos = mp_list[ 1 ]
                           
                        mismatch_count = ( indel[ type ][ key ][ '-' ] + indel[ type ][ key ][ '+' ])
                        mismatch_rate = (float(mismatch_count) / float(depth))

                        if mismatch_rate >= max_mismatch_rate:
                            start_pos = int(start_pos)
                            end_pos   = int(start_pos)

                            if (type == '-'):
                                start_pos = int(start_pos) + 1
                                end_pos = int(start_pos) + len((key.split('\t'))[3]) - 1 

                            if ((start_pos - self.neighbor <= int(start) + 1 and int(start) + 1 <= self.neighbor + end_pos) 
                              or(start_pos - self.neighbor <= int(end)      and  int(end)       <= self.neighbor + end_pos)): 

                                max_mismatch_count = mismatch_count
                                max_mismatch_rate  = mismatch_rate
                                
        pileup.stdout.close()
        pileup.wait()
        FNULL.close()
        return (max_mismatch_count, max_mismatch_rate)


    ############################################################
    def filter_main_anno(self, in_mutation_file, in_bam, output, thread_idx, is_multi_thread):

        srcfile = open(in_mutation_file,'r')
        if is_multi_thread:
            hResult = open(output,'w')
        else:
            hResult = open(output,'a')

        if self.header_flag and thread_idx == 1:
            # skip header
            srcfile.readline()

        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
    
            # input file is annovar format (not zero-based number)
            chr, start, end, ref, alt  = (itemlist[0], (int(itemlist[1]) - 1), int(itemlist[2]), itemlist[3], itemlist[4])

            max_mismatch_count, max_mismatch_rate = self.filter_main(chr, start, end, ref, alt, in_bam)            
            
            ####
            if(max_mismatch_count <= self.min_mismatch or max_mismatch_rate <= self.af_thres):
                print(line +"\t"+ str(max_mismatch_count) +"\t"+ str('{0:.3f}'.format(float(max_mismatch_rate))), file=hResult) 

        ####
        hResult.close()
        srcfile.close()


    ############################################################
    def Print_header(self, in_mutation_file, hResult):
        with open(in_mutation_file,'r') as srcfile:
            header = srcfile.readline().rstrip('\n')  
            newheader = ("indel_mismatch_count\tindel_mismatch_rate")
            print(header +"\t"+ newheader, file=hResult)
            

    ############################################################
    def filter(self, in_mutation_file, in_bam, output):
    
        with open(output, 'w') as w:
            if self.header_flag:
                self.Print_header(in_mutation_file, w)
        self.filter_main_anno(in_mutation_file, in_bam, output, 1, False)
        
        
        for idx in range(1, 2): 
            if os.path.exists(in_mutation_file +"."+str(idx)): os.unlink(in_mutation_file +"."+str(idx))
            if os.path.exists(output +"."+str(idx)): os.unlink(output +"."+str(idx))
        