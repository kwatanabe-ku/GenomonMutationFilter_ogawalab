import sys
import os
import re
import logging
import pysam

#
# Class definitions
#
class Simple_repeat_filter:


    def __init__(self, simple_repeat_db, header_flag):
        self.simple_repeat_db = simple_repeat_db
        self.header_flag = header_flag


    def write_result_file(self, line, file_handle, db_pos, db_seq):
        print(line + "\t" +db_pos+ "\t" +db_seq, file=file_handle)


    def filter_main(self, chr, start, end, tb):

        chridx = chr.find('chr')
        if chridx < 0:
            chr = 'chr' + str(chr)
        # tabix databese is a zero-based number 
        position_array = []
        sequence_array = []
        try:
            records = tb.fetch(chr, start, end)
            for record_line in records:
                record = record_line.split('\t')
                position_array.append(record[0] +":"+ record[1] +"-"+ record[2])
                sequence_array.append(record[15])
                    
        except Exception:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            logging.error( ("{0}: {1}:{2}".format( exc_type, fname, exc_tb.tb_lineno) ) )
            return ("","")
          
        ####
        db_pos = ";".join(map(str,position_array[0:5]))
        db_seq = ";".join(map(str,sequence_array[0:5]))
        return (db_pos, db_seq)


    ############################################################
    def filter(self, in_mutation_file, output):
    
        tb = pysam.TabixFile(self.simple_repeat_db)

        # tabix open
        srcfile = open(in_mutation_file,'r')
        hResult = open(output,'w')
        if self.header_flag:
            header = srcfile.readline().rstrip('\n')  
            newheader = "simple_repeat_pos\tsimple_repeat_seq"
            print(header +"\t"+ newheader, file=hResult)

        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
    
            # input file is annovar format (not zero-based number)
            chr = itemlist[0]
            start = (int(itemlist[1]) - 1)
            end = int(itemlist[2])

            db_pos, dbseq = self.filter_main(chr,start,end,tb)
            if db_pos == "": db_pos = "---"
            if dbseq == "": dbseq = "---"
            self.write_result_file(line, hResult, db_pos, dbseq)

        ####
        hResult.close()
        srcfile.close()
