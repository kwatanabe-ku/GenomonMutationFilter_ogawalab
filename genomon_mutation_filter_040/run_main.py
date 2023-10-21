# -*- coding: utf-8 -*-

import sys
import os
import math
import argparse
import logging
from realignment_filter import realignment_filter
from indel_filter import Indel_filter
from breakpoint_filter import Breakpoint_filter
from simple_repeat_filter import Simple_repeat_filter


#
# Main
#
def run_realignment_filter(arg):

    logging.info( 'realignment filter start')
    #####Fix:v0.3.0 20221114: removedup_byumi, removedup_thresを追加
    #####Fix:v0.3.2 20230104: debugを追加
    realignf = realignment_filter(
                    arg.ref_genome,
                    arg.tumor_min_mismatch,
                    arg.normal_max_mismatch,
                    arg.window_size,
                    arg.score_difference,
                    arg.blat_path,
                    arg.header_flag,
                    arg.max_depth,
                    arg.exclude_sam_flags,
                    arg.mapping_quality,
                    arg.base_quality,
                    arg.mutation_cluster,
                    arg.removedup_byumi,
                    arg.removedup_duplex,
                    arg.removedup_thres,
                    arg.cruciform_filter,
                    arg.cruciform_search_window,
                    arg.cruciform_min_matchrate,
                    arg.search_length_indel,
                    arg.debug
                )
    realignf.filter(arg.bam1, arg.bam2, arg.output, arg.target_mutation_file)
    logging.info( 'realignment filter end')


def run_indel_filter(arg):

    logging.info( 'indel filter start')
    indelf = Indel_filter(
                    arg.search_length,
                    arg.min_depth,
                    arg.min_mismatch,
                    arg.af_thres,
                    arg.neighbor,
                    arg.header_flag,
                    arg.samtools_path,
                    arg.samtools_params
                )
    indelf.filter(arg.target_mutation_file, arg.bam2, arg.output)
    logging.info( 'indel filter end')


def run_breakpoint_filter(arg):

    logging.info( 'breakpoint filter start')
    bpf = Breakpoint_filter(
                    arg.max_depth,
                    arg.min_clip_size,
                    arg.junc_num_thres,
                    arg.mapq_thres,
                    arg.header_flag,
                    arg.exclude_sam_flags
                )
    bpf.filter(arg.target_mutation_file, arg.bam2, arg.output)
    logging.info( 'breakpoint filter end')


def run_simple_repeat_filter(arg):

    logging.info( 'simple repeat filter start')
    simplef = Simple_repeat_filter(
                    arg.simple_repeat_db,
                    arg.header_flag
                )
    simplef.filter(arg.target_mutation_file, arg.output)
    logging.info( 'simple repeat filter end')
