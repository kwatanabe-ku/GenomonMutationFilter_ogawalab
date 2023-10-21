# GenomonMutationFilter

GenomonMutationFilter is a software package for filtering false poistive somatic mutations from cancer genome sequencing data.
This version is modified for use in ogawa-lab to remove artifacts caused by cruciform DNA and to improve the sensitivity of base insertion mutation detection.

## Dependency

Python (>= 3.7), pysam packages

You can choose realignment using blat.

* [blat](http://genome.ucsc.edu/)

## Preparation

  **target somatic mutation candidats**: the somatic mutation candidates (should be .tsv format).  
  **target tumor bam**: the indexed bam file of the target tumor sample.  
  **target normal bam**: the indexed bam file of the target normal sample.  


## Run

```
usage: mutfilter_0.4.0 realignment [-h] -t TARGET_MUTATION_FILE -1 BAM1 [-2 BAM2]
                             -o OUTPUT -r REF_GENOME
                             [-b BLAT_PATH] [-m tumor_min_mismatch]
                             [-M normal_max_mismatch] [-s score_difference]
                             [-w window_size] [-d max_depth]
                             [-F exclude_sam_flags] [--header]
                             [--mapping_quality mapping_quality]
                             [--base_quality base_quality]
                             [--mutation_cluster mutation_cluster]
                             [--removedup_byumi] [--removedup_duplex]
							 [--removedup_thres removedup_thres]
                             [--cruciform_filter]
							 [--cruciform_search_window cruciform_search_window]
							 [--search_length_indel search_length_indel]
							 [--debug]
```

```
usage: mutfilter_0.4.0 indel [-h] -t TARGET_MUTATION_FILE -2 BAM2 [-A SAMPLE1]
                       [-B SAMPLE2] -o OUTPUT [-l search_length] [-n neighbor]
                       [-d min_depth] [-m min_mismatch]
                       [-a allele_frequency_thres] [--header] -s SAMTOOLS_PATH
                       [-S SAMTOOLS_PARAMS] [-O {vcf,anno}] [-r REF_GENOME]
```

```
usage: mutfilter_0.4.0 breakpoint [-h] -t TARGET_MUTATION_FILE -2 BAM2 [-A SAMPLE1]
                            [-B SAMPLE2] -o OUTPUT [-d max_depth]
                            [-c min_clip_size] [-j junc_num_thres]
                            [-m mapping_quality_thres] [-F exclude_sam_flags]
                            [--header] [-O {vcf,anno}] [-r REF_GENOME]
```

```
usage: mutfilter_0.4.0 simplerepeat [-h] -t TARGET_MUTATION_FILE -o OUTPUT -S
                              SIMPLE_REPEAT_DB [--header] [-O {vcf,anno}]
```
