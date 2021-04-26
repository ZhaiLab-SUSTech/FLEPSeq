# FLEPSeq
The analysis pipeline for FLEP-seq

This toolset can be used to analyze the sequencing data generated either by Nanopore or by PacBio using the FLEP-seq method.

The example scripts for ploting: https://nbviewer.jupyter.org/github/ZhaiLab-SUSTech/FLEPSeq/blob/master/script/FLEP.figure_example.ipynb

Reference: https://www.nature.com/articles/s41477-020-0688-1

Jia, J., Long, Y., Zhang, H. et al. Post-transcriptional splicing of nascent RNA contributes to widespread intron retention in plants. Nat. Plants 6, 780–788 (2020). https://doi.org/10.1038/s41477-020-0688-1

# Software and package requirements

* MinKNOW (MinION software) (https://community.nanoporetech.com/downloads) (required for Nanopore sequencing)
* Guppy v4.0.11 or above (https://community.nanoporetech.com/downloads) (required for Nanopore data basecalling )
* CCS (https://github.com/PacificBiosciences/ccs) (required for PacBio data analysis)
* Lima (https://github.com/PacificBiosciences/barcoding) (required for PacBio data analysis)
* Minimap2 (https://github.com/lh3/minimap2)
* SAMtools (http://www.htslib.org/)
* BLAST+ (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (require for Nanopore data analysis)
* Python 3.7 or above, and following packages:
  * Pysam (https://github.com/pysam-developers/pysam)
  * ont_fast5_api (https://github.com/nanoporetech/ont_fast5_api)
  * pandas (https://pandas.pydata.org/)
  * NumPy (https://numpy.org/)
  * Matplotlib (https://matplotlib.org/)
  * Joblib (https://github.com/joblib/joblib)
  * click (https://click.palletsprojects.com/en/7.x/)
* R 3.5.2 or above, and following packages:
  * Tidyverse (https://www.tidyverse.org/)
  * optparse (https://cran.r-project.org/web/packages/optparse/index.html)

# Step by step workflow

All scripts are deposied in `scripts` directory, and can use `--help` to view the help information.

## Nanopore data analysis pipeline

1.	Nanopore basecalling

You can use MinKNOW to perform real-time basecalling while sequencing, or use the GPU version of Guppy to speed up basecalling after sequencing. Both MinKNOW and Guppy are available via Nanopore community site (https://community.nanoporetech.com). Command-line for running Guppy basecalling is as follow:

```
$ guppy_basecaller -i raw_fast5_dir -s out_fast5_dir -c dna_r9.4.1_450bps_hac.cfg --recursive --fast5_out --disable_pings --qscore_filtering --device "cuda:all:100%"
```

2.	Convert FASTQ files to FASTA format

```
$ python fastqdir2fasta.py --indir path/to/fastq_pass --out all.fasta
```

3.	Use minimap2 to map reads to reference genome

```
$ minimap2 -ax splice --secondary=no genome.fasta all.fasta > tmp.sam
```

CAUTION: For the organisms with short introns, such as Arabidopsis, it might be better to use the parameter “-G” to set the max intron length, for example, “-G 12000”. You also can set “-t number_of_threads” to use more threads to speed up.

```
$ samtools sort -o mapped.bam tmp.sam
$ samtools index mapped.bam
$ rm tmp.sam
```

4.	(Optional) Remove rRNA and tRNA derived reads

```
$ python filter_rRNA_bam.py --inbam mapped.bam --inbed rRNAtRNAetc.bed --out clean.bam
$ samtools index clean.bam
```

5.	Find 3’ adapter in reads

```
$ python adapterFinder.py --inbam clean.bam --inseq all.fasta --out adapter.result.txt --threads 36
```

6.	Identify polyA tail and estimate its length

```
$ python PolyACaller.py --inadapter adapter.result.txt --summary sequencing_summary.txt --fast5dir fast5_dir --out polyA_tail.result.txt --threads 36
```

7.	Extract read information

This pipeline will produce a table containing intron retention information and Pol II position.

```
$ python extract_read_info.py --inbam clean.bam --inbed lib/exon_intron_pos.repr.bed --out read_info.result.txt
```

8.	Merge the above analysis results

```
$ Rscript merge_read_info.R --type Nanopore --inreadinfo read_info.result.txt --inadapter adapter.result.txt --inpolya polyA_tail.result.txt --out read.info.txt
```

9.	Analyze splicing kinetics

```
$ python prepare_data_for_splice_kinetics.py --inreadinfo read.info.txt --inbed lib/exon_intron_pos.repr.txt --out read.intron.pos.splicing.txt
$ Rscript plot_intron_splicing_kinetics.R --inrelpos read.intron.pos.splicing.txt --inreadinfo read.info.txt --inintron lib/select_introns.txt --out read.splicing_kinetics.txt --pdf read.splicing_kinetics.pdf 
```

10.	Calculate intron retention ratio of polyadenylated transcripts

```
$ Rscript cal_polya_transcript_ir.R --inrelpos read.intron.pos.splicing.txt --inreadinfo read.info.txt --outrna mRNA.incompletely_spliced_ratio.txt --outintron intron.unspliced_ratio.txt
```

## PacBio data analysis pipeline
1.	Generate highly accurate single-molecule consensus reads using ccs 

It is the most time-consuming step in PacBio data analysis. For the data we used (~12M subreads), it requires ~30 h to generate consensus sequences (Hifi-reads) on Intel Xeon 6140 CPU at 2.3GHz.

```
$ ccs --num-threads 36 --min-rq 0.9 --report-file ccs.report input.subreads.bam ccs.bam
```

2.	Remove adapter using lima

```
$ echo '>primer_3p\nAAGCAGTGGTATCAACGCAGAGTACATTGATGGTGCCTACAG\n>primer_5p\nAAGCAGTGGTATCAACGCAGAGTACATGGG\n' > primer.fasta
$ lima -j 36 ccs.bam primer.fasta lima.bam --isoseq --peek-guess
$ python lima_bam2fasta.py --infile lima.primer_3p--primer_5p.bam --out all.fasta
```

3.	Use minimap2 to map reads to reference genome, and remove rRNA and tRNA derived reads. 

The same as Step 3-4 in Nanopore pipeline.

4.	Identify polyA tail and estimate its length

```
$ python pacbio_find_polyA.py --inbam clean.bam --inseq all.fasta --out polyA_tail.result.txt
```

5.	Extract read information.

The same as Step 7 in Nanopore pipeline.

6.	Merge the above analysis results

```
$ Rscript merge_read_info.R --type PacBio --inreadinfo read_info.result.txt --inpolya polyA_tail.result.txt --out read.info.txt
```

7.	Analyze splicing kinetics, and calculate intron retention ratio of polyadenylated transcripts.

The same as Step 9-10 in Nanopore pipeline

## Illumina data analysis pipeline

1. Use hisat2 to map reads to reference genome.

```
$ hisat2 -x /lib_path/ath_hisat2_genome -p 20 --min-intronlen 20 --max-intronlen 12000  --dta  --time  -1 R1.fastq.gz -2 R2.fastq.gz  -S hisat2.sam
samtools sort -o hisat2.bam hisat2.sam
samtools index hisat2.bam
rm hisat2.sam
```

2. Use picard to remove PCR duplication reads.

```
$ java -jar /soft_path/picard.jar MarkDuplicates REMOVE_DUPLICATES=true SORTING_COLLECTION_SIZE_RATIO=0.01 I=hisat2.bam O=hisat2.markdump.bam M=hisat2.markdump.txt
samtools index hisat2.markdump.bam
```

3. Calculate the unspliced and spliced ratio of introns.

```
$ python ASCaller.py -i hisat2.markdump.bam -o irratio.txt --file_intron_pos lib/intron_pos.repr.txt --strand_flag 0 --min_overlap 6
```

For strand-specific library, you can set “--strand_flag 1” when the right-most end of the fragment is the first sequenced and set “--strand_flag 2” when the left-most end of fragment is the first sequenced.


# Output

The most important output is read.info.txt generated by `merge_read_info.R`. It is a plain file seprated by Tab including these important columns: 

`read_core_id`, `mRNA`, `retention_introns`, `polya_length`, `type`.

type value:
```
elongating
splicing_intermediate
elongating_3_mapping_low_accuracy
elongating_5lost
polya
polya_3_not_in_last_exon
polya_5lost
```


PS:
All columns in the read.info.txt file (most of columns are only used for script developing)
```
read_core_id            de4d183f-9187-43e5-ab03-56ce63e7de88,chr1,7030,8176

#generate by extract_read_info.py

chr                     chr1         #chromatin name 
read_start              7030         #alignment start 1-based
read_end                8176         #alignment end
read_strand             +            #alignment strand
read_exon_total_num     4            #total number of read exon
mRNA                    AT1G01020.1  #mRNA name       
mRNA_intron_num         8            #total intron number of mRNA
mRNA_start              6788                
mRNA_end                9130                
mRNA_strand             -                   
mRNA_length             2343                
mRNA_pos5               954          #the 5' end position of read relative to the 5' end of mRNA (0-based) 
mRNA_pos3               -242         #the 3' end poistion of read relative to the 3' end of mRNA (0-based)
rel_mRNA_pos5           954          #see below. indicate wrong annotation of mapping if rel_mRNA_pos5 != mRNA_pos5
rel_mRNA_pos3           -242         #see below. indicate wrong annotation of mapping if rel_mRNA_pos3 != mRNA_pos3
total_coverage          797          #total coverage length with mRNA
f_read_exon_num         4            #see below
f_feature_type          intron
f_feature_num           3
f_feature_length        248
f_pos5                  59
f_pos3                  226
l_read_exon_num         1
l_feature_type          exon
l_feature_num           9
l_feature_length        282
l_pos5                  0
l_pos3                  -242
end_type                2            #see below
retention_introns       3:4:6        #retained intron, seprated by ":", the intron number is ordered from 5' mRNA to 3' mRNA
retention_intron_num    3            #the total number of retained intron
span_intron_num         6            #the total number of introns the read overlapping (3-8). The first span intron is f_feature_num, the last span intron is l_feature_num - 1.

#generate by PolyACaller.py or pacbio_find_polyA.py
#polya_start_raw, polya_end_raw, polya_score is only generated for Nanopore data
polya_start_raw     	1185         #The raw signal event index of potential polyA region
polya_end_raw       	1924	     #The raw signal event index of potential polyA region
polya_start_base    	73
polya_end_base      	77
polya_length        	83.34
polya_score         	740
polya_start_raw						
polya_end_raw

#generate by adapterFinder.py (only for Nanopore data)
read_align_strand         +            
rna_strand                -            
read_length               839          # read total length
primer_type	              R-F          
genome_align_start        71
genome_align_end          820
primer_score              1
polyA_type                T            # not shown for now
f_primer_type             R
f_primer_start            1
f_align_end               70
r_primer_type             F
r_primer_start            1
r_align_start             821

#generated by merge_read_info.R
end3_alignment_score      -1           # only for Nanopore data
end_polyA_type			  3	
end5ss_type				  0
type					  elongating
```

some informations used for script developing:

1. read_core_id

Due to that one read may be mapped to multiple positions, we merge read id and mapping information to
represent each alignment:

```
read_id:      7fab6cbf-a9ed-4427-951f-741515ddba0b
read_core_id: 7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687
```

Note: Logically, the read_core_id may be also not uniquely, thus need to remove duplication. (At merge step)

2. read exon

The exon information of each read from bam file based on the 6th column (CIGAR string).
Only consider N (represent intron junction) in CIGAR string, not consider indel.
Split read  to exon by junction.

```
exon1   intron1    exon2   intron2     exon3   intron3    exon4
>>>>> - - - - - - >>>>>>> - - - - - - >>>>>>>> - - - - - >>>>>>>> gene  (>>>>> for exon, - - - - for intron)
    >>>>>>>>>>>>>>>>>>>>>- - - - - -  >>>>>>>>>>>>>>>>>           read  (CIGAR: 110M60N70M)
      read_exon1(110bp)      (60bp)    read_exon2(70bp)
```

3. mRNA_pos5, mRNA_pos3, rel_mRNA_pos5, and rel_mRNA_pos3

We don't care about the alignment direction of read, because we haven't identify the 3' adapter position,  thus we don't know the orignal mRNA strand of this read.

```
gene (length: l)  ----->>>>>>>>>>>>>>>>>>>>>>---------      >>>> repesent gene
                  |    |   |              | |    |
                 -4    0   4             -2 0    5
                  (mRNA pos5)            (mRNA pos3)

    gene          -----<<<<<<<<<<<<<<<<<<<<<<---------
                       |                    |
                  mRNA pos3             mRNA_pos5
```

rel_mRNA_pos3, rel_mRNA_pos5:

```
 >>>> exon,  .... intron
mRNA            >>>>>>>>>>>.......>>>>>>>>>         
read >>>......>>>>>>>>>>>>>.......>>>>>>>>>>>>......>>>>>>>>>>
     |        |                              |               |
     |        |                              |               |
  mRNA_pos5   |                              |          mRNA_pos3
          rel_mRNA_pos5                 rel_mRNA_pos3
```

The inconsistent results (<1% reads in Arabidopsis.) between rel_mRNA_pos and mRNA_pos indicate the wrong mapping or wrong transcription annotation, this pipeline only consider the represent transcripts, it should be removed (in merge step). This is mainly due to alternative splicing or mismapping of polyA You can try set less max intron length during mapping step, but it is not suitable for genome with long introns. And also we can remove the two side mapping read exon based on polyA information.

4. f_feature, l_feature

```
                  exon1        intron1          exon2        intron2       exon3
mRNA           >>>>>>>>>>>>...............>>>>>>>>>>>>>>>...............>>>>>>>>>>>
read                >>>>>>>...............>>>>>>>>>>>>>>>>>>>>>>>
```

f_feature: the 5' end feature overlapping with this read, in this example, it is the first exon: exon1, thus `f_feature_type` is "exon", `f_feature_num` is 1, it is the first read exon of the read overlapping with it, thus `f_read_exon_num` is 1. f_pos5 and f_pos3 is the position of the read exon 1 relative the exon1, the value is like 
mRNA_pos5 and mRNA_pos3.

l_feature: the 3' end feature overlapping with this read, in this example, it is the second intron: intron2, thus `l_feature_type` is "intron", `l_feature_num` is 2, it is the second read exon of the read overlapping with it,  thus `l_read_exon_num` is 2.

5. intron retention

for each intron overlapping with this read:

(a) if it is the f_feature: it would be retented. but this may be affected by 3' SS.
		
(b) if it is the l_feature: the transcription hasn't been finished. not retented.

(c) in the inner: if coverage_ratio >= 0.8: retention.

6. end_type:

```
mRNA                >>>>>>>>>>>.......>>>>>>>>>.........>>>>>>>>>>      
read  end_type: 0                  >>>>>>>>>>>>>>>
      end_type: 1     >>>>>>>>>>>>>>>>>>>>>>
      end_type: 2                  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      end_type: 3   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    #end5: read overlapped with the the first exon
    #end3: read overlapped with the the last exon
    #	     noend5 end5
    #noend3      0    1
    #end3        2    3
```

6. primer/adapter

```
primer_5p = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
primer_3p = "AAGCAGTGGTATCAACGCAGAGTACATTGATGGTGCCTACAG"
```

primer_type: `F` means `primer_5p`, `R` means `primer_3p`

f_ indicate the 5' read, r_ indicate the 3' read

Note: f may be F or R, r may be F or R.

```
genome 5’--------------------------------------------------------3'  
read1      f_primer--------------------------->r_primer
                   ||                          ||
                   1|                          |2
                    3                          4

\b
read2      r_primer<---------------------------f_primer
                   ||                          ||
                   2|                          |1
                    4                          3    
						
1: f_align_end, 2: r_align_start, 3: genome_align_start, 4: genome_align_end 
```                  


coordinate of read: f_align_end, r_align_start, genome_align_start, genome_align_end

coordinate of adapter: f_primer_start, r_primer_start

```
if r_primer_start if 4, the first two base of r_primer is not mapped to read

for example:
                                             AAGGGAGAGAGAG
                                                ||||||||||
read1      f_primer--------------------------->GGAGAGAGAG
```

one end primer type

(a) no alignment
       ["N", ["", "", -1, -1, -1, -1, ""]]
	   
(b) If "F" primer, if primer_start > primer_f_max_start: UF

(c) if "R" primer, if primer_start > primer_r_max_start2: UUR

    if primer_r_max_start2 >= primer_start > primer_r_max_start1: UR

read adapter type based on the primer types of both ends. for example: F-R indicate the 5' primer is F, the 3' primer is F. Each kind of primer type is assigned a score, and also can predict the rna_strand by primer type, The score 1 or 2 is reliable.

7. end_type and end_polyA_type:

```
end_type:
    mRNA                >>>>>>>>>>>.......>>>>>>>>>.........>>>>>>>>>>      
    read  end_type: 0                  >>>>>>>>>>>>>>>
          end_type: 1     >>>>>>>>>>>>>>>>>>>>>>
          end_type: 2                  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
          end_type: 3   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        #end5: read overlapped with the the first exon
        #end3: read overlapped with the the last exon
        #	     noend5 end5
        #noend3      0    1
        #end3        2    3
polyA_type:
      0:  no polyA. polyA_length < 15  
 	   1:  has polyA. polyA_length >= 15
end_polyA_type:
	              end_type=0    end_type=1    end_type=2    end_type=3
    polyA_type = 0         0             1             2             3
    polyA_type = 1         4             5             6             7
```

# generated by IRCaller.py
Output: tab-sperated file. including clolumns:
intron_id	            AT1G01020.1_intron4
chr_name	            chr1
intron_start	        7836      # 1-based
intron_end	            7941
intron_strand	        -
intron_coverage_ratio	1
a	                    126
b	                    231
ab	                    103
c	                    215
o	                    20
t	                    695
iratio	                0.545014521  # unspliced ratio of intron
sratio	                0.416263311  # spliced ratio of intron
oratio	                0.038722168  # the ratio of alternative splicing
o1ratio	                0.013552759
o2ratio	                0.011616651
o1_type	                o||1|:7829-7941
o1_count	            7
o2_type	                o|1|0|:7902-8235
o2_count	            6
other_o	                o|1|0|:7907-8235=3;o|1|1|:7650-8235=2;o||1|:7832-7941=1;o|1||:7836-8235=1


a, b, ab, c, o columns store the count of read with specific read type:
--------...............--------- gene structure
      --...............---       c
      -----............---       o  #alternative splicing
      --------                   a
                   ----------    b
      -------------------------  ab
            -------              i
            
t = a + b + c + o + ab
iratio = (a + b + 2*ab)/(a + b + 2*ab + 2*o + 2*c)
sratio = 2*c/(a + b + 2*ab + 2*o + 2*c)
oratio = 2*o/(a + b + 2*ab + 2*o + 2*c)

o: alternative splicing
The format of alternative splicing type:
    o|x1|x2|x3:start1-end1:< start2-end2 >
    x1, x2 can be one value of "" or "0" or "1"
        "" means splicing at the annotated splicing site
        "0" means splicing downstream the 5' splicing site or upstream the 3' splicing site
        "1" means splicing upstream the 3' splicing site or downstream the 5' splicing site
    x1 for 5' splicing site
    x2 for 3' splicing site
    x3 for exon skipping, can be one value of "" or "e1" or "e2" or ...
    x3 only care about the exon skipping occuring inner the annotated intron. the nubmer after
        "e" means the number of skipped exons.
    start1-end1 means the junction position of the read overlapped with the annotated intron.
    If one read has multiple junction positions, they are sperated by ":".
Considering that some introns have multiple alternative splicing types. The most frequent 
type is marked as o1, the second frequent type is marked as o2. Their type, count and ratio are 
recored in `o1ratio`, `o2_ratio`, `o1_count`, `o2_count`, `o1_type`, `o2_count` columns. The 
inforamtion of all other types are recored in `other_o` column in type3=count3;type4=cout4 format.

Some exampels of alternative splicing type:

--------...............--------- gene structure
    ------........----------     o|0|0|         5' inner, 3' inner
      ------------.....--------- o|0||          5' inner, 3' right
      ------------.........----  o|0|1|         5' inner, 3' outer
 ----......................----  o|1|1|         5' outer, 3' outer

--------...............--------- gene structure
  ------.....--------------      o||0|          5' right, 3' inner
  ------.....-----.....------    o|||e1         exon skipping
  ------.....-----.......-----   o||1|e1        5' right, 3' inner, exon skipping
  ------.....-----               o||0|

  ---........-------------       o|1|0|         5' outer, 3' inner
  ---........-----.....------    o|1||e1
  ---........-----.......-----   o|1|1|e1       5' outer, 3' outer, exon skipping
  ---........-----               o|1|0| 

             -----.....------    o|0|| 
             -----........---    o|0|1|                        

