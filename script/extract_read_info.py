import pysam
import sys
import collections
import click

"""
wirte by Jia Jinbu 2020.09.12

Assign the read to mRNAs via comparing the read exon position with the feature (exon/intron)
position of represent transcripts. And extract the first (5') and last (3') feature of mRNA
overlappign with the read, and the intron retention inforamtion.

See more by --help.
"""

@click.command()
@click.option('-i', '--inbam', help='Input bam file', 
                    required=True, type=click.Path(exists=True))
@click.option('-b', '--inbed', help="""Input bed file. The positions of exon and intron of represent transcripts.
                                       For example:


                        \b
                         chr1    3630    3913    AT1G01010.1_exon1       .       +
                         chr1    3913    3995    AT1G01010.1_intron1     .       +
                         chr1    3995    4276    AT1G01010.1_exon2       .       +
                         chr1    4276    4485    AT1G01010.1_intron2     .       +
                         chr1    4485    4605    AT1G01010.1_exon3       .       +
 
                         Note: the script will use "AT1G01010.1_exon1".split("_")[0] to extract the mRNA 
                         name, and use "AT1G01010.1_exon1".split("_")[1].startwith("e) check whether this 
                         feature is exon or intron. The number 1, 2, is the ordered from 5' mRNA to 3' mRNA.
                                     """, 
                    required=True, type=click.Path(exists=True))                    
@click.option('-o', '--out', help="""Output read info file. Tab-seperated.

                    \b
                    read_core_id            de4d183f-9187-43e5-ab03-56ce63e7de88,chr1,7030,8176
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
                    """, 
                    required=True)
def main(inbam, inbed, out):
    """
    Assign the read to mRNAs via comparing the read exon position with the feature (exon/intron)
    position of represent transcripts. And extract the first (5') and last (3') feature of mRNA
    overlappign with the read, and the intron retention inforamtion.
    
    1. Frist extract read exon information:
    Extract the exon information of each read from bam file based on the 6th column (CIGAR string).
    Only consider N (represent intron junction) in CIGAR string, not consider indel.
    Split read  to exon by junction.

    \b
    exon1   intron1    exon2   intron2     exon3   intron3    exon4
    >>>>> - - - - - - >>>>>>> - - - - - - >>>>>>>> - - - - - >>>>>>>> gene  (>>>>> for exon, - - - - for intron)
        >>>>>>>>>>>>>>>>>>>>>- - - - - -  >>>>>>>>>>>>>>>>>           read  (CIGAR: 110M60N70M)
          read_exon1(110bp)      (60bp)    read_exon2(70bp)

    Output:
    1-based.
    
    Due to that one read may be mapped to multiple positions, we merge read id and mapping information to
    represent each alignment: 
    
    \b
    read_core_id: 7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687 #7fab6cbf-a9ed-4427-951f-741515ddba0b is read_id
    read_name: 7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687,2  #2 is total read exon number of this read
    read_exon_name: 7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687,2,exon1
    read_exon_name: 7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687,2,exon2
    exon1 and exon2 ordered by chromation position. not consider alignment direction.

    Note: Logically, the read_core_id may be also not uniquely, thus need to remove duplication.

    2. Compare read exon position with intron exon position of represent transcripts.

    3. Calculate the coverage of a mRNA (included introns), if a read overlapped with multiple mRNAs, select the 
    mRNA with maxmium coverage with the read. 
    \b
    4. Extract the frist (5') and last (3') feature (exon/intron) of mRNA overlapping with the read. Extract 
    intron retention information of read. 
    
    mRNA_pos5, mRNA_pos3:
        We don't care about the alignment direction of read, because we haven't identify the 3' adapter position,  thus
        we don't know the orignal mRNA strand of this read.
    
    \b
        gene (length: l)  ----->>>>>>>>>>>>>>>>>>>>>>---------      >>>> repesent gene
                          |    |   |              | |    |
                         -4    0   4             -2 0    5
                          (mRNA pos5)            (mRNA pos3)

            gene          -----<<<<<<<<<<<<<<<<<<<<<<---------
                               |                    |
                          mRNA pos3             mRNA_pos5
    
    rel_mRNA_pos3, rel_mRNA_pos5:
        The inconsistent results (<1% reads in Arabidopsis.) between rel_mRNA_pos and mRNA_pos indicate the wrong mapping 
        or wrong transcription annotation, this pipeline only consider the represent transcripts, it should be removed 
        (in merge step). This is mainly due to alternative splicing or mismapping of polyA You can try set less max intron 
        length during mapping step, but it is not suitable for genome with long introns. And also we can 
        remove the two side mapping read exon based on polyA information.
    
    \b
         >>>> exon,  .... intron
        mRNA            >>>>>>>>>>>.......>>>>>>>>>         
        read >>>......>>>>>>>>>>>>>.......>>>>>>>>>>>>......>>>>>>>>>>
             |        |                              |               |
             |        |                              |               |
          mRNA_pos5   |                              |          mRNA_pos3
                  rel_mRNA_pos5                 rel_mRNA_pos3

    \b
    f_feature, l_feature: 
                          exon1        intron1          exon2        intron2       exon3
        mRNA           >>>>>>>>>>>>...............>>>>>>>>>>>>>>>...............>>>>>>>>>>>
        read                >>>>>>>...............>>>>>>>>>>>>>>>>>>>>>>>
        f_feature: the 5' end feature overlapping with this read, in this example, it is 
        the first exon: exon1, thus `f_feature_type` is "exon", `f_feature_num` is 1, it is the 
        first read exon of the read overlapping with it, thus `f_read_exon_num` is 1.
        f_pos5 and f_pos3 is the position of the read exon 1 relative the exon1, the value is like 
        mRNA_pos5 and mRNA_pos3.
        l_feature: the 3' end feature overlapping with this read, in this example, it is the 
        second intron: intron2, thus `l_feature_type` is "intron", `l_feature_num` is 2, it is the
        second read exon of the read overlapping with it,  thus `l_read_exon_num` is 2.
        

    \b
    intron retention:
        for each intron overlapping with this read:
        (a) if it is the f_feature: it would be retented. but this may be affected by 
                3' SS.
        (b) if it is the l_feature: the transcription hasn't been finished. not retented.
        (c) in the inner: if coverage_ratio >= 0.8: retention.
    
    \b
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
    """
    with open(out, 'w') as o:
        header = ("read_core_id\tchr\tread_start\tread_end\tread_strand\t"
                  "read_exon_total_num\tmRNA\tmRNA_intron_num\tmRNA_start\t"
                  "mRNA_end\tmRNA_strand\tmRNA_length\tmRNA_pos5\tmRNA_pos3\t"
                  "rel_mRNA_pos5\trel_mRNA_pos3\ttotal_coverage\tf_read_exon_num\tf_feature_type\t"
                  "f_feature_num\t f_feature_length\tf_pos5\tf_pos3\tl_read_exon_num\t"
                  "l_feature_type\tl_feature_num\tl_feature_length\tl_pos5\tl_pos3\tend_type\t"
                  "retention_introns\tretention_intron_num\tspan_intron_num\n")
        o.write(header)
        for d in extract_read_information(inbam, inbed):
            output = "\t".join([str(i) for i in d]) + "\n"
            o.write(output)

def revcom(seq):
    """
    !!!The function in included in both adapterFinder.py and 
    pacbio_find_polyA.py and extract_read_info.py. 
    They are same function, but haven't be put in a module to 
    keep each script can be run independently. If you want to 
    modify one of them, please modify them at the same time.
    
    Return the reverse complement sequence of origin sequence.
    
    The origin sequence can be in uppercase and lowercase letters.
    But the output is in uppercase letters.
    All letters not in `ATCGatcg` will be converted to `N`. 
    """
    def complement(seq):
        seq = seq.upper()
        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','N':'N'}
        def _com(base):
            try:
                return basecomplement[base]
            except:
                return "N"
        letters = list(seq)
        letters = [_com(base) for base in letters]
        return ''.join(letters)
            
    return complement(seq[::-1])
    
def get_read_blocks(read):
    '''
    Input:
    read: the alignment object generated by pysam
    
    Return:
    a list of read exon poistion [[start, end], [start, end]], 1-based
    ordered by read exon start position (on chroamtin).
    
    Skip `I` (indel) and `D` (deletion).
    Only consider `N` (skipped region) to split read to read exon.
    
    M	BAM_CMATCH	0
    I	BAM_CINS	1
    D	BAM_CDEL	2
    N	BAM_CREF_SKIP	3
    S	BAM_CSOFT_CLIP	4
    H	BAM_CHARD_CLIP	5
    P	BAM_CPAD	6
    =	BAM_CEQUAL	7
    X	BAM_CDIFF	8
    B	BAM_CBACK	9
    '''
    blocks = []
    start = read.reference_start + 1
    end = start - 1
    for (flag, length) in read.cigartuples:
        if flag == 4:
            continue
        elif flag == 5:
            continue
        elif flag == 0:
            end += length
        elif flag == 1:
            continue
        elif flag == 2:
            end += length
        elif flag == 3:
            blocks.append([start, end])
            start = end + length + 1
            end = start - 1
    blocks.append([start, end])
    return blocks

def iter_all_reads_blocks(file_bam):
    """
    Input:
    file_bam 
    
    Return:
    a iterator, each element is a tuple (read_core_id, read_exon_blocks):
    read_exon_blocks is the list of the read exons of this read, the element is a list:
    [read_exon_name, chromatin_name, start, end, read_alignment_strand]
    
    read_core_id: 3217606a-2c2c-4598-9929-8ff20c1bf938,chr1,3618,4603
    read_name: 3217606a-2c2c-4598-9929-8ff20c1bf938,chr1,3618,4603,+,3
    read_exon_name: 3217606a-2c2c-4598-9929-8ff20c1bf938,chr1,3618,4603,+,3,1
    start, end is 1-based
    
    see detail in `get_read_blocks`
    """
    bam_obj = pysam.AlignmentFile(file_bam, "rb")
    have_extract_read_core_ids = set()
    for read in bam_obj.fetch():
        #reference_id, reference_end, reference_start, 
        read_strand = "-" if read.is_reverse else "+"        
        read_blocks = get_read_blocks(read)
        total_block_num = len(read_blocks)
        read_core_id = ",".join([read.query_name,
                              read.reference_name, 
                              str(read.reference_start + 1), 
                              str(read.reference_end)])
        if read_core_id in have_extract_read_core_ids:
            continue
        have_extract_read_core_ids.add(read_core_id)
        read_name = ",".join([read_core_id,
                              read_strand,
                              str(total_block_num)])
        this_read_exons = []
        for i, (start, end) in enumerate(read_blocks):
            #exon_num = i + 1 if read_strand == "+" else total_block_num - i
            exon_num = i + 1
            block_name = read_name + "," + str(exon_num)
            this_read_exons.append([block_name, read.reference_name, start, end, read_strand])
        if len(read_blocks):
            yield([read_core_id, this_read_exons])
            
def get_all_reads_blocks(file_bam):
    """
    Input:
    file_bame
    
    Return:
    a dictory: key:   chromatin_name
               value: a list of all read_exons in this chromatin, 
                      each element is :
                      [read_exon_name, chromatin_name, start, end, read_alignment_strand]
    
    read_exon_name: 3217606a-2c2c-4598-9929-8ff20c1bf938,chr1,3618,4603,+,3,1
    start, end is 1-based
                      
    see detail in `iter_all_reads_blocks` and `get_read_blocks`
    """
    pos_dict = collections.defaultdict(list)
    for read_core_id, this_read_exons in iter_all_reads_blocks(file_bam):
        for block_name, chr_name, start, end, read_strand in this_read_exons:
            pos_dict[chr_name].append([block_name, start, end, read_strand])
    for chr_name, chr_data in pos_dict.items():
        chr_data.sort(key=lambda x: x[1])
    return pos_dict
    
def read_bed_file(filein, header = False, sorted = False):
    '''
    Input:
    filein: bed file
          column name: chromosome, start, end, id, na, direction, ...
          0-based
    
    Return:
    a dictory: key:   chromatin_name
               value: a list of all read_exons in this chromatin, 
                      each element is :
                      [exon_intron_name, start, end, read_alignment_strand]
               exon_intron_name: AT1G01010.1_exon1   
               1-based
    '''
    
    pos_dict = collections.defaultdict(list)
    with open(filein) as f:
        if header: f.readline()
        num = 0
        for raw_line in f:
            d = raw_line.rstrip("\n").split()
            chr_, start, end = d[:3]
            start, end = int(start) + 1, int(end)
            try:
                id_ = d[3]
            except:
                num += 1
                id_ = str(num)
            try:
                strand = d[5]
            except:
                strand = "+"
            pos_dict[chr_].append([id_, start, end, strand])
            
    if sorted:
        #sort by start position. 1 is the index of start.
        for chr_, chr_data in pos_dict.items():
            chr_data.sort(key=lambda x: x[1])

    return(pos_dict)
                    
def blocks_in_blocks(block_xs, block_ys, up=0, down=0):
    
    """
    Input:
    Both block_xs and block_ys are dict:
        key:   chromatin_name
        value: a list of all read_exons in this chromatin, 
               each element is :
               [block_name, start, end, read_alignment_strand]
    
    Output:
    A iterator. Each element is a list represent a overlap relationship: 
    [x_id, y_id, to_block_dir, type5, pos5, type3, pos3, 
        chr_name, x_start, x_end, x_dir, y_start, y_end, y_dir]
    x_id is block_name in block_xs
    y_id is block_name in block_ys
    to_block_dir "+":  the strand of x_id and x_id is same, "-": not same
    
    0-based
    
    type5, pos5: the 5' indicate the 5' positon of x_id relative to the 5' of y_id
    type3, pos3: the 3' indicate the 3' positon of x_id relative to the 5' of y_id
    
    
    y_id (length: l)  ----->>>>>>>>>>>>>>>>>>>>>>---------      >>>> repesent y_id
                      |    |   |                |    |
                      |  on5(0)|            on3(l-1) |
                    up(-4)  in(4)                    |
                                                down (l+4)
    
    x_id                    
    type5: on5 type3: on3  >>>>>>>>>>>>>>>>>>>>>>
    type5: on3 type3: on5  <<<<<<<<<<<<<<<<<<<<<<
    """
    
    def pos_2_block(pos, block_start, block_end, block_dir="+", onisin = False):
        # give a block 
        type_, relative_position = "", 0
        #type list : "down","on5","in","on3","up"
        if block_dir == "+":
            relative_position = pos - block_start  # 0 for the same position
            if pos < block_start:
                type_ = "up"
            elif pos == block_start:
                type_ = "on5"
            elif pos < block_end:
                type_ = "in"
            elif pos == block_end:
                type_ = "on3"
            else:
                type_ = "down"
        else:
            relative_position = block_end - pos
            if pos > block_end:
                type_ = "up"
            elif pos == block_end:
                type_ = "on5"
            elif pos > block_start:
                type_ = "in"
            elif pos == block_start:
                type_ = "on3"
            else:
                type_ = "down"
        if onisin:
            if type_ in ["on5", "on3"]:
                type_ = "in"
        return((type_, relative_position))
        
    def blocks_in_blocks_one_chr(x, y, up=0, down=0, chr_name=""):
        '''
        Find overlapped genome features (blocks) when give two
        kind of genome features. 
        x, y are lists of two kind of genome features sorted by start position.
        Blocks are also list: [id_, start, end, dir_].
        The start position is smaller than the end position no matther 
        what the direction is.
        '''
        
        results = []
        for data in y:
            start,end = data[1:3]
            up_pos = start - up
            if up_pos < 1: up_pos = 1
            data.append(up_pos)
            data.append(end+down)
        
        for x_id,x_start,x_end,x_dir in x:
            del_pos = []
            for i, (y_id,y_start,y_end,y_dir,y_up,y_down) in enumerate(y):

                if x_end < y_up: 
                    break
                elif x_start > y_down:
                    del_pos.append(i)
                    continue

                if x_dir == "+": 
                    x_5_terminal, x_3_terminal = x_start, x_end
                else:
                    x_5_terminal, x_3_terminal = x_end, x_start
            
                if x_dir == y_dir:
                    to_block_dir = "+"
                else:
                    to_block_dir = "-"

                type5, pos5 = pos_2_block(x_5_terminal, y_start, y_end, y_dir)
                type3, pos3 = pos_2_block(x_3_terminal, y_start, y_end, y_dir)

                results.append([x_id, y_id, to_block_dir, 
                                type5, pos5, type3, pos3, 
                                chr_name, x_start, x_end, x_dir,
                                y_start, y_end, y_dir])
            del_pos.reverse()
            for i in del_pos:
                del y[i]
        return results
        
    for chr_, this_chr_block_xs in block_xs.items():
        if chr_ not in block_ys:
            continue
        this_chr_block_ys = block_ys[chr_]
        results = blocks_in_blocks_one_chr(this_chr_block_xs, this_chr_block_ys, 
                                  up, down, chr_)
        for r in results:
            yield(r)

def compare_bam_read_to_exon_intron(file_bam, file_exon_intron):
    """
    Extract read exon features and compare their position to those of the exon and 
    intron in represent transcripts.
    
    Input: 
    file_bam
    file_exon_intron exon_intron bed file, 
            chr1    3630    3913    AT1G01010.1_exon1       .       +
            chr1    3913    3995    AT1G01010.1_intron1     .       +
    
    Output:
    A tuple of two element.
    The first element is a dict: key: read_core_id
                                 value: a list of read exon ordered by number
    The second element is a list of overlap relationship. Each element is a list 
    represent a overlap relationship: 
    [x_id, y_id, to_block_dir, type5, pos5, type3, pos3, 
        chr_name, x_start, x_end, x_dir, y_start, y_end, y_dir]
    x_id: read_exon_num: 7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687,2,exon1
    y_id: AT1G01010.1_exon1, AT1G01010.1_intron1
    
    see detail in `blocks_in_blocks`
    """
    
    def sorted_results_by_read(results):
        data = []
        id_orders = {}
        id_num = 0
        for d in results:
            #not filter strand
            #if d[2] == "-": continue
            data.append(d)
            id_ = d[0][:-len(d[0].split(",")[-1])]
            if id_ not in id_orders:
                id_num += 1
                id_orders[id_] = id_num
        data.sort(key = lambda x: [id_orders[x[0][:-len(x[0].split(",")[-1])]], x[8]])
        return(data)
    
    block_xs = get_all_reads_blocks(file_bam)
    block_ys = read_bed_file(file_exon_intron, header = False, sorted=True)
    
    results = blocks_in_blocks(block_xs, block_ys)
    results = sorted_results_by_read(results)
    return results
    
def read_mRNA_exon_num_and_mRNA_pos(fileexon_bed):
    """
    Extract exon information from fileexon_bed.
    Input:
    file_exon_intron exon_intron bed file (0-based), 
            chr1    3630    3913    AT1G01010.1_exon1       .       +
            chr1    3913    3995    AT1G01010.1_intron1     .       +
    
    Output:
    A dict: key:  mRNA name, such as AT1G01010.1. exon_name.split("_")[0]
            value:[chr_name, mRNA_start, mRNA_end, strand, total_exon_num]
            1-based
    """
    mRNA2exons = collections.defaultdict(list) 
    for l in open(fileexon_bed):
        d = l.rstrip("\n").split("\t")
        chr_name, start, end, exon_or_intron_name, na, strand = d
        start, end = int(start) + 1, int(end)
        mRNA, exon_or_intron_id = exon_or_intron_name.split("_")
        if exon_or_intron_id.startswith("i"):
            continue
        mRNA2exons[mRNA].append((chr_name, start, end, strand))
    mRNA_pos_exon_nums = {}
    for mRNA, mRNA_exons in mRNA2exons.items():
        total_exon_num = len(mRNA_exons)
        chr_name = mRNA_exons[0][0]
        strand = mRNA_exons[0][3]
        mRNA_start = min([i[1] for i in mRNA_exons])
        mRNA_end = max([i[2] for i in mRNA_exons])
        mRNA_pos_exon_nums[mRNA] = [chr_name, mRNA_start, mRNA_end, strand, total_exon_num]
    return mRNA_pos_exon_nums
    
def extract_read_information(file_bam, file_exon_intron):
    
    '''
    Input: 
    file_bam
    file_exon_intron exon_intron bed file, 
            chr1    3630    3913    AT1G01010.1_exon1       .       +
            chr1    3913    3995    AT1G01010.1_intron1     .       +
    '''
    
    def get_retention_introns(mRNA_data, f_feature_num, f_feature_type, l_feature_num):
        """
        Output: a set. retention_introns.
        f_feature_num, l_feature_num, and the values in retention_introns is integer.
        
        Only: for feature_type == "intron":
        (a) if it is the f_feature: it would be retented. but this may be affected by 
                3' SS.
        (b) if it is the l_feature: the transcription hasn't been finished. not retented.
        (c) in the inner: if coverage_ratio >= 0.8: retention.
        """
        #f_feature_num和l_feature_num在read_read_exon_info里已转变为
        #int格式(feature_num)，因此可以直接用。
        #只考虑feature_type == "intron"的情况
        #此时，如果feature_num == l_feature_num，也就是该feature是最后一个feature，
        #那就是转录还没转录到该内含子后，所以该内含子肯定不是retention intron。
        #注意如果l_feature_type，feature_num不可能等于l_feature_num。
        #如果feature_num == f_feature_num，且f_feature_type == "intron"
        #也就是该片段5’端断在内含子上，这个内含子应该还没剪切。
        #注意，这种情况可信性可能不高。譬如可能受可变3'SS
        #的影响，也可能受比对的影响。
        retention_introns = set()
        for (read_exon_num, feature_type, feature_num, 
                  type5, pos5, type3, pos3, 
                  feature_start, feature_end, feature_strand, feature_length,
                  read_exon_start, read_exon_end, coverage) in mRNA_data:
            if feature_type == "intron" and feature_num < l_feature_num:
                #feature_num < l_feature_num can remove the last intron which transcription 
                #hasn't been finished.
                if feature_num == f_feature_num and (f_feature_type == "intron"):
                    #if f_feature_type is intron, it would be retentaed
                    retention_introns.add(feature_num)
                elif coverage * 1.0 / feature_length >= 0.8:
                    retention_introns.add(feature_num)
        return retention_introns

    def read_read_exon_info(read_blocks_compare_exon_intron_results):
        """
        Input: 
        the result generated by `compare_bam_read_to_exon_intron`, see detail in this function.
        
        Output:
        A dicat: key:   read_core_id  (a7ac3172-124d-48bd-bb38-b458a994f21c,chr1,3648,5792)
                 value: a dict:
                        {
                        "info": [chr_, read_start, read_end, read_strand, read_exon_total_num]
                        "mRNA": { mRNA_name1: [record1, reacord2, record3] }
                        }
                 record is the overlapped infromation between read_core_id and a exon or intron in the mRNA_name1
                 record is a list: 
                 1-based [read_exon_num, feature_type, feature_num, 
                                type5, pos5, type3, pos3, 
                                feature_start, feature_end, feature_strand, feature_length,
                                read_exon_start, read_exon_end, coverage]
        
        feature_name is exon or intron
        AT1G01010.1_exon1 or AT1G01010.1_intron1  #orderred by the exon or intron order in mRNA
        mRNA_name is feature_name.split("_")[0]
        feature_type: exon or intron
        feature_num:  1        
        
        1. We don't care about the direction of read exon. Thus, we redefined the 
        means of type5 and type3. The 5 is the less index of gene_exon_intron 
        the 3 is the more. reverse them if pos3 < pos5. Thus now pos5 <= pos3.
        gene_exon            >>>>>>>>>>>>>>>>>>>>>>>>>>>>
        original: 
            read_exon1       >>>>>>>>>>>>>>>>>>>>>>>>>>>>
                             |                          |
                           type5                      type3
            read_exon2       <<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             |                          |
                           type3                      type5
        Converted:
            read_exon1       >>>>>>>>>>>>>>>>>>>>>>>>>>>>
                             |                          |
                           type5                      type3
            read_exon2       <<<<<<<<<<<<<<<<<<<<<<<<<<<<
                             |                          |
                           type5                      type3
        
        2. caluculate overlapped length: converage
        
        3. extract read mappign and read exon information from read exon num
        read exon name:
        a7ac3172-124d-48bd-bb38-b458a994f21c,chr1,3648,5792,-,6,1
        read core id:
        a7ac3172-124d-48bd-bb38-b458a994f21c,chr1,3648,5792
        6 means there are six read exons in this read
        1 means this read exon is the first read exon (ordered by chromatin position)
        
        feature_name is exon or intron
        AT1G01010.1_exon1 or AT1G01010.1_intron1
        mRNA_name is feature_name.split("_")[0]
        
        #4. record read inforamtion to the dict data
        """
                
        data = {}
        for (read_exon_name, feature_name, strand, 
                type5, pos5, type3, pos3, 
                chr_, read_exon_start, read_exon_end, read_strand, 
                feature_start, feature_end, feature_strand) in read_blocks_compare_exon_intron_results:
            
            read_exon_start, read_exon_end = int(read_exon_start), int(read_exon_end)
            feature_start, feature_end = int(feature_start), int(feature_end)
        
            #1. reverse type5 and type3, and reverse pos5, pos3 if pos3 < pos5
            pos5, pos3 = int(pos5), int(pos3)
            if pos3 < pos5:
                pos5, pos3 = pos3, pos5
                type5, type3 = type3, type5
            feature_length = feature_end - feature_start + 1
            pos3 = pos3 + 1 - feature_length 

            #2. caluculate coverage
            coverage = feature_length
            if pos5 > 0: coverage -= pos5
            if pos3 < 0: coverage += pos3
            if coverage < 0: raise 1/0
        
            #3. extract read mappign and read exon information from read exon num
            _ids = read_exon_name.split(",")
            read_core_id = ",".join(_ids[:4]) 
            read_start, read_end = int(_ids[2]), int(_ids[3])
            read_exon_total_num = int(_ids[5])
            read_exon_num = int(_ids[6])
            
            mRNA, feature_id = feature_name.split("_")
            if feature_id.startswith("e"):
                feature_type = "exon"
                feature_num = int(feature_id[4:])
            else:
                feature_type = "intron"
                feature_num = int(feature_id[6:])
        
            #4. record read inforamtion to the dict data
            read_info = [chr_, read_start, read_end, read_strand, 
                            read_exon_total_num]
            record = [read_exon_num, feature_type, feature_num, 
                      type5, pos5, type3, pos3, 
                      feature_start, feature_end, feature_strand, feature_length,
                      read_exon_start, read_exon_end, coverage]
    
            if read_core_id not in data:
                data[read_core_id] = {"info": read_info,
                                 "mRNA": { mRNA: [record] }
                                }
            elif mRNA not in data[read_core_id]["mRNA"]:
                data[read_core_id]["mRNA"][mRNA] = [record]
            else:
                data[read_core_id]["mRNA"][mRNA].append(record)
            
        return data
    
    mRNA_pos_exon_nums = read_mRNA_exon_num_and_mRNA_pos(file_exon_intron)
    results = compare_bam_read_to_exon_intron(file_bam, file_exon_intron)
    data = read_read_exon_info(results)
    
    for read_core_id, read_data in data.items():
        chr_, read_start, read_end, read_strand, read_exon_total_num = read_data["info"]
        mRNA_datas = read_data["mRNA"]
        
        #1. Only keep the mRNA with longest overlapped region with the read_core_id 
        #sorted by feature_num, reverse = True
        mRNA_coverages = []
        for mRNA, mRNA_data in mRNA_datas.items():
            ##i[-1] is coverage, the length of overlapped region
            total_coverage = sum([i[-1] for i in mRNA_data]) 
            mRNA_coverages.append([mRNA, mRNA_data, total_coverage])
        mRNA_coverages.sort(key = lambda x: x[2], reverse = True)
        mRNA, mRNA_data, total_coverage = mRNA_coverages[0]
        
        #2. sorted mRNA_data by the order in mRNA.
        #[read_exon_num, feature_type, feature_num,type5, pos5, 
        #  type3, pos3, feature_start, feature_end, feature_strand,
        #  feature_length, read_exon_start, read_exon_end, coverage]
        #Note: one feature may overlapped with two read exon, thus 
        #sort by (x[7], x[11]), not only x[7]
        mRNA_strand = mRNA_data[0][9]
        if mRNA_strand == "+":
            mRNA_data.sort(key = lambda x: (x[7], x[11])) #7 is feature_start, 11 is read_exon_start
        else:
            mRNA_data.sort(key = lambda x: (x[7], x[11]), reverse = True)
        chr_, mRNA_start, mRNA_end, mRNA_strand, mRNA_exon_num = mRNA_pos_exon_nums[mRNA]
                        
        if mRNA_strand == "+":
            mRNA_pos5 = read_start - mRNA_start
            mRNA_pos3 = read_end - mRNA_end
        else:
            mRNA_pos5 = mRNA_end - read_end
            mRNA_pos3 = mRNA_start - read_start
        mRNA_length = mRNA_end - mRNA_start + 1

        #3. care about the position of 5' end and 3' end of read
        (f_read_exon_num, f_feature_type, f_feature_num,
            f_type5, f_pos5, f_type3, f_pos3, 
            f_feature_start, f_feature_end, f_feature_strand, f_feature_length,
            f_read_exon_start, f_read_exon_end, f_coverage) = mRNA_data[0]
        
        (l_read_exon_num, l_feature_type, l_feature_num,
            l_type5, l_pos5, l_type3, l_pos3, 
            l_feature_start, l_feature_end, l_feature_strand, l_feature_length,
            l_read_exon_start, l_read_exon_end, l_coverage) = mRNA_data[-1]
        
        #4. Calculate rel_mRNA_pos5, rel_mRNA_pos3, and filter not consitent
        if mRNA_strand == "+":
            rel_mRNA_pos5 = f_read_exon_start - mRNA_start
            rel_mRNA_pos3 = l_read_exon_end - mRNA_end
        else:
            rel_mRNA_pos5 = mRNA_end - f_read_exon_end
            rel_mRNA_pos3 = mRNA_start - l_read_exon_start
                
        #5. extract retention_introns
        retention_introns = get_retention_introns(mRNA_data, f_feature_num, f_feature_type, l_feature_num)
        retention_intron_num = len(retention_introns)
        str_retention_introns = ":".join([str(i) for i in sorted(list(retention_introns))])
        
        #6. calcualte end type
        #read overlapped with the the first exon is end5.
        #read overlapped with the the last exon is end3.
        #	     noend5 end5
        #noend3      0    1
        #end3        2    3
        end5_type = 1 if f_feature_type == "exon" and f_feature_num == 1 else 0
        end3_type = 1 if l_feature_type == "exon" and l_feature_num == mRNA_exon_num else 0
        end_type = end5_type + end3_type * 2
        
        #7. span intron number
        #the first intron number is `f_feature_num` no matter whether the f_feature_type is exon or intron
        #the last intron number is `l_feature_num - 1` no matter whether the f_feature_type is exon or intron
        span_intron_num = l_feature_num - f_feature_num        
        mRNA_intron_num = mRNA_exon_num - 1
        
        d = [
            read_core_id, chr_, read_start, read_end, read_strand, read_exon_total_num, 
            mRNA, mRNA_intron_num, mRNA_start, mRNA_end, mRNA_strand, mRNA_length, mRNA_pos5, mRNA_pos3,
            rel_mRNA_pos5, rel_mRNA_pos3, 
            total_coverage, f_read_exon_num, f_feature_type, f_feature_num, f_feature_length, f_pos5, f_pos3, 
            l_read_exon_num, l_feature_type, l_feature_num, l_feature_length, l_pos5, 
            l_pos3, end_type, str_retention_introns, retention_intron_num, span_intron_num
        ]
        yield(d)
        
if __name__ == "__main__":
    main()
    
