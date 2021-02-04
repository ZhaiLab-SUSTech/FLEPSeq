import pysam
import sys
import click
import collections

"""

By Jia Jinbu
2021.02 version 0.3
2019.12 version 0.2
2018.06 version 0.1

Usage:
python ASCaller.py -i file_bam -o fileout --file_intron_pos intron.bed --strand_flag 0 --min_overlap 6 --debug_file read_intron_type.txt
"""


@click.command()
@click.option('-i', '--inbam', help=('Input the bam file.'), 
                    required=True, type=click.Path(exists=True))
@click.option('-l', '--file_intron_pos', help=('Input intron pos file'), 
                    required=True, type=click.Path(exists=True))
@click.option('-o', '--fileout', help=("""Output File of unspliced ratio intron results.
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
                    iratio	                0.545014521
                    sratio	                0.416263311
                    oratio	                0.038722168
                    o1ratio	                0.013552759
                    o2ratio	                0.011616651
                    o1_type	                o||1|:7829-7941
                    o1_count	            7
                    o2_type	                o|1|0|:7902-8235
                    o2_count	            6
                    other_o	                o|1|0|:7907-8235=3;o|1|1|:7650-8235=2;o||1|:7832-7941=1;o|1||:7836-8235=1
                    
                    iratio: unspliced ratio of intron
                    sratio: spliced ratio of intron
                    oratio: the ratio of alternative splicing
                    
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
                        o|x1|x2|x3:start1-end1:<start2-end2>
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
                                 
                                 #readname2flag is used for pair-end read.
                                 #if one is startswith("p"), we just filter it.
                                 #if flags of two pair-end reads are consistent, only count once.
                                 #if not consistent, justify as fellow:
                                 readname2flag = {}
                                 for flag, read_name, str_blocks in this_intron_data:
                                     if flag.startswith("p"): continue
                                     if read_name in readname2flag:
                                         #potential flag value:
                                         #a b c ab i
                                         #o|x|x|:x-x  #such as: o|0|1|:214456-214607
                                         #pa pb pc pab pi
                                         #po|x|x|:x-x
                                         #potential cases when the flags of pair-end reads are not consistent.
                                         #a/b ab   --> ab
                                         #o   ab   --> ab  #wrong mapping?
                                         #a/b c    --> c   #wrong mapping?
                                         #o   c    --> c   #wrong mapping?
                                         #a/b o    --> o   #a may be o
                                         #a/b a/b    --> ab
                                         #o   o    --> o   #wrong mapping? select o with shorter junction
                                         #i a/b/ab/o/c
                                 
                                 
                                 
                    For pair-end reads, only count once. If the types of pair-end reads are different, the type is set 
                    as the order:
                    one is ab    --> ab
                    one a, one b --> ab
                    one is c     --> c #wrong mapping?
                    one is o     --> o #diffrent alternative splicing type means wrong mapping, especially for genes with similar sequence
                    """), required=True)
@click.option('-s', '--strand_flag', required=False, default=0, help=("""Strand-specific format. 0: no. 1: RF. 2: FR. Default: 0."""))
@click.option('-m', '--min_overlap', required=False, default=6, help=("""Min overlap. Default: 6."""))
@click.option('-d', '--debug_file', required=False, default="", help='Write junction read information to this file. Default: Not output.')
    
def main(inbam, file_intron_pos, fileout, strand_flag, min_overlap, debug_file):
    with open(fileout, 'w') as OUT:
        OUT.write("intron_id\tchr_name\tintron_start\tintron_end\tintron_strand\tintron_coverage_ratio\ta\tb\tab\tc\to\tt\tiratio\tsratio\toratio\to1ratio\to2ratio\to1_type\to1_count\to2_type\to2_count\tother_o\n")
        for this_intron_stat_result in iter_stat_intron_read_type(inbam, file_intron_pos, min_overlap, strand_flag, debug_file):
            OUT.write("\t".join([str(i) for i in this_intron_stat_result]) + "\n")
                            
def read_bed_file(filein, header = False, will_sort = False):
    '''
    Input:
    filein: bed file
          column name: chromosome, start, end, id, na, direction, ...
          0-based
    
    Return:
    a dictory: key:   chromatin_name
               value: a list of all features in this chromatin, 
                      each element is :
                      [exon_intron_name, start, end, strand]
               exon_intron_name: AT1G01010.1_exon1   
               1-based
    Also used in extract_read_info.py
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
            
    if will_sort:
        #sort by start position. 1 is the index of start.
        for chr_, chr_data in pos_dict.items():
            chr_data.sort(key=lambda x: x[1])

    return(pos_dict)

def convert_pos_feature(d):
    
    """
    Flatten a dict which element is a list of features.
    Input: features. such as the result of `read_bed_file` function.
           A dictory: key is chromatin_name, value is a list of all features in this chromatin.
                      each element is: [feature_id, start, end, strand]
    Output: a iterator of all features. each element is: [[feature_id, chr_name, start, end, strand]
    """
    for chr_, chr_features in d.items():
        for (id_, start, end, strand) in chr_features:
            yield([id_, chr_, start, end, strand])                


def filter_bam_aligns(aligns, filter_tag=["mapped"]):
    
    """
    Filter read objects and convert read_obj to [chr_name, start, end, strand, read_obj]
    Parameters
    ----------
    aligns: a iterator of read objects generated by pysam.
    filter_tag: a list of str, such as ["unmapped"]:
            "mapped": removed unmapped reads
            "noduplicated": removed duplicated reads
            "uniquely_mapped": removed uniquley mapped reads of hisat2 output based NH tag.
    Output
    ------
        a iterator of [chr_name, start, end, strand, read_obj]
    """
    
    def _generate_bam_filter_func(filter_tag=["mapped"]):    
        tag2func = {
            "mapped": lambda read: read.is_unmapped,
            "noduplicated": lambda read: read.is_duplicate,
            "uniquely_mapped": lambda read: read.get_tag("NH") != 1 #for hisat2
        }
    
        funcs = [tag2func[tag] for tag in filter_tag]
    
        def _filter(read):
            for _func in funcs:
                if _func(read):
                    return True
            return False
    
        return _filter
    
    filter_read_func = _generate_bam_filter_func(filter_tag)
    for read in aligns:                
        if filter_read_func(read): continue
        chr_ = read.reference_name
        #id_ += 0
        start = read.reference_start + 1
        strand = "-" if read.is_reverse else "+" 
        end = read.reference_end
        yield([ chr_, start, end, strand, read]) #str(id_),

def iter_pos_feature(aligns, features, updown=0, is_point=0, simple_chr=False):

    """    
    Parameters
    ----------
    aligns : iterator. each element is a feature (`align_data`):
            if is_point == 0: [chr_, start, end, ....] 
            if is_point == 0: [chr_, start, ...]
            Nomatter whether the strand is "+", start is less than end,
             and the elemetns are sorted by chr_ start position.
            strand not used in this function.
            start, end: 1-based
            strand not used     
    features : dict. key is `chr_`, value is the list of feature 
            (`feature_data`) sorted by start position. each feature is represented by a list. 
            the function only use the second and the third element as start 
            and end position, a typical `feature_data` is [id, start, end, strand].
            The start is less than end no matther what the strand is.
            feature also can be bed file.
    
    updown : Upstream or Downstream bases of feature. nt, default 0.
    
    down : Downstream bases of feature. nt, default 0.
    
    pos_format: if is_point == 0: [chr_, start, end, ....] 
                if is_point == 0: [chr_, start, ...]
        
    Return values
    ----------
    Iterator. [chr_, feature_data]. Each element represents a feature and aligns overlapped with the feature.
    `feature_data` is a copy of `feature_data` stored in `features`, but append a additional element which 
    stores the list of `align_data`. In fact, we append a `align_data` in orignal `features`, and copy `features`
    to a new list, and delete `align_data` from orignal `features`.
    """
    
    if isinstance(features, str):
        features = read_bed_file(features, will_sort = True)
    
    this_chr = ""
    this_feature_aligns = []
    scaned_chr = set()
    scan_this_chr_feature = 0
    ls_test_i = 0
    for align_data in aligns:
        ls_test_i += 1
        if is_point:
            align_chr, align_start = align_data[:2]
            align_end = align_start
        else:
            align_chr, align_start, align_end = align_data[:3]
        if align_chr != this_chr:
            scaned_chr.add(align_chr)
            if this_chr:
                #After scan a chr, yield the features haven't been oupput.
                for d in this_feature_data[feature_indicator:]:
                    yield([this_chr, d.copy()])
                    d.pop()
            this_chr = align_chr
            this_feature_data = features[this_chr]
            for d in this_feature_data:
                d.append([])
            feature_indicator = 0
            scan_this_chr_feature = 0
            feature_num = len(this_feature_data)
        if scan_this_chr_feature: continue
        pad_align_start = align_start - updown #down
        pad_align_end = align_end + updown #up
        for i in range(feature_indicator, feature_num):
            feature_info = this_feature_data[i]
            if pad_align_end < feature_info[1]:
                break
            elif pad_align_start > feature_info[2]:
                if i == feature_indicator:
                    feature_indicator = i + 1
                    yield([this_chr, feature_info.copy()])
                    this_feature_data[i].pop()
                    if feature_indicator == feature_num:
                        scan_this_chr_feature = 1
                    continue
            else:
                feature_info[-1].append(align_data)
    if this_chr:
        #After scan the last chr, yield the features haven't been oupput.
        for d in this_feature_data[feature_indicator:]:
            yield([this_chr, d.copy()])
            d.pop()
    for chr_, feature_data in features.items():
        #After scanning, yield the chrs haven't been scanned. 
        if chr_ not in scaned_chr:
            for d in feature_data:
                d1 = d.copy()
                d1.append([])
                yield([chr_, d1])
                
def iter_bam_feature(file_bam, file_bed, filter_bam_func=filter_bam_aligns, 
                    filter_bam_params={"filter_tag":["mapped"]}, updown=0, method=1):
    
    with pysam.AlignmentFile(file_bam, "rb") as bam_obj:
        if method == 1:
            features = read_bed_file(file_bed, will_sort = True)
            features = convert_pos_feature(features)
            for feature in features:
                id_, chr_, start, end, strand = feature
                start = start - 1 - updown
                end = end + updown
                if start <0:
                    start = 0 
                aligns = bam_obj.fetch(contig=chr_, start=start, end=end)
                aligns = filter_bam_func(aligns, **filter_bam_params)
                yield([feature, aligns])
        else:
            features = read_bed_file(file_bed, will_sort = True)
            aligns = bam_obj.fetch()
            aligns = filter_bam_func(aligns, **filter_bam_params)       
            for chr_name, (gene_id, start, end, block_strand, this_aligns) in iter_pos_feature(aligns, features, updown):
                yield([[gene_id, chr_name, start, end, block_strand], this_aligns])

def get_read_blocks(read):
    
    """
    Input: read objects generated by pysam.
    Output: a list of mapped blocks (genome position, 1-based) of read.
         Each element is a list, [start, end].
         Don't care small InDel and mismatches.
         Sorted by genome position.
    """
    #don't use read.get_blocks
    #insert also seprate block
    #bam_obj.find_introns !!!
    blocks = []
    start = read.reference_start + 1
    end = start - 1
    for (flag, length) in read.cigartuples:
        if flag == 4 or flag == 5 or flag == 1: continue
        if flag == 0 or flag == 2:
            end += length
        if flag == 3:
            blocks.append([start, end])
            start = end + length + 1
            end = start - 1
    blocks.append([start, end])
    return blocks

def trim_read_blocks(read_blocks, start, end):
    """
    Parameters
    ----------
    read_blocks
    start
    end
    ----------
    cut read bocks to three parts generated by `get_read_blocks` based on the position 
    relative to [start, end] region. Not care strand.
    Each part is also a list, and element is [block_start, block_end],
    genome position, 1-based.
    The first part: The blocks locating upstream of [start, end] region.
    The second part: The blocks overlapped with [start, end] region.
    The third part: The blocks locating downstream of [start, end] region.
    """
    
    block_total_num = len(read_blocks)
    #in_flag: 0: left, 1: overlap, 2: right
    in_flags = []
    for s, e in read_blocks:
        if e < start: 
            in_flag = 0
        elif s > end:
            in_flag = 2
        else:
            in_flag = 1
        in_flags.append(in_flag)
    
    #find inner continul blocks in [start, end]
    blocks_left_end_index = -1
    blocks_start_index = block_total_num
    blocks_end_index = -1
    blocks_right_start_index = block_total_num
    for i, in_flag in enumerate(in_flags):
        if in_flag == 0:
            blocks_left_end_index = i
        elif in_flag == 1:
            if blocks_start_index == block_total_num: #indicate that the frist is 1
                blocks_end_index = blocks_start_index = blocks_left_end_index + 1
            else:
                blocks_end_index = i
        elif in_flag == 2:
            blocks_right_start_index = i
            break
    
    left_blocks = read_blocks[:(blocks_left_end_index+1)] if blocks_left_end_index != -1 else []
    inner_blocks = read_blocks[blocks_start_index:(blocks_end_index+1)] if blocks_start_index != block_total_num else []
    right_blocks = read_blocks[blocks_right_start_index:] if blocks_right_start_index != block_total_num else []
    return([left_blocks, inner_blocks, right_blocks])
    
def iter_intron_read_type(file_bam, file_intron_pos, min_overlap=6, strand_flag=0, debug_file=""):
    
    filter_tag=["mapped", "noduplicated"]
    feature_aligns = iter_bam_feature(file_bam, file_intron_pos, filter_bam_params={"filter_tag":filter_tag}, method=2)
    if debug_file:
        READ_INFO_OUT = open(debug_file, 'w')
        READ_INFO_OUT.write("intron_id\tchr_name\tintron_start\tintron_end\tintron_strand\tflag\tstr_blocks\tis_read1\tstrand\tread_name\n")
    
    for [intron_id, chr_name, intron_start, intron_end, intron_strand], aligns in feature_aligns:
        
        #border position (include)
        #The first 1 indicate the left of intorn, 2 indicate the right of intron
        #The second 1 indicate the left region, 2 indicate the right region
        xover11 = intron_start - min_overlap
        xover12 = intron_start + min_overlap - 1
        xover21 = intron_end - min_overlap + 1
        xover22 = intron_end + min_overlap
        
        intron_reads = [] #Supporting intron retention reads, including inner reads
        this_intron_data = []
        for na, start, end, strand, read in aligns:
            
            #1. For strand specific library, removed antisense reads.
            if strand_flag == 1:
                if read.is_read1:
                    if strand == intron_strand: continue
                else:
                    if strand != intron_strand: continue
            elif strand_flag == 2:
                if read.is_read1:
                    if strand != intron_strand: continue
                else:
                    if strand == intron_strand: continue
            
            #2. Cut read blocks into three part based on the [intron_start, intron_end] region.
            left_blocks, blocks, right_blocks = trim_read_blocks(get_read_blocks(read), intron_start, intron_end)
                            
            #3. Extract junction information
            junction = ""
            if left_blocks:
                junction_start = left_blocks[-1][1] + 1
                if blocks:
                    if blocks[0][0] >= intron_start:
                        junction = str(junction_start) + "-" + str(blocks[0][0] - 1)
                else:
                    junction = str(junction_start) + "-" + str(right_blocks[0][0] - 1)
            if len(blocks) >= 2:
                for i in range(len(blocks)-1):
                    if junction: junction += ":"
                    junction += str(blocks[i][1] + 1) + "-" + str(blocks[i+1][0] - 1)
            if right_blocks:
                junction_end = right_blocks[0][0] - 1
                if blocks and (blocks[-1][1] <= intron_end):
                    if junction: junction += ":"
                    junction += str(blocks[-1][1] + 1) + "-" + str(junction_end)                        
            
            #4. Calcualte read type
            if blocks:
                block_start = blocks[0][0]
                block_end = blocks[-1][1]
                
            if not blocks:
                #When no blocks, must have left_blocks and right_blocks
                #Must have a junction span the intron (is the intron, or longer than the intron)
                #--------...............--------- gene structure
                #  ------...............-----     ["c"]
                #  ------.................---     ["o", "", "1"]
                #  ----.................-----     ["o", "1"]
                #  ----...................---     ["o", "1", "1"]
                if junction_start == intron_start:
                    if junction_end == intron_end:
                        flag = ["c"]   #spliced as the annotated intron
                    else:
                        flag = ["o", "", "1"] #3' outer spliced, may be exon skipping
                else:
                    if junction_end == intron_end:
                        flag = ["o", "1"] #5’ outer spliced, may be exon skipping
                    else:
                        #5' and 3' outer spliced, may be exon skipping
                        flag = ["o", "1", "1"]
            elif len(blocks) == 1:
                if block_start < intron_start:
                    #--------...............--------- gene structure
                    #      -------------------------  ["ab"]
                    #      ------------.....--------- ["o", "0"]  5' inner spliced
                    #      ------------.........----  ["o", "0", "1"] 5' inner 3' outer
                    #      --------                   ["a"]
                    if block_end > intron_end:
                        flag = ["ab"]
                    else:
                        if right_blocks:
                            #junction_end must >= intron_end
                            if junction_end == intron_end:
                                flag = ["o", "0"]
                            else:
                                flag = ["o", "0", "1"]
                        else:
                            flag = ["a"]
                else:
                    #o||0u: Most are o||0,  but maybe o|||e1,  o||1|e1  #intron 3' not coveraged
                    #o|1|0u:Most are o|1|0, but maybe o|1||e1, o|1|1|e1 #intron 3' not coveraged
                    #o|0u|: Most are o|0,   but maybe o|||e1,  o|1||e1  #intron 5' not coveraged
                    #o|0u|1:Most are o|0,   but maybe o||1|e1, o|1|1|e1 #intron 5' not coveraged
                    #i multple possibles                         #intron 5' and 3' not coveraged
                    #--------...............--------- gene structure
                    #  ------.....--------------      ["o","", "0"] 3' inner
                    #  ------.....-----               ["o","", "0u"]
                    #  ------.....-----.....------    ["o","","","e1"] ExonSkipping
                    #  ------.....-----.......-----   ["o","", "1", "e1"] 
                    #
                    #  ---........-------------       ["o","1", "0"] 5' outer 3' inner
                    #  ---........-----.....------    ["o","1","", "e1"]
                    #  ---........-----.......-----   ["o","1", "1", "e1"]
                    #  ---........-----               ["o","1", "0u"]
                    #
                    #             ---------------     ["b"]
                    #             -----.....------    ["o","0u"]
                    #             -----........---    ["o","0u", "1"]
                    #             ------              ["i"]
                    if left_blocks:
                        #junction_start must <= intron_start
                        if junction_start == intron_start:
                            if block_end > intron_end:
                                flag = ["o","", "0"]
                            else:
                                if right_blocks:
                                    if junction_end == intron_end:
                                        flag = ["o","","","e1"]
                                    else:
                                        flag = ["o","", "1", "e1"]
                                else:
                                    flag = ["o","", "0u"]                                   
                        else:
                            if block_end > intron_end:
                                flag = ["o","1", "0"]
                            else:
                                if right_blocks:
                                    if junction_end == intron_end:
                                        flag = ["o","1","", "e1"]
                                    else:
                                        flag = ["o","1", "1", "e1"]
                                else:
                                    flag = ["o","1", "0u"]
                    else:
                        if block_end > intron_end:
                            flag = ["b"]
                        else:
                            if right_blocks:
                                #junction_end must >= intron_end
                                if junction_end == intron_end:
                                    flag = ["o","0u"]
                                else:
                                    flag = ["o","0u", "1"]
                            else:
                                flag = ["i"]
            elif len(blocks) >= 2: 
                #--------...............--------- gene structure
                #  ------..---...---....---       ["o", "", "", "e2"] #n=block_num
                #     ------.......--------       ["o", "0", "0"]  #n=block_num - 2 
                #     ------...--..--.-----       ["o", "0", "0", "e2"] #n=block_num
                #  ----....---...---....---       ["o", "1", "", "e2"] #n=block_num
                #  ------..---...---......-       ["o", "", "1", "e2"] #n=block_num
                #  ----....---........-----       ["o", "1", "0", "e1"] #n=block_num - 1
                #  --------..---..........-       ["o", "0", "1", "e2"] #n=block_num - 1
                #          ---...---....---       ["o", "0u", "", "e1"] #n=block_num - 2
                #          ---........-----       ["o", "0u", "0"] #n=block_num - 1
                #          ---...---......-       ["o", "0u", "1", "e1"] #n=block_num - 1
                #   -----...---.---               ["o", "", "0u", "e2"] #n=block_num - 1
                #   --------....---               ["o", "0","0u"] #n=block_num - 2
                #   ---.....---.---               ["o", "1","0u", "e1"] #n=block_num - 1
                e_num = len(blocks)
                flag = ["o", "", ""]
                if block_start < intron_start:
                    e_num -= 1
                    flag[1] = "0"
                else:
                    if left_blocks:
                        if junction_start < intron_start:
                            flag[1] = "1"
                        else:
                            flag[1] = ""
                    else:
                        flag[1] = "0u"
                        e_num -= 1
                if block_end > intron_end:
                    e_num -= 1
                    flag[2] = "0"
                else:
                    if right_blocks:
                        if junction_end > intron_end:
                            flag[2] = "1"
                        else:
                            flag[2] = ""
                    else:
                        flag[2] = "0u"
                        e_num -= 1
                if e_num:
                    flag.append("e" + str(e_num))

            #5. Calculate whether the overlap of both ends is ok. flag_overlap == True: ok
            #If have left block or right block, don't care the overlap.
            flag_overlap = True
            if flag[0] == "a":
                flag_left_block_overlap = (len(left_blocks) > 0) or (blocks[0][0] <= xover11)
                flag_right_block_overlap = blocks[0][1] >= xover12
                flag_overlap = flag_left_block_overlap and flag_right_block_overlap
            elif flag[0] == "b":
                flag_left_block_overlap =  blocks[0][0] <= xover21
                flag_right_block_overlap = (len(right_blocks) > 0) or (blocks[0][1] >= xover22)
                flag_overlap = flag_left_block_overlap and flag_right_block_overlap
            elif flag[0] == "ab":
                flag_left_block_overlap = (len(left_blocks) > 0) or (blocks[0][0] <= xover11)
                flag_right_block_overlap = (len(right_blocks) > 0) or (blocks[0][1] >= xover22)
                if flag_left_block_overlap:
                    if not flag_right_block_overlap:
                        flag[0] = "a"
                    else:
                        flag[0] = "ab"
                else:
                    if not flag_right_block_overlap:
                        flag[0] = "ab"
                        flag_overlap = False
                    else:
                        flag[0] = "b"
            elif flag[0] == "i":
                flag_overlap = True
            else:
                flag_left_block_overlap = False
                flag_right_block_overlap = False
                if left_blocks:
                    flag_left_block_overlap = (len(left_blocks)>1) or (left_blocks[0][1] - left_blocks[0][0] + 1 >= min_overlap)
                else:
                    if block_start < intron_start:
                        if block_start <= xover11:
                            flag_left_block_overlap = True
                    else:
                        if (blocks[0][1] - blocks[0][0] + 1) >= min_overlap:
                            flag_left_block_overlap = True
                if right_blocks:
                    flag_right_block_overlap = (len(right_blocks)>1) or (right_blocks[0][1] - right_blocks[0][0] + 1 >= min_overlap)
                else:
                    if block_end > intron_end:
                        if block_end >= xover22:
                            flag_right_block_overlap = True
                    else:
                        if (blocks[-1][1] - blocks[-1][0] + 1) >= min_overlap:
                            flag_right_block_overlap = True
                flag_overlap = flag_left_block_overlap and flag_right_block_overlap
            
            #6. Extend flags start with "o" to 4 elements by "".
            if flag[0] == "o":
                flag = flag + [""]*(4-len(flag))
            
            #6. If intron_strand is "-", modify flag
            if intron_strand == "-":
                if flag[0] == "a":
                    flag = ["b"]
                elif flag[0] == "b":
                    flag = ["a"]
                elif flag[0] == "ab":
                    pass
                elif flag[0] == "c":
                    pass
                elif flag[0] == "i":
                    pass
                else: #for "o"
                    ls5, ls3 = flag[2], flag[1]
                    flag[1] = ls5
                    flag[2] = ls3
            
            #7. Convert flag list to str, don't care for now. Delete u.          
            flag = "|".join(flag)
            flag = flag.replace("u", "")
            
            #8. Add junction to flag
            if junction and (flag[0] != "c"):
                flag = flag + ":" + junction
            
            #9. recored intron_reads, use for calculate intron coverage
            if flag == "a" or flag == "b" or flag == "ab" or flag == "i":
                intron_reads.append([block_start, block_end])
            
            #10. If overlap in both end is not ok, add 'p' (pass) before flag
            if not flag_overlap:
                flag = "p" + flag
            # str(start), str(end), strand, 
            str_blocks = ":".join([",".join([str(s) + "-" + str(e) for s, e in left_blocks]),
                                   ",".join([str(s) + "-" + str(e) for s, e in blocks]),
                                   ",".join([str(s) + "-" + str(e) for s, e in right_blocks])])
            if debug_file:
                READ_INFO_OUT.write("\t".join([intron_id, chr_name, str(intron_start), str(intron_end), intron_strand, 
                                    flag, str_blocks, str(read.is_read1), strand, read.query_name]) + "\n") 
            this_intron_data.append([flag, read.query_name, str_blocks])
            
        yield([intron_id, chr_name, intron_start, intron_end, intron_strand], this_intron_data, intron_reads)
    try:
        READ_INFO_OUT.close()
    except:
        pass

def stat_intron_coverage(intron_reads, intron_start, intron_end):
    if not intron_reads: return 0
    intron_length = intron_end - intron_start + 1
    intron_reads.sort(key=lambda x: x[0])
    start, end = intron_reads[0]
    if start < intron_start: start = intron_start
    if end >= intron_end:
        return((intron_end - start + 1)/intron_length)
    total_coverage = 0
    if len(intron_reads) > 1:
        for s, e in intron_reads[1:]:
            if s <= end:
                if e > end:
                    if e >= intron_end:
                        return((total_coverage + (intron_end - start +1))/intron_length)
                    end = e
            else:
                total_coverage += (end - start + 1)
                start = s
                if e >= intron_end:
                    return((total_coverage + (intron_end - start +1))/intron_length)
                else:
                    end = e
    total_coverage += end - start + 1
    return(total_coverage/intron_length)

def cal_junction_length(flag):
    d = flag.split(":")
    if len(d) == 1: return 0
    d = d[1:]
    length = 0
    for j in d:
        s, e = j.split("-")
        s, e = int(s), int(e)
        length += abs(e - s) + 1
    return length

def iter_stat_intron_read_type(inbam, file_intron_pos, min_overlap=6, strand_flag=0, debug_file=""):
    for [intron_id, chr_name, intron_start, intron_end, intron_strand], this_intron_data, intron_reads in iter_intron_read_type(inbam, file_intron_pos, min_overlap, strand_flag, debug_file):
        intron_coverage_ratio = stat_intron_coverage(intron_reads, intron_start, intron_end)
        #readname2flag is used for pair-end read.
        #if one is startswith("p"), we just filter it.
        #if flags of two pair-end reads are consistent, only count once.
        #if not consistent, justify as fellow:
        readname2flag = {}
        for flag, read_name, str_blocks in this_intron_data:
            if flag.startswith("p"): continue
            if read_name in readname2flag:
                #potential flag value:
                #a b c ab i
                #o|x|x|:x-x  #such as: o|0|1|:214456-214607
                #pa pb pc pab pi
                #po|x|x|:x-x
                #potential cases when the flags of pair-end reads are not consistent.
                #a/b ab   --> ab
                #o   ab   --> ab  #wrong mapping?
                #a/b c    --> c   #wrong mapping?
                #o   c    --> c   #wrong mapping?
                #a/b o    --> o   #a may be o
                #a/b a/b    --> ab
                #o   o    --> o   #wrong mapping? select o with shorter junction
                #i a/b/ab/o/c
                f = readname2flag[read_name]
                if f == flag:
                    pass
                elif f.startswith("ab") or flag.startswith("ab"):
                    #ab a/b/c/o/i, a/b/c/o/i ab
                    flag = "ab"
                elif f.startswith("c") or flag.startswith("c"):
                    #c a/b/o/i, a/b/o/i c
                    flag = "c"
                elif f.startswith("o") or flag.startswith("o"):
                    #About thousands reads in cyto4 sample is o o, or o c.
                    #many from AT1G07930 and AT1G07920, due the high similar of these two genes.
                    #one read were spliced into two genes.
                    if (f[0] == "o" and flag[0] == "o"):        
                        #o o                
                        junction_length1 = cal_junction_length(f)
                        junction_length2 = cal_junction_length(flag)
                        if junction_length2 > junction_length1:
                            flag = f
                    else: #most are a, b, or i
                        #o a/b/i, a/b/i o
                        if f[0] == "o": 
                            flag = f
                elif f == "a":
                    if flag == "b":
                        flag = "ab"
                    elif flag == "i":
                        flag = "a"
                    else:
                        1/0
                elif f == "b":
                    if flag == "a":
                        flag = "ab"
                    elif flag == "i":
                        flag = "b"
                    else:
                        1/0
                elif f == "i":
                    if flag == "a":
                        flag = "a"
                    elif flag == "b":
                        flag = "b"
                    elif flag == "i":
                        flag = "i"
                    else:
                        1/0
                else:
                    1/0
                #elif (f[0]=="a" and flag[0]=="b") or (f[0]=="b" and flag[0]=="a"):
                #    flag = "ab"
                readname2flag[read_name] = flag
            else:
                readname2flag[read_name] = flag
        flag_stat = {"a": 0, "b": 0, "ab": 0, "c": 0, "o": 0, "i": 0}
        flag_o_stat = {}
        for flag in readname2flag.values():
            if flag[0] == "o":
                flag_stat["o"] += 1
                if flag not in flag_o_stat:
                    flag_o_stat[flag] = 1
                else:
                    flag_o_stat[flag] += 1
            else:
                flag_stat[flag] +=1
        #intron_id chr_name start end strand a b ab c o o1_type o1_count o2_type o2_count  other_o(type=count;)
        o_stat_list = sorted(flag_o_stat.items(), key=lambda x: x[1], reverse=True)
        if len(o_stat_list) >= 1:
            o1_type, o1_count = o_stat_list[0]
        else:
            o1_type, o1_count = "", 0
        if len(o_stat_list) >= 2:
            o2_type, o2_count = o_stat_list[1]
        else:
            o2_type, o2_count = "", 0
        if len(o_stat_list) > 2:
            other_o = ";".join([s + "=" + str(v) for s, v in o_stat_list[2:]])
        else:
            other_o = ""
        a = flag_stat["a"]
        b = flag_stat["b"]
        ab = flag_stat["ab"]
        c = flag_stat["c"]
        o = flag_stat["o"]
        t = a + b + ab + c + o
        if t == 0: continue
        ls_t = a + b + 2*(ab + c + o)
        iratio = (a + b + 2*ab)/ls_t
        sratio = 2*c/ls_t
        oratio = 2*o/ls_t
        o1ratio = 2*o1_count/ls_t
        o2ratio = 2*o2_count/ls_t
        stat_result = [intron_id, chr_name, intron_start, intron_end, intron_strand,
                       intron_coverage_ratio,
                       a, b, ab, c, o, t, iratio, sratio, oratio, o1ratio, o2ratio,
                       o1_type, o1_count, o2_type, o2_count, other_o]
        yield(stat_result)

if __name__ == "__main__":
    main()
                           
"""
PS: Few introns are very short (e.g. 4bp). For example: specific introns in gene 
#AT1G02305
#AT1G63290
#AT1G63770
#AT2G39780
#AT2G41520
#AT3G01850
#AT4G01610
#AT5G51700
only two gene has `os`, although no genome annotation, but reads supported the alternative splicing
#AT2G35020, AT3G47560
"""
