import subprocess
import pysam
import sys
import time
import os
import tempfile
import joblib
import click
import multiprocessing

"""
Data : August 29, 2020

Author: Jia Jinbu

see more by --help.
"""

@click.command()
@click.option('-i', '--inbam', help='Input bam file', 
                    required=True, type=click.Path(exists=True))
@click.option('-s', '--inseq', help='Input read fastq file', 
                    required=True, type=click.Path(exists=True))                    
@click.option('-o', '--out', help='Output adapter information of each read', required=True)
@click.option('-t', '--threads', required=False, default=10, help='Number of threads to use. (default: 10)')
def main(inbam, inseq, out, threads):
    """
    \b
    Require: pysam, joblib pacage
    Require: blastn makeblast
    
    \b
    Output:
    read_core_id              7fab6cbf-a9ed-4427-951f-741515ddba0b,chr1,7026,8687
    read_align_strand         +            
    rna_strand                -            
    read_length               839          # read total length
    primer_type	              R-F          
    genome_align_start        71
    genome_align_end          820
    primer_score              1
    polyA_type                T
    f_primer_type             R
    f_primer_start            1
    f_align_end               70
    r_primer_type             F
    r_primer_start            1
    r_align_start             821

    \b
    primer means adapter.
    primer_5p = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
    primer_3p = "AAGCAGTGGTATCAACGCAGAGTACATTGATGGTGCCTACAG"
    primer_type: `F` means `primer_5p`, `R` means `primer_3p`
    
    
    f_ indicate the 5' read, r_ indicate the 3' read
    Note: f may be F or R, r may be F or R.
    
    \b
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

    coordinate of read: f_align_end, r_align_start,
           genome_align_start, genome_align_end

    #coordinate of adapter: f_primer_start, r_primer_start
    #if r_primer_start if 4, the first two base of r_primer is not mapped to read
    
    \b
    #for example:
                                                 AAGGGAGAGAGAG
                                                    ||||||||||
    #read1      f_primer--------------------------->GGAGAGAGAG

    \b
    one end primer type
    (a) no alignment
           ["N", ["", "", -1, -1, -1, -1, ""]]
    (b) if "F" primer, if primer_start > primer_f_max_start: UF
    (c) if "R" primer, if primer_start > primer_r_max_start2: UUR
           if primer_r_max_start2 >= primer_start > primer_r_max_start1: UR

    read adapter type based on the primer types of both ends:
    for example: F-R indicate the 5' primer is F, the 3' primer is F
    Each kind of primer type is assigned a score, and also can predict the 
    rna_strand by primer type, see `type2strand` variable value: The score 1 
    or 2 is reliable.
    """
    
    max_threads = multiprocessing.cpu_count()
    if threads > max_threads:
        max_threads = max_threads
    
    pad_length = 20
    min_primer_align_length = 15
    primer_r_max_start1 = 8
    primer_r_max_start2 = 12
    primer_f_max_start = 1
    match = 1
    mis = -1.5
    primer_5p = "AAGCAGTGGTATCAACGCAGAGTACATGGG"
    primer_3p = "AAGCAGTGGTATCAACGCAGAGTACATTGATGGTGCCTACAG"
        
    seq_lists = extract_bam_clip_fasta_seq_split(inbam, inseq, split=threads)
    threads = len(seq_lists)
    
    with joblib.Parallel(threads) as pool:  
        
        def _core(seqs):
            blastoutput = blast_nanopore_adapter(seqs, primer_5p, primer_3p)
            return extract_read_primer_type_from_balst(blastoutput, 
                       None, 
                       False,
                       min_primer_align_length,
                       primer_r_max_start1,
                       primer_r_max_start2,
                       primer_f_max_start,
                       match,
                       mis)
                        
        res = pool(
                joblib.delayed(_core)(seqs)
                for seqs in seq_lists
            )
            
    with open(out, 'w') as o:
        header = ("read_core_id\tread_align_strand\trna_strand\tread_length\t"
                  "primer_type\tgenome_align_start\tgenome_align_end\tprimer_score\t"
                  "polyA_type\tf_primer_type\tf_primer_start\tf_align_end\t"
                  "r_primer_type\tr_primer_start\tr_align_start\n")
        o.write(header)
        for r in res:
            o.write(r)
    
######### 1. clip sequence function ###############
def iter_group_by_id(d, id_index=0, sort_flag=False):
    '''
    Split a list of lists based on id.

    Input: d is a list of lists
    [["id1", 1, 2], ["id2", 2, 3], ["id2", 4, 5], ["id3", 4, 5]]

    Note: The core function is _iter_id,
    two `id2` is continous, or they would not be grouped, like the 
    characterize of `uniq` command. In this case, you can first sort the list.
    _iter_id(sorted(d, key=lambda x: x[0]), 0). Thus, iter_id provide 
    sort_flag option. You can set sort_flag as True if you want to sort the 
    list by id first. Default: sort_flag is False.

    id_index is the index of id in each child list.

    Output:
    ["id1", [["id1", 1, 2]]],
    ["id2", [["id2", 2, 3], ["id2", 4, 5]]],
    ["id3", [["id3", 4, 5]]]
    '''
    
    def _iter_id(d, id_index=0):
        '''
        
        '''
        record_id = ""
        record_d = []
        for l in d:
            if record_id != l[id_index]:
                if record_d:
                    yield([record_id, record_d])
                record_d = []
                record_id = l[id_index]
            record_d.append(l)
        yield([record_id, record_d])
        
    if sort_flag: d = sorted(d, key=lambda x: x[id_index])
    return _iter_id(d, id_index)

def read_fasta_to_dict(filein):
    """
    !!!The function in included in both adapterFinder.py and 
    pacbio_find_polyA.py. They are same function, but haven't
    be put in a module to keep each script can be run independently.
    If you want to modify one of them, please modify them at the 
    same time.
    
    Input: fasta files
    
    Output: dict
    key: sequence name
         #>seq1 npr1
         #run: l[1:].split()[0]
         #will get `seq1`
    value: sequence
    """
    
    id2seq = {}
    
    seq_id, seq = "", ""
    for line_num, l in enumerate(open(filein)):
        l = l.strip()
        if l.startswith(">"):
            if seq_id:
                id2seq[seq_id] = seq
            seq = ""
            seq_id = l[1:].split()[0]
        else:
            seq += l
    if seq_id:
        id2seq[seq_id] = seq
        
    return(id2seq)

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
    
def iter_bam_clip_seq(filein_bam, filein_seq, pad_length=20):

    '''
    !!!The function in included in both adapterFinder.py and 
    pacbio_find_polyA.py. They are same function, but haven't
    be put in a module to keep each script can be run independently.
    If you want to modify one of them, please modify them at the 
    same time. 
    
    genome 5'------------------------------------------------3'
    mapping region     |||||||||||||||||||||||||||
    read         ---------------------------------------
      5' clip    <<<<<<---                     --->>>>>> 3' clip
                       pad                     pad
    The 5' or 3' of clip sequence is based on the genome fasta 
    direction, not the original read direction.

    Extract 5' and 3' clip sequence, due to the clip sequence is removed 
    from bam file when alignment contain hard clip, thus you also need 
    to provide original sequence file.

    Input: 
    filein_bam: aligned file
    filein_seq: The origin seq file used to generate filein_bam.

    Output:
    iterator, each element is :
    [read, [read_strand, seq_length, read_name, \
            left_clip_name, left_clip_seq, \
            right_clip_name, right_clip_seq]]        
    
    read_id:        00000449-e191-4eda-92e7-c4b2843daba2
    read_core_id:   00000449-e191-4eda-92e7-c4b2843daba2,chr4,8523083,8523682
    read_name:      00000449-e191-4eda-92e7-c4b2843daba2,chr4,8523083,8523682,+,660,58,55,20
    left_clip_name: 00000449-e191-4eda-92e7-c4b2843daba2,chr4,8523083,8523682,+,660,58,55,20,5
    right_clip_name:00000449-e191-4eda-92e7-c4b2843daba2,chr4,8523083,8523682,+,660,58,55,20,3
    

    read:           the object of each line in bam file, generated by pysam
    read_id:        query_name in bam file, e.g.
    read_core_id:
        Due to one read may be mapped to multiple positions, thus only read
        name cannot uniquly represent the alignment, so we combine the read
        name (the first column in bam file) and the mapping position as the 
        read_core_id by `,`: `read_name,chr_name,start,end`.

        Note: Logically, the read_core_id may be also not uniquely because: 
        read contain duplicated regions, all regions can mapped to same 
        positions in genome. And `read_name,chr_name,start,end` didn't 
        consider strand information. The output didn't remove duplicated 
        read_core_id row, thusbefore generating the final result in this 
        pipeline, we need to removed them in specific step (now is in merge 
        step).
    read_strand:    read alignment direction `+` or `-`
    seq_length:     original sequence length
    read_name:      final read_name
        `read_name,chr_name,start,end,strand,seq_length,left_clip_length,
        right_clip_length,pad_length`.

        Note: This is used to record the clip information, but the final 
        read_name has been used in the downstream pipeline, thus you can 
        just replace it by read_core_id, which is exactly used in the 
        downstream pipeline.
    left_clip_name:      `read_name,5`
    left_clip_seq:  5' clip sequence contain adjecent pad_length mapping bases
    right_clip_name:     `read_name,3`
    right_clip_seq: 3' clip sequence contain adjecent pad_length mapping bases

    For developer:
    Note: The tenth column of bam file store query sequence, you can use
    read.query_sequence (pysam package) to extract it. In minimap2 alignment
    file, this column store the reverse complement of read sequence if 
    the alignment direction is minus (read.is_reverse is True). That means
    the sequence is in same direction with genome sequence, but it is not 
    genome sequence if mismatch or indel exist. If you want to extract genome 
    sequence, you need to use original genome fasta files. If the alignment 
    contain soft cliped bases, but didn't contain hard cliped bases, thus if 
    the alignment contain hard cliped bases (extract from read.cigartuples), 
    you need to extract original read sequence from fastq files, and reverse 
    complement it if the alignment direction is minus.
    
    left_clip_seq and right_clip_seq may be empty string.
    '''

    def get_clip_seq(read, origin_seqs, pad_length=20):
        #supposed that the first and the last type in cigartuples is in [0, 4, 5]
        read_strand = "-" if read.is_reverse else "+"
        left_clip_type, left_clip_length = read.cigartuples[0]
        right_clip_type, right_clip_length = read.cigartuples[-1]
        if left_clip_type == 0:
            left_clip_length = 0
        if right_clip_type == 0:
            right_clip_length = 0
        left_clip_length += pad_length
        right_clip_length += pad_length
    
        seq = read.query_sequence
        if left_clip_type == 5 or right_clip_type == 5:
            #if Hard clip, need origin_seqs to extract origin seq
            #In this case, if alignment direction is minus, 
            #need to reverse the origin seq
            seq = origin_seqs[read.query_name]
            if read_strand == "-": seq = revcom(seq)
        seq_length = len(seq)

        left_clip_seq = seq[:left_clip_length]
        left_clip_seq = revcom(left_clip_seq)
        right_clip_seq = seq[(-right_clip_length):]
    
        read_name = ",".join([read.query_name, 
                              read.reference_name,
                              str(read.reference_start + 1), 
                              str(read.reference_end),
                              read_strand,
                              str(seq_length),
                              str(left_clip_length),
                              str(right_clip_length),
                              str(pad_length)])
        left_clip_name, right_clip_name = [read_name + ",5", read_name + ",3"]
        return [read_strand, seq_length, read_name, left_clip_name, left_clip_seq, right_clip_name, right_clip_seq]

    origin_seqs = read_fasta_to_dict(filein_seq)

    with pysam.AlignmentFile(filein_bam, "rb") as bam_obj:
        for read in bam_obj.fetch():
            #remove unampped read
            if read.is_unmapped:
                continue
            yield([read, get_clip_seq(read, origin_seqs, pad_length)])
            
def extract_bam_clip_fasta_seq(filein_bam, filein_seq, pad_length=20):
    seqs = ""
    for read, (read_strand, seq_length, read_name, left_name, left_clip_seq, right_name, right_clip_seq) \
                in iter_bam_clip_seq(filein_bam, filein_seq, pad_length):
        if left_clip_seq:
            seqs += f">{left_name}\n{left_clip_seq}\n"
        if right_clip_seq:
            seqs += f">{right_name}\n{right_clip_seq}\n"
    return(seqs)
    
def cut_lists(data, split_number=10):
    
    #data is a list, split data into split_number.
    results = []
    if not data:
        results.append([])
    else:
        each_part_element_number = (len(data) - 1)//split_number + 1
        ls_d = []
        for i, d in enumerate(data):
            ls_d.append(d)
            if i % each_part_element_number == each_part_element_number - 1:
                results.append(ls_d)
                ls_d = []
        if ls_d:
            results.append(ls_d)
    return results
      
def extract_bam_clip_fasta_seq_split(filein_bam, filein_seq, pad_length=20, split=10):
    seq_list = []
    for read, (read_strand, seq_length, read_name, left_name, left_clip_seq, right_name, right_clip_seq) \
                in iter_bam_clip_seq(filein_bam, filein_seq, pad_length):
        seqs = ""
        if left_clip_seq:
            seqs += f">{left_name}\n{left_clip_seq}\n"
        if right_clip_seq:
            seqs += f">{right_name}\n{right_clip_seq}\n"
        seq_list.append(seqs)
    seqs_list = []
    for d in cut_lists(seq_list, split):
        seqs_list.append("".join(d))
    return seqs_list

####### 2. find adapter by blast from clip sequence ##################
type2strand = {
            #key: primer_type  
            #value: [rna_strand, score]
            #score 1 or 2 is reliable.
            #1 reliable F, R:  F-R, R-F
            #2 reliable R, but R not reliable
            #3 not reliable R
            #4 two side is contradictory
            #5 contain N: F-N, N-N
        
            "F-R": ["+", 1],
            "F-UR": ["+", 2],
            "F-UUR": ["+", 3],
            "F-UF": ["+", 4],
            "F-F": ["+", 5],
            "F-N": ["+", 5],
        
            "R-F": ["-", 1],
            "R-R": ["-", 5],
            "R-N": ["-", 2],
            "R-UF": ["-", 2],
            "R-UUR": ["-", 2], #2 or 3
            "R-UR": ["-", 3], #2 or 3
        
            "N-F": ["-", 5],
            "N-R": ["+", 2],
            "N-N": ["+", 5], #5 or 6
            "N-UF": ["-", 5],
            "N-UR": ["+", 3],
            "N-UUR": ["+", 5],
        
            "UF-F": ["-", 4],
            "UF-R": ["+", 2],
            "UF-N": ["+", 5],
            "UF-UF": ["+", 4],
            "UF-UR": ["+", 3],
            "UF-UUR": ["+", 4],
        
            "UR-F": ["-", 2],
            "UR-R": ["+", 3],
            "UR-N": ["-", 3],
            "UR-UF": ["-", 3],
            "UR-UR": ["-", 4],
            "UR-UUR": ["-", 3],
        
            "UUR-F": ["-", 3],
            "UUR-R": ["+", 2],
            "UUR-N": ["-", 5],
            "UUR-UF": ["-", 4],
            "UUR-UR": ["+", 3],
            "UUR-UUR": ["-", 4]
}

def blast(seqs, db_path=None, lib_seqs=None):
    """
    Run blast and return result in outfmt 6 in string.
    
    Input:
    seqs:  string, fasta format
    you can provide a blast db path via set db_path, 
    or you can directly provide a string contain lib seqeunces in fasta 
    format. If you provide db_path, it didn't perform makeblastdb.
    
    Output:
    String. blast results in outfmt 6.
    """
    def _blast(seqs, db_path):
        proc = subprocess.run(
                ['blastn', '-query', "-", '-db', db_path,
                 '-word_size', '6', '-outfmt', "6"],
                input=str.encode(seqs), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return(proc.stdout.decode())
    
    def _makeblastdb(fasta_file):
        proc = subprocess.run(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl'])
    
    if db_path:
        return _blast(seqs, db_path)
    else:
        output = ""
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_fasta_file = os.path.join(tmp_dir, "lib.fa")
            with open(tmp_fasta_file, 'w') as o:
                o.write(lib_seqs)
            _makeblastdb(tmp_fasta_file)
            output = _blast(seqs, tmp_fasta_file)
        return output
    
def blast_nanopore_adapter(seqs, 
                           primer_5p="AAGCAGTGGTATCAACGCAGAGTACATGGG", 
                           primer_3p="AAGCAGTGGTATCAACGCAGAGTACATTGATGGTGCCTACAG"):
    """
    Input:
    seqs is the output of extract_bam_clip_fasta_seq function,
        the usage in main function.
        the seqs is fasta seqs of clip sequences. 
        format in extract_bam_clip_fasta_seq function.
    
    Output: string, blast result in outfmt 6.
            The primer name is F and R.
    """
    primer_5p_rev = revcom(primer_5p)
    primer_3p_rev = revcom(primer_3p)
    primer_seqs = f">F\n{primer_5p_rev}\n>R\n{primer_3p_rev}\n"
    return blast(seqs, lib_seqs=primer_seqs)

def iter_blast_by_read(filein, is_file=True):
    '''
    Input:
    filein: can be blast output file (outfmt=6) when `is_file` is set as True (
            default). When `is_file` is set as False, the filein is string 
            represent the content store in blast output file. This usually be 
            generated by subprocess which run blast and capture the output (see
            `blastn` function in this script).
    
    Output:
    each element represent one read's information.
    #both 5' clip seq and 3' clip seq can blast to adapter
    [read_name, [[left_clip_name, read_clip_name_split, blast_aligns], 
                 [right_clip_name, read_clip_name_split, blast_results]]] or 
    ##only 5' clip seq or 3' clip seq can blast to adapter
    [read_name, [[left_clip_name or right_clip_name, read_clip_name_split, blast_results]]]
    
    read_name: 00000449-e191-4eda-92e7-c4b2843daba2,chr4,8523083,8523682,+,660,58,55,20
    
    origin query name (read_clip_name) in filein: 
       5 clip: 00000449-e191-4eda-92e7-c4b2843daba2,chr4,8523083,8523682,+,660,58,55,20,5
       3 clip: 00000449-e191-4eda-92e7-c4b2843daba2,chr4,8523083,8523682,+,660,58,55,20,3
    
    read_clip_name_split: = read_clip_name.split(",")
    
    blast_aligns is a list of the blast result, one element represent one line in blast 
        result file: 
        [read_clip_name, primer_name, start, end, primer_start, primer_end, strand]
    start, end is the position of clip seq, primer_start, primer_end is the position of primer
    primer_name: F or R
    '''
    
    def iter_blast_by_query(filein, is_file=True):
        
        def _iterline_blast_file(filein):
            if is_file:
                IN = open(filein)
            else:
                IN = filein.splitlines()
            for l in IN:
                d = l.rstrip("\n").split("\t")
                query = d[0]
                primer_name = d[1]
                start, end = int(d[6]), int(d[7])
                primer_start, primer_end = int(d[8]), int(d[9]) 
                primer_align_length = primer_end - primer_start + 1
                strand = "+"
                if primer_start > primer_end:
                    strand = "-"
                    primer_start, primer_end = primer_end, primer_start
                yield [query, primer_name, start, end, primer_start, primer_end, strand]
            if is_file:
                IN.close()
    
        return iter_group_by_id(_iterline_blast_file(filein))
        
    def _iterread_balst_file(filein, is_file=True):
        for read_clip_name, blast_aligns in iter_blast_by_query(filein, is_file):
            read_clip_name_split = read_clip_name.split(",")
            read_name = ",".join(read_clip_name_split[:-1])
            yield([read_name, read_clip_name, read_clip_name_split, blast_aligns])
                
    return iter_group_by_id(_iterread_balst_file(filein, is_file))

def extract_end_primer(blast_aligns, 
                       min_primer_align_length = 15,
                       primer_r_max_start1 = 8,
                       primer_r_max_start2 = 12,
                       primer_f_max_start = 1):
    """
    Input:
    blast_aligns: the list represent blast results of one side clip sequence
    each element is [read_clip_name, primer_name, start, end, primer_start, primer_end, strand]
    
    Output:
    [primer_type, [read_clip_name, primer_name, start, end, primer_start, primer_end, strand]]
    After filtering, the first alignment result based on the start position of clip sequence
    was return.
    After filtering, if no alignment exist, return:
    ["N", ["", "", -1, -1, -1, -1, ""]]
                       
    1. filering out:
    (a) strand (alignment direction) == "-"
    (b) alignment length (in primer) < min_primer_align_length(default 18)
    2. sorted by start position of clip sequence
    3. Calcute primer type, only care about the first alignment
    (a) no alignment after filtering, return 
           ["N", ["", "", -1, -1, -1, -1, ""]]
    (b) if "F" primer, if primer_start > primer_f_max_start: UF
    (c) if "R" primer, if primer_start > primer_r_max_start2: UUR
           if primer_r_max_start2 >= primer_start > primer_r_max_start1: UR
                       
    Note: the second element of output is a list which is also a element of the input 
          blast_aligns. They are share the same memory address, so change one, another 
          is also changed.
    """
    
    #1. filering out:
    filtered_blast_aligns = []
    for d in blast_aligns:
        na, primer_name, start, end, primer_start, primer_end, strand = d
        if strand == "-": continue
        if primer_end - primer_start + 1 < min_primer_align_length: continue
        filtered_blast_aligns.append(d)
        
    #2. sorted by start position of clip sequence
    filtered_blast_aligns.sort(key=lambda x: x[2])
    
    #3. Calcute primer type
    if len(filtered_blast_aligns) == 0:
        primer_d = ["", "", -1, -1, -1, -1, ""]
        primer_type = "N"
    else:
        primer_d = filtered_blast_aligns[0]
        primer_type = primer_d[1]
        primer_start = primer_d[4]
        if primer_type == "F" and primer_start > primer_f_max_start:
            primer_type = "UF"
        if primer_type == "R" and primer_start > primer_r_max_start1:
            if primer_start > primer_r_max_start2:
                primer_type = "UUR"
            else:
                primer_type = "UR"
    return [primer_type, primer_d]

def extract_read_primer_type_from_balst(filein_blast, 
                       fileout=None, 
                       input_is_file=True,
                       min_primer_align_length = 15,
                       primer_r_max_start1 = 8,
                       primer_r_max_start2 = 12,
                       primer_f_max_start = 1,
                       match = 1,
                       mis = -1.5):
    """
    Input: 
    filein_blast: result file of clip sequence blast aganist adapter sequence.
                  you can provide string if set input_is_file = False
    Output:
    fileout: if None, not write to fileout, but return string (no header line)
             if fileout, write to fileout, have header line
    columns:
       header = ("read_core_id\tread_align_strand\trna_strand\tread_length\t"
                 "primer_type\tgenome_align_start\tgenome_align_end\tprimer_score\t"
                 "polyA_type\tf_primer_type\tf_primer_start\tf_align_end\t"
                 "r_primer_type\tr_primer_start\tr_align_start\n")
    
    see extract_end_primer for detail.
    """
    def _core():
        for read_name, pair_end in iter_blast_by_read(filein_blast, input_is_file):
            read_core_id = ",".join(read_name.split(",")[:4])
            #please see the data structure of pair_end in `iter_blast_by_read`
            #[[left_clip_name, read_clip_name_split, blast_aligns], 
            #     [right_clip_name, read_clip_name_split, blast_results]]
            #one element represent 5' clip alignment information
            #one element represent 3' clip alignment inforamtion.
            #if only 3' clip or 5' clip blast to adapter, only one element
            
            #1. get basic alignment information of read from the first element
            #pair_end[0][2] is read_clip_name.split(",")
            #00000449-e191-4eda-92e7-c4b2843daba2,chr4,8523083,8523682,\
            #+,660,58,55,20,5   --> 4:10
            read_align_strand, read_length, \
                left_length, right_length, \
                clip_inner_length, first_end_type = pair_end[0][2][4:10]
            read_length = int(read_length)
            left_length = int(left_length)
            right_length = int(right_length)
            clip_inner_length = int(clip_inner_length)
    
            #2. extract adapter information based on blast alignments
            left_primer_info = extract_end_primer(pair_end[0][3], \
                                    min_primer_align_length, \
                                    primer_r_max_start1, \
                                    primer_r_max_start2, \
                                    primer_f_max_start)
            if len(pair_end) == 2:
                right_primer_info = extract_end_primer(pair_end[1][3], \
                                        min_primer_align_length, \
                                        primer_r_max_start1, \
                                        primer_r_max_start2, \
                                        primer_f_max_start)
            else:
                right_primer_info = ["N", ["", "", 0, 0, 0, 0, ""]]
            
            #although the 5' clip should be in the first element, we also 
            #check it, so you don't need to care about it.
            if first_end_type == "3":
                left_primer_info, right_primer_info = right_primer_info, left_primer_info
    
            #3. Based on the adapter inforamtion in both side, 
            #calculate adapter type of read, whether is reliable, strand
            left_primer_type = left_primer_info[0]
            left_primer_start = left_primer_info[1][4]
            left_align_start = left_primer_info[1][2]
            
            right_primer_type = right_primer_info[0]
            right_primer_start = right_primer_info[1][4]
            right_align_start = right_primer_info[1][2]
            #F, R, N, UF, UR, UUR
            rna_strand, primer_score = type2strand[left_primer_type + "-" + right_primer_type]
            #if strand is `-`, reverse left and right
            #genome 5’--------------------------------------------------------3'
            #read1   + + A   5'adapter-------------------AAAA>3'adapter  
            #read2   + - T   3'adapterTTTT------------------->5'adapter 
            #read3   - + T   5'adapter<-------------------TTTT3'adapter
            #read4   - - A   3'adapter<AAAA-------------------5'adapter
            # (read_align_strand rna_strand polyT/A) --> alignment direction
            #5' adapter in left, rna_strand is "+", or rna_strand is "-"
            if read_align_strand == rna_strand:
                polyA_type = "A"
            else:
                polyA_type = "T"
            
            #f_ indicate the 5' read, r_ indicate the 3' read
            #left_ indicate the 5' genome, right_ indicate the 3' genome
            #genome 5’--------------------------------------------------------3'
            #read1      left_primer--------------------------->right_primer
            #           left_length                            right_length
            #              f_primer                            r_primer
            #                 genome_align_start   genome_align_end
            #read2      left_primer<---------------------------right_primer
            #           left_length                            right_length
            #              r_primer                            f_primer
            #                 genome_align_end     genome_align_start
            #   --> indicate alignment direction (read_align_strand)
            #    left_length and right_length include clip_inner_length
            
            #coordinate of adapter: f_primer_start, r_primer_start
            #coordinate of read: f_align_end, r_align_start,
            #       genome_align_start, genome_align_end
            if read_align_strand == "+":
                f_primer_type = left_primer_type
                f_primer_start = left_primer_start
                r_primer_type = right_primer_type
                r_primer_start = right_primer_start
                genome_align_start = left_length - clip_inner_length + 1
                genome_align_end = read_length  - right_length + clip_inner_length
                f_align_end = left_length - left_align_start + 1
                r_align_start = read_length - right_length + right_align_start
            else:
                f_primer_type = right_primer_type
                f_primer_start = right_primer_start 
                r_primer_type = left_primer_type
                r_primer_start = left_primer_start
                genome_align_start = right_length - clip_inner_length + 1
                genome_align_end = read_length - left_length + clip_inner_length
                f_align_end = right_length - right_align_start + 1
                r_align_start = read_length - left_length + left_align_start
            
            primer_type = f_primer_type + "-" + r_primer_type
            
            #4. Output
            r = [read_core_id, read_align_strand, rna_strand,  str(read_length), primer_type, 
                 str(genome_align_start), str(genome_align_end),
                 str(primer_score), polyA_type, f_primer_type, str(f_primer_start), str(f_align_end), 
                 r_primer_type, str(r_primer_start), str(r_align_start)]
            yield("\t".join(r) + "\n")
    
    header = ("read_core_id\tread_align_strand\trna_strand\tread_length\t"
              "primer_type\tgenome_align_start\tgenome_align_end\tprimer_score\t"
              "polyA_type\tf_primer_type\tf_primer_start\tf_align_end\t"
              "r_primer_type\tr_primer_start\tr_align_start\n")
    
    if fileout:
        with open(fileout, 'w') as o:
            o.write(header)
            for l in _core():
                o.write(l)
    else:
        return("".join(list(_core())))

if __name__ == '__main__':
    main()

