import pysam 
import sys
import click

"""
Write by Jia Jinbu 2020.08.31

Remove the reads overlapped with specific features recorded in a bed file
from original bam file. You need use samtools index to index the resulted 
bam file.

See more by --help.
"""

@click.command()
@click.option('-i', '--inbam', help='Input bam file', 
                    required=True, type=click.Path(exists=True))
@click.option('-f', '--inbed', help='Input bed file', 
                    required=True, type=click.Path(exists=True))                    
@click.option('-o', '--out', help='Output Bam file', required=True)
def main(inbam, inbed, out):
    """
    Remove the reads overlapped with specific features recorded in a bed file
    from original bam file. You need use samtools index to index the resulted 
    bam file.
    """
    remove_read_ids = get_read_id_by_bed(inbam, inbed)
    filter_read_id_bam(inbam, out, remove_read_ids)

def get_read_id_by_bed(file_bam, file_bed):
    #0-based position
    def read_bed_file(filein):
        with open(filein) as f:
            for raw_line in f:
                d = raw_line.rstrip("\n").split()
                chr_, start, end = d[:3]
                start, end = int(start), int(end)
                yield((chr_, start, end))
    
    read_ids = set()
    regions = read_bed_file(file_bed)
    origin_bam_obj = pysam.AlignmentFile(file_bam, "rb")
    for chr_name, start, end in read_bed_file(file_bed):
        for read in origin_bam_obj.fetch(chr_name, start, end):
            read_ids.add(read.query_name)
    origin_bam_obj.close()
    return read_ids

def filter_read_id_bam(file_origin_bam, file_out_bam, ids, extract=0):
    origin_bam_obj = pysam.AlignmentFile(file_origin_bam, "rb")
    out_bam_obj = pysam.AlignmentFile(file_out_bam, "wb", template=origin_bam_obj)
    for line_data in origin_bam_obj:
        id_ = line_data.qname
        if id_ in ids:
            if not extract:
                continue
        else:
            if extract: continue
        out_bam_obj.write(line_data)
    origin_bam_obj.close()
    out_bam_obj.close()    

if __name__ == "__main__":
    main()
