import os
import sys
import click

"""
Data : August 27, 2020

Author: Jia Jinbu

Convert fastq files with suffirx ".fastq" in specific directory to one fasta file.

see more in help.
"""

@click.command()
@click.option('-i', '--indir', help='Input directory containg fastq files with suffirx ".fastq"', 
                    required=True, type=click.Path(exists=True))
@click.option('-o', '--out', help='Output fasta file', required=True)
def main(indir, out):
    """
    Convert fastq files with suffirx ".fastq" in specific directory to one fasta file.
    """
    FASTQ_SUFFIX = ".fastq"
    fastq_dir_to_fasta(indir, out, FASTQ_SUFFIX)
             
def list_all_fastq(in_dir, suffix=".fastq", add_path=1):
    files = [filename for filename in os.listdir(in_dir) if filename.endswith(suffix)]
    if add_path:
        files = [os.path.join(in_dir, filename) for filename in files]
    return files
    
def iter_fastq(fastq_file):
    with open(fastq_file) as f:
        while 1: 
            seq_name = f.readline().rstrip()
            if not seq_name: break
            seq = f.readline().rstrip()
            nm = f.readline().rstrip()
            quality = f.readline().rstrip()
            yield((seq_name, seq, nm, quality))
        
def fastqs_to_fasta(file_fastqs, file_fasta):
    with open(file_fasta, 'w') as o:
        for file_fastq in file_fastqs:
            for seq_name, seq, nm, quality in iter_fastq(file_fastq):
                seq_name = seq_name[1:].split()[0]
                o.write(">%s\n%s\n" % (seq_name, seq))

def fastq_dir_to_fasta(fastq_dir, file_fasta, suffix=".fastq"):
    file_fastqs = list_all_fastq(fastq_dir, suffix)
    fastqs_to_fasta(file_fastqs, file_fasta)


if __name__ == '__main__':
    main()
    
    






        
        