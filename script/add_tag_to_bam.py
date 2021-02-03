#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
@Author       : windz
@Date         : 2020-04-14 10:46:37
LastEditTime : 2020-12-22 11:02:11
@Description  : add polya length and remove wrong strandness reads
                tags:
                    'pa': polya tail length
                    'gi': gene_id
                    'mi': mRNA intron count
                    'ir': retention introns
'''


import pysam
import pandas as pd
import numpy as np
import click


import logging
logging.basicConfig(level=logging.DEBUG,  
                    format='%(asctime)s %(filename)s: %(message)s',  
                    datefmt='%m/%d/%Y %I:%M:%S %p',
                    )


def load_read_info(infile):
    '''
    load read info from FLEP-seq read_info result

    Return:
        '55c7bb93-ebb1-479f-a8ad-73168cbd51ef': {
            'read_strand': '+',
            'mRNA': 'AT1G01010.1',
            'mRNA_intron_num': 5,
            'mRNA_strand': '+',
            'retention_introns': None,
            'retention_intron_num': 0,
            'polya_length': 0.84,
            'type': 'splicing_intermediate'
            }, ...
    '''
    read_info_df = pd.read_csv(
    infile, sep='\t',
    usecols=[
        'read_core_id', 'read_strand', 'mRNA', 'mRNA_intron_num',
        'mRNA_strand', 'retention_introns', 'retention_intron_num',
        'polya_length', 'type'
    ],
    dtype={'retention_introns': str}
)
    # only use column: 
    # read_core_id, read_strand, mRNA, mRNA_intron_num, mRNA_strand, 
    # retention_introns, retention_intron_num, polya_length, type

    read_info_df['retention_introns'] = read_info_df['retention_introns'].replace({np.nan: None})
    read_info_df['read_id'] = read_info_df['read_core_id'].map(lambda x: x.split(',')[0])
    read_info_df.drop_duplicates(subset=['read_id'], keep=False, inplace=True)
    read_info_df.set_index('read_id', inplace=True)
    read_info_dict = read_info_df.iloc[:, 1:].to_dict(orient='index')
    
    return read_info_dict


@click.command()
@click.option('-i', '--infile', required=True)
@click.option('-o', '--outfile', required=True)
@click.option('--read_info', required=True)  # 提供polya长度和3'adapter信息
def main(infile, outfile, read_info):
    logging.info('Load read_info')
    read_info_dict = load_read_info(read_info)
    logging.info('Load read_info done.')

    logging.info('Start main function')
    inbam = pysam.AlignmentFile(infile, 'rb') 
    outbam = pysam.AlignmentFile(outfile, 'wb', template=inbam)

    for read in inbam:
        if not read.is_unmapped and read.query_name in read_info_dict:
            read_info_dict_ = read_info_dict[read.query_name]
            read_is_reverse = False if read_info_dict_['mRNA_strand'] == '+' else True
            if read.is_reverse is not read_is_reverse:
                if read.is_reverse:
                    read.flag += -16
                else:
                    read.flag += 16
            
            retention_introns = read_info_dict_['retention_introns'] if read_info_dict_['retention_introns'] is not None else 'None'
            
            read.set_tag('gi', read_info_dict_['mRNA'])  # gi: gene id
            read.set_tag('mi', read_info_dict_['mRNA_intron_num'])  # mi: mRNA intron count
            read.set_tag('pa', read_info_dict_['polya_length'])  # pa: polya length
            read.set_tag('ri', retention_introns)  # ir: retention introns
            read.set_tag('ty', read_info_dict_['type'])
            outbam.write(read)

    
    inbam.close()
    outbam.close()
    # build bam index
    pysam.index(outfile, '-@ 60')

if __name__ == "__main__":
    main()
    