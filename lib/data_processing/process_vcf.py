from argparse import ArgumentParser
import os
import sys
import re
from tqdm import tqdm
from lib.utils import general as g
from lib.config.config_class import VCFConfigClass

vcf_config = VCFConfigClass()

def check_vcf(vcf_path):
    with g.writing(vcf_path) as vf:
        vcf_dict = None
        for line in tqdm(vf,desc='find header variants'):
            if line.startswith('#CHROM'):
                items = g.line_to_items(line,'\t')
                vcf_dict = g.list_to_dict(items)
                break
        assert type(vcf_dict) is not None, "vcf file must have header row"
        for line in tqdm(vf,desc='convert to true format'):
            items = g.line_to_items(line,'\t')
            chrom = items[vcf_dict[vcf_config.chrom_col]]
            if chrom.startswith('chr') == False:
                chrom = 'chr'+chrom
                vf.rea()
        pass