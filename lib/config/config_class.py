import numpy as np

class VCFConfigClass():
    def __init__(self):
        self.vcf_header_line_startswith = '#CHROM'  # header line of data
        self.vcf_header_sample_startindex = 9
        self.vcf_data_split_params = '\t'
        self.header_dict_sample_key = 'samples'
        self.chrom_col = '#CHROM'
        self.position_col = 'POS'
        self.snp_col = 'ID'
        self.ref_col = 'REF'
        self.alt_col = 'ALT'
        self.qual_col = 'QUAL'
        self.filter_col = 'FILTER'
        self.info_col = 'INFO'
        self.format_col = 'FORMAT'
        self.chrom_values = list(np.hstack([np.arange(1,23).astype(str),np.array(['X','Y'])]))

