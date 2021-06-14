import numpy as np
from ..utils import general as g

class ManiConfig():
    def __init__(self) -> None:
        self.start_header_end_line = '[Assay]'
        self.IlmnStrand_col = 'IlmnStrand'
        self.Chr_col = 'Chr'
        self.MapInfo_col = 'MapInfo'
        self.SNP_col = 'SNP'
        self.SourceStrand_col = 'SourceStrand'
        self.SourceSeq_col = 'SourceSeq'
        self.RefStrand_col = 'RefStrand'
        self.start_data_end_line = '[Controls]'
        self.data_line_split_params = ','
        self.alleles_split_params = '/'
        self.alleles_position_0 = ['I','D']
        self.alleles_missing = ['-']
        self.alleles_character_replace_missing = ''
        self.sourceSeq_split_params = '[\[\/\]]'        

mani_config = ManiConfig()

class VCFInfoConfig():
    def __init__(self) -> None:
        self.af = 'AF'

class VCFConfig():
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
        self.genotype_split_params = '|'
        self.info_key = VCFInfoConfig()
        self.chrom_values = list(np.hstack([np.arange(1,23).astype(str),np.array(['X','Y'])]))

vcf_config = VCFConfig()

class LegendConfig():
    def __init__(self) -> None:
        self.unobserve = '0'
        self.observe = '1'
        self.hap_tail = '.hap.gz'
        self.legend_tail = '.legend.gz'
        #legend
        self.snp = 'id'
        self.chrom = 'chr'
        self.position = 'position'
        self.ref = 'ref'
        self.alt = 'alt'
        self.af = 'af'
        self.maf = 'maf'
        self.array_marker_flag = 'array_marker_flag'
        self.legend_header = [self.snp,self.chrom,self.position,self.ref,self.alt,self.af,self.maf,self.array_marker_flag]
        self.legend_split_params = ' '
        #hap
        self.hap_split_params= ' '
        #Sample
        self.sample_tail = '.sample'        
        self.sample_first_line = '0 0 0'
        self.sample_id_01 = 'ID_1'
        self.sample_id_02 = 'ID_2'
        self.sample_split_params = ' '
        self.sample_missing = 'missing'
        self.sample_header = [self.sample_id_01,self.sample_id_02,self.sample_missing]
    
    @property
    def legend_header_line(self)->str:
        return self.legend_split_params.join(self.legend_header)

    @property
    def sample_header_line(self)->str:
        return self.sample_split_params.join(self.sample_header)

    def get_hap_file_name(self,output_prefix:str) -> str:
        return output_prefix+self.hap_tail
    
    def get_legend_file_name(self,output_prefix:str) ->str:
        return output_prefix+self.legend_tail
    
    def get_legend_dict(self,snp, chrom, position, ref, alt, af, flag, nb_alt:int)->dict:
        reDict = {}
        reDict[self.snp] = snp
        reDict[self.chrom]=chrom
        reDict[self.position] = position
        reDict[self.ref] = ref
        reDict[self.alt] = alt
        reDict[self.af] = af
        reDict[self.array_marker_flag] = flag
        if reDict[self.snp] == '.':
            reDict[self.snp] = '{:s}:{:s}:{:s}:{:s}'.format(chrom, position, ref, alt)
        elif nb_alt >= 2:
            reDict[self.snp] += ':{:s}:{:s}'.format(ref, alt)
        return reDict
    
    def get_legend_value_line(self,legend_dict:dict)->str:
        values = list(map(legend_dict.get,self.legend_header))
        return self.legend_split_params.join(values)

    def get_sample_file_name(self,output_prefix:str) ->str:
        return output_prefix+self.sample_tail
    
    def get_sample_dict(self,id_1:str,id_2:str,missing:str) -> dict:
        reDict = {}
        reDict[self.sample_id_01] = id_1
        reDict[self.sample_id_02] = id_2
        reDict[self.sample_missing] = missing

        return reDict
    
    def get_sample_value_line(self,sample_dict:dict)->str:
        values = list(map(sample_dict.get,self.sample_header))
        return self.sample_split_params.join(values)

legend_config = LegendConfig()