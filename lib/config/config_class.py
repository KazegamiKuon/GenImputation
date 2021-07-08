import os
import json
import numpy as np
import pandas as pd
from ..utils import general as g
from multimethod import multimethod

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
        self.vcf_tails = ['.vcf','vcf.gz','vcf.tar.gz']

vcf_config = VCFConfig()

class VCFpageNLP:
    def __init__(self) -> None:
        self.page = "page"
        self.variant = "variant"
        self.info = "info"
        self.token = 'token'
        self.masked = 'masked'

        self.__file_tyeps = [self.page,self.variant,self.info,self.token,self.masked]

        self.page_split_params = ' '
        self.tail = dict({
            self.page: '.page.gz',
            self.variant: '.variant.gz',
            self.info: '.info',
            self.token: '.token.json.gz',
            self.masked:'.masked.json.gz'
        })
        pass
    
    def get_file_path(self,output_prefix:str,file_type:str)->str:
        assert file_type in self.__file_tyeps, 'type must be in this list [{}]'.format(', '.join(self.__file_tyeps))
        return output_prefix+self.tail.get(file_type)
    
    def get_file_path_from_page(self,path:str,file_type:str,detail_name:str='')->str:
        assert file_type in self.__file_tyeps, 'type must be in this list [{}]'.format(', '.join(self.__file_tyeps))
        assert path.endswith(self.tail[self.page]), "Must be page type data"
        return path.replace(self.tail[self.page],detail_name+ self.tail[file_type])
    
    def token_masked_to_json(self,data:dict,path:str)->None:
        with g.writing(path) as wjson:
            json.dump(data,wjson)


page_config = VCFpageNLP()    

class RegionConfig():
    @multimethod
    def __init__(self,num_inputs:int,fw:list,bw:list,config_json:str) -> None:
        with g.reading(config_json) as default_config:
            self.__dict__ = json.load(default_config)
        self.num_inputs = num_inputs
        self.output_points_bw = bw
        self.output_points_fw = fw
        self.num_outputs = len(bw)
    @multimethod
    def __init__(self,config_json:str) -> None:
        with g.reading(config_json) as config:
            self.__dict__ = json.load(config)
    
    def to_json(self,json_path):
        with g.writing(json_path) as wjson:
            json.dump(self.__dict__,wjson)

class LegendConfig():
    def __init__(self) -> None:
        self.unobserve = '0'
        self.observe = '1'
        self.hap_tail = '.hap.gz'
        self.gtrue_hap_tail = '_gtrue'+self.hap_tail
        self.legend_tail = '.legend.gz'
        self.gtrue_legend_tail = '_gtrue'+self.legend_tail
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
        # region
        self.region_folder ='region'
        # config
        self.config_tail = '.config.json'
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
        reDict[self.maf] = str(1-float(af)) if float(af) > 0.5 else af
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

    def make_region_dir(self,file_path:str)->None:
        dirname = file_path
        if os.path.isfile(file_path):
            dirname = os.path.dirname(file_path)
        dirname = os.path.join(dirname,self.region_folder)
        os.makedirs(dirname,exist_ok=True)
        return dirname

    def get_legend_region_file_name(self,file_path:str,bin:int,nb_character:int,output_folder = None)->str:
        dirname, basename = g.get_dir_and_base_name(file_path)
        if output_folder is not None and os.path.isdir(output_folder):
            dirname = output_folder
        region_tail = '_{:0'+str(nb_character)+'d}'+self.legend_tail
        basename = basename.replace(self.legend_tail,region_tail.format(bin))
        region_path = os.path.join(dirname,basename)
        return region_path

    def get_hap_region_file_name(self,file_path:str,bin:int,nb_character:int,output_folder = None)->str:
        dirname, basename = g.get_dir_and_base_name(file_path)
        if output_folder is not None and os.path.isdir(output_folder):
            dirname = output_folder
        region_tail = '_{:0'+str(nb_character)+'d}'+self.hap_tail
        basename = basename.replace(self.hap_tail,region_tail.format(bin))
        region_path = os.path.join(dirname,basename)
        return region_path
    
    def get_legend_gtrue_file(self,legend_file:str)->str:
        return legend_file.replace(self.legend_tail,self.gtrue_legend_tail)
    
    def get_hap_gtrue_file(self,hap_file:str)->str:
        return hap_file.replace(self.hap_tail,self.gtrue_hap_tail)
    
    def legend_dataframe_to_csv(self,path,df:pd.DataFrame)->None:
        df.to_csv(path,sep=self.legend_split_params,index=False)
    
    def legend_to_json_path(self,legend_path:str)->str:
        return legend_path.replace(self.legend_tail,self.config_tail)

    def to_region_config(self,num_inputs:int,fw:list,bw:list,legend_path:str,default_config_path:str):
        json_path = self.legend_to_json_path(legend_path)
        region_config = RegionConfig(num_inputs,fw,bw,default_config_path)
        region_config.to_json(json_path)

legend_config = LegendConfig()