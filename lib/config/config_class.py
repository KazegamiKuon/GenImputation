import os
import json
import typing
from marshmallow.fields import Function, Number
import numpy as np
import pandas as pd
from transformers.file_utils import ExplicitEnum
from transformers import ElectraConfig
from ..utils import general as g
from multimethod import multimethod
from dataclasses_json import dataclass_json
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
from ..model.overwriter import OTrainingArguments
from copy import deepcopy, copy


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
        self.alleles_position_0 = ['I', 'D']
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
        self.chrom_values = list(
            np.hstack([np.arange(1, 23).astype(str), np.array(['X', 'Y'])]))
        self.vcf_tails = ['.vcf', 'vcf.gz', 'vcf.tar.gz']


vcf_config = VCFConfig()


class VCFpageNLP:
    def __init__(self) -> None:
        self.page = "page"
        self.variant = "variant"
        self.info = "info"
        self.token = 'token'
        self.masked = 'masked'

        self.file_types = [self.page, self.variant,
                           self.info, self.token, self.masked]

        self.page_split_params = ' '
        self.tail = dict({
            self.page: '.page.gz',
            self.variant: '.variant.gz',
            self.info: '.info.csv',
            self.token: '.token.json.gz',
            self.masked: '.masked.json.gz'
        })

        self.model_args = "model_args"
        self.train_args = "train_args"
        self.file_train_prefix = "file_train_prefix"
        self.file_test_prefix = "file_test_prefix"
        self.vocab_file = "vocab_file"
        self.save_dir = "save_dir"
        pass

    def get_file_path(self, output_prefix: str, file_type: str, region: int = 0, batch: int = 0) -> str:
        assert file_type in self.file_types, 'type must be in this list [{}]'.format(
            ', '.join(self.file_types))
        path = output_prefix+'.r{0:04d}.b{1:04d}'.format(region, batch) + self.tail.get(file_type)
        return os.path.abspath(path)

    def get_file_path_from_page(self, path: str, file_type: str, detail_name: str = '') -> str:
        assert file_type in self.file_types, 'type must be in this list [{}]'.format(
            ', '.join(self.file_types))
        assert path.endswith(self.tail[self.page]), "Must be page type data"
        if file_type == self.variant:
            temps = path.split('.')
            temps[-3] = 'b0000'
            path = '.'.join(temps)
        return path.replace(self.tail[self.page], detail_name + self.tail[file_type])

    def get_file_paths(self, prefix: str, file_type: str, regions: typing.List[int], batchs: typing.List[int]) -> list:
        assert file_type in self.file_types, 'type must be in this list [{}]'.format(
            ', '.join(self.file_types))
        region_file_paths = [[self.get_file_path(
            prefix, file_type, region, batch) for batch in batchs] for region in regions]
        return region_file_paths

    def get_file_paths_in_dir(self, dir: str, file_type: str):
        assert file_type in self.file_types, 'type must be in this list [{}]'.format(
            ', '.join(self.file_types))
        assert os.path.isdir(dir), "{} not found!".format(dir)
        return [os.path.join(dir, file_name) for file_name in os.listdir(dir) if file_name.endswith(page_config.tail[file_type])]

    def token_masked_to_json(self, data: dict, path: str) -> None:
        with g.writing(path) as wjson:
            json.dump(data, wjson)


page_config = VCFpageNLP()


class RegionConfig():
    @multimethod
    def __init__(self, num_inputs: int, fw: list, bw: list, config_json: str) -> None:
        with g.reading(config_json) as default_config:
            self.__dict__ = json.load(default_config)
        self.num_inputs = num_inputs
        self.output_points_bw = bw
        self.output_points_fw = fw
        self.num_outputs = len(bw)

    @multimethod
    def __init__(self, config_json: str) -> None:
        with g.reading(config_json) as config:
            self.__dict__ = json.load(config)

    def to_json(self, json_path):
        with g.writing(json_path) as wjson:
            json.dump(self.__dict__, wjson)


class LegendConfig():
    def __init__(self) -> None:
        self.unobserve = '0'
        self.observe = '1'
        self.hap_tail = '.hap.gz'
        self.gtrue_hap_tail = '_gtrue'+self.hap_tail
        self.legend_tail = '.legend.gz'
        self.gtrue_legend_tail = '_gtrue'+self.legend_tail
        # legend
        self.snp = 'id'
        self.chrom = 'chr'
        self.position = 'position'
        self.ref = 'ref'
        self.alt = 'alt'
        self.af = 'af'
        self.maf = 'maf'
        self.array_marker_flag = 'array_marker_flag'
        self.legend_header = [self.snp, self.chrom, self.position,
                              self.ref, self.alt, self.af, self.maf, self.array_marker_flag]
        self.legend_split_params = ' '
        # hap
        self.hap_split_params = ' '
        # Sample
        self.sample_tail = '.sample'
        self.sample_first_line = '0 0 0'
        self.sample_id_01 = 'ID_1'
        self.sample_id_02 = 'ID_2'
        self.sample_split_params = ' '
        self.sample_missing = 'missing'
        self.sample_header = [self.sample_id_01,
                              self.sample_id_02, self.sample_missing]
        # region
        self.region_folder = 'region'
        # config
        self.config_tail = '.config.json'

    @property
    def legend_header_line(self) -> str:
        return self.legend_split_params.join(self.legend_header)

    @property
    def sample_header_line(self) -> str:
        return self.sample_split_params.join(self.sample_header)

    def get_hap_file_name(self, output_prefix: str) -> str:
        return output_prefix+self.hap_tail

    def get_legend_file_name(self, output_prefix: str) -> str:
        return output_prefix+self.legend_tail

    def get_legend_dict(self, snp, chrom, position, ref, alt, af, flag, nb_alt: int) -> dict:
        reDict = {}
        reDict[self.snp] = snp
        reDict[self.chrom] = chrom
        reDict[self.position] = position
        reDict[self.ref] = ref
        reDict[self.alt] = alt
        reDict[self.af] = af
        reDict[self.maf] = str(1-float(af)) if float(af) > 0.5 else af
        reDict[self.array_marker_flag] = flag
        if reDict[self.snp] == '.':
            reDict[self.snp] = '{:s}:{:s}:{:s}:{:s}'.format(
                chrom, position, ref, alt)
        elif nb_alt >= 2:
            reDict[self.snp] += ':{:s}:{:s}'.format(ref, alt)
        return reDict

    def get_legend_value_line(self, legend_dict: dict) -> str:
        values = list(map(legend_dict.get, self.legend_header))
        return self.legend_split_params.join(values)

    def get_sample_file_name(self, output_prefix: str) -> str:
        return output_prefix+self.sample_tail

    def get_sample_dict(self, id_1: str, id_2: str, missing: str) -> dict:
        reDict = {}
        reDict[self.sample_id_01] = id_1
        reDict[self.sample_id_02] = id_2
        reDict[self.sample_missing] = missing

        return reDict

    def get_sample_value_line(self, sample_dict: dict) -> str:
        values = list(map(sample_dict.get, self.sample_header))
        return self.sample_split_params.join(values)

    def make_region_dir(self, file_path: str) -> None:
        dirname = file_path
        if os.path.isfile(file_path):
            dirname = os.path.dirname(file_path)
        dirname = os.path.join(dirname, self.region_folder)
        os.makedirs(dirname, exist_ok=True)
        return dirname

    def get_legend_region_file_name(self, file_path: str, bin: int, nb_character: int, output_folder=None) -> str:
        dirname, basename = g.get_dir_and_base_name(file_path)
        if output_folder is not None and os.path.isdir(output_folder):
            dirname = output_folder
        region_tail = '_{:0'+str(nb_character)+'d}'+self.legend_tail
        basename = basename.replace(self.legend_tail, region_tail.format(bin))
        region_path = os.path.join(dirname, basename)
        return region_path

    def get_hap_region_file_name(self, file_path: str, bin: int, nb_character: int, output_folder=None) -> str:
        dirname, basename = g.get_dir_and_base_name(file_path)
        if output_folder is not None and os.path.isdir(output_folder):
            dirname = output_folder
        region_tail = '_{:0'+str(nb_character)+'d}'+self.hap_tail
        basename = basename.replace(self.hap_tail, region_tail.format(bin))
        region_path = os.path.join(dirname, basename)
        return region_path

    def get_legend_gtrue_file(self, legend_file: str) -> str:
        return legend_file.replace(self.legend_tail, self.gtrue_legend_tail)

    def get_hap_gtrue_file(self, hap_file: str) -> str:
        return hap_file.replace(self.hap_tail, self.gtrue_hap_tail)

    def legend_dataframe_to_csv(self, path, df: pd.DataFrame) -> None:
        df.to_csv(path, sep=self.legend_split_params, index=False)

    def legend_to_json_path(self, legend_path: str) -> str:
        return legend_path.replace(self.legend_tail, self.config_tail)

    def to_region_config(self, num_inputs: int, fw: list, bw: list, legend_path: str, default_config_path: str):
        json_path = self.legend_to_json_path(legend_path)
        region_config = RegionConfig(num_inputs, fw, bw, default_config_path)
        region_config.to_json(json_path)


legend_config = LegendConfig()


class TrainModeType(ExplicitEnum):
    NOPRETRAIN = "nopretrain"
    PRETRAIN = "pretrain"


class MaskedModeType(ExplicitEnum):
    RANDOM = "random"
    NORMAL = "normal"


class DirType(ExplicitEnum):
    OUTPUT = "output"
    LOGGING = "logging"
    CHECKPOINT = "checkpoint"


@dataclass_json
@dataclass
class TrainElectraArguments:
    regions: List[int] = field(
        metadata={"help": "The region will be trained"}
    )
    batchs: List[int] = field(
        metadata={"help": "The region's batchs will be trained"}
    )
    save_dir_format: str = field(
        metadata={"help": "Save dir format path use to save data."}
    )
    vocab_file: str = field(
        metadata={"help": "vocab file passing to tokenizer."}
    )
    file_train_prefix: str = field(
        metadata={
            "help": "the file train prefix passing to page_config.get_file_paths()"}
    )
    file_test_prefix: str = field(
        metadata={
            "help": "the file train prefix passing to page_config.get_file_paths()"}
    )
    output_dir_format: str = field(
        metadata={
            "help": "Output dir format path use to create output_dir in trainning args."}
    )
    logging_dir_format: str = field(
        metadata={
            "help": "logging dir format path use to create output_dir in logging args."}
    )
    train_args: dict = field(
        metadata={"help": "the trainning arguments"}
    )
    model_args: dict = field(
        metadata={"help": "the electra config"}
    )
    mode: TrainModeType = field(
        default="nopretrain",
        metadata={"help": "the training mode to run"}
    )
    model_name: str = field(
        default="",
        metadata={"help": "prefix will adding to output"}
    )
    masked_mode: MaskedModeType = field(
        default=MaskedModeType.NORMAL,
        metadata={
            "help": "masked mode use to dataset param when it masked by flag is true"}
    )
    seed: int = field(
        default=42,
        metadata={"help": "seed use for train and orther"}
    )
    resume_from_checkpoint_format: str = field(
        default=None,
        metadata={
            "help": "resume checkpoint format path which you will use for get checkpoint path"}
    )
    pretrain_path_format: str = field(
        default=None,
        metadata={
            "help": "pretrain path format path which you will use for pretrain path"}
    )

    def get_vocab_file(self):
        return os.path.abspath(self.vocab_file)
        
    def get_pretrain_path_format(self):
        if self.pretrain_path_format is None:
            return self.save_dir_format
        else:
            return self.pretrain_path_format

    def get_trainning_args(self, func_format: Function):
        trainning_args = deepcopy(self.train_args)
        trainning_args['output_dir'] = func_format(
            DirType.OUTPUT, self.output_dir_format)
        trainning_args['logging_dir'] = func_format(
            DirType.LOGGING, self.logging_dir_format)
        trainning_args['seed'] = self.seed
        return OTrainingArguments(**trainning_args)

    def get_model_args(self):
        model_args = deepcopy(self.model_args)
        return ElectraConfig(**model_args)

    def get_resume_from_checkpoint(self, func_format: Function):
        resume_from_checkpoint = None
        if self.resume_from_checkpoint_format is not None:
            resume_from_checkpoint = func_format(
                DirType.CHECKPOINT, self.resume_from_checkpoint_format)
        return resume_from_checkpoint

    def save_to_json(self, path):
        with g.writing(path) as jf:
            json.dump(self.to_dict(), jf)


class TrainElectraConfig():
    def __init__(self) -> None:
        self.learning_rate = 'learning_rate'
        self.loss = "loss"
        self.epoch = "epoch"
        self.eval_R2_score_SV = "eval_R2 score SV"
        self.eval_R2_score_VS = "eval_R2 score VS"
        self.eval_R2_score_pred = "eval_R2 score pred"
        self.eval_R2_score_sum = "eval_R2 score sum"
        self.eval_loss = "eval_loss"
        self.log_history = "log_history"
        self.history_file = "trainer_state.json"
        pass

    def get_history_data(self, checkpoint_dir: str, data_name: str):
        paths = os.listdir(checkpoint_dir)
        paths.sort()
        latest_path = os.path.join(
            checkpoint_dir, paths[-1], self.history_file)
        data = []
        with g.reading(latest_path) as hf:
            log_historys = json.load(hf)[self.log_history]
            data = [log_history[data_name]
                    for log_history in log_historys if data_name in log_history]
        return data

    def load_from_json(self, path: str) -> TrainElectraArguments:
        data = None
        with g.reading(path) as jf:
            data = json.load(jf)
        assert data is not None, "data from json can't none"
        return TrainElectraArguments.from_dict(data)

    def get_func_format(self, region: int, detail: str) -> Function:
        def func_format(type: DirType, path_format: str) -> str:
            return path_format.format(region)+detail
        return func_format


train_electra_config = TrainElectraConfig()

@dataclass_json
@dataclass
class ILogHistory:
    epoch: float = field(
        metadata={"help": "epoch log at"}
    )
    step: float = field(
        metadata={"help": "step log at"}
    )
    def get_score(self,name:str):
        pass

    def check_key(self,key)->bool:
        pass
    
    def is_this_class(data: dict)->bool:
        "overwrite this"
        pass

@dataclass_json
@dataclass
class TrainLogHistory(ILogHistory):
    learning_rate: float = field(
        metadata={"help": "lr log at"}
    )
    loss: float = field(
        metadata={"help": "loss log at"}
    )
    def check_key(self, key) -> bool:
        return key in self.__dict__
    
    def get_score(self,name:str):
        return self.__getattribute__(name)

    def is_this_class(data: dict) -> bool:
        if "learning_rate" in data:
            return True
        return False

class EvalLogHistory(ILogHistory):
    def __init__(self,data:dict) -> None:
        self.epoch = data['epoch']
        self.step = data["step"]
        self.eval_loss = data["eval_loss"]
        self.eval_runtime = data["eval_runtime"]
        self.eval_samples_per_second = data["eval_samples_per_second"]
        self.eval_steps_per_second = data["eval_steps_per_second"]
        # overwrite data from IEvalMetrics
        self.__eval_name = 'eval_{}'
        self.__other_data = data
    
    def __check_data(self,key)->bool:
        if key in self.__other_data:
            return True
        return False
    
    def get_keys(self)->List:
        return self.__other_data.keys()

    def check_key(self, key) -> bool:
        return key in self.__other_data
    
    def get_score(self,name:str):
        if self.__check_data(name):
            return self.__other_data[name]
        return None

    def is_this_class(data: dict) -> bool:
        if "eval_runtime" in data:
            return True
        return False

@dataclass_json
@dataclass
class CheckpointLogConfig:
    best_metric: float = field(
        metadata={"help": "best metric"}
    )
    best_model_checkpoint: str = field(
        metadata={"help": "best model checkpoint"}
    )
    epoch: float = field(
        metadata={"help": "epoch"}
    )
    global_step: float = field(
        metadata={"help": "global step"}
    )
    is_hyper_param_search: bool = field(
        metadata={"help": "is hyper param search"}
    )
    is_local_process_zero: bool = field(
        metadata={"help": "is local process zero"}
    )
    is_world_process_zero: bool = field(
        metadata={"help": "is world process zero"}
    )
    log_history: List[dict] = field(
        metadata={"help": "log history"}
    )
    def get_train_log(self):
        data = [TrainLogHistory.from_dict(logh) for logh in self.log_history if TrainLogHistory.is_this_class(logh)]
        return data
    def get_eval_log(self):
        data = [EvalLogHistory(logh) for logh in self.log_history if EvalLogHistory.is_this_class(logh)]
        return data

class CheckpointConfig:
    def __init__(self) -> None:
        self.log_data = 'trainer_state.json'
        pass
    
    def get_checkpoint_name(self,step:Number)->str:
        return 'checkpoint-{}'.format(step)

    def get_last_checkpoint(self,checkpoint_dir:str)->str:
        assert os.path.isdir(checkpoint_dir), "{} doesnt exist".format(checkpoint_dir)
        checkpoints = os.listdir(checkpoint_dir)
        checkpoints.sort(key=lambda x: int(x.split('-')[-1]))
        last_checpoint = os.path.join(checkpoint_dir,checkpoints[-1])
        return last_checpoint
    
    def get_checkpoint_data(self,checkpoint_dir:str,is_checpoint_dir:bool=True)->CheckpointLogConfig:
        last_checpoint = self.get_last_checkpoint(checkpoint_dir)
        if not is_checpoint_dir:
            last_checpoint = checkpoint_dir
        data_path = os.path.join(last_checpoint,self.log_data)
        with g.reading(data_path) as cf:
            data = json.load(cf)
            return CheckpointLogConfig.from_dict(data)
    
    def get_best_model_from_checpoint(self,checkpoint_dir:str,key:str)->Tuple[str,EvalLogHistory]:
        assert os.path.isdir(checkpoint_dir), "{} doesnt exist".format(checkpoint_dir)
        data_log = self.get_checkpoint_data(checkpoint_dir)
        eval_logs = data_log.get_eval_log()
        scores = [eval_log.get_score(key) for eval_log in eval_logs]
        scores = list(map(lambda x: -np.inf if np.isnan(x) else x,scores))
        index = np.argmax(scores)
        best_eval_log = eval_logs[index]
        checkpoint_name = self.get_checkpoint_name(best_eval_log.step)
        return os.path.join(checkpoint_dir,checkpoint_name), best_eval_log
    
    def get_best_score_each_cycle(self,checpoint_dir:str,key:str,nb_cycle:int)->List[EvalLogHistory]:
        data_log = self.get_checkpoint_data(checpoint_dir)
        eval_logs = data_log.get_eval_log()
        data = [{'epoch': eval_log.epoch,'score':eval_log.get_score(key),'index': i} for i, eval_log in enumerate(eval_logs)]
        data = pd.DataFrame(data)
        nbdata_cycle = max(data['epoch'])/nb_cycle
        data['cycle'] = data['epoch']//nbdata_cycle
        idx = data.groupby(['cycle'])['score'].transform(max) == data['score']
        maxdata = data[idx]
        return list(np.array(eval_logs)[maxdata['index']])
        pass

cpconfig = CheckpointConfig()