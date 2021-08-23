import sys
sys.path.insert(0,'/client/user1/cuongdev/GenImputation/')
from lib.config.config_class import cpconfig
from transformers import ElectraForMaskedLM, ElectraTokenizer, EvalPrediction
import os
from lib.config.config_class import train_electra_config as teconfig
from lib.config.config_class import TrainModeType, page_config
from lib.data_processing import GenNLPMaskedDataset
from transformers import ElectraForMaskedLM, ElectraTokenizer
from lib.utils import general as g
from lib.model.overwriter import OTrainer
from lib.utils.metrics import eval_metrics
import json
from IPython.display import clear_output

config_path = '/client/user1/cuongdev/GenImputation/data/train/electra_G1K_22_hs37d5/config_pretrain.json'
config = teconfig.load_from_json(config_path)
regions = config.regions
batchs = config.batchs
detail = config.model_name
mode = config.mode

masked_mode = config.masked_mode

test_region_paths = page_config.get_file_paths(config.file_test_prefix,page_config.page,regions,[0])
vocab_file = config.vocab_file
save_dir_format = config.save_dir_format+detail
pretrain_path_format = config.pretrain_path_format if config.pretrain_path_format is not None else save_dir_format
tokenizer = ElectraTokenizer(vocab_file=vocab_file)
seed = config.seed

key = eval_metrics.get_score_key(eval_metrics.pearsonr_pred)
scores = []
for i, region in enumerate(regions):
    if region != 1:
        continue
    clear_output(wait=True)
    print('Region {} viewing...'.format(region))
    # train arg and model arg
    func_format = teconfig.get_func_format(region,detail)
    training_args = config.get_trainning_args(func_format)
    checkpoint_dir = training_args.output_dir
    ## test data
    test_batch_paths = test_region_paths[i]
    test_dataset = GenNLPMaskedDataset(test_batch_paths,tokenizer,seed=seed,masked_by_flag=True,only_input=True,force_create=True)
    ## model
    checkpoint_path, eval_log = cpconfig.get_best_model_from_checpoint(checkpoint_dir,key = key)
    electra_model = ElectraForMaskedLM.from_pretrained(checkpoint_path)
    input_ids=test_dataset.maskeds['input_ids']
    output_model = electra_model(input_ids=input_ids)
    output_model = EvalPrediction(predictions=output_model.logits,label_ids=test_dataset.labels)
    score = eval_metrics.get_score_transformers(output_model,input_ids,tokenizer,True)
    scores.append(score)
    pass