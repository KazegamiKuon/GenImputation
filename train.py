import sys
import typing
sys.path.insert(0,'./')
from lib.data_processing import GenNLPMaskedDataset
from transformers import ElectraForMaskedLM, ElectraTokenizer
from lib.utils import general as g
from lib.config.config_class import TrainModeType, page_config
from lib.config.config_class import train_electra_config as teconfig
from lib.model.overwriter import OTrainer
from lib.utils.metrics import eval_metrics
import json
import os
from IPython.display import clear_output

config_path = '/client/user1/cuongdev/GenImputation/data/train/electra_G1K_22_hs37d5/config_pretrain_v2.json'

config = teconfig.load_from_json(config_path)

regions = config.regions
batchs = config.batchs
detail = config.model_name
mode = config.mode

masked_mode = config.masked_mode

train_region_paths = page_config.get_file_paths(config.file_train_prefix,page_config.page,regions,batchs)
test_region_paths = page_config.get_file_paths(config.file_test_prefix,page_config.page,regions,[0])
vocab_file = config.vocab_file
save_dir_format = config.save_dir_format+detail
pretrain_path_format = config.pretrain_path_format if config.pretrain_path_format is not None else save_dir_format
tokenizer = ElectraTokenizer(vocab_file=vocab_file)
seed = config.seed

for i, region in enumerate(regions):
    clear_output(wait=True)
    print('Region {} trainning...'.format(region))
    print('Prevert region {} trainning...'.format(region-1))
    #get path and params
    save_path = save_dir_format.format(region)
    prevert_path = pretrain_path_format.format(region-1)
    func_format = teconfig.get_func_format(region,'')
    resume_from_checkpoint_path = config.get_resume_from_checkpoint(func_format)
    # train arg and model arg
    func_format = teconfig.get_func_format(region,detail)
    training_args = config.get_trainning_args(func_format)
    model_config = config.get_model_args()
    ## Train and eval data
    train_batch_paths = train_region_paths[i]
    train_dataset = GenNLPMaskedDataset(train_batch_paths, tokenizer, seed=seed, masked_by_flag=True, only_input=True,force_create=True,masked_mode = masked_mode)
    # eval_dataset = GenNLPMaskedDataset(train_batch_paths[-1:], tokenizer, seed=seed, masked_by_flag=True, only_input=True)
    ## test data
    test_batch_paths = test_region_paths[i]
    test_dataset = GenNLPMaskedDataset(test_batch_paths,tokenizer,seed=seed,masked_by_flag=True,only_input=True)
    ## model
    model_config.vocab_size = tokenizer.vocab_size
    if mode == TrainModeType.PRETRAIN:
        model_config.max_position_embeddings = 1300
    elif mode == TrainModeType.NOPRETRAIN:
        model_config.max_position_embeddings = train_dataset.max_position_embeddings()
    else:
        model_config.max_position_embeddings = train_dataset.max_position_embeddings()
    electra_model = ElectraForMaskedLM(model_config)
    if os.path.isdir(prevert_path):
        electra_model = ElectraForMaskedLM.from_pretrained(prevert_path,config=model_config)
    trainer = OTrainer(
        model = electra_model,
        args=training_args,
        train_dataset = train_dataset,
        eval_dataset = test_dataset,
        compute_metrics = lambda eval_prediction: eval_metrics.get_score_transformers(eval_prediction,test_dataset.maskeds['input_ids'],tokenizer),
    )
    trainer.train(resume_from_checkpoint = resume_from_checkpoint_path)
    trainer.save_model(save_path)
    output_test = trainer.predict(test_dataset)
    metrics = output_test.metrics
    test_result_path = os.path.join(save_path,'test_result.json')
    with g.writing(test_result_path) as trf:
        json.dump(metrics,trf)