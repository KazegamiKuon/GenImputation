from sklearn.metrics import r2_score, matthews_corrcoef
from scipy.stats import pearsonr
from transformers import EvalPrediction
import torch
from torch.nn import functional as F
from transformers.file_utils import ExplicitEnum
from transformers.models.electra.tokenization_electra import ElectraTokenizer
import numpy as np
from typing import List

def evalpred_to_word(eval_prediction: EvalPrediction):
    label_ids = eval_prediction.label_ids
    logits = eval_prediction.predictions
    softmax = F.softmax(torch.as_tensor(logits),dim=-1)
    top_words = torch.argmax(softmax,dim=-1)
    return label_ids, top_words

def to_zero_one(data:List[List[int]],tokenizer:ElectraTokenizer)->List[List[int]]:
    one_value = tokenizer('1')
    one_value = one_value['input_ids'][1]
    return list(map(lambda x: list(map(lambda y: 1 if y == one_value else 0,x)),data))

class ScoreTypeEnum(ExplicitEnum):
    EVAL = "eval_"
    TEST = "test_"

class EvalMetrics:
    def __init__(self) -> None:
        self.r2_vs = "R2 score VS"
        self.r2_sv = "R2 score SV"
        self.r2_sum = "R2 score sum"
        self.r2_pred = "R2 score pred"
        self.pearsonr_vs = "Pearsonr VS"
        self.pearsonr_pred = "Pearsonr pred"
        self.__score_keys = [self.r2_vs,self.r2_sv,self.r2_sum,self.r2_pred,self.pearsonr_vs,self.pearsonr_pred]
        pass

    def get_score_key(self,key:str,ktype:ScoreTypeEnum=ScoreTypeEnum.EVAL)->str:
        assert key in self.__score_keys, "key must be in [{}]".format(", ".join(self.__score_keys))
        return ktype.value+key

    def get_score_transformers(self,eval_prediction: EvalPrediction,input_ids:torch.Tensor, tokenizer:ElectraTokenizer,zero_one:bool=False)->dict:
        label_ids, top_words = evalpred_to_word(eval_prediction)
        # convert data to 0 and 1
        if zero_one:
            label_ids = to_zero_one(label_ids,tokenizer)
            label_ids = np.array(label_ids)
            top_words = to_zero_one(top_words,tokenizer)
            top_words = np.array(top_words)
            pass
        data = {}
        data[self.r2_vs] = r2_score(label_ids.T,top_words.T)
        data[self.r2_sv] = r2_score(label_ids,top_words)
        data[self.r2_sum] = data[self.r2_vs] + data[self.r2_sv]
        r2_score_pred_only = []
        maskeds = []
        for i, input_id in enumerate(input_ids):
            label_id, top_word = label_ids[i], top_words[i]
            masked = input_id == tokenizer.mask_token_id
            label_id = label_id[masked]
            top_word = top_word[masked]
            masked = np.array(masked)
            maskeds.append(masked)
            r2_score_pred_only.append(r2_score(label_id,top_word))
        data[self.r2_pred] = np.mean(r2_score_pred_only)
        maskeds = np.array(maskeds)
        maskeds = maskeds.flatten()
        label_ids = label_ids.flatten()
        top_words = top_words.flatten()
        data[self.pearsonr_vs] = pearsonr(label_ids,top_words)[0]**2
        data[self.pearsonr_pred] = pearsonr(label_ids[maskeds],top_words[maskeds])[0]**2
        return data

eval_metrics = EvalMetrics()