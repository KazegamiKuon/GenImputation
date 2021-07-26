from sklearn.metrics import r2_score
from transformers import EvalPrediction
import torch
from torch.nn import functional as F

def evalpred_to_word(eval_prediction: EvalPrediction):
    label_ids = eval_prediction.label_ids
    logits = eval_prediction.predictions
    softmax = F.softmax(torch.as_tensor(logits),dim=-1)
    top_word = torch.argmax(softmax,dim=-1)
    return label_ids, top_word

def r2_score_transformers(eval_prediction: EvalPrediction)->dict:
    label_ids, top_word = evalpred_to_word(eval_prediction)
    r2_score_vs = r2_score(label_ids.T,top_word.T)
    r2_score_sv = r2_score(label_ids,top_word)
    r2_score_sum = r2_score_vs+r2_score_sv
    return {
        'R2 score VS':r2_score_vs,
        'R2 score SV':r2_score_sv,
        'R2 score sum':r2_score_sum
    }