# from torch.utils.data.dataloader import DataLoader
# from torch.utils.data.dataset import IterableDataset
# import numpy as np
# import torch
# from transformers.trainer_pt_utils import (
#     IterableDatasetShard,
#     find_batch_size,
#     nested_concat,
#     nested_numpify,
#     nested_truncate,
# )
# from transformers.file_utils import (
#     is_torch_tpu_available,
# )
# from transformers.trainer_utils import (
#     EvalLoopOutput,
#     EvalPrediction,
#     denumpify_detensorize,
# )
# import collections
# from transformers.utils import logging
# logger = logging.get_logger(__name__)
# from transformers.deepspeed import deepspeed_init
import os
from typing import List, Optional
from dataclasses import dataclass, field
from transformers import Trainer, TrainingArguments, get_cosine_with_hard_restarts_schedule_with_warmup

@dataclass
class OTrainingArguments(TrainingArguments):
    num_cycles: Optional[int] = field(
        default=1,
        metadata={"help": "num cycles at hard restarts schedule"},
    )


# if is_torch_tpu_available():
#     import torch_xla.core.xla_model as xm
#     import torch_xla.debug.metrics as met
#     import torch_xla.distributed.parallel_loader as pl

class OTrainer(Trainer):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def create_scheduler(self, num_training_steps: int):
        """
        Setup the scheduler. The optimizer of the trainer must have been set up before this method is called.

        Args:
            num_training_steps (int): The number of training steps to do.
        """
        if self.lr_scheduler is None:
            self.lr_scheduler = get_cosine_with_hard_restarts_schedule_with_warmup(
                self.optimizer,
                num_warmup_steps=self.args.warmup_steps,
                num_training_steps=num_training_steps,
                num_cycles=self.args.num_cycles
            )
