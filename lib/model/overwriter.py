from dataclasses import dataclass, field
from typing import Optional
from transformers import Trainer, TrainingArguments, get_cosine_with_hard_restarts_schedule_with_warmup, get_scheduler

@dataclass
class OTrainingArguments(TrainingArguments):
    def __init_subclass__(cls) -> None:
        return super().__init_subclass__()
    
    def __post_init__(self):
        return super().__post_init__()
    
    num_cycles: Optional[int] = field(
        default=1,
        metadata={"help": "num cycles at hard restarts schedule"},
    )
    
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
                num_cycles = self.args.num_cycles
            )