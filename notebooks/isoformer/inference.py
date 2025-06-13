from datasets import load_dataset
from transformers import AutoTokenizer, AutoModelForMaskedLM
import numpy as np
import torch

# Load the dataset
transcript_expression_dataset = load_dataset(
    "InstaDeepAI/multi_omics_transcript_expression",
    task_name="transcript_expression_expression",
    sequence_length=12288,
    filter_out_sequence_length=196608,
    revision="un-standardized-values",
    split="test"
)

# Import the tokenizer and the model
tokenizer = AutoTokenizer.from_pretrained("InstaDeepAI/isoformer", trust_remote_code=True)
model = AutoModelForMaskedLM.from_pretrained("InstaDeepAI/isoformer",trust_remote_code=True)

# Sample data
sample_data = transcript_expression_dataset[0]
protein_sequences = [sample_data["Protein"]]
rna_sequences = [sample_data["RNA"]]
dna_sequences = [sample_data["DNA"]]
sequence_length = 196_608
rng = np.random.default_rng(seed=0)

torch_tokens = tokenizer(
    dna_input=dna_sequences, rna_input=rna_sequences, protein_input=protein_sequences
)
dna_torch_tokens = torch.tensor(torch_tokens[0]["input_ids"])
rna_torch_tokens = torch.tensor(torch_tokens[1]["input_ids"])
protein_torch_tokens = torch.tensor(torch_tokens[2]["input_ids"])

torch_output = model.forward(
    tensor_dna=dna_torch_tokens,
    tensor_rna=rna_torch_tokens,
    tensor_protein=protein_torch_tokens,
    attention_mask_rna=rna_torch_tokens != 1,
    attention_mask_protein=protein_torch_tokens != 1,
)

print(f"Gene expression predictions: {torch_output['gene_expression_predictions']}")
print(f"Final DNA embedding: {torch_output['final_dna_embeddings']}")



