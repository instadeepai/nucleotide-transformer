# Nucleotide Transformers and SegmentNT

[![Python Version](https://img.shields.io/badge/python-3.8-blue.svg)](https://docs.python.org/3.8/library/index.html)
[![Jax Version](https://img.shields.io/badge/jax-0.3.25-informational)](https://jax.readthedocs.io/en/latest/)
[![license](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-blue.svg)](LICENSE)

Welcome to this InstaDeep Github repository, where are featured:
1. A collection of transformer based genomic language models from both of our research works, [The Nucleotide
Transformer ](https://www.biorxiv.org/content/10.1101/2023.01.11.523679v3) and [Agro Nucleotide Transformer](https://www.biorxiv.org/content/10.1101/2023.10.24.563624v2).
2. A collection of segmentation models using the Nucleotide Transformers as a backbone, allowing segmentation of a DNA sequence's
genomic elements at single-nucleotide resolution: the [SegmentNT](https://www.biorxiv.org/content/10.1101/2024.03.14.584712v1.full.pdf) models.
3. Similarly to the SegmentNT models, SegmentEnformer and SegmentBorzoi, allowing segmentation of a DNA sequence's genomic elements at single-nucleotide resolution, using respectively [Enformer](https://www.nature.com/articles/s41592-021-01252-x) and [Borzoi](https://www.nature.com/articles/s41588-024-02053-6)

We are thrilled to open-source these works and provide the community with access to the
code and pre-trained weights for these nine genomics language models and 2 segmentation models. Models from [The Nucleotide Transformer
](https://www.biorxiv.org/content/10.1101/2023.01.11.523679v3) project were
developed in collaboration with Nvidia and TUM, and the models were trained on DGX
A100 nodes on Cambridge-1. The model from the [Agro
Nucleotide Transformer](https://www.biorxiv.org/content/10.1101/2023.10.24.563624v1)
project was develop in collaboration with Google, and the model trained on TPU-v4
accelerators.

Overall, our works provides novel insights related to the pretraining and application
of language foundational models, as well as the training of models using them as
a backbone encoder, to genomics with ample opportunities of their applications in the field.

In this repository, you will find the following:

- Inference code for our models
- Pre-trained weights for all 9 NT models and 2 SegmentNT models
- Instructions for using the code and pre-trained models

---
## The Nucleotide Transformer Models
Compared to other approaches, our models do not only integrate information from single reference genomes,
but leverage DNA sequences from over 3,200 diverse human genomes, as well as 850 genomes from a wide range of species,
including model and non-model organisms. Through robust and extensive evaluation,
we show that these large models provide extremely accurate molecular phenotype prediction compared to existing methods.

<img src="imgs/nt_results_rebuttal_2.png" alt= "Performance on downstream tasks" width="800" height="800">

*Fig. 1: The Nucleotide Transformer model accurately predicts diverse genomics tasks
after fine-tuning. We show the performance results across downstream tasks for fine-tuned transformer models. Error bars represent 2 SDs derived from 10-fold cross-validation.*


## Agro Nucleotide Transformer Model
In this work we present a novel foundational large language model trained
on reference genomes from 48 plant species with a predominant focus on crop
species. We assessed the performance of AgroNT across several prediction tasks
ranging from regulatory features, RNA processing, and gene expression, and show that
AgroNT can obtain state-of-the art performance.

<img src="imgs/Agro_NT_Gene_Expression.png" alt="AgroNT Performance on Gene Expression">

*Fig. 2: AgroNT provides gene expression prediction across different plant species.
Gene expression prediction on holdout genes across all tissues are correlated with
observed gene expression levels. The coefficient of determination (R<sup>2</sup>) from a linear model
and associated P -values between predicted and observed values are shown.*


#### Get started 🚀

To use the code and pre-trained models, simply:

1. Clone the repository to your local machine.
2. Install the package by running `pip install .`.

You can then download and do the inference with any of our nine models in only a few
lines of codes:
```python
import haiku as hk
import jax
import jax.numpy as jnp
from nucleotide_transformer.pretrained import get_pretrained_model

# Get pretrained model
parameters, forward_fn, tokenizer, config = get_pretrained_model(
    model_name="500M_human_ref",
    embeddings_layers_to_save=(20,),
    max_positions=32,
)
forward_fn = hk.transform(forward_fn)

# Get data and tokenize it
sequences = ["ATTCCGATTCCGATTCCG", "ATTTCTCTCTCTCTCTGAGATCGATCGATCGAT"]
tokens_ids = [b[1] for b in tokenizer.batch_tokenize(sequences)]
tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)

# Initialize random key
random_key = jax.random.PRNGKey(0)

# Infer
outs = forward_fn.apply(parameters, random_key, tokens)

# Get embeddings at layer 20
print(outs["embeddings_20"].shape)
```
Supported model names are:
- **500M_human_ref**
- **500M_1000G**
- **2B5_1000G**
- **2B5_multi_species**
- **50M_multi_species_v2**
- **100M_multi_species_v2**
- **250M_multi_species_v2**
- **500M_multi_species_v2**
- **1B_agro_nt**

You can also run our models and find more example code in google colab [![Open All Collab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/instadeepai/nucleotide-transformer/blob/main/examples/inference.ipynb)

The code runs both on GPU and TPU thanks to Jax!

#### Nucleotide Transformers v2 models
Our second version Nucleotide Transformer v2 models include a series of architectural changes that proved more efficient: instead of using learned positional embeddings, we use Rotary Embeddings that are used at each attention layer and Gated Linear Units with swish activations without bias. These improved models also accept sequences of up to 2,048 tokens leading to a longer context window of 12kbp.
Inspired by Chinchilla scaling laws, we also trained our NT-v2 models on our multi-species dataset for longer duration (300B tokens for the 50M and 100M models; 1T tokens for the 250M and 500M model) compared to the v1 models (300B tokens for all four models).


#### Embeddings retrieval
The transformer layers are 1-indexed, which means that calling `get_pretrained_model` with the arguments `model_name="500M_human_ref"` and `embeddings_layers_to_save=(1, 20,)` will result in extracting embeddings after the first and 20-th transformer layer. For transformers using the Roberta LM head, it is common practice to extract the final embeddings after the first layer norm of the LM head rather than after the last transformer block. Therefore, if `get_pretrained_model` is called with the following arguments `embeddings_layers_to_save=(24,)`, the embeddings will not be extracted after the final transformer layer but rather after the first layer norm of the LM head.

---

## The SegmentNT Models

SegmentNT models leverage a Nucleotide Transformer (NT) transformer from which we removed the language model head and replaced by a 1-dimensional U-Net segmentation head to predict the location of several types of genomics elements in a sequence at a single nucleotide resolution. We present two different model variants on 14 different classes of human genomics elements in input sequences up to 30kb. These include gene (protein-coding genes, lncRNAs, 5’UTR, 3’UTR, exon, intron, splice acceptor and donor sites) and regulatory (polyA signal, tissue-invariant and tissue-specific promoters and enhancers, and CTCF-
bound sites) elements. SegmentNT achieves superior performance over the state-of-the-art U-Net segmentation architecture, benefiting from the pre-trained weights of NT, and demonstrates zero-shot generalization up to 50kbp.

<img src="imgs/segment_nt_panel1_screen.png" alt= "Performance on downstream tasks" >

*Fig. 1: SegmentNT localizes genomics elements at nucleotide resolution.*

#### Get started 🚀

To use the code and pre-trained models, simply:

1. Clone the repository to your local machine.
2. Install the package by running `pip install .`.

You can then download and infer on a sequence with any of our models in only a few
lines of codes:

⚠️ The SegmentNT models have been trained on a sequences of 30,000 nucleotides, or 5001 tokens (accounting for the CLS token). However, SegmentNT has been shown to generalize up to sequences of 50,000 bp. For training on 30,000 bps, which is a length
superior than the maximum length of 2048 6-mers tokens that the nucleotide transformer can handle, Yarn rescaling is employed.
By default, the `rescaling factor` is set to the one used during the training. In case you need to infer on sequences between 30kbp and 50kbp, make sure to pass the `rescaling_factor` argument in the `get_pretrained_segment_nt_model` function with
the value `rescaling_factor = max_num_nucleotides / max_num_tokens_nt` where `num_dna_tokens_inference` is the number of tokens at inference (i.e 6669 for a sequence of 40008 base pairs) and `max_num_tokens_nt` is the max number of tokens on which the backbone nucleotide-transformer was trained on, i.e `2048`.

🔍 The notebook `examples/inference_segment_nt.ipynb` showcases how to infer on a 50kb sequence and plot the probabilities to reproduce the Fig.3 of the paper.

🚧 The SegmentNT models do not handle any "N" in the input sequence because each nucleotides need to be tokenized as 6-mers, which can not be the case when using sequences containing one or multiple "N" base pairs.

```python
import haiku as hk
import jax
import jax.numpy as jnp
from nucleotide_transformer.pretrained import get_pretrained_segment_nt_model

# Initialize CPU as default JAX device. This makes the code robust to memory leakage on
# the devices.
jax.config.update("jax_platform_name", "cpu")

backend = "cpu"
devices = jax.devices(backend)
num_devices = len(devices)
print(f"Devices found: {devices}")

# The number of DNA tokens (excluding the CLS token prepended) needs to be dividible by
# 2 to the power of the number of downsampling block, i.e 4.
max_num_nucleotides = 8

assert max_num_nucleotides % 4 == 0, (
    "The number of DNA tokens (excluding the CLS token prepended) needs to be dividible by"
     "2 to the power of the number of downsampling block, i.e 4.")

parameters, forward_fn, tokenizer, config = get_pretrained_segment_nt_model(
    model_name="segment_nt",
    embeddings_layers_to_save=(29,),
    attention_maps_to_save=((1, 4), (7, 10)),
    max_positions=max_num_nucleotides + 1,
)
forward_fn = hk.transform(forward_fn)
apply_fn = jax.pmap(forward_fn.apply, devices=devices, donate_argnums=(0,))


# Get data and tokenize it
sequences = ["ATTCCGATTCCGATTCCAACGGATTATTCCGATTAACCGATTCCAATT", "ATTTCTCTCTCTCTCTGAGATCGATGATTTCTCTCTCATCGAACTATG"]
tokens_ids = [b[1] for b in tokenizer.batch_tokenize(sequences)]
tokens = jnp.asarray(tokens_ids, dtype=jnp.int32)

random_key = jax.random.PRNGKey(seed=0)
keys = jax.device_put_replicated(random_key, devices=devices)
parameters = jax.device_put_replicated(parameters, devices=devices)
tokens = jax.device_put_replicated(tokens, devices=devices)

# Infer on the sequence
outs = apply_fn(parameters, keys, tokens)
# Obtain the logits over the genomic features
logits = outs["logits"]
# Transform them in probabilities
probabilities = jnp.asarray(jax.nn.softmax(logits, axis=-1))[...,-1]
print(f"Probabilities shape: {probabilities.shape}")

print(f"Features inferred: {config.features}")

# Get probabilities associated with intron
idx_intron = config.features.index("intron")
probabilities_intron = probabilities[..., idx_intron]
print(f"Intron probabilities shape: {probabilities_intron.shape}")
```

Supported model names are:
- **segment_nt**
- **segment_nt_multi_species**

The code runs both on GPU and TPU thanks to Jax!

---
## Tokenization :abc:

The models are trained on sequences of length up to 1000 tokens, including the
\<CLS> token prepended automatically to the beginning of the sequence. The tokenizer
starts tokenizing from left to right by grouping the letters "A", "C", "G" and "T" in
6-mers. The "N" letter is chosen not to be grouped inside the k-mers, therefore
whenever the tokenizer encounters a "N", or if the number of nucleotides in the sequence
is not a multiple of 6, it will tokenize the nucleotides without grouping them. Examples
are given below:

```python
dna_sequence_1 = "ACGTGTACGTGCACGGACGACTAGTCAGCA"
tokenized_dna_sequence_1 = [<CLS>,<ACGTGT>,<ACGTGC>,<ACGGAC>,<GACTAG>,<TCAGCA>]

dna_sequence_2 = "ACGTGTACNTGCACGGANCGACTAGTCTGA"
tokenized_dna_sequence_2 = [<CLS>,<ACGTGT>,<A>,<C>,<N>,<TGCACG>,<G>,<A>,<N>,<CGACTA>,<GTCTGA>]
```

All the v1 and v2 transformers can therefore take sequences of up to 5994 and 12282 nucleotides respectively if there are
no "N" inside.

---

## SegmentEnformer

SegmentEnformer leverages [Enformer](https://www.nature.com/articles/s41592-021-01252-x) by removing the prediction head and replacing it by a 1-dimensional U-Net segmentation head to predict the location of several types of genomics elements in a sequence at a single nucleotide resolution.


#### Get started 🚀

To use the code and pre-trained models, simply:

1. Clone the repository to your local machine.
2. Install the package by running `pip install .`.

You can then download and infer on a sequence with any of our models in only a few
lines of codes:

🔍 The notebook `examples/inference_segment_enformer.ipynb` showcases how to infer on a 196608bp sequence and plot the probabilities.

```python
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np

from enformer.pretrained import get_pretrained_segment_enformer_model
from enformer.features import FEATURES

# Initialize CPU as default JAX device. This makes the code robust to memory leakage on
# the devices.
jax.config.update("jax_platform_name", "cpu")

backend = "cpu"
devices = jax.devices(backend)
num_devices = len(devices)

# Load model
parameters, state, forward_fn, tokenizer, config = get_pretrained_segment_enformer_model()
forward_fn = hk.transform_with_state(forward_fn)

apply_fn = jax.pmap(forward_fn.apply, devices=devices, donate_argnums=(0,))
random_key = jax.random.PRNGKey(seed=0)

# Replicate over devices
keys = jax.device_put_replicated(random_key, devices=devices)
parameters = jax.device_put_replicated(parameters, devices=devices)
state = jax.device_put_replicated(state, devices=devices)

# Get data and tokenize it
sequences = ["A" * 196_608]
tokens_ids = [b[1] for b in tokenizer.batch_tokenize(sequences)]
tokens = jnp.stack([jnp.asarray(tokens_ids, dtype=jnp.int32)] * num_devices, axis=0)

# Infer
outs, state = apply_fn(parameters, state, keys, tokens)

# Obtain the logits over the genomic features
logits = outs["logits"]
# Transform them on probabilities
probabilities = np.asarray(jax.nn.softmax(logits, axis=-1))[..., -1]

# Get probabilities associated with intron
idx_intron = FEATURES.index("intron")
probabilities_intron = probabilities[..., idx_intron]
print(f"Intron probabilities shape: {probabilities_intron.shape}")
```

## SegmentBorzoi
SegmentBorzoi leverages [Borzoi](https://www.nature.com/articles/s41588-024-02053-6) by removing the prediction head and replacing it by a 1-dimensional U-Net segmentation head to predict the location of several types of genomics elements in a sequence.

#### Get started 🚀

To use the code and pre-trained models, simply:

1. Clone the repository to your local machine.
2. Install the package by running `pip install .`.

You can then download and infer on a sequence with any of our models in only a few
lines of codes:

🔍 The notebook `examples/inference_segment_borzoi.ipynb` showcases how to infer on a 196608bp sequence and plot the probabilities.

```python
import haiku as hk
import jax
import jax.numpy as jnp
import numpy as np

from borzoi.pretrained import get_pretrained_segment_borzoi_model
from enformer.features import FEATURES

# Initialize CPU as default JAX device. This makes the code robust to memory leakage on
# the devices.
jax.config.update("jax_platform_name", "cpu")

backend = "cpu"
devices = jax.devices(backend)
num_devices = len(devices)

# Load model
parameters, state, forward_fn, tokenizer, config = get_pretrained_segment_borzoi_model()
forward_fn = hk.transform_with_state(forward_fn)
apply_fn = jax.pmap(forward_fn.apply, devices=devices, donate_argnums=(0,))
random_key = jax.random.PRNGKey(seed=0)

# Replicate over devices
keys = jax.device_put_replicated(random_key, devices=devices)
parameters = jax.device_put_replicated(parameters, devices=devices)
state = jax.device_put_replicated(state, devices=devices)

# Get data and tokenize it
sequences = ["A" * 524_288]
tokens_ids = [b[1] for b in tokenizer.batch_tokenize(sequences)]
tokens = jnp.stack([jnp.asarray(tokens_ids, dtype=jnp.int32)] * num_devices, axis=0)

# Infer
outs, state = apply_fn(parameters, state, keys, tokens)

# Obtain the logits over the genomic features
logits = outs["logits"]
# Transform them on probabilities
probabilities = np.asarray(jax.nn.softmax(logits, axis=-1))[..., -1]

# Get probabilities associated with intron
idx_intron = FEATURES.index("intron")
probabilities_intron = probabilities[..., idx_intron]
print(f"Intron probabilities shape: {probabilities_intron.shape}")

```

---

## HuggingFace 🤗

The collection of models presented in this repository are available on Instadeep's
huggingface spaces here: [The Nucleotide Transformers space](https://huggingface.co/collections/InstaDeepAI/nucleotide-transformer-65099cdde13ff96230f2e592)
and [Agro Nucleotide Transformer space](https://huggingface.co/collections/InstaDeepAI/agro-nucleotide-transformer-65b25c077cd0069ad6f6d344)!

- **Nucleotide Transformer**: Two
example notebooks showing how to finetune any of the models [with regular finetuning](https://github.com/huggingface/notebooks/blob/main/examples/nucleotide_transformer_dna_sequence_modelling.ipynb)
and [with LoRA](https://github.com/huggingface/notebooks/blob/main/examples/nucleotide_transformer_dna_sequence_modelling_with_peft.ipynb) on any of the Nucleotide Transformer tasks are also available in HuggingFace example notebooks.
- **SegmentNT**:
An [inference notebook](https://colab.research.google.com/#fileId=https%3A//huggingface.co/InstaDeepAI/segment_nt/blob/main/inference_segment_nt.ipynb) shows how to use the torch [SegmentNT model](https://huggingface.co/InstaDeepAI/segment_nt)
to infer on a given 50kb sequence.

---

## Acknowledgments 🙏

We thank Maša Roller, as well as members of the Rostlab, particularly Tobias Olenyi, Ivan Koludarov,
and Burkhard Rost for constructive discussions that helped identify interesting research directions.
Furthermore, we extend gratitude to all those who deposit experimental data in public databases, to
those who maintain these databases, and those who make analytical and predictive methods freely
available. We also thank the Jax development team.

## Citing our work 📚

If you find this repository useful in your work, please add a relevant citation to
either of our associated papers:

[The Nucleotide Transformer paper](https://www.biorxiv.org/content/10.1101/2023.01.11.523679v2):
```bibtex
@article{dalla2023nucleotide,
  title={The Nucleotide Transformer: Building and Evaluating Robust Foundation Models for Human Genomics},
  author={Dalla-Torre, Hugo and Gonzalez, Liam and Mendoza Revilla, Javier and Lopez Carranza, Nicolas and Henryk Grywaczewski, Adam and Oteri, Francesco and Dallago, Christian and Trop, Evan and Sirelkhatim, Hassan and Richard, Guillaume and others},
  journal={bioRxiv},
  pages={2023--01},
  year={2023},
  publisher={Cold Spring Harbor Laboratory}
}
```


[Agro Nucleotide Transformer paper](https://www.nature.com/articles/s42003-024-06465-2):
```bibtex
@article{mendoza2024foundational,
  title={A foundational large language model for edible plant genomes},
  author={Mendoza-Revilla, Javier and Trop, Evan and Gonzalez, Liam and Roller, Ma{\v{s}}a and Dalla-Torre, Hugo and de Almeida, Bernardo P and Richard, Guillaume and Caton, Jonathan and Lopez Carranza, Nicolas and Skwark, Marcin and others},
  journal={Communications Biology},
  volume={7},
  number={1},
  pages={835},
  year={2024},
  publisher={Nature Publishing Group UK London}
}
```

[SegmentNT paper](https://www.biorxiv.org/content/biorxiv/early/2024/03/15/2024.03.14.584712.full.pdf)
```bibtex
@article{de2024segmentnt,
  title={SegmentNT: annotating the genome at single-nucleotide resolution with DNA foundation models},
  author={de Almeida, Bernardo P and Dalla-Torre, Hugo and Richard, Guillaume and Blum, Christopher and Hexemer, Lorenz and Gelard, Maxence and Pandey, Priyanka and Laurent, Stefan and Laterre, Alexandre and Lang, Maren and others},
  journal={bioRxiv},
  pages={2024--03},
  year={2024},
  publisher={Cold Spring Harbor Laboratory}
}
```

If you have any questions or feedback on the code and models, please feel free to reach out to us.

Thank you for your interest in our work!
