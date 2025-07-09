# MOJO

MOJO (MultiOmics JOint representation learning) is a multimodal model designed learn embeddings of multi-omics data. It integrates bulk RNA-seq and DNA methylation data to generate powerful joint representations. These representations are tailored for improving predictive performance on downstream tasks like cancer-type classification and survival analysis.

* 📜 **[Read the Paper (ICML Workshop on Generative AI and Biology 2025)](https://www.biorxiv.org/content/10.1101/2025.06.25.661237v1)**
* 🤗 **[Hugging Face Link](https://huggingface.co/InstaDeepAI/MOJO)**

## Training data

MOJO is trained on multi-modal data from The Cancer Genome Atlas (TCGA) dataset. The training data consists of paired bulk RNA-seq and DNA methylation profiles, enabling the model to learn the complex interplay between transcription and epigenetic regulation.
The model takes as input sequences of 17116 gene expression and DNA methylation. The gene ids that must be used and in which order they should appear in the sequence as provided in `../notebooks/data/mojo_gene_names.txt`

## Training procedure

MOJO uses a bimodal masked language modeling objective. It is trained to simultaneously predict masked values in both the RNA-seq and DNA methylation data modalities. This process forces the model to learn the intricate cross-modal relationships between gene expression and epigenetic regulation, leading to robust, integrated representations.

## How to use 🚀

We make MOJO available in Jax in this repository and in PyTorch on HuggingFace. Examples on how to use it at:
- Jax: `../notebooks/mojo/inference_mojo_jax_example.ipynb`.
- PyTorch: `../notebooks/mojo/inference_mojo_pytorch_example.ipynb`.

## Citing our work 📚

You can cite our model at:

```bibtex
@article {G{\'e}lard2025.06.25.661237,
	author = {G{\'e}lard, Maxence and Benkirane, Hakim and Pierrot, Thomas and Richard, Guillaume and Courn{\`e}de, Paul-Henry},
	title = {Bimodal masked language modeling for bulk RNA-seq and DNA methylation representation learning},
	elocation-id = {2025.06.25.661237},
	year = {2025},
	doi = {10.1101/2025.06.25.661237},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/06/27/2025.06.25.661237},
	journal = {bioRxiv}
}
```
