# sCT

sCT (single-Cell Transformer) is our foundational transformer model for single-cell and spatial transcriptomics data. sCT aims to learn rich representations from complex, high-dimensional single-cell datasets to improve various downstream analytical tasks.

sCT processes raw gene expression profiles across multiple cells to predict discretized gene expression levels for unseen cells without retraining. 
The model can handle up to 20,000 protein-coding genes and a bag of 50 cells in the same sample. 
This ability (around a million-gene expressions tokens) allows it to learn cross-cell relationships and 
capture long-range dependencies in gene expression data, and to mitigate the sparsity typical in single-cell datasets.

sCT is trained on a large dataset of single-cell RNA-seq and finetuned on spatial transcriptomics data. Evaluation tasks include zero-shot imputation of masked gene expression, and zero-shot prediction of cell types.

* 📜 **[Read the Paper (OpenReview preprint)](https://openreview.net/forum?id=VdX9tL3VXH)**
* 🤗 **[Hugging Face Link](https://huggingface.co/InstaDeepAI/sCellTransformer)**

## Training data

The model was trained following a two-step procedure: pre-training on single-cell data, then finetuning on spatial transcriptomics data. The single-cell data used for pre-training, comes from the Cellxgene Census collection datasets used to train the scGPT models. It consists of around 50 millions cells and approximately 60,000 genes. The spatial data comes from both the human breast cell atlas and the human heart atlas.

## Training procedure

As detailed in the paper, the gene expressions are first binned into a pre-defined number of bins. This allows the model to better learn the distribution of the gene expressions through sparsity mitigation, noise reduction, and extreme-values handling. Then, the training objective is to predict the masked gene expressions in a cell, following a BERT-like style training.

## How to use 🚀

We make sCT available in Jax in this repository and in PyTorch on HuggingFace. Examples on how to use it at:
- Jax: `../notebooks/sct/inference_sCT_jax_example.ipynb`
- PyTorch: `../notebooks/sct/inference_sCT_pytorch_example.ipynb`.

## Citing our work 📚

You can cite our model at:

```bibtex
@misc{joshi2025a,
    title={A long range foundation model for zero-shot predictions in single-cell and
    spatial transcriptomics data},
    author={Ameya Joshi and Raphael Boige and Lee Zamparo and Ugo Tanielian and Juan Jose
    Garau-Luis and Michail Chatzianastasis and Priyanka Pandey and Janik Sielemann and
    Alexander Seifert and Martin Brand and Maren Lang and Karim Beguir and Thomas PIERROT},
    year={2025},
    url={https://openreview.net/forum?id=VdX9tL3VXH}
}
```