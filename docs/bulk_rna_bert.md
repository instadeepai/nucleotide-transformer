# BulkRNABert

BulkRNABert is a transformer-based, encoder-only foundation model designed for bulk RNA-seq data. It learns biologically meaningful representations from large-scale transcriptomic profiles. Once pre-trained, BulkRNABert can be fine-tuned for various cancer-related downstream tasks, such as cancer type classification or survival analysis, by using its learned embeddings.

* 📜 **[Read the Paper (Machine Learning for Health 2024)](https://proceedings.mlr.press/v259/gelard25a.html)**
* 🤗 **[Hugging Face Link](https://huggingface.co/InstaDeepAI/BulkRNABert)**

## Training data

The model was pre-trained on a large dataset of bulk RNA-seq profiles from The Cancer Genome Atlas (TCGA) dataset.
The model takes as input sequences of 19062 gene expression. The gene ids that must be used and in which order they should appear in the sequence as provided in `../notebooks/data/bulkrnabert_gene_ids.txt`



## Training procedure

Following the original BERT framework, BulkRNABert uses a self-supervised, masked language modeling objective. During pre-training, gene expression values are randomly masked, and the model is tasked with reconstructing these values from their surrounding genomic context. This process enables the model to learn rich, contextual representations of transcriptomic profiles.

## How to use 🚀

We make BulkRNABert available in Jax in this repository and in PyTorch on HuggingFace. Examples on how to use it at:
- Jax: `../notebooks/bulk_rna_bert/inference_bulkrnabert_jax_example.ipynb`.
- PyTorch: `../notebooks/bulk_rna_bert/inference_bulkrnabert_pytorch_example.ipynb`.

## Citing our work 📚

You can cite our model at:

```bibtex
@InProceedings{pmlr-v259-gelard25a,
  title = 	 {BulkRNABert: Cancer prognosis from bulk RNA-seq based language models},
  author =       {G{\'{e}}lard, Maxence and Richard, Guillaume and Pierrot, Thomas and Courn{\`{e}}de, Paul-Henry},
  booktitle = 	 {Proceedings of the 4th Machine Learning for Health Symposium},
  pages = 	 {384--400},
  year = 	 {2025},
  editor = 	 {Hegselmann, Stefan and Zhou, Helen and Healey, Elizabeth and Chang, Trenton and Ellington, Caleb and Mhasawade, Vishwali and Tonekaboni, Sana and Argaw, Peniel and Zhang, Haoran},
  volume = 	 {259},
  series = 	 {Proceedings of Machine Learning Research},
  month = 	 {15--16 Dec},
  publisher =    {PMLR},
  url = 	 {https://proceedings.mlr.press/v259/gelard25a.html},
}
```
