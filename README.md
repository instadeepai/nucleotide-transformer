<p align="center">
  <img src="imgs/instadeep_logo.png" alt="InstaDeep AI for Genomics Logo" width="200"/>
</p>

<h1 align="center">AI Foundation Models for Genomics</h1>

<p align="center">
  <strong>A hub for InstaDeep's cutting-edge deep learning models and research for genomics, originating from the Nucleotide Transformer and its evolutions.</strong>
</p>

<p align="center">
  <a href="./LICENSE"> <img src="https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-blue.svg" alt="License: CC BY-NC-SA 4.0">
  </a>
  <a href="https://docs.python.org/3.8/library/index.html">
    <img src="https://img.shields.io/badge/python-3.11-blue.svg" alt="Python 3.8">
  </a>
  <a href="https://jax.readthedocs.io/en/latest/"> <img src="https://img.shields.io/badge/jax-0.3.25+-informational" alt="Jax 0.3.25+">
  </a>
  <a href="https://huggingface.co/InstaDeepAI"> <img src="https://img.shields.io/badge/%F0%9F%A4%97%20Hugging%20Face-InstaDeepAI-orange" alt="Hugging Face Models">
  </a>
</p>

---

## 🎯 Our Focus: Advancing Genomics with AI

Welcome to the InstaDeep AI for Genomics repository! This is where we feature our collection of transformer-based genomic language models and innovative downstream applications. Our work in the genomics space began with **The Nucleotide Transformer**, developed in collaboration with Nvidia and TUM and trained on Cambridge-1, and has expanded to include projects like the **Agro Nucleotide Transformer** (in collaboration with Google, trained on TPU-v4 accelerators), **SegmentNT**, and **ChatNT**.

Our mission is to provide the scientific community with powerful, reproducible, and accessible tools to unlock new insights from biological sequences. This repository serves as the central place for sharing our models, inference code, pre-trained weights, and research contributions in the genomics domain, with explorations into future areas like single-cell transcriptomics.

We are thrilled to open-source these works and provide the community with access to the code and pre-trained weights for our diverse set of genomics language models and segmentation models.

## ✨ Featured Models & Research Evolutions

This section highlights the key models and research directions from our team. Each entry provides a brief overview and links to detailed documentation, publications, and resources. *(Detailed code examples, setup for specific models, and in-depth figures are now located in their respective documentation pages within the `./docs` folder.)*

---

### 🧬 The Nucleotide Transformer (NT)

Our foundational language models leverage DNA sequences from over 3,200 diverse human genomes and 850 genomes from a wide range of species. These models provide extremely accurate molecular phenotype prediction compared to existing methods. *This family includes multiple variants (e.g., 500M_human_ref, 2B5_1000G, NT-v2 series) which are detailed further in the specific documentation.*

* **Keywords:** Foundational Model, Genomics, DNA/RNA, Pre-trained, Sequence Embeddings, Phenotype Prediction
* ➡️ **[Model Details, Variants & Usage](./docs/nucleotide_transformer.md)**
* 📜 **[Read the Paper (Nature Methods 2025)](https://www.nature.com/articles/s41592-024-02523-z)**
* 🤗 **[Hugging Face Collection](https://huggingface.co/collections/InstaDeepAI/nucleotide-transformer-65099cdde13ff96230f2e592)**
* 🚀 **Fine-tuning Notebooks (HF): ([LoRA](https://github.com/huggingface/notebooks/blob/main/examples/nucleotide_transformer_dna_sequence_modelling_with_peft.ipynb) and [regular](https://github.com/huggingface/notebooks/blob/main/examples/nucleotide_transformer_dna_sequence_modelling.ipynb))**

---

### 🌾 Agro Nucleotide Transformer (AgroNT)

A novel foundational large language model trained on reference genomes from 48 plant species, with a predominant focus on crop species. AgroNT demonstrates state-of-the-art performance across several prediction tasks ranging from regulatory features, RNA processing, and gene expression in plants.

* **Keywords:** Plant Genomics, Foundational Model, Crop Science, Gene Expression, Agriculture AI
* ➡️ **[Model Details & Usage](./docs/agro_nucleotide_transformer.md)**
* 📜 **[Read the Paper (Communications Biology 2024)](https://www.nature.com/articles/s42003-024-06465-2)**
* 🤗 **[Hugging Face Collection](https://huggingface.co/collections/InstaDeepAI/agro-nucleotide-transformer-65b25c077cd0069ad6f6d344)**

---

### 🧩 SegmentNT (& family: SegmentEnformer, SegmentBorzoi)

Segmentation models using transformer backbones (Nucleotide Transformers, Enformer, Borzoi) for predicting genomic elements at single-nucleotide resolution. SegmentNT, for instance, predicts 14 different classes of human genomic elements in sequences up to 30kb (generalizing to 50kbp) and demonstrates superior performance.

* **Keywords:** Genome Segmentation, Single-Nucleotide Resolution, Genomic Elements, U-Net, Enformer, Borzoi
* ➡️ **[Model Details & Usage](./docs/segment_nt.md)** (Covers SegmentNT, SegmentEnformer, SegmentBorzoi)
* 📜 **[Read the Paper (bioRxiv preprint)](https://www.biorxiv.org/content/10.1101/2024.03.14.584712v1)**
* 🤗 **[Hugging Face Collection](https://huggingface.co/collections/InstaDeepAI/segmentnt-65eb4941c57808b4a3fe1319)**
* 🚀 **[SegmentNT Inference Notebook (HF)](https://colab.research.google.com/#fileId=https%3A//huggingface.co/InstaDeepAI/segment_nt/blob/main/inference_segment_nt.ipynb)**

---

### 💬 ChatNT

A multimodal conversational agent designed with a deep understanding of DNA biological sequences, enabling interactive exploration and analysis of genomic data through natural language.

* **Keywords:** Conversational AI, Multimodal, DNA Analysis, Genomics Chatbot, Interactive Biology
* ➡️ **[Model Details & Usage](./docs/chat_nt.md)**
* 📜 **[Read the Paper (Nature Machine Intelligence 2025)](https://www.nature.com/articles/s42256-025-01047-1)**
* 🤗 **[ChatNT on Hugging Face](https://huggingface.co/InstaDeepAI/ChatNT)**
* 🚀 **[ChatNT Inference Notebook (Jax)](./notebooks/chat_nt/inference.ipynb)**

---

### 3️⃣ Codon-NT (Exploring 3-mer Tokenization)

A Nucleotide Transformer model variant trained on 3-mers (codons). This work investigates alternative tokenization strategies for genomic language models and their impact on downstream performance and interpretability.

* **Keywords:** Genomics, Language Model, Codon, Tokenization, 3-mers, Nucleotide Transformer Variant
* ➡️ **[Model Details & Usage](./docs/codon_nt.md)**
* 📜 **[Read the Paper (Bioinformatics 2024)](https://academic.oup.com/bioinformatics/article/40/9/btae529/7745814)**
* 🤗 **[Hugging Face Link](https://huggingface.co/InstaDeepAI/nucleotide-transformer-v2-50m-3mer-multi-species)**

---

### 🧬 Isoformer

A model designed for learning isoform-aware embeddings directly from RNA-seq data, enabling a deeper understanding of transcript-specific expression and regulation.

* **Keywords:** RNA-seq, Transcriptomics, Isoforms, Gene Expression, Embeddings
* ➡️ **[Model Details & Usage](./docs/isoformer.md)**
* 📜 **[Read the Paper (NeurIPS 2024)](https://papers.nips.cc/paper_files/paper/2024/file/8f6b3692297e49e5d5c91ba00281379c-Paper-Conference.pdf)**
* 🤗 **[Hugging Face Link](https://huggingface.co/InstaDeepAI/isoformer)**
* 🚀 **[Isoformer Inference Notebook (HF)](./notebooks/isoformer/inference.ipynb)**

---

### 🔬 sCT (single-Cell Transformer)

Our foundational transformer model for single-cell and spatial transcriptomics data. sCT aims to learn rich representations from complex, high-dimensional single-cell datasets to improve various downstream analytical tasks.

* **Keywords:** Single-cell RNA-seq, Spatial Transcriptomics, Foundational Model, Transformer, Gene Expression
* ➡️ **[Model Details & Usage](./docs/sct.md)**
* 📜 **[Read the Paper (OpenReview preprint)](https://openreview.net/forum?id=VdX9tL3VXH)**
* 🤗 **[Hugging Face Link](https://huggingface.co/InstaDeepAI/sCellTransformer)**
* 🚀 **[sCT Inference Notebook (HF)](./notebooks/sct/inference_sCT_pytorch_example.ipynb)**

---

### 🧪 BulkRNABert

BulkRNABert is a transformer-based, encoder-only foundation model designed for bulk RNA-seq data. It learns biologically meaningful representations from large-scale transcriptomic profiles.

* **Keywords:** Bulk RNA-seq, Foundational Model, Transformer, Cancer prognosis
* ➡️ **[Model Details & Usage](./docs/bulk_rna_bert.md)**
* 📜 **[Read the Paper (Machine Learning for Health 2024)](https://proceedings.mlr.press/v259/gelard25a.html)**
* 🤗 **[Hugging Face Link](https://huggingface.co/InstaDeepAI/BulkRNABert)**
* 🚀 **[BulkRNABert Inference Notebook (HF)](notebooks/bulk_rna_bert/inference_bulkrnabert_pytorch_example.ipynb)**

---

### 🔗 MOJO (<u>M</u>ulti-<u>O</u>mics <u>JO</u>int representation)

MOJO is a multimodal model designed learn embeddings of multi-omics data. It integrates bulk RNA-seq and DNA methylation data to generate powerful joint representations tailored for cancer-type classification and survival analysis.

* **Keywords:** Bulk RNA-seq, DNA Methylation, Foundational Model, Transformer, Multimodal, Cancer prognosis
* ➡️ **[Model Details & Usage](./docs/mojo.md)**
* 📜 **[Read the Paper (ICML Workshop on Generative AI and Biology 2025)](https://www.biorxiv.org/content/10.1101/2025.06.25.661237v1)**
* 🤗 **[Hugging Face Link](https://huggingface.co/InstaDeepAI/MOJO)**
* 🚀 **[MOJO Inference Notebook (HF)](./notebooks/mojo/inference_mojo_pytorch_example.ipynb)**

---

## 💡 Why Choose InstaDeep's Genomic Models?

* **Built on Strong Foundations:** Leveraging large-scale pre-training and diverse genomic datasets.
* **Cutting-Edge Research:** Incorporating the latest advancements in deep learning for biological sequence analysis.
* **High Performance:** Designed and validated to achieve state-of-the-art results on challenging genomic tasks.
* **Open and Accessible:** We provide pre-trained weights, usage examples, and aim for easy integration into research workflows.
* **Collaborative Spirit:** Developed with leading academic and industry partners.
* **Focused Expertise:** Created by a dedicated team specializing in AI for genomics at InstaDeep.

## 🚀 Getting Started

To begin using models from this repository:

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/instadeepai/nucleotide-transformer.git
    cd nucleotide-transformer
    ```
2.  **Set up your environment (virtual environment recommended):**
    ```bash
    python -m venv .venv
    source .venv/bin/activate # On Windows use `source .venv\Scripts\activate`
    ```
3.  **Install the package and dependencies:**
    ```bash
    pip install . # Installs the local package
    # Or, for a general requirements file if you have one:
    # pip install -r requirements.txt
    ```

For detailed instructions on individual models, including specific dependencies, downloading pre-trained weights, and Python usage examples, please refer to their dedicated documentation pages linked in the "Featured Models & Research Evolutions" section above (e.g., `./docs/nucleotide_transformer.md`).

## 🤝 Community & Support

* **Questions & Bug Reports:** Please use the [GitHub Issues](https://github.com/instadeepai/nucleotide-transformer/issues) page.
* **Discussions:** For broader discussions or questions, please use the [GitHub Discussions](https://github.com/instadeepai/nucleotide-transformer/discussions) tab (if enabled).
* **Stay Updated:** Follow InstaDeep's official channels for announcements on new model releases and research updates.
