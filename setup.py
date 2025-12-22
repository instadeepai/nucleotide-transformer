import os

from setuptools import find_packages, setup

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(CURRENT_DIR, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="nucleotide_transformer",
    version="0.0.1",
    packages=find_packages(),
    url="https://github.com/instadeepai/nucleotide-transformer",
    license="CC BY-NC-SA 4.0",
    author="InstaDeep Ltd",
    python_requires=">=3.9",
    description="The Nucleotide Transformer: Building and Evaluating "
    "Robust Foundation Models for Human Genomics ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        "absl-py>=1.0.0",
        "jax>=0.3.25",
        "jaxlib>=0.3.25",
        "dm-haiku>=0.0.9",
        "numpy>=1.23.5,<2.0.0",
        "typing_extensions>=3.10.0",
        "joblib>=1.2.0",
        "tqdm>=4.56.0",
        "regex>=2022.1.18",
        "huggingface-hub>=0.23.0",
        "dm-haiku>=0.0.14",
        "anndata>=0.11.0",
        "jax>=0.6.0",
        "scanpy>=1.11.2",
        "cellxgene_census>=1.17.0",
        "transformers>=4.52.4",
        "torch>=2.7.1",
        "einops>=0.8.1",
        "pydantic==1.10.13",
        "flax==0.10.4",
        "ipykernel>=6.0.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "pyfaidx>=0.7.0",
        "requests>=2.25.0",
    ],
    dependency_links=[
        "https://storage.googleapis.com/jax-releases/jax_releases.html",
    ],
    keywords=["Genomics", "Language Model", "Deep Learning", "JAX"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
)
