import os

from setuptools import find_packages, setup

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(CURRENT_DIR, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

# Read rather than repeated, so these pins cannot drift from the ones CI
# installs. See the file for why the JAX stack has to move as one unit.
with open(os.path.join(CURRENT_DIR, "requirements-jax.txt"), encoding="utf-8") as f:
    jax_requirements = [
        line.strip()
        for line in f
        if line.strip() and not line.startswith("#")
    ]

# PyPI's jaxlib is CPU-only, and nothing fails when it is used on a GPU host:
# JAX falls back silently and ntv3-annotate takes ~7.5 days over a genome
# instead of ~45 minutes. The `gpu` extra adds CUDA over the SAME bounds,
# derived rather than retyped so the two cannot disagree about which jax is
# meant -- and raising here rather than shipping an empty extra, which would
# reintroduce the silent CPU fallback this exists to prevent.
gpu_requirements = [
    req.replace("jax>", "jax[cuda12]>")
    for req in jax_requirements
    if req.startswith("jax>")
]
if not gpu_requirements:
    raise RuntimeError(
        "requirements-jax.txt no longer has a 'jax>' line, so the gpu extra "
        "would install nothing and a GPU host would quietly run on the CPU."
    )

setup(
    name="nucleotide_transformer",
    version="0.0.1",
    packages=find_packages(),
    package_data={"nucleotide_transformer_v3.gene_annotation": ["model_cfg/*.yaml"]},
    extras_require={"gpu": gpu_requirements},
    entry_points={
        "console_scripts": [
            "ntv3-annotate=nucleotide_transformer_v3.gene_annotation.cli:main",
        ],
    },
    url="https://github.com/instadeepai/nucleotide-transformer",
    license="CC BY-NC-SA 4.0",
    author="InstaDeep Ltd",
    python_requires=">=3.10",
    description="The Nucleotide Transformer: Building and Evaluating "
    "Robust Foundation Models for Human Genomics ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=[
        *jax_requirements,
        "absl-py>=1.0.0",
        "dm-haiku>=0.0.9",
        "numpy>=1.23.5",
        "typing_extensions>=3.10.0",
        "joblib>=1.2.0",
        "tqdm>=4.56.0",
        "regex>=2022.1.18",
        "huggingface-hub>=0.23.0",
        "dm-haiku>=0.0.14",
        "anndata>=0.11.0",
        "scanpy>=1.11.2",
        "cellxgene_census>=1.17.0",
        "transformers>=4.52.4",
        "einops>=0.8.1",
        "pydantic==1.10.13",
        "ipykernel>=6.0.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "pyfaidx>=0.7.0",
        # ntv3 predict reads ENCODE ground-truth BigWigs with it, but only from a
        # build with libcurl (pyBigWig.remote). With the PyPI wheel, predict
        # saves the tracks and records why ground truth was not compared.
        "pyBigWig>=0.3.22",
        "requests>=2.25.0",
        "pyyaml>=6.0",
    ],
    dependency_links=[
        "https://storage.googleapis.com/jax-releases/jax_releases.html",
    ],
    keywords=["Genomics", "Language Model", "Deep Learning", "JAX"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
    ],
)
