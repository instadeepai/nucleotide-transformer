"""Small run directories the report tests build on. No model, no network."""

import json
from pathlib import Path

import numpy as np

# Production ntv3-annotate output: gene, CDS and intron rows, 1-based closed, with
# the ##ntv3- provenance it writes. gene_1_p starts before the HBE1 window.
ANNOTATION_GFF3 = """\
##gff-version 3
##sequence-region chr11 1 135086622
##ntv3-checkpoint /cache/NTv3_annot_mammals/jax_model
##ntv3-model-name new_ntv3_700M_annotate_mammals_1mb
##ntv3-sequence-length 1048576
##ntv3-keep-target-center-fraction 0.75
##ntv3-decoder hmm_constrained
##ntv3-apply-biological-filters true
##ntv3-min-cds-length 201
##ntv3-use-bf16-compute true
##ntv3-backend gpu
chr11\tNTv3_HMM\tgene\t5266001\t5267000\t.\t+\t.\tID=gene_1_p
chr11\tNTv3_HMM\tCDS\t5266001\t5266200\t.\t+\t0\tID=gene_1_p_CDS_5266001;Parent=gene_1_p
chr11\tNTv3_HMM\tintron\t5266201\t5266700\t.\t+\t.\tID=gene_1_p_intron_5266201;Parent=gene_1_p
chr11\tNTv3_HMM\tCDS\t5266701\t5267000\t.\t+\t2\tID=gene_1_p_CDS_5266701;Parent=gene_1_p
chr11\tNTv3_HMM\tgene\t5268001\t5269500\t.\t-\t.\tID=gene_2_m
chr11\tNTv3_HMM\tCDS\t5268001\t5268300\t.\t-\t1\tID=gene_2_m_CDS_5268001;Parent=gene_2_m
chr11\tNTv3_HMM\tintron\t5268301\t5269200\t.\t-\t.\tID=gene_2_m_intron_5268301;Parent=gene_2_m
chr11\tNTv3_HMM\tCDS\t5269201\t5269500\t.\t-\t0\tID=gene_2_m_CDS_5269201;Parent=gene_2_m
"""

# What run_annotate records for a mammals ask.
ANNOTATION_PLAN = {
    "clade": "mammals",
    "model": "InstaDeepAI/NTv3_annot_mammals",
    "decoder": "hmm_constrained",
    "fasta": "/data/hbe1_locus.fa",
    "locus": None,
    "records": [{"seqid": "chr11", "length": 135_086_622}],
    "tiles": 172,
    "window_bp": 1_048_576,
    "keep_fraction": 0.75,
    "min_cds_length": 201,
    "apply_biological_filters": True,
    "output": "annotation.gff3",
    "want_critic": True,
    "approved": True,
}


def write_annotation(run_dir: Path, gff3: str = ANNOTATION_GFF3) -> None:
    """Write annotation_plan.json and a finished annotation.gff3."""
    (run_dir / "annotation_plan.json").write_text(json.dumps(ANNOTATION_PLAN))
    (run_dir / "annotation.gff3").write_text(gff3)


def predicted_run(run_dir: Path) -> None:
    """Write the plan, lookup, and prediction artifacts a report reads."""
    Path(run_dir, "plan.yaml").write_text(
        "species: human\nassembly: hg38\ngene: HBE1\n"
        "chrom: chr11\nstart: 5266144\nend: 5270240\n"
        "tissue: K562\ntracks: [CAGE]\n"
        "device: cpu\nwindow_bp: 4096\n"
        "want_report: true\nwant_critic: false\nwant_notebook: false\n"
        "analyses: [predict]\nucsc_tracks: [knownGene, clinvarMain]\n"
        "model: NTv3_tracks\n"
    )
    Path(run_dir, "track_resolve.json").write_text(
        json.dumps(
            {
                "species": "human",
                "assembly": "hg38",
                "tissue": "K562",
                "tracks": [
                    {
                        "query": "CAGE",
                        "gt_accession": "ENCSR000CIO",
                        "model_head": "CNhs11250_P",
                    }
                ],
            }
        )
    )
    tracks = Path(run_dir, "tracks")
    tracks.mkdir()
    np.save(tracks / "tracks.npy", np.array([[0.2, 1.5, 0.4]]))
    (tracks / "accessions.json").write_text('["CNhs11250_P"]')
    (tracks / "metadata.json").write_text(
        '{"chrom":"chr11","start":5267424,"end":5268960,"bin_size":512}'
    )
    ucsc = Path(run_dir, "ucsc")
    ucsc.mkdir()
    (ucsc / "knownGene.json").write_text(
        json.dumps(
            {
                "genome": "hg38",
                "chrom": "chr11",
                "start": 5266144,
                "end": 5270240,
                "knownGene": [
                    {
                        "chrom": "chr11",
                        "chromStart": 5268344,
                        "chromEnd": 5269945,
                        "geneName": "HBE1",
                        "strand": "-",
                        "exonStarts": "5268344,5269128,5269813,",
                        "exonEnds": "5268692,5269350,5269945,",
                    }
                ],
                "_request_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=knownGene",
                "_label": "knownGene",
            }
        )
    )
    (ucsc / "clinvarMain.json").write_text(
        json.dumps(
            {
                "genome": "hg38",
                "chrom": "chr11",
                "clinvarMain": [
                    {
                        "chrom": "chr11",
                        "chromStart": 5268500,
                        "chromEnd": 5268501,
                        "name": "NM_005330.5(HBE1):c.1A>G",
                    }
                ],
                "_request_url": "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=clinvarMain",
                "_label": "clinvarMain",
            }
        )
    )


def write_gene_locus(run_dir: Path) -> None:
    """Write locus.json as `ntv3 plan` does for a gene lookup (text source)."""
    Path(run_dir, "locus.json").write_text(
        json.dumps(
            {
                "gene": "HBE1",
                "assembly": "hg38",
                "chrom": "chr11",
                "gene_start": 5268344,
                "gene_end": 5269945,
                "start": 5266144,
                "end": 5270240,
                "source": "UCSC /search (hgnc) and knownGene (ENST00000292896.3)",
                "tss": 5269944,
                "strand": "-",
                "transcript": "ENST00000292896.3",
            }
        )
    )


def add_saliency_and_compare(run_dir: Path) -> None:
    """Approve saliency and one SNV in the plan and write what predict writes."""
    plan = Path(run_dir, "plan.yaml")
    plan.write_text(
        plan.read_text().replace(
            "analyses: [predict]\n",
            "analyses: [predict, saliency, compare]\n"
            "variant: {chrom: chr11, position: 5268501, ref: A, alt: G}\n",
        )
    )
    figures = Path(run_dir, "figures")
    figures.mkdir(exist_ok=True)
    for name in ("saliency_CNhs11250_P", "disrupt_log2fc"):
        (figures / f"{name}.png").write_bytes(b"png")
        (figures / f"{name}.pdf").write_bytes(b"pdf")
    interpret = Path(run_dir, "interpret")
    interpret.mkdir()
    scores = np.zeros((16, 4), dtype=np.float16)
    scores[3, 2] = 0.5
    scores[9, 0] = -0.25
    np.save(interpret / "CNhs11250_P.npy", scores)
    (interpret / "CNhs11250_P.json").write_text(
        json.dumps(
            {
                "accession": "CNhs11250_P",
                "start": 5268400,
                "end": 5268416,
                "seq": "ACGTACGTACGTACGT",
                "logo": {
                    "selection": "strongest_peak",
                    "center_bp": 5268403,
                    "start": 5268400,
                    "end": 5268416,
                    "png": "figures/saliency_CNhs11250_P.png",
                    "pdf": "figures/saliency_CNhs11250_P.pdf",
                },
            }
        )
    )
    disrupt = Path(run_dir, "disrupt")
    disrupt.mkdir()
    (disrupt / "variant.json").write_text(
        json.dumps(
            {
                "chrom": "chr11",
                "position": 5268501,
                "ref": "A",
                "alt": "G",
                "checkpoint": "NTv3_tracks",
                "window": {"chrom": "chr11", "start": 5266144, "end": 5270240},
                "prediction": {
                    "chrom": "chr11",
                    "start": 5267424,
                    "end": 5268960,
                    "bin_size": 512,
                },
            }
        )
    )
    (disrupt / "effect_percentile.json").write_text(
        json.dumps(
            {
                "signed_effect": -0.3,
                "absolute_effect": 0.3,
                "percentile": None,
                "null_n": 0,
                "seed": 0,
                "per_track": [
                    {
                        "track": "CAGE",
                        "model_head": "CNhs11250_P",
                        "signed_effect": -0.3,
                        "percentile": None,
                    }
                ],
                "method": {"track_reduction": "signed abs-max"},
            }
        )
    )
    np.save(disrupt / "delta_tracks.npy", np.array([[0.0, -0.3, 0.1]]))
