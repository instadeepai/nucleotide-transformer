# Copyright 2022 InstaDeep Ltd
#
# Licensed under the Creative Commons BY-NC-SA 4.0 License (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      https://creativecommons.org/licenses/by-nc-sa/4.0/
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Tests that every stitched seam lies on the grid ``distance_to_seam`` assumes."""

import pytest

from nucleotide_transformer_v3.gene_annotation.check_equivalence import distance_to_seam
from nucleotide_transformer_v3.gene_annotation.utils import plan_contig_windows

# (sequence_length, keep_target_center_fraction): the shipped 1 Mb geometry and
# the golden tests' 4 kb one.
_GEOMETRIES = [(1_048_576, 0.75), (4_096, 0.75)]


def _retained_regions(
    sequence_length: int, fraction: float, contig_length: int
) -> tuple[list[tuple[int, int]], int]:
    """Each window's stitched ``[start, end)`` slice of the contig, and the stride."""
    windows, crop = plan_contig_windows(
        sequence_length=sequence_length,
        keep_target_center_fraction=fraction,
        contig_length=contig_length,
        chrom="chr1",
    )
    regions = [
        (w.keep_start, w.keep_start + w.slice_stop - w.slice_start) for w in windows
    ]
    return regions, sequence_length - 2 * crop


@pytest.mark.parametrize(("sequence_length", "fraction"), _GEOMETRIES)
@pytest.mark.parametrize(
    "contig_strides",
    # Contig length in strides: below one, exactly one, between, an exact
    # multiple, and long contigs whose tail window is shifted back a long way.
    [0.4, 1.0, 1.3, 3.0, 6.36, 53.4],
)
def test_every_seam_of_a_contig_is_on_the_stride_grid(
    sequence_length: int, fraction: float, contig_strides: float
) -> None:
    """The shifted tail window's overlap is dropped, so no seam leaves the grid."""
    _, stride = _retained_regions(sequence_length, fraction, sequence_length)
    contig_length = round(contig_strides * stride)
    regions, stride = _retained_regions(sequence_length, fraction, contig_length)

    assert regions[0][0] == 0 and regions[-1][1] == contig_length
    seams = [end for (_, end), (start, _) in zip(regions, regions[1:])]
    assert seams == [start for start, _ in regions[1:]], "the tiles do not abut"
    assert all(seam % stride == 0 for seam in seams)
    assert all(distance_to_seam(seam + 1, stride) == 0 for seam in seams)


def test_distance_to_seam_measures_to_the_nearest_grid_point() -> None:
    """Distances are taken from the 1-based start to the nearest multiple."""
    assert distance_to_seam(1, 100) == 0
    assert distance_to_seam(101, 100) == 0
    assert distance_to_seam(131, 100) == 30
    assert distance_to_seam(171, 100) == 30
