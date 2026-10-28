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

"""Class layout and per-base codes for the 15-state annotation task.
"""

# The 15 states the annotation head predicts, in index order. The order is the
# checkpoint's: it is what the head's output axis means, so it is not free.
GENOME_ANNOTATION_CLASSES: tuple[str, ...] = (
    "intergenic",  # 0
    "intron_p0",  # 1
    "intron_p1",  # 2
    "intron_p2",  # 3
    "CDS_p0",  # 4
    "CDS_p1",  # 5
    "CDS_p2",  # 6
    "splice_acceptor_p0",  # 7
    "splice_acceptor_p1",  # 8
    "splice_acceptor_p2",  # 9
    "splice_donor_p0",  # 10
    "splice_donor_p1",  # 11
    "splice_donor_p2",  # 12
    "start_codon",  # 13
    "stop_codon",  # 14
)

CLASS_NAME_TO_INDEX: dict[str, int] = {
    label: i for i, label in enumerate(GENOME_ANNOTATION_CLASSES)
}

NUM_GENOME_ANNOTATION_CLASSES: int = len(GENOME_ANNOTATION_CLASSES)

STRAND_NAMES: list[str] = ["plus", "minus"]

# A=0, C=1, G=2, T=3, N=4. reverse-complement relies on T=3 / N=4:
# complement(x) = NUC_T - x for x in {A, C, G, T}, and N maps to itself.
NUC_A, NUC_C, NUC_G, NUC_T, NUC_N = 0, 1, 2, 3, 4
