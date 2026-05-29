"""OmniSyn — núcleo de análise de sequências (Biopython)."""

from analyzer.core import (
    LIMITATIONS,
    AnalysisResult,
    analyze_sequence,
    find_orfs,
    pairwise_align,
    parse_sequences,
)

__all__ = [
    "LIMITATIONS",
    "AnalysisResult",
    "analyze_sequence",
    "find_orfs",
    "pairwise_align",
    "parse_sequences",
]
