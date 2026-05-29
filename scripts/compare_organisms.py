"""Compara GC%, ORFs e uso de códon entre amostras reais em sample_data/."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Bio.SeqIO import parse

from analyzer.core import analyze_sequence, find_orfs

SAMPLE_FILES = [
    ("solo (bactéria)", ROOT / "sample_data" / "bacillus_subtilis_soil.fasta"),
    ("marinho (bactéria)", ROOT / "sample_data" / "vibrio_cholerae_marine.fasta"),
    ("archaea (halófila)", ROOT / "sample_data" / "haloferax_volcanii_archaea.fasta"),
]

MIN_ORF_BP = 90
TOP_CODONS = 10


def top_codons(codon_usage: dict[str, int], n: int = TOP_CODONS) -> list[tuple[str, int]]:
    items = sorted(codon_usage.items(), key=lambda x: -x[1])
    return items[:n]


def compare() -> str:
    lines = [
        "OmniSyn — comparação entre organismos (sample_data/)",
        "Métricas: GC%, ORFs (min {} bp), top {} códons.".format(MIN_ORF_BP, TOP_CODONS),
        "",
        "# ORF prediction is simplistic and may not fully represent fragmented metagenomic assemblies.",
        "",
    ]

    for label, path in SAMPLE_FILES:
        if not path.exists():
            lines.append(f"[ERRO] Arquivo ausente: {path}")
            continue

        record = next(parse(path, "fasta"))
        result = analyze_sequence(record)
        orfs = find_orfs(record, min_length=MIN_ORF_BP)

        lines.extend(
            [
                f"## {label}",
                f"ID: {result.id}",
                f"Comprimento: {result.length:,} bp",
                f"GC%: {result.gc_percent}",
                f"AT%: {result.at_percent}",
                f"ORFs encontrados (>={MIN_ORF_BP} bp): {len(orfs)}",
            ]
        )

        if orfs:
            longest = orfs[0]
            lines.append(
                f"Maior ORF: {longest.length} bp (frame {longest.strand}{longest.frame}, "
                f"pos {longest.start}-{longest.end})"
            )

        lines.append("Top códons (contagem):")
        for codon, count in top_codons(result.codon_usage):
            lines.append(f"  {codon}: {count}")
        lines.append("")

    return "\n".join(lines)


def main() -> None:
    report = compare()
    out_path = ROOT / "sample_data" / "comparison_report.txt"
    out_path.write_text(report, encoding="utf-8")
    print(report)
    print(f"\nRelatório salvo em: {out_path}")


if __name__ == "__main__":
    main()
