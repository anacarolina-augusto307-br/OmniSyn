"""OmniSyn — análise de sequências com Biopython."""

from __future__ import annotations

# Limitações explícitas (ver também LIMITATIONS.md):
# - ORF prediction is simplistic and may not fully represent fragmented metagenomic assemblies.
# - Codon usage e GC são descritivos; não substituem anotação profissional ou QC de sequenciamento.

LIMITATIONS: tuple[str, ...] = (
    "ORF prediction is simplistic and may not fully represent fragmented metagenomic assemblies.",
    "Codon usage counts raw triplets on the forward frame; not expression-weighted.",
    "Global pairwise alignment is not suitable for distant or rearranged regions.",
)

# Biopython: tradução, alinhamento, GC e manipulação de sequências.
from dataclasses import dataclass, field
from io import StringIO
from typing import Iterator

from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqIO import parse as seqio_parse
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis


@dataclass
class AnalysisResult:
    """Summary statistics for a single DNA/RNA sequence."""

    id: str
    description: str
    length: int
    gc_percent: float
    at_percent: float
    nucleotide_counts: dict[str, int]
    nucleotide_percent: dict[str, float]
    is_rna: bool
    molecular_weight: float
    reverse_complement: str
    mrna: str
    protein: str
    codon_usage: dict[str, int] = field(default_factory=dict)
    reading_frames: dict[str, str] = field(default_factory=dict)

# Guardar as OFRs encontradas durante a varredura e depois salvar as proteínas traduzidas com os códons de entrada e depois de parada
# A interface vai fazer o cofere por lá 
@dataclass
class ORF:
    """Open reading frame hit."""

    frame: int
    strand: str
    start: int
    end: int
    length: int
    sequence: str
    protein: str

# Resultados finais e % de identidade, source e números de matches que aparecem na interface
@dataclass
class AlignmentResult:
    """Pairwise alignment between two sequences."""

    seq_a_id: str
    seq_b_id: str
    aligned_a: str
    aligned_b: str
    score: float
    identity_percent: float
    matches: int
    length: int


def _clean_sequence(text: str) -> str:
    return "".join(text.upper().split())

# O usuário pode colar uma sequência pura ou enviar conteúdo em FASTA/GenBank
def parse_sequences(
    raw: str,
    fmt: str = "fasta",
) -> list[SeqRecord]:
    """Parse sequences from FASTA, GenBank, or plain text."""
    raw = raw.strip()
    if not raw:
        return []
        
  # Se o usuário colou apenas a sequência, crio manualmente um SeqRecord, para manter o restante do pipeline funcionando do mesmo jeito
    if fmt == "plain":
        seq = _clean_sequence(raw)
        if not seq:
            return []
        return [SeqRecord(Seq(seq), id="sequence_1", description="Pasted sequence")]
        
    # Deixar o Biopython fazer o trabalho pesado, para FASTA e GenBank
    handle = StringIO(raw)
    records = list(seqio_parse(handle, fmt))
    return records


def _is_rna(seq: str) -> bool:
    return "U" in seq and "T" not in seq


def _dna_from_record(record: SeqRecord) -> Seq:
    seq = str(record.seq).upper()
    # Padronizar tudo como DNA para evitar duplicação de lógica
    if _is_rna(seq):
        return Seq(seq.replace("U", "T"))
    return Seq(seq)


def _nucleotide_composition(seq: str) -> tuple[dict[str, int], dict[str, float]]:
    bases = ("A", "T", "G", "C")
    counts = {b: seq.count(b) for b in bases}
    total = sum(counts.values()) or 1
    percents = {b: round(100 * counts[b] / total, 2) for b in bases}
    return counts, percents

# Presente: identificar reconrrência dos códons
# Futuro:comparar padrões de uso de códons entre organismos diferentes
def _codon_usage(seq: str) -> dict[str, int]:
    usage = {}
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3 and "N" not in codon:
            if codon in usage:
                usage[codon] += 1
            else:
                usage[codon] = 1
    return dict(sorted(usage.items(), key=lambda x: x[1], reverse=True))


def _reading_frames(seq: Seq) -> dict[str, str]:
    """Calcula as proteínas para as 6 janelas de leitura (três forward, três reverse)."""
    frames = {}
    
    # Sentido Direto (+)
    for frame in range(3):
        offset = frame
        sub = str(seq[offset:])
        protein = str(Seq(sub).translate(table=1, to_stop=False))
        frames[f"+{frame + 1}"] = protein[:500]
        
    # Sentido Reverso Complementar (-)
    rev = seq.reverse_complement()
    for frame in range(3):
        offset = frame
        sub = str(rev[offset:])
        protein = str(Seq(sub).translate(table=1, to_stop=False))
        frames[f"-{frame + 1}"] = protein[:500]
        
    return frames


def analyze_sequence(record: SeqRecord) -> AnalysisResult:
    """Run full analysis on one sequence record."""
    dna = _dna_from_record(record)
    seq_str = str(dna)

    counts, percents = _nucleotide_composition(seq_str)
    gc = round(gc_fraction(dna, ambiguous="ignore") * 100, 2)
    at = round(100 - gc, 2)

    mrna = seq_str.replace("T", "U")
    protein = str(dna.translate(table=1, to_stop=True))

    try:
        mw = round(dna.molecular_weight(), 2)
    except Exception:
        mw = 0.0

    return AnalysisResult(
        id=record.id or "unknown",
        description=(record.description or "").strip(),
        length=len(dna),
        gc_percent=gc,
        at_percent=at,
        nucleotide_counts=counts,
        nucleotide_percent=percents,
        is_rna=_is_rna(str(record.seq).upper()),
        molecular_weight=mw,
        reverse_complement=str(dna.reverse_complement()),
        mrna=mrna,
        protein=protein,
        codon_usage=_codon_usage(seq_str),
        reading_frames=_reading_frames(dna),
    )

#ORF finder propositalmente simples
    #fornecer uma exploração rápida de regiões codificadoras potenciais dentro do Omnisyn
    #não substitui um pipeline
def find_orfs(
    record: SeqRecord,
    min_length: int = 90,
    start_codons: tuple[str, ...] = ("ATG",),
    stop_codons: tuple[str, ...] = ("TAA", "TAG", "TGA"),
) -> list[ORF]:
    """Find open reading frames on both strands."""
    dna = _dna_from_record(record)
    orfs: list[ORF] = []

    def scan(seq: Seq, strand: str, frame_offset: int) -> Iterator[ORF]:
        seq_str = str(seq)
        frame_num = frame_offset + 1
        strand_label = "+" if strand == "forward" else "-"

        for i in range(frame_offset, len(seq_str) - 2, 3):
            codon = seq_str[i : i + 3]
            if codon in start_codons:
                for j in range(i + 3, len(seq_str) - 2, 3):
                    stop = seq_str[j : j + 3]
                    if stop in stop_codons:
                        orf_seq = seq_str[i:j + 3]
                        if len(orf_seq) >= min_length:
                            prot = str(Seq(orf_seq).translate(table=1, to_stop=True))
                            yield ORF(
                                frame=frame_num,
                                strand=strand_label,
                                start=i + 1,
                                end=j + 3,
                                length=len(orf_seq),
                                sequence=orf_seq,
                                protein=prot,
                            )
                        break

    for frame in range(3):
        orfs.extend(scan(dna, "forward", frame))
    rev = dna.reverse_complement()
    for frame in range(3):
        orfs.extend(scan(rev, "reverse", frame))

    orfs.sort(key=lambda o: -o.length)
    return orfs


def protein_properties(protein: str) -> dict[str, float | str | int]:
    """Compute physicochemical properties for a protein sequence."""
    clean = "".join(c for c in protein.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
    if len(clean) < 1:
        return {}
    analysis = ProteinAnalysis(clean)
    return {
        "length": len(clean),
        "molecular_weight": round(analysis.molecular_weight(), 2),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "aromaticity": round(analysis.aromaticity(), 4),
        "instability_index": round(analysis.instability_index(), 2),
        "gravy": round(analysis.gravy(), 4),
        "secondary_structure": analysis.secondary_structure_fraction(),
    }

# Sem criação de algoritmos de alinhamento para evitar lentidão da Omnisyn
# Usar biopython para deixar no pique de caroo de fórmula 1
def pairwise_align(
    record_a: SeqRecord,
    record_b: SeqRecord,
    match: int = 2,
    mismatch: int = -1,
    gap_open: int = -2,
    gap_extend: int = -0.5,
) -> AlignmentResult | None:
    a = str(_dna_from_record(record_a))
    b = str(_dna_from_record(record_b))

    if not a or not b:
        return None
# Extração de estátistica de identidade
    alignments = pairwise2.align.globalms(
        a, b, match, mismatch, gap_open, gap_extend, one_alignment_only=True
    )
    if not alignments:
        return None

    aln = alignments[0]
    aligned_a, aligned_b = aln.seqA, aln.seqB
    matches = sum(1 for x, y in zip(aligned_a, aligned_b) if x == y and x != "-")
    length = len(aligned_a)
    identity = round(100 * matches / length, 2) if length else 0.0

    return AlignmentResult(
        seq_a_id=record_a.id or "A",
        seq_b_id=record_b.id or "B",
        aligned_a=aligned_a,
        aligned_b=aligned_b,
        score=aln.score,
        identity_percent=identity,
        matches=matches,
        length=length,
    )

#Uso das "sliding window" (50bp)
    # Evitar erros ao calcular sequências de 5000 caracteres/digitos/codons 
    # Fazer a média global do GC para coletar áreas de baixa e alta presença de GC
        # Lembrar que GC ajuda a achar os genes ativos do genoma 
def sliding_gc(seq: str, window: int = 50) -> list[tuple[int, float]]:
     """Compute local GC content using overlapping windows."""
    seq = _clean_sequence(seq)
    if len(seq) < window:
        gc = round(gc_fraction(Seq(seq), ambiguous="ignore") * 100, 2) if seq else 0
        return [(1, gc)]

    points: list[tuple[int, float]] = []
    for i in range(0, len(seq) - window + 1, max(1, window // 10)):
        chunk = seq[i : i + window]
        gc = round(gc_fraction(Seq(chunk), ambiguous="ignore") * 100, 2)
        points.append((i + 1, gc))
    return points
