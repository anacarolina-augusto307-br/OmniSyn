#!/usr/bin/env python3
""" Interface de linha de comando - Omnisyn """

import argparse
import sys
from pathlib import Path

#Mantendo lógica do meu núcleo
from analyzer.core import (
  analyze_sequence,
  find_orfs,
  pairwise-align,
  parse_sequences, 
)

# Para executar analise individual da sequência
def handle_analyze(args):
    if not Path(args.input).exists():
        print(f"error: file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)
        
    raw_content = Path(args.input).read_text(encoding="utf-8")
    records = parse_sequences(raw_content, fmt=args.format)
    
    if not records:
        print("No valid sequences were found", file=sys.stderr)
        return

    for record in records:
        res = analyze_sequence(record)
        print("\n" + "="*50)
        print(f"sequence id: {res.id}")
        print(f"description: {res.description}")
        print("-" * 50)
        print(f"length: {res.length:,} bp")
        print(f"gc content: {res.gc_percent}%  |  at content: {res.at_percent}%")
        print(f"molecular weight: {res.molecular_weight:,.2f} Da")
        print(f"type: {'RNA' if res.is_rna else 'DNA'}")
        print("\nnucleotide composition:")
        for base, count in res.nucleotide_counts.items():
            pct = res.nucleotide_percent[base]
            print(f"  [{base}]: {count:<8} ({pct:.2f}%)")
        print("="*50)

def handle_orfs(args):
    """Find ORFs."""
    if not Path(args.input).exists():
        print(f"error: file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)
        
    raw_content = Path(args.input).read_text(encoding="utf-8")
    records = parse_sequences(raw_content, fmt=args.format)
    
    for record in records:
        print(f"\nsearching ORFs in '{record.id}' (min: {args.min_len} bp)")
        orfs = find_orfs(record, min_length=args.min_len)
        
        if not orfs:
            print("no ORFs found")
            continue
            
        print(f"Found {len(orfs)} ORF(s). Showing the largest ones:")
        print(f"{'Strand':<5} | {'Frame':<5} | {'Start':<8} | {'End':<8} | {'Length':<12}")
        print("-" * 48)
        for idx, o in enumerate(orfs[:args.limit]):
            print(f"{o.strand:<5} | {o.frame:<5} | {o.start:<8} | {o.end:<8} | {o.length:<12}")

def handle_align(args):
    """Align two sequences."""
    if not Path(args.file1).exists():
        print(f"error: file '{args.file1}' not found.", file=sys.stderr)
        sys.exit(1)
        
    content1 = Path(args.file1).read_text(encoding="utf-8")
    records_a = parse_sequences(content1, fmt=args.format1)
    
    if args.file2:
        if not Path(args.file2).exists():
            print(f"error: file '{args.file2}' not found.", file=sys.stderr)
            sys.exit(1)
        content2 = Path(args.file2).read_text(encoding="utf-8")
        records_b = parse_sequences(content2, fmt=args.format2)
    else:
        # Se o segundo arquivo não foi passado, assume-se que as duas sequências estão dentro do primeiro arquivo
        records_b = records_a

    if len(records_a) < 1 or len(records_b) < (2 if args.file2 is None else 1):
        print("❌ Erro: Sequências insuficientes para realizar o alinhamento comparativo.", file=sys.stderr)
        sys.exit(1)
        
    rec_a = records_a[0]
    rec_b = records_b[1] if args.file2 is None else records_b[0]
    
    print(f"aligning '{rec_a.id}' against '{rec_b.id}'...")
    aln = pairwise_align(rec_a, rec_b)
    
    if aln is None:
        print("Alignment failed.")
        return
        
    print("\n" + "="*50)
    print("GLOBAL alignment")
    print("-" * 50)
    print(f"score: {aln.score:.1f}")
    print(f"identity: {aln.identity_percent}%")
    print(f"matches: {aln.matches}/{aln.length} posições")
    print("="*50)
    
    # Gerando as barrinhas visuais clássicas do terminal
    match_line = "".join("|" if a == b and a != "-" else " " for a, b in zip(aln.aligned_a, aln.aligned_b))
    
    print(f"\n{aln.seq_a_id}:\n{aln.aligned_a}")
    print(f"{match_line}")
    print(f"{aln.aligned_b}\n{aln.seq_b_id}\n")

def main():
    parser = argparse.ArgumentParser(
        description="OmniSyn CLI — Ferramenta Computacional/Bioinformática via Terminal"
    )
    subparsers = parser.add_subparsers(dest="command", required=True, help="Commands")

    # Comando: analyze
    p_analyze = subparsers.add_parser("analyze", help="Análise molecular e composição/nucleotídeos")
    p_analyze.add_argument("input", type=str, help="Input sequence file")
    p_analyze.add_argument("-f", "--format", type=str, default="fasta", choices=["fasta", "genbank"], help="Input format (default: fasta)")
    p_analyze.set_defaults(func=handle_analyze)

    # Comando: orfs
    p_orfs = subparsers.add_parser("orfs", help="Predição e varredura/Open Reading Frames (ORFs)")
    p_orfs.add_argument("input", type=str, help="Input sequence file")
    p_orfs.add_argument("-f", "--format", type=str, default="fasta", choices=["fasta", "genbank"], help="Input format")
    p_orfs.add_argument("-m", "--min-len", type=int, default=90, help="Minimum ORF length")
    p_orfs.add_argument("-l", "--limit", type=int, default=10, help="Limite/exibição na tabela/resultados")
    p_orfs.set_defaults(func=handle_orfs)

    # Comando: align
    p_align = subparsers.add_parser("align", help="pairwise alignment")
    p_align.add_argument("file1", type=str, help="First input file")
    p_align.add_argument("file2", type=str, nargs="?", default=None, help="Second input file (optional)")
    p_align.add_argument("--format1", type=str, default="fasta", help="First file format")
    p_align.add_argument("--format2", type=str, default="fasta", help="Second file format")
    p_align.set_defaults(func=handle_align)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
  
        

  
