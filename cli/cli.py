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
        print(f"Erro: file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)

    raw_content = Path(args.input).read_text(encoding="utf-8")
    records = parse_sequences(raw_content, fmt=args.format)

    if not records:
      print("No valid sequence found with the chosen format.")

  for record in records:
      res = analyze_sequences(records):
      print("\n" + "="*50)
      print(f"SEQUENCE ID: {res.id}")
      print(f"Description: {res.description}")
      print(f"-" * 50)
      print(f"Length: {res.length:,} bp")
      print(f"GC percent: {res.gc_percent}%  | AT percent: {res.at_percent}%")
      print(f"Molecular weight: {res.molecular_weight:,.2f} Da")
      print(f"Tipo: {'RNA' if res.is_rna else 'DNA'}")
      print(f"\n📊 Nucleotide Composition:")
      for base, count in res.base_nucleotide_counts_items():
        pct = res.nucleotide_percent[base]
        print(f"  [{base}]: {count:<8} ({pct:.2f}%)")
      print("="*50)
    
#Para buscar os orfs
def handle_orfs(args):
  if not Path(args.input).exists():
    print(f"Erro: file '{args.input} not found.", file=sys.stderr)
    sys.exit(1)

  raw_content = Path(args.input).read_text(encoding="utf-8")
  records = parse_sequences(raw_content, fmt=args.format)

  for record in records:
    print(f"\n Searching for ORFs in '{record.id}' (Minimum: {min_args_len} bp)... ")
    orfs = find_orfs(record, min_length=args.min_len)

      if not orfs:
        print(f"No ORF found for the denifed criteria.")
        continue

      print(f"  Found {len(orfs)} ORF(s). Displaying the largest ones:")
      print(f"{'Strand':<5} | {'Frame':<5} | {'Star':<8} | {'End':<8} | {'Lenght (bp)':<12}")
      print("-" * 48)
      for idx, o in enumerate(orfs[:args.limit]):
          print(f"{o.strand:<5} | {o.frame:<5} | {o.start:<8} | {o.end:<8} | {o.length:<12}")

#Referente ao alinhamento de duas ou mais sequências contidas em um ou dois arquivos
def handle_align(args):
  if not Path(args.input).exists():
    print(f"Erro: file '{args.input}' not fond.", file=sys.stderr)
    sys.exit(1)

raw_content = Path(args.input).read_text(enconding="utf-8")
records = parse_sequences(raw_content

    
  
        

  
