# Dados em `sample_data`

## Organismos reais (NCBI GenBank)

Fragmentos de **25 kb** do início de cada cromossomo/região de referência. Baixados via NCBI E-utilities (`efetch`), maio/2026.

| Arquivo | Organismo | Contexto | Acesso NCBI |
|---------|-----------|----------|-------------|
| `bacillus_subtilis_soil.fasta` | *Bacillus subtilis* 168 | Bactéria de solo, modelo laboratorial | [NC_000964.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000964.3) (nt 1–25000) |
| `vibrio_cholerae_marine.fasta` | *Vibrio cholerae* O1 biovar El Tor | Bactéria marinha / estuarina | [NC_002505.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_002505.1) (nt 1–25000) |
| `haloferax_volcanii_archaea.fasta` | *Haloferax volcanii* DS2 | Archaea halófila (ambientes salinos) | [NC_013967.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_013967.1) (nt 1–25000) |

Fontes alternativas para amostras metagenômicas: [EBI Metagenomics](https://www.ebi.ac.uk/metagenomics), [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/).

## Demo curta

- `example.fasta` — fragmentos artificiais/pequenos para testar a interface web rapidamente (não são genomas completos).

## Comparação

```bash
python scripts/compare_organisms.py
```

Gera `sample_data/comparison_report.txt` com GC%, contagem de ORFs e uso de códon (top 10) por organismo.
