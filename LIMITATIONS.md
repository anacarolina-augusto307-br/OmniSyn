# Limitações conhecidas (OmniSyn)

Este projeto é um exercício de **microbiologia computacional em nível de graduação**, não um pipeline de produção certificado.

## Análise de sequências (`analyzer/`)

- **ORF prediction is simplistic and may not fully represent fragmented metagenomic assemblies.** O buscador usa start/stop canônicos em janelas fixas; não há modelagem de genes partidos, shift de frame nem anotação estrutural.
- **Codon usage** é contagem bruta em leitura +1; não corrige viés de composição, codon adaptation index (CAI) nem expressão.
- **GC content** em janela deslizante usa passo fixo; não substitui perfis de qualidade de sequenciamento.
- **Alinhamento global** (Needleman–Wunsch) é adequado para trechos parecidos; distorce comparações entre genomas distantes ou regiões rearranjadas.
- **Tradução** assume código genético padrão (bacterial/archaeal table 11 na maioria dos casos, mas não é validada por organismo).

## Pipeline FASTQ (`omnisyn_pipeline.py`)

- Parâmetros de trimming derivados do FastQC são **heurísticos** (WARN/FAIL → flags fastp), não substituem revisão manual.
- Requer `fastqc` e `fastp` instalados fora do Python; versões diferentes podem alterar resultados.
- Mock data em `data/raw/` não representa complexidade real de metagenomas.

## Dados em `sample_data/`

- Fragmentos de **25 kb** de referências GenBank; não são genomas completos nem MAGs curados.
- Comparações entre solo, mar e archaea são **ilustrativas**, não conclusões filogenéticas ou ecológicas.

## Interface web

- Streamlit roda análises na memória; arquivos muito grandes podem travar a sessão.
- Não há autenticação, fila de jobs nem rastreabilidade tipo LIMS.
