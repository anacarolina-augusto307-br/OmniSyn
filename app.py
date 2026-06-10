'''OmniSyn — interface web (Streamlit) para análise de sequências.'''

from __future__ import annotations

from pathlib import Path
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
import requests  # Importado explicitamente para a comunicação web básica
from Bio.SeqRecord import SeqRecord

from omnisyn_meta import (
    PROJECT_MEANING,
    PROJECT_NAME,
    PROJECT_SUBTITLE,
)

# 1. CONFIGURAÇÃO INICIAL DA PÁGINA (Comando obrigatório do Streamlit)
st.set_page_config(
    page_title=PROJECT_NAME,
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# 2. RENDERIZAÇÃO DO LOGOTIPO CENTRALIZADO
coluna_esquerda, coluna_central, coluna_direita = st.columns([1.1, 1, 1])

with coluna_central:
    st.image(
        "assets/logo_omnisyn_sf.png.png",
        width=240
    )

# 3. TÍTULOS DA INTERFACE (Via blocos simples de Markdown)
st.markdown(
    """
    <h1 style="text-align:center; margin-top:-35px; margin-bottom:0px; font-weight:700;">
        OmniSyn
    </h1>
    """,
    unsafe_allow_html=True
)

st.markdown(
    f"""
    <p style="text-align:center; font-size:18px; color:#a1a1aa;">
        {PROJECT_SUBTITLE}
    </p>
    """,
    unsafe_allow_html=True
)

# 4. IMPORTAÇÃO DAS FUNÇÕES LÓGICAS DO CORES (Mantidas conforme o laboratório estruturou)
from analyzer.core import (
    LIMITATIONS,
    analyze_sequence,
    find_orfs,
    pairwise_align,
    parse_sequences,
    protein_properties,
    sliding_gc,
)

# Dicionários e strings estáticas de exemplo
ORGANISM_SAMPLES = {
    "Bacillus subtilis (soil)": "sample_data/bacillus_subtilis_soil.fasta",
    "Vibrio cholerae (marine)": "sample_data/vibrio_cholerae_marine.fasta",
    "Haloferax volcanii (archaea)": "sample_data/haloferax_volcanii_archaea.fasta",
    "Short Demo": "sample_data/example.fasta",
}

SAMPLE_FASTA = """>sample_gene
ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTA
AAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTA
AAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTA
AAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTA
AAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACATG
"""

# Estilos CSS básicos da página
st.markdown(
    """
    <style>
    .main-header {
        font-size: 2.4rem;
        font-weight: 700;
        background: linear-gradient(90deg, #0f766e, #0891b2);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin-bottom: 0.25rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


# 5. FUNÇÕES INDIVIDUAIS PARA CONSTRUÇÃO DE CADA GRÁFICO (Lógica básica de listas)

def plot_nucleotide_bar(percents: dict[str, float], title: str) -> go.Figure:
    # Cria listas manuais a partir do dicionário para montar o DataFrame de forma visível
    lista_bases = list(percents.keys())
    lista_valores = list(percents.values())
    
    dados_grafico = {
        "Base": lista_bases,
        "Percent": lista_valores
    }
    df = pd.DataFrame(dados_grafico)
    
    colors = {"A": "#ef4444", "T": "#3b82f6", "G": "#22c55e", "C": "#eab308"}
    fig = px.bar(
        df,
        x="Base",
        y="Percent",
        title=title,
        color="Base",
        color_discrete_map=colors,
        text="Percent",
    )
    fig.update_traces(texttemplate="%{text:.1f}%", textposition="outside")
    fig.update_layout(showlegend=False, height=380, template="plotly_dark", margin=dict(t=50, b=40))
    return fig


def plot_gc_pie(gc: float, at: float) -> go.Figure:
    lista_rotulos = ["GC", "AT"]
    lista_proporcoes = [gc, at]
    
    fig = go.Figure(
        data=[
            go.Pie(
                labels=lista_rotulos,
                values=lista_proporcoes,
                hole=0.45,
                marker_colors=["#14b8a6", "#6366f1"],
            )
        ]
    )
    fig.update_layout(title="GC / AT Ratio", height=380, template="plotly_dark", margin=dict(t=50, b=40))
    return fig


def plot_sliding_gc(points: list[tuple[int, float]], seq_id: str) -> go.Figure:
    df = pd.DataFrame(points, columns=["Position", "GC %"])
    fig = px.line(
        df,
        x="Position",
        y="GC %",
        title="Sliding-Window GC Content — " + str(seq_id),
        markers=True,
    )
    fig.update_layout(height=400, template="plotly_dark", margin=dict(t=50, b=40))
    fig.update_traces(line_color="#0d9488")
    return fig


def plot_codon_heatmap(codon_usage: dict[str, int]) -> go.Figure:
    if len(codon_usage) == 0:
        return go.Figure()

    # Separação manual dos top 20 códons mais utilizados
    todos_os_codons = list(codon_usage.items())
    top_codons = todos_os_codons[:20]

    df = pd.DataFrame(top_codons, columns=["Codon", "Count"])

    fig = px.bar(
        df,
        x="Codon",
        y="Count",
        title="Top Codon Usage",
        color="Count",
        color_continuous_scale="Teal",
    )
    fig.update_layout(height=400, template="plotly_dark", showlegend=False, margin=dict(t=50, b=40))
    return fig


# 6. CONFIGURAÇÃO DOS MENUS NA BARRA LATERAL (SIDEBAR)

with st.sidebar:
    st.header("Data Source")
    
    input_format = st.selectbox(
        "Sequence Format",
        ["fasta", "plain", "genbank"],
        format_func=lambda x: {"fasta": "FASTA", "plain": "Plain text", "genbank": "GenBank"}[x],
    )
    
    uploaded = st.file_uploader(
        "Upload Sequence File",
        type=["fasta", "fa", "fna", "gb", "gbk", "txt"],
        help="Supported formats: FASTA, GenBank, or raw text.",
    )
    
    organism_sample = st.selectbox(
        "Reference Sample Genomic Data",
        ["(None)"] + list(ORGANISM_SAMPLES.keys()),
        help="Select verified reference fragments from genomic data repositories.",
    )
    
    use_sample = st.checkbox("Enable short sandbox demo", value=False)

    st.divider()
    st.markdown("**About " + str(PROJECT_NAME) + "**")
    st.caption(PROJECT_MEANING)
    
    with st.expander("Analytical Framework Constraints"):
        for note in LIMITATIONS:
            st.markdown("- " + str(note))
        st.caption("See `LIMITATIONS.md` for full details.")

    st.divider()
    
    # Executa o botão de comparação gerando o relatório local de forma explícita
    clicou_comparar = st.button("Execute Cross-Organism Comparison")
    if clicou_comparar == True:
        report_path = Path("sample_data/comparison_report.txt")
        try:
            from scripts.compare_organisms import compare
            conteudo_relatorio = compare()
            report_path.write_text(conteudo_relatorio, encoding="utf-8")
            st.success("Report successfully written: " + str(report_path))
        except Exception as exc:
            st.error("Failed to generate automated comparison matrix: " + str(exc))


# 7. CAPTURA E LEITURA DA SEQUÊNCIA (Lógica condicional linear sem atalhos)

raw = ""
fmt = input_format

if uploaded is not None:
    raw = uploaded.read().decode("utf-8", errors="replace")
    fmt = input_format
    nome_arquivo_minusculo = uploaded.name.lower()
    
    if nome_arquivo_minusculo.endswith(".gb") or nome_arquivo_minusculo.endswith(".gbk"):
        fmt = "genbank"
    elif nome_arquivo_minusculo.endswith(".fa") or nome_arquivo_minusculo.endswith(".fasta") or nome_arquivo_minusculo.endswith(".fna"):
        fmt = "fasta"
        
elif organism_sample != "(None)" and organism_sample in ORGANISM_SAMPLES:
    sample_path = Path(ORGANISM_SAMPLES[organism_sample])
    if sample_path.exists() == False:
        st.error("File system path error: Sample file could not be resolved.")
        st.stop()
    raw = sample_path.read_text(encoding="utf-8")
    fmt = "fasta"
    
elif use_sample == True:
    raw = SAMPLE_FASTA
    fmt = "fasta"
    
else:
    # Se nenhuma opção de arquivo foi marcada, captura o input digitado na caixa de texto
    raw = st.text_area(
        "Paste Sequence Data Workspace",
        height=200,
        placeholder=">seq_id_1\nATGAAACGC...\n\nOr paste raw text: ATGAAACGC...",
    )
    fmt = input_format

# Validação se a string de entrada está vazia
if raw == "" or raw.strip() == "":
    st.info("Awaiting structural or genetic input. Upload a file, choose a sample reference, or enter data in the workspace text field.")
    st.stop()


# 8. PROCESSAMENTO DO PARSER DA SEQUÊNCIA VIA BIOPYTHON

try:
    records = parse_sequences(raw, fmt=fmt)
except Exception as exc:
    st.error("Methodological parsing failure: " + str(exc))
    st.stop()

if len(records) == 0:
    st.warning("No validated records matched the selected structural parser settings.")
    st.stop()

st.success("Successfully loaded and cached " + str(len(records)) + " sequence entity record(s).")


# 9. SELETOR DE PÁGINAS DO INTERFACE
page = st.radio(
    "View Mode Selector",
    ["Single sequence", "Compare sequences"],
    horizontal=True,
    label_visibility="collapsed",
)


# ==================== MODO 1: SINGLE SEQUENCE ====================
if page == "Single sequence":
    
    # Se houver mais de uma sequência carregada, exibe o selectbox de escolha
    if len(records) > 1:
        indices_disponiveis = range(len(records))
        
        def formatar_opcao_seletor(opcao_i):
            sequencia_opcao = records[opcao_i]
            retorno_texto = str(sequencia_opcao.id) + " (" + str(len(sequencia_opcao.seq)) + " bp)"
            return retorno_texto
            
        idx = st.selectbox(
            "Select active sequence entity",
            indices_disponiveis,
            format_func=formatar_opcao_seletor,
        )
        record = records[idx]
    else:
        record = records[0]

    # Processamento dos cálculos estruturais no núcleo lógico
    result = analyze_sequence(record)

    st.subheader("Sequence: " + str(result.id))
    if result.description and result.description != result.id:
        st.caption(result.description)

    # Criação das colunas numéricas de metadados
    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Length", f"{result.length:,} bp")
    m2.metric("GC%", str(result.gc_percent) + "%")
    m3.metric("AT%", str(result.at_percent) + "%")
    m4.metric("Mol. Weight", f"{result.molecular_weight:,.0f} Da")
    
    if result.is_rna == True:
        m5.metric("Type", "RNA")
    else:
        m5.metric("Type", "DNA")

    # Renderização dos gráficos principais de nucleotídeos
    c1, c2 = st.columns(2)
    with c1:
        grafico_barras_bases = plot_nucleotide_bar(result.nucleotide_percent, "Nucleotide Composition")
        st.plotly_chart(grafico_barras_bases, use_container_width=True)
    with c2:
        grafico_pizza_gc = plot_gc_pie(result.gc_percent, result.at_percent)
        st.plotly_chart(grafico_pizza_gc, use_container_width=True)

    # Gráfico de Janela Deslizante de GC
    sequencia_limpa_gc = str(record.seq).upper().replace("U", "T")
    gc_points = sliding_gc(sequencia_limpa_gc)
    grafico_linha_deslizante = plot_sliding_gc(gc_points, result.id)
    st.plotly_chart(grafico_linha_deslizante, use_container_width=True)

    # Criação das abas de detalhamento técnico da sequência
    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        ["Sequences", "Translation", "ORFs", "Codons", "Protein Properties"]
    )

    with tab1:
        st.text_area("Original Sequence", str(record.seq), height=120, disabled=True)
        st.text_area("Reverse Complement", result.reverse_complement, height=120, disabled=True)

    with tab2:
        st.text_area("mRNA (T→U)", result.mrna[:2000], height=100, disabled=True)
        st.text_area("Protein (stop at first stop codon)", result.protein[:2000], height=100, disabled=True)
        
        with st.expander("Reading Frames (Preview)"):
            for frame, prot in result.reading_frames.items():
                st.markdown("**Frame " + str(frame) + "**")
                texto_proteina_exibicao = prot[:300]
                if len(prot) > 300:
                    texto_proteina_exibicao = texto_proteina_exibicao + "..."
                st.code(texto_proteina_exibicao, language=None)

    with tab3:
        min_len = st.slider("Minimum ORF length (bp)", 60, 300, 90, key="orf_slider_" + str(result.id))
        orfs = find_orfs(record, min_length=min_len)
        
        if len(orfs) == 0:
            st.info("No ORFs found at this threshold.")
        else:
            lista_orfs_tabela = []
            for o in orfs[:50]:
                item_dicionario = {
                    "Frame": o.frame,
                    "Strand": o.strand,
                    "Start": o.start,
                    "End": o.end,
                    "Length (bp)": o.length,
                    "Protein Length": len(o.protein),
                }
                lista_orfs_tabela.append(item_dicionario)
                
            orf_df = pd.DataFrame(lista_orfs_tabela)
            st.dataframe(orf_df, use_container_width=True, hide_index=True)
            
            def formatar_opcao_orf(indice_orf):
                orf_especifica = orfs[indice_orf]
                tamanho_orf_bp = orf_especifica.length
                sentido_cadeia = str(orf_especifica.strand) + str(orf_especifica.frame)
                return "ORF " + str(indice_orf + 1) + " — " + str(tamanho_orf_bp) + " bp (frame " + sentido_cadeia + ")"
                
            total_itens_caixa = min(10, len(orfs))
            sel = st.selectbox(
                "Inspect ORF",
                range(total_itens_caixa),
                format_func=formatar_opcao_orf,
                key="orf_sel_box_" + str(result.id),
            )
            
            orf_selecionada_exibir = orfs[sel]
            st.text_area("ORF DNA Sequence", orf_selecionada_exibir.sequence, height=80, disabled=True)
            st.text_area("ORF Protein Sequence", orf_selecionada_exibir.protein, height=80, disabled=True)

    with tab4:
        grafico_contagem_codons = plot_codon_heatmap(result.codon_usage)
        st.plotly_chart(grafico_contagem_codons, use_container_width=True)
        
        if len(result.codon_usage) > 0:
            lista_tabela_codons = []
            for k, v in result.codon_usage.items():
                lista_tabela_codons.append({"Codon": k, "Count": v})
                
            df_tabela_codons = pd.DataFrame(lista_tabela_codons)
            st.dataframe(df_tabela_codons, use_container_width=True, hide_index=True)

    with tab5:
        props = protein_properties(result.protein)
        if len(props) == 0:
            st.info("No protein sequence available for physicochemical analysis.")
        else:
            p1, p2, p3, p4 = st.columns(4)
            p1.metric("Protein Length", props.get("length", 0))
            p2.metric("Isoelectric Point", props.get("isoelectric_point", 0))
            p3.metric("GRAVY Index", props.get("gravy", 0))
            p4.metric("Instability Index", props.get("instability_index", 0))
            
            st.divider()
            
            ss = props.get("secondary_structure")
            if ss is not None:
                col_chart, col_desc = st.columns([3, 2])
                
                with col_chart:
                    labels = ["Helix", "Turn", "Sheet"]
                    values = [ss[0], ss[1], ss[2]]
                    
                    fig = go.Figure(data=[
                        go.Bar(
                            x=labels, 
                            y=values,
                            text=[f"{v*100:.1f}%" for v in values],
                            textposition='auto',
                            marker=dict(
                                color=["#14b8a6", "#f59e0b", "#6366f1"],
                                line=dict(color='#000', width=1)
                            )
                        )
                    ])
                    
                    fig.update_layout(
                        title="Predicted Secondary Structure Fraction",
                        template="plotly_dark",
                        height=380,
                        margin=dict(t=50, b=40),
                        yaxis=dict(title="Fraction (0 to 1.0)", range=[0, 1])
                    )
                    st.plotly_chart(fig, use_container_width=True)
                    
                    # 10. INTEGRAÇÃO DA IA DE PREDIÇÃO MOLECULAR 3D REAIS (ESMFold da Meta AI)
                    st.write("#### 🕶️ AI-Predicted 3D Macromolecular Structure (ESMFold)")
                    sequencia_aminoacidos_limpa = result.protein.replace("*", "").strip()
                    
                    if len(sequencia_aminoacidos_limpa) > 5:
                        # Executa o envio do pedido POST para a API pública do ESMFold
                        link_api_esmfold = "https://api.esmatlas.com/foldSequence/v1/pdb/"
                        try:
                            resposta_servidor = requests.post(link_api_esmfold, data=sequencia_aminoacidos_limpa, timeout=30)
                            
                            # Se a API retornou o arquivo PDB predito com sucesso
                            if resposta_servidor.status_code == 200:
                                dados_coordenadas_pdb = resposta_servidor.text
                                
                                # Carrega as bibliotecas interativas de modelagem de forma explícita
                                import py3Dmol
                                from stmol import showmol
                                
                                visualizador_3d = py3Dmol.view(width=400, height=350)
                                visualizador_3d.addModel(dados_coordenadas_pdb, "pdb")
                                visualizador_3d.setStyle({'cartoon': {'color': 'spectrum'}})
                                visualizador_3d.addSurface(py3Dmol.SAS, {'opacity': 0.2, 'color': '#0d9488'})
                                visualizador_3d.zoomTo()
                                
                                showmol(visualizador_3d, height=350, width=500)
                            else:
                                st.error("⚠️ ESMFold API returned an error status code: " + str(resposta_servidor.status_code))
                        except Exception as erro_conexao:
                            st.error("⚠️ Connection failure with ESMFold AI API: " + str(erro_conexao))
                    else:
                        st.warning("Protein sequence too short for high-fidelity 3D structural prediction.")

                with col_desc:
                    st.markdown("### 🧬 Structural Distribution Overview")
                    st.caption("Predicted secondary structure fractions for the analyzed protein module.")
                    
                    h_pct = ss[0] * 100
                    t_pct = ss[1] * 100
                    s_pct = ss[2] * 100
                    
                    st.markdown(f"""
                    * **Helix (Alpha-Helix):** Represents **{h_pct:.1f}%** of the conformation. Stabilized by regular intrachain hydrogen bonds, highly prevalent in transmembrane domains and globular structural folds.
                    * **Sheet (Beta-Sheet):** Represents **{s_pct:.1f}%** of the structural composition. Formed by lateral hydrogen bonds between extended peptide chains, providing central stability.
                    * **Turn (Beta-Turn):** Represents **{t_pct:.1f}%** of the backbone. These short loop regions invert the three-dimensional direction of the polypeptide chain, crucial for folding.
                    """)
                    
                    if h_pct > t_pct and h_pct > s_pct:
                        st.info("💡 **Structure Insight:** The high alpha-helix content suggests a structural profile typical of globular receptors or transmembrane channels.")
                    elif s_pct > h_pct:
                        st.info("💡 **Structure Insight:** The predominance of beta-sheets suggests a structural fold often found in stable fibrous proteins or core beta-barrels.")


# ==================== MODO 2: COMPARE SEQUENCES ====================
else:
    if len(records) < 2:
        st.warning("Comparative algorithms require at least two distinct sequence records.")
        st.stop()
        
    i, j = st.columns(2)
    with i:
        idx_a = st.selectbox("Sequence Dataset A", range(len(records)), format_func=lambda x: records[x].id)
    with j:
        idx_b = st.selectbox(
            "Sequence Dataset B",
            range(len(records)),
            index=min(1, len(records) - 1),
            format_func=lambda x: records[x].id,
        )
        
    if idx_a == idx_b:
        st.warning("Comparative matrix requires distinct datasets. Select different entities.")
        st.stop()

    col_a, col_b = st.columns(2)
    with col_a:
        ra = analyze_sequence(records[idx_a])
        grafico_comparativo_a = plot_nucleotide_bar(ra.nucleotide_percent, str(ra.id) + " Composition Profile")
        st.plotly_chart(grafico_comparativo_a, use_container_width=True)
    with col_b:
        rb = analyze_sequence(records[idx_b])
        grafico_comparativo_b = plot_nucleotide_bar(rb.nucleotide_percent, str(rb.id) + " Composition Profile")
        st.plotly_chart(grafico_comparativo_b, use_container_width=True)

    st.subheader("Global Pairwise Sequence Alignment")
    aln = pairwise_align(records[idx_a], records[idx_b])
    
    if aln is None:
        st.warning("Could not align the selected sequences.")
    else:
        col_res1, col_res2, col_res3, col_res4 = st.columns(4)
        col_res1.metric("Alignment Score", f"{aln.score:.1f}")
        col_res2.metric("Identity", str(aln.identity_percent) + "%")
        col_res3.metric("Matches", aln.matches)
        col_res4.metric("Length", aln.length)

        match_line = []
        for char_a, char_b in zip(aln.aligned_a, aln.aligned_b):
            if char_a == char_b and char_a != "-":
                match_line.append("|")
            else:
                match_line.append(" ")
        linha_barras_casamento = "".join(match_line)

        texto_alinhamento_final = (
            str(aln.seq_a_id) + "\n" + 
            str(aln.aligned_a) + "\n" + 
            linha_barras_casamento + "\n" + 
            str(aln.aligned_b) + "\n" + 
            str(aln.seq_b_id)
        )
        st.code(texto_alinhamento_final, language=None)

    st.subheader("Dynamic GC Divergence Metric")
    fig = go.Figure()
    
    # Executa o loop explicitamente para criar as duas linhas sobrepostas no gráfico
    for rec in (records[idx_a], records[idx_b]):
        sequencia_rec_limpa = str(rec.seq).upper().replace("U", "T")
        pts = sliding_gc(sequencia_rec_limpa)
        
        df_linha_gc = pd.DataFrame(pts, columns=["Position", "GC %"])
        fig.add_trace(
            go.Scatter(x=df_linha_gc["Position"], y=df_linha_gc["GC %"], mode="lines", name=rec.id)
        )
        
    fig.update_layout(height=400, template="plotly_dark", title="Comparative Sliding-Window GC Overlay")
    st.plotly_chart(fig, use_container_width=True)
