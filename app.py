'''OmniSyn — interface web (Streamlit) para análise de sequências.'''

from __future__ import annotations

from pathlib import Path

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from Bio.SeqRecord import SeqRecord

from omnisyn_meta import (
    PROJECT_MEANING,
    PROJECT_NAME,
    PROJECT_SUBTITLE,
)

# Configuração do código para vizualização molecular 3D -> teste 4
    # o py3Dmol e o stmol foram adicionados no arquivo requirements desse github
    # O teste 12 apresentou menor inconsistência no kaggle, mas devido a erros de escrita falhou na importação. LEMBRAR QUE PYTHON É EXTREMAMENTE SENSÍVEL
    # O teste 19 não apareceu mais erros de escrita até então - eu acho -, mas tive problema com o efeito cascata. ATENÇÃO AOS ESPAÇOS E PARÁGRAFOS DO PYTHON
    # Tentativa 27 com sucesso
try:
    import py3Dmol
    from stmol import showmol
    STMOL_AVAILABLE = True

except Exception as e:
    STMOL_AVAILABLE = False
    st.error(f"3D Import Error: {e}")

# Configuração da página (deve ser o primeiro comando Streamlit)
    # Substituir em breve o símbolo de DNA da página para o da omnisyn
    #

st.set_page_config(
    page_title=PROJECT_NAME,
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Logo da omnisyn
    # Melhorar a qualidade da logo em breve, atualmente está em 400x400 - feita pelo remove.bg gratuito
    # Atualmente o "right = st.columns([1.1, 1, 1])" funciona pq a aba está evidente na tela, se ela fechar a logo fica torta
left, center_col, right = st.columns([1.1, 1, 1])

with center_col:
    st.image(
        "assets/logo_omnisyn_sf.png.png", #Sim, o nome da logo está estranho assim porque me confundi na hora de salvar o arquivo e o png duplicou
        width=240
    )

# Título da Omnisyn
    # Para modificar a frase do título - acho improvável, mas enfim - trocar no arquivo omnisyn_meta.py desse github
st.markdown(
    """
    <h1 style="
        text-align:center;
        margin-top:-35px;
        margin-bottom:0px;
        font-weight:700;
    ">
        OmniSyn
    </h1>
    """,
    unsafe_allow_html=True
)

# Subtítulo da Omnisyn
    # Para modificar a frase do subtítulo trocar no arquivo omnisyn_meta.py desse github
st.markdown(
    f"""
    <p style="text-align:center; font-size:18px; color:#a1a1aa;">
        {PROJECT_SUBTITLE}
    </p>
    """,
    unsafe_allow_html=True
)

# Importações das funções lógicas do núcleo do analisador
from analyzer.core import (
    LIMITATIONS,
    analyze_sequence,
    find_orfs,
    pairwise_align,
    parse_sequences,
    protein_properties,
    sliding_gc,
)

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


def apply_styles() -> None:
    #Adicionando estilo custumizado na header do projeto + CSS

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


def plot_nucleotide_bar(percents: dict[str, float], title: str) -> go.Figure:
    df = pd.DataFrame({"Base": list(percents.keys()), "Percent": list(percents.values())})
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
    fig = go.Figure(
        data=[
            go.Pie(
                labels=["GC", "AT"],
                values=[gc, at],
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
        title=f"Sliding-Window GC Content — {seq_id}",
        markers=True,
    )
    fig.update_layout(height=400, template="plotly_dark", margin=dict(t=50, b=40))
    fig.update_traces(line_color="#0d9488")
    return fig

#Receptor do dicionário de códons, seguindo para a sua verificação
def plot_codon_heatmap(codon_usage: dict[str, int]) -> go.Figure:
    if len(codon_usage) == 0:
        return go.Figure()

    top_codons = list(codon_usage.items())[:20]

    df = pd.DataFrame(top_codons, columns=["Codon", "Count"])

    fig = px.bar(
        df,
        x="Codon",
        y="Count",
        title="Top Codon Usage",
        color="Count",#Colore as barras e deixa o fundo azul escuro
        color_continuous_scale="Teal",
    )

    fig.update_layout(
        height=400,
        template="plotly_dark",
        showlegend=False,
        margin=dict(t=50, b=40),
    )

    return fig

def render_alignment(aln) -> None:
    if aln is None:
        st.warning("Could not align the selected sequences.")
        return

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Alignment Score", f"{aln.score:.1f}")
    col2.metric("Identity", f"{aln.identity_percent}%")
    col3.metric("Matches", aln.matches)
    col4.metric("Length", aln.length)

    match_line = []
    for a, b in zip(aln.aligned_a, aln.aligned_b):
        match_line.append("|" if a == b and a != "-" else " ")

    st.code(
        f"{aln.seq_a_id}\n{aln.aligned_a}\n{''.join(match_line)}\n{aln.aligned_b}\n{aln.seq_b_id}",
        language=None,
    )


def build_3d_structure_preview(style_type: str = "helix") -> py3Dmol.view:
    """Gera coordenadas moleculares geométricas reais em 3D para renderização nativa."""
    view = py3Dmol.view(width=400, height=350)
    
    # Mapeamento de coordenadas PDB fixas para demonstração estrutural
    pdb_samples = {
        "sheet": (
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\n"
            "ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C\n"
            "ATOM      3  C   ALA A   1       2.010   1.420   0.000  1.00 20.00           C\n"
            "ATOM      4  O   ALA A   1       1.230   2.370   0.000  1.00 20.00           O\n"
            "ATOM      5  N   ALA A   2       3.330   1.570   0.000  1.00 20.00           N\n"
            "ATOM      6  CA  ALA A   2       4.010   2.870   0.000  1.00 20.00           C\n"
            "ATOM      7  C   ALA A   2       5.510   2.740   0.000  1.00 20.00           C\n"
            "ATOM      8  O   ALA A   2       6.230   3.730   0.000  1.00 20.00           O"
        ),
        "helix": (
            "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   GLY A   1       2.009   1.427   0.000  1.00  0.00           C\n"
            "ATOM      4  O   GLY A   1       1.221   2.375   0.000  1.00  0.00           O\n"
            "ATOM      5  N   GLY A   2       3.328   1.576   0.000  1.00  0.00           N\n"
            "ATOM      6  CA  GLY A   2       3.987   2.879   0.000  1.00  0.00           C\n"
            "ATOM      7  C   GLY A   2       5.467   2.693   0.000  1.00  0.00           C\n"
            "ATOM      8  O   GLY A   2       6.230   3.663   0.000  1.00  0.00           O"
        )
    }

    # Define o dado PDB com base no tipo (fallback para helix caso não seja sheet)
    pdb_data = pdb_samples.get(style_type, pdb_samples["helix"])
    color = "#6366f1" if style_type == "sheet" else "#14b8a6"

    # Renderização e configuração do modelo 3D
    view.addModel(pdb_data, "pdb")
    view.setStyle({'cartoon': {'color': color}})
    view.addSurface(py3Dmol.SAS, {'opacity': 0.4, 'color': '#0d9488'})
    view.zoomTo()
    
    return view

def render_sequence_analysis(record: SeqRecord) -> None:
    result = analyze_sequence(record)

    st.subheader(f"Sequence: `{result.id}`")
    if result.description and result.description != result.id:
        st.caption(result.description)

    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("Length", f"{result.length:,} bp")
    m2.metric("GC%", f"{result.gc_percent}%")
    m3.metric("AT%", f"{result.at_percent}%")
    m4.metric("Mol. Weight", f"{result.molecular_weight:,.0f} Da")
    m5.metric("Type", "RNA" if result.is_rna else "DNA")

    c1, c2 = st.columns(2)
    with c1:
        st.plotly_chart(
            plot_nucleotide_bar(result.nucleotide_percent, "Nucleotide Composition"),
            use_container_width=True,
        )
    with c2:
        st.plotly_chart(plot_gc_pie(result.gc_percent, result.at_percent), use_container_width=True)

    gc_points = sliding_gc(str(record.seq).upper().replace("U", "T"))
    st.plotly_chart(plot_sliding_gc(gc_points, result.id), use_container_width=True)

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
                st.markdown(f"**Frame {frame}**")
                st.code(prot[:300] + ("..." if len(prot) > 300 else ""), language=None)

    with tab3:
        min_len = st.slider("Minimum ORF length (bp)", 60, 300, 90, key=f"orf_{result.id}")
        orfs = find_orfs(record, min_length=min_len)
        if not orfs:
            st.info("No ORFs found at this threshold.")
        else:
            orf_df = pd.DataFrame(
                [
                    {
                        "Frame": o.frame,
                        "Strand": o.strand,
                        "Start": o.start,
                        "End": o.end,
                        "Length (bp)": o.length,
                        "Protein Length": len(o.protein),
                    }
                    for o in orfs[:50]
                ]
            )
            st.dataframe(orf_df, use_container_width=True, hide_index=True)
            sel = st.selectbox(
                "Inspect ORF",
                range(min(10, len(orfs))),
                format_func=lambda i: f"ORF {i+1} — {orfs[i].length} bp (frame {orfs[i].strand}{orfs[i].frame})",
                key=f"orf_sel_{result.id}",
            )
            st.text_area("ORF DNA Sequence", orfs[sel].sequence, height=80, disabled=True)
            st.text_area("ORF Protein Sequence", orfs[sel].protein, height=80, disabled=True)

    with tab4:
        st.plotly_chart(plot_codon_heatmap(result.codon_usage), use_container_width=True)
        if result.codon_usage:
            st.dataframe(
                pd.DataFrame(
                    [{"Codon": k, "Count": v} for k, v in result.codon_usage.items()]
                ),
                use_container_width=True,
                hide_index=True,
            )

    with tab5:
        props = protein_properties(result.protein)
        if not props:
            st.info("No protein sequence available for physicochemical analysis.")
        else:
            p1, p2, p3, p4 = st.columns(4)
            p1.metric("Protein Length", props.get("length", 0))
            p2.metric("Isoelectric Point", props.get("isoelectric_point", 0))
            p3.metric("GRAVY Index", props.get("gravy", 0))
            p4.metric("Instability Index", props.get("instability_index", 0))
            
            st.divider()
            
            ss = props.get("secondary_structure")  # Returns tuple: (helix, turn, sheet)
            if ss:
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
                    
                    # Visualização Tridimensional Geométrica aproximada
                    if STMOL_AVAILABLE:
                        st.write("#### 🕶️ Interactive 3D Structural View")
                        style = "sheet" if ss[2] > ss[0] else "helix"
                        mol_view = build_3d_structure_preview(style)
                        showmol(mol_view, height=350, width=500)
                    else:
                        st.info("To enable 3D molecular structures, run: pip install stmol py3Dmol")

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


def main() -> None:
    apply_styles()
    
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
        st.markdown(f"**About {PROJECT_NAME}**")
        st.caption(PROJECT_MEANING)
        with st.expander("Analytical Framework Constraints"):
            for note in LIMITATIONS:
                st.markdown(f"- {note}")
            st.caption("See `LIMITATIONS.md` for full details.")

        st.divider()
        if st.button("Execute Cross-Organism Comparison"):
            report_path = Path("sample_data/comparison_report.txt")
            try:
                from scripts.compare_organisms import compare
                report_path.write_text(compare(), encoding="utf-8")
                st.success(f"Report successfully written: `{report_path}`")
            except Exception as exc:
                st.error(f"Failed to generate automated comparison matrix: {exc}")

    # Fluxo condicional de resolução de entrada de dados de sequência -> sequência - fluxo - resolução <-
    if uploaded:
        raw = uploaded.read().decode("utf-8", errors="replace")
        fmt = input_format
        if uploaded.name.lower().endswith((".gb", ".gbk")):
            fmt = "genbank"
        elif uploaded.name.lower().endswith((".fa", ".fasta", ".fna")):
            fmt = "fasta"
    elif organism_sample != "(nenhuma)" and organism_sample in ORGANISM_SAMPLES:
        sample_path = Path(ORGANISM_SAMPLES[organism_sample])
        if not sample_path.exists():
            st.error(f"File system path error: {sample_path} could not be resolved.")
            st.stop()
        raw = sample_path.read_text(encoding="utf-8")
        fmt = "fasta"
    elif use_sample:
        raw = SAMPLE_FASTA
        fmt = "fasta"
    else:
        raw = st.text_area(
            "Paste Sequence Data Workspace",
            height=200,
            placeholder=">seq_1\nATGAAACGC...\n\nOr paste raw text: ATGAAACGC...",
        )
        fmt = input_format

    if not raw or not raw.strip():
        st.info("Awaiting structural or genetic input. Upload a file, choose a sample reference, or enter data in the workspace text field.")
        st.stop()

    try:
        records = parse_sequences(raw, fmt=fmt)
    except Exception as exc:
        st.error(f"Methodological parsing failure: {exc}")
        st.stop()

    if not records:
        st.warning("No validated records matched the selected structural parser settings.")
        st.stop()

    st.success(f"Successfully loaded and cached **{len(records)}** sequence entity record(s).")

    page = st.radio(
        "View Mode Selector",
        ["Single sequence", "Compare sequences"],
        horizontal=True,
        label_visibility="collapsed",
    )

    if page == "Single sequence":
        if len(records) > 1:
            idx = st.selectbox(
                "Select active sequence entity",
                range(len(records)),
                format_func=lambda i: f"{records[i].id} ({len(records[i].seq)} bp)",
            )
            record = records[idx]
        else:
            record = records[0]
        render_sequence_analysis(record)
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
            st.plotly_chart(plot_nucleotide_bar(ra.nucleotide_percent, f"{ra.id} Composition Profile"), use_container_width=True)
        with col_b:
            rb = analyze_sequence(records[idx_b])
            st.plotly_chart(plot_nucleotide_bar(rb.nucleotide_percent, f"{rb.id} Composition Profile"), use_container_width=True)

        st.subheader("Global Pairwise Sequence Alignment")
        aln = pairwise_align(records[idx_a], records[idx_b])
        render_alignment(aln)

        st.subheader("Dynamic GC Divergence Metric")
        fig = go.Figure()
        for rec in (records[idx_a], records[idx_b]):
            pts = sliding_gc(str(rec.seq).upper().replace("U", "T"))
            df = pd.DataFrame(pts, columns=["Position", "GC %"])
            fig.add_trace(
                go.Scatter(x=df["Position"], y=df["GC %"], mode="lines", name=rec.id)
            )
        fig.update_layout(height=400, template="plotly_dark", title="Comparative Sliding-Window GC Overlay")
        st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
