"""
MicrobiomeDash — Introduction / landing page.

Static informational page (no callbacks). Provides an overview of the
application, a quick-start workflow guide, setup instructions, and a
summary of supported analysis methods.
"""
import dash_bootstrap_components as dbc
from dash import html, dcc


# ── Helper builders ──────────────────────────────────────────────────────────

def _tool_card(title, subtitle, items):
    """A card describing one of the three core tools."""
    return dbc.Col(
        dbc.Card(
            dbc.CardBody([
                html.H5(title, className="card-title"),
                html.P(subtitle, className="card-subtitle text-muted mb-3"),
                html.Ul([html.Li(i) for i in items], className="mb-0"),
            ]),
            className="h-100",
        ),
        md=4,
        className="mb-3",
    )


def _step(number, title, description, link_href, link_label):
    """A numbered workflow step with a page link."""
    return dbc.ListGroupItem([
        html.Div([
            dbc.Badge(str(number), color="primary", className="me-3",
                      style={"fontSize": "1rem", "padding": "0.5rem 0.75rem"}),
            html.Div([
                html.Strong(title),
                html.Br(),
                html.Span(description, className="text-muted"),
                html.Span(" "),
                dcc.Link(link_label, href=link_href, className="ms-1"),
            ]),
        ], className="d-flex align-items-start"),
    ])


# ── Methods summary table ────────────────────────────────────────────────────

methods_table = dbc.Table(
    [
        html.Thead(html.Tr([
            html.Th("Category"),
            html.Th("Methods"),
        ])),
        html.Tbody([
            html.Tr([
                html.Td("Alpha Diversity"),
                html.Td("Shannon, Simpson, Observed OTUs, Chao1, Pielou Evenness"),
            ]),
            html.Tr([
                html.Td("Beta Diversity"),
                html.Td("Bray-Curtis, Jaccard; PCoA, NMDS; PERMANOVA"),
            ]),
            html.Tr([
                html.Td("Differential Abundance"),
                html.Td("ALDEx2, ANCOM-BC, DESeq2, LinDA, MaAsLin2"),
            ]),
            html.Tr([
                html.Td("Taxonomy"),
                html.Td("SILVA 138.2 classifier; heatmap at all ranks (Phylum to Species)"),
            ]),
            html.Tr([
                html.Td("Pathway Analysis"),
                html.Td("PICRUSt2 functional prediction; MetaCyc / KEGG pathway comparison"),
            ]),
        ]),
    ],
    bordered=True,
    hover=True,
    responsive=True,
    className="mt-3",
)


# ── Page layout ──────────────────────────────────────────────────────────────

layout = dbc.Container([

    # Header
    html.Div([
        html.H2("16S Pipeline", className="mb-2"),
        html.P(
            "An end-to-end 16S rRNA amplicon analysis platform — from raw FASTQ "
            "files to publication-ready results.",
            className="lead text-muted",
        ),
    ], className="mt-4 mb-4"),

    html.Hr(),

    # ── About ────────────────────────────────────────────────────────────────
    html.H4("About", className="mt-4 mb-3"),
    html.P(
        "16S Analyzer integrates quality control, amplicon sequence variant (ASV) "
        "inference, taxonomic classification, data curation, and statistical "
        "analysis into a single browser-based interface. The platform is built "
        "around three core tools:"
    ),

    dbc.Alert([
        html.Strong("Supported input: "),
        "Illumina paired-end or single-end amplicon FASTQ files targeting specific "
        "16S variable regions (V1-V2, V3-V4, V4, V4-V5, V5-V6). "
        "Full-length 16S long reads (PacBio HiFi, Nanopore) also supported (auto-detected).",
    ], color="info", className="mb-3"),

    dbc.Row([
        _tool_card(
            "1. Pipeline Engine",
            "FASTQ to Dataset",
            [
                "FastQC quality reports",
                "Cutadapt primer trimming",
                "DADA2 ASV inference & chimera removal",
                "SILVA 138.2 taxonomic classification",
                "Phylogenetic tree construction",
                "PICRUSt2 functional prediction",
            ],
        ),
        _tool_card(
            "2. Data Management",
            "Curate & Transform",
            [
                "BIOM file browser & inspector",
                "Rare ASV removal (prevalence / abundance)",
                "Rarefaction & sample filtering",
                "Dataset combination (merge studies)",
                "V-region extraction",
                "BIOM & MOTHUR format conversion",
            ],
        ),
        _tool_card(
            "3. Analysis & Visualization",
            "Statistics & Figures",
            [
                "Alpha diversity (5 metrics, group comparisons)",
                "Beta diversity (PCoA, NMDS, PERMANOVA)",
                "Taxonomy heatmap plots",
                "Differential abundance (5 methods, volcano plots)",
                "PICRUSt2 pathway comparison",
                "KEGG pathway map viewer",
            ],
        ),
    ], className="mt-3"),

    html.Hr(),

    # ── Quick Start ──────────────────────────────────────────────────────────
    html.H4("Quick Start", className="mt-4 mb-3"),
    html.P("Follow these steps to go from raw sequences to analysis results:"),

    dbc.ListGroup([
        _step(1, "Register FASTQ files",
              "Upload or link your paired-end FASTQ files.",
              "/files", "File Registration"),
        _step(2, "Run the DADA2 pipeline",
              "Set primers, trimming parameters, and launch the pipeline.",
              "/pipeline", "Pipeline"),
        _step(3, "Curate your data",
              "Remove rare ASVs, rarefy, filter samples, or combine datasets.",
              "/biom-browser", "Data Management"),
        _step(4, "Analyze & visualize",
              "Explore diversity, taxonomy, and differential abundance.",
              "/alpha", "Analysis"),
        _step(5, "Pathway analysis",
              "Run PICRUSt2 functional prediction and compare pathways.",
              "/picrust2", "Pathways"),
    ], flush=True, className="mb-4"),

    html.Hr(),

    # ── Setup ────────────────────────────────────────────────────────────────
    html.H4("Setup", className="mt-4 mb-3"),

    html.H6("Prerequisites"),
    html.Ul([
        html.Li("Linux or WSL2 on Windows"),
        html.Li(["conda or mamba (", html.A("Miniforge",
                 href="https://github.com/conda-forge/miniforge",
                 target="_blank"), " recommended)"]),
    ]),

    html.H6("Installation", className="mt-3"),
    dbc.Card(
        dbc.CardBody(
            html.Pre(
                "git clone https://github.com/tatsu1207/16S-Pipeline.git && cd 16S-Pipeline\n"
                "bash setup_ubuntu.sh      # creates conda envs, installs R packages, downloads SILVA\n"
                "conda activate microbiome_16S\n"
                "./run.sh                   # starts the server on port 7000 + UID",
                className="mb-0",
                style={"whiteSpace": "pre-wrap"},
            ),
        ),
        className="bg-dark border-secondary",
    ),

    html.P(
        "The setup script creates four conda environments: microbiome_16S (web app), "
        "dada2_16S (R + DADA2), analysis_16S (R + DA tools), and picrust2_16S "
        "(PICRUSt2). It also downloads the SILVA 138.1 reference database.",
        className="text-muted mt-2",
    ),

    html.Hr(),

    # ── Supported Methods ────────────────────────────────────────────────────
    html.H4("Supported Analysis Methods", className="mt-4 mb-2"),
    methods_table,

], fluid=True, className="pb-5")
