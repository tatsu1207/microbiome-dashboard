"""
16S Pipeline — Introduction / landing page.

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
                html.Td("Shannon, Simpson, Observed ASVs, Chao1, ACE, Pielou Evenness"),
            ]),
            html.Tr([
                html.Td("Beta Diversity"),
                html.Td("Bray-Curtis, Jaccard; PCoA, NMDS; PERMANOVA"),
            ]),
            html.Tr([
                html.Td("Differential Abundance"),
                html.Td("ALDEx2, ANCOM-BC2, DESeq2, LinDA, MaAsLin2 (multi-method consensus)"),
            ]),
            html.Tr([
                html.Td("Taxonomy"),
                html.Td("SILVA v138.1 classifier; stacked bar plots at all ranks (Phylum to Species)"),
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
            "An end-to-end 16S rRNA amplicon analysis platform \u2014 from raw FASTQ "
            "files to publication-ready results.",
            className="lead text-muted",
        ),
    ], className="mt-4 mb-4"),

    html.Hr(),

    # ── About ────────────────────────────────────────────────────────────────
    html.H4("About", className="mt-4 mb-3"),
    html.P(
        "16S Pipeline integrates quality control, amplicon sequence variant (ASV) "
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
                "Auto-detection of sequencing type, region & platform",
                "Automated quality parameter selection",
                "DADA2 ASV inference & chimera removal",
                "SILVA v138.1 taxonomic classification",
                "Phylogenetic tree construction (MAFFT + FastTree)",
                "PICRUSt2 functional prediction (optional)",
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
                "SRA download (fetch public datasets by accession)",
                "SRA submission helper (generate metadata spreadsheet)",
            ],
        ),
        _tool_card(
            "3. Analysis & Visualization",
            "Statistics & Figures",
            [
                "Alpha diversity (6 metrics, group comparisons)",
                "Beta diversity (PCoA, NMDS, PERMANOVA)",
                "Taxonomy stacked bar plots (all ranks)",
                "Differential abundance (5 methods + consensus)",
                "PICRUSt2 pathway comparison",
                "KEGG pathway map viewer",
                "Analysis report PDF generation",
            ],
        ),
    ], className="mt-3"),

    html.Hr(),

    # ── Quick Start ──────────────────────────────────────────────────────────
    html.H4("Quick Start", className="mt-4 mb-3"),
    html.P("Follow these steps to go from raw sequences to analysis results:"),

    dbc.ListGroup([
        _step(1, "Upload FASTQ files",
              "Upload your FASTQ files or download from NCBI SRA.",
              "/files", "File Manager"),
        _step(2, "Attach metadata",
              "Upload a CSV/TSV with sample grouping information.",
              "/files", "File Manager"),
        _step(3, "Run the pipeline",
              "Auto-detects parameters and runs DADA2 denoising + taxonomy.",
              "/pipeline", "Pipeline"),
        _step(4, "Curate your data",
              "Remove rare ASVs, rarefy, filter samples, or combine datasets.",
              "/biom-browser", "Data Management"),
        _step(5, "Analyze & visualize",
              "Explore diversity, taxonomy, differential abundance, and pathways.",
              "/alpha", "Analysis"),
        _step(6, "Generate report",
              "Create a comprehensive PDF report with all analysis results.",
              "/report", "Report"),
    ], flush=True, className="mb-4"),

    html.Hr(),

    # ── Setup ────────────────────────────────────────────────────────────────
    html.H4("Setup", className="mt-4 mb-3"),

    html.H6("Supported Platforms"),
    html.Ul([
        html.Li("Linux (Ubuntu/Debian)"),
        html.Li("Windows (via WSL2)"),
        html.Li("macOS (Apple Silicon)"),
    ]),

    html.H6("Prerequisites", className="mt-3"),
    html.Ul([
        html.Li(["conda or mamba (", html.A("Miniforge",
                 href="https://github.com/conda-forge/miniforge",
                 target="_blank"), " recommended)"]),
        html.Li("8 GB RAM minimum (16 GB for PICRUSt2)"),
    ]),

    html.H6("Installation", className="mt-3"),
    dbc.Card(
        dbc.CardBody(
            html.Pre(
                "git clone https://github.com/tatsu1207/16S-Pipeline.git && cd 16S-Pipeline\n"
                "\n"
                "# Native install (Linux / WSL2)\n"
                "bash setup_ubuntu.sh\n"
                "\n"
                "# Or use Docker (any OS)\n"
                "docker compose up -d\n"
                "\n"
                "# Start the server (native install)\n"
                "conda activate microbiome_16S\n"
                "./run.sh",
                className="mb-0",
                style={"whiteSpace": "pre-wrap"},
            ),
        ),
        className="bg-dark border-secondary",
    ),

    html.P(
        "The setup script creates five conda environments: microbiome_16S (web app + CLI tools), "
        "dada2_16S (R + DADA2), analysis_16S (R + ALDEx2, DESeq2, ANCOM-BC2), "
        "maaslin2_16S (R + MaAsLin2, LinDA, vegan), and picrust2_16S (PICRUSt2). "
        "It also downloads the SILVA v138.1 reference database.",
        className="text-muted mt-2",
    ),

    html.Hr(),

    # ── Supported Methods ────────────────────────────────────────────────────
    html.H4("Supported Analysis Methods", className="mt-4 mb-2"),
    methods_table,

], fluid=True, className="pb-5")
