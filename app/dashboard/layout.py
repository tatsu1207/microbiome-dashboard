"""
MicrobiomeDash — Main layout with sidebar navigation and page routing.
"""
import dash_bootstrap_components as dbc
from dash import Input, Output, dcc, html

from app.dashboard.app import app as dash_app


def _nav_link(label, href, disabled=False):
    return dbc.NavLink(label, href=href, active="exact", disabled=disabled)


sidebar = html.Div(
    [
        html.H4("16S Pipeline", className="text-center my-3"),
        html.Hr(),
        dbc.Nav(
            [
                _nav_link("Introduction", "/"),
                html.Hr(className="my-2"),
                # File Manager
                html.H6("FILE MANAGER", className="text-muted mt-3 mb-2 px-3"),
                _nav_link("File Registration", "/files"),
                # Pipeline
                html.H6("PIPELINE", className="text-muted mt-4 mb-2 px-3"),
                _nav_link("DADA2", "/pipeline"),
                # Data Management
                html.H6("DATA MANAGEMENT", className="text-muted mt-4 mb-2 px-3"),
                _nav_link("BIOM Browser", "/biom-browser"),
                _nav_link("Rare ASV Removal", "/rare-asv"),
                _nav_link("Rarefy & Filter", "/subsample"),
                _nav_link("Combine Biom Files", "/combine"),
                _nav_link("V-Region Extraction", "/datasets"),
                _nav_link("BIOM & MOTHUR", "/mothur"),
                # Tool 3: Analysis
                html.H6("ANALYSIS", className="text-muted mt-4 mb-2 px-3"),
                _nav_link("Alpha Diversity", "/alpha"),
                _nav_link("Beta Diversity", "/beta"),
                _nav_link("Taxonomy", "/taxonomy"),
                _nav_link("Diff. Abundance", "/diff-abundance"),
                # Pathway Analysis
                html.H6("PATHWAY ANALYSIS", className="text-muted mt-4 mb-2 px-3"),
                _nav_link("PICRUSt2", "/picrust2"),
                _nav_link("Comparison", "/pathways"),
                _nav_link("KEGG Map", "/kegg-map"),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    className="bg-dark vh-100 position-fixed p-3",
    style={"width": "260px", "overflowY": "auto"},
)

SIDEBAR_WIDTH = "260px"


def create_layout():
    return html.Div(
        [
            dcc.Location(id="url", refresh=False),
            sidebar,
            html.Div(
                id="page-content",
                className="p-4",
                style={"marginLeft": SIDEBAR_WIDTH},
            ),
        ]
    )


# ── Page routing callback ────────────────────────────────────────────────────


@dash_app.callback(Output("page-content", "children"), Input("url", "pathname"))
def render_page(pathname):
    if pathname == "/":
        from app.dashboard.pages.intro_page import layout as intro_layout

        return intro_layout

    if pathname in ("/upload", "/pipeline"):
        from app.dashboard.pages.pipeline_status import layout as pipeline_layout

        return pipeline_layout

    if pathname == "/files":
        from app.dashboard.pages.file_manager import layout as fm_layout

        return fm_layout

    if pathname == "/picrust2":
        from app.dashboard.pages.picrust2_page import get_layout as picrust2_layout

        return picrust2_layout()

    if pathname == "/datasets":
        from app.dashboard.pages.datasets_page import get_layout as datasets_layout

        return datasets_layout()

    if pathname == "/combine":
        from app.dashboard.pages.combine_page import get_layout as combine_layout

        return combine_layout()

    if pathname == "/mothur":
        from app.dashboard.pages.mothur_page import get_layout as mothur_layout

        return mothur_layout()

    if pathname == "/subsample":
        from app.dashboard.pages.subsampling_page import get_layout as subsample_layout

        return subsample_layout()

    if pathname == "/rare-asv":
        from app.dashboard.pages.rare_asv_page import get_layout as rare_asv_layout

        return rare_asv_layout()

    if pathname == "/biom-browser":
        from app.dashboard.pages.biom_browser_page import get_layout as bb_layout

        return bb_layout()

    if pathname == "/alpha":
        from app.dashboard.pages.alpha_page import get_layout as alpha_layout

        return alpha_layout()

    if pathname == "/beta":
        from app.dashboard.pages.beta_page import get_layout as beta_layout

        return beta_layout()

    if pathname == "/taxonomy":
        from app.dashboard.pages.taxonomy_page import get_layout as taxonomy_layout

        return taxonomy_layout()

    if pathname == "/diff-abundance":
        from app.dashboard.pages.diff_abundance_page import get_layout as da_layout

        return da_layout()

    if pathname == "/pathways":
        from app.dashboard.pages.pathways_page import get_layout as pathways_layout

        return pathways_layout()

    if pathname == "/kegg-map":
        from app.dashboard.pages.kegg_map_page import get_layout as kegg_map_layout

        return kegg_map_layout()

    return dbc.Container(
        [
            html.H3("Coming Soon", className="text-muted mt-5"),
            html.P(
                "This page will be available in a future update.",
                className="text-muted",
            ),
        ]
    )
