"""
16S Pipeline — Analysis Report page.

Generate a comprehensive PDF report for a completed pipeline dataset.
"""
import logging

import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, no_update

from app.analysis.shared import get_dataset_metadata_df, get_pipeline_biom_options
from app.dashboard.app import app as dash_app

logger = logging.getLogger(__name__)

ALL_SECTIONS = [
    {"label": "Dataset Summary", "value": "summary"},
    {"label": "Materials & Methods", "value": "methods"},
    {"label": "Alpha Diversity", "value": "alpha"},
    {"label": "Beta Diversity (PCoA + PERMANOVA)", "value": "beta"},
    {"label": "Taxonomy Composition", "value": "taxonomy"},
    {"label": "Read Tracking", "value": "reads"},
]

DEFAULT_SECTIONS = [s["value"] for s in ALL_SECTIONS]


def get_layout():
    return dbc.Container(
        [
            html.H3("Analysis Report", className="mb-2"),
            html.P(
                "Generate a comprehensive PDF report combining key analyses for a completed dataset.",
                className="text-muted mb-4",
            ),
            dcc.Store(id="rpt-biom-path"),
            dcc.Store(id="rpt-dataset-id"),
            dcc.Download(id="rpt-download"),

            dbc.Row(
                [
                    # Left column: settings
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5("Settings", className="mb-3"),

                                    # Dataset selector
                                    dbc.Label("Dataset"),
                                    dcc.Dropdown(
                                        id="rpt-dataset-dd",
                                        placeholder="Select a completed dataset...",
                                        className="mb-3",
                                    ),

                                    # Group column
                                    dbc.Label("Group Column"),
                                    dcc.Dropdown(
                                        id="rpt-group-dd",
                                        placeholder="Select grouping variable...",
                                        className="mb-1",
                                    ),
                                    html.Small(
                                        "Required for alpha/beta diversity sections.",
                                        className="text-muted d-block mb-3",
                                    ),

                                    # Sections checklist
                                    dbc.Label("Sections to Include"),
                                    dbc.Checklist(
                                        id="rpt-sections",
                                        options=ALL_SECTIONS,
                                        value=DEFAULT_SECTIONS,
                                        className="mb-4",
                                    ),

                                    # Generate button
                                    dbc.Button(
                                        "Generate PDF Report",
                                        id="rpt-generate-btn",
                                        color="primary",
                                        className="w-100",
                                        disabled=True,
                                    ),
                                ]
                            )
                        ),
                        md=5,
                    ),

                    # Right column: status
                    dbc.Col(
                        dbc.Card(
                            dbc.CardBody(
                                [
                                    html.H5("Report Status", className="mb-3"),
                                    html.Div(id="rpt-status", children=[
                                        html.P(
                                            "Select a dataset and click Generate.",
                                            className="text-muted",
                                        ),
                                    ]),
                                ]
                            )
                        ),
                        md=7,
                    ),
                ],
                className="g-3",
            ),
        ],
        fluid=True,
    )


# ── Callbacks ────────────────────────────────────────────────────────────────

@dash_app.callback(
    Output("rpt-dataset-dd", "options"),
    Input("rpt-dataset-dd", "id"),
)
def populate_datasets(_):
    """Load completed datasets into the dropdown."""
    options = get_pipeline_biom_options()
    return options


@dash_app.callback(
    [Output("rpt-biom-path", "data"),
     Output("rpt-dataset-id", "data"),
     Output("rpt-group-dd", "options"),
     Output("rpt-group-dd", "value")],
    Input("rpt-dataset-dd", "value"),
)
def on_dataset_selected(biom_path):
    """When a dataset is selected, load its metadata columns."""
    if not biom_path:
        return None, None, [], None

    # Extract dataset ID from dropdown option
    from app.db.database import SessionLocal
    from app.db.models import Dataset

    db = SessionLocal()
    try:
        ds = db.query(Dataset).filter(Dataset.asv_table_path == biom_path).first()
        ds_id = ds.id if ds else None
    finally:
        db.close()

    meta_df, sid_col = get_dataset_metadata_df(biom_path)
    if meta_df is None:
        return biom_path, ds_id, [], None

    group_cols = [c for c in meta_df.columns if c != sid_col]
    options = [{"label": c, "value": c} for c in group_cols]
    default = group_cols[0] if group_cols else None
    return biom_path, ds_id, options, default


@dash_app.callback(
    Output("rpt-generate-btn", "disabled"),
    [Input("rpt-dataset-dd", "value"),
     Input("rpt-sections", "value")],
)
def toggle_generate_btn(biom_path, sections):
    """Enable generate button when a dataset and at least one section are selected."""
    return not biom_path or not sections


@dash_app.callback(
    [Output("rpt-status", "children"),
     Output("rpt-download", "data")],
    Input("rpt-generate-btn", "n_clicks"),
    [State("rpt-dataset-id", "data"),
     State("rpt-group-dd", "value"),
     State("rpt-sections", "value")],
    prevent_initial_call=True,
)
def on_generate(n_clicks, dataset_id, group_col, sections):
    """Generate the PDF report."""
    if not n_clicks or not dataset_id:
        return no_update, no_update

    # Warn if alpha/beta selected without group column
    needs_group = {"alpha", "beta"}
    selected_needs = needs_group & set(sections)
    if selected_needs and not group_col:
        sections = [s for s in sections if s not in needs_group]
        if not sections:
            return (
                dbc.Alert(
                    "Alpha and Beta diversity require a grouping column. "
                    "Please select a group column or uncheck those sections.",
                    color="warning",
                ),
                no_update,
            )

    try:
        from app.report.report_generator import generate_report_pdf

        pdf_path = generate_report_pdf(
            dataset_id=dataset_id,
            group_col=group_col,
            sections=sections,
        )

        status = dbc.Alert(
            [
                html.Strong("Report generated successfully!"),
                html.Br(),
                html.Small(f"Saved to: {pdf_path}"),
            ],
            color="success",
        )
        return status, dcc.send_file(pdf_path)

    except Exception as e:
        logger.exception("Report generation failed")
        return (
            dbc.Alert(f"Report generation failed: {e}", color="danger"),
            no_update,
        )
