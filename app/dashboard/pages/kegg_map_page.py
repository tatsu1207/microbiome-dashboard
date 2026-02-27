"""
MicrobiomeDash — KEGG Map page for targeted pathway inspection.
"""
import io
import traceback

import dash_bootstrap_components as dbc
import pandas as pd
import plotly.graph_objects as go
from dash import Input, Output, State, callback_context, dcc, html, no_update

from app.analysis.kegg_map import (
    build_kegg_color_url,
    compute_pathway_activity,
    cross_reference_ids,
    fetch_pathway_ecs,
    fetch_pathway_kos,
    fetch_pathway_name,
    ko_coverage_by_group,
    normalize_pathway_id,
)
from app.analysis.pathways import load_picrust2_table, load_prediction_file, parse_uploaded_prediction_file, run_pathway_da
from app.analysis.shared import (
    find_metadata_for_samples,
    get_group_columns,
    get_picrust2_run_options,
    parse_uploaded_metadata,
    validate_metadata_vs_biom,
)
from app.dashboard.app import app as dash_app

HIDDEN = {"display": "none"}

# Preset common KEGG pathways
COMMON_PATHWAYS = [
    {"label": "-- None --", "value": ""},
    {"label": "Methane metabolism (00680)", "value": "00680"},
    {"label": "Nitrogen metabolism (00910)", "value": "00910"},
    {"label": "Sulfur metabolism (00920)", "value": "00920"},
    {"label": "Carbon fixation - Calvin cycle (00710)", "value": "00710"},
    {"label": "Carbon fixation - prokaryotes (00720)", "value": "00720"},
    {"label": "Glycolysis / Gluconeogenesis (00010)", "value": "00010"},
    {"label": "TCA cycle (00020)", "value": "00020"},
    {"label": "Pentose phosphate pathway (00030)", "value": "00030"},
    {"label": "Oxidative phosphorylation (00190)", "value": "00190"},
    {"label": "Purine metabolism (00230)", "value": "00230"},
    {"label": "Amino sugar & nucleotide sugar (00520)", "value": "00520"},
    {"label": "Fatty acid biosynthesis (00061)", "value": "00061"},
    {"label": "Butanoate metabolism (00650)", "value": "00650"},
    {"label": "Propanoate metabolism (00640)", "value": "00640"},
    {"label": "Pyruvate metabolism (00620)", "value": "00620"},
    {"label": "Benzoate degradation (00362)", "value": "00362"},
    {"label": "Flagellar assembly (02040)", "value": "02040"},
    {"label": "ABC transporters (02010)", "value": "02010"},
    {"label": "Two-component system (02020)", "value": "02020"},
    {"label": "Quorum sensing (02024)", "value": "02024"},
    {"label": "Biofilm formation (02025)", "value": "02025"},
]


def _try_auto_load_meta_for_run(run_id: str):
    """Try to auto-load metadata from the pipeline dataset that sourced this PICRUSt2 run."""
    from app.analysis.shared import get_dataset_metadata_df
    from app.db.database import SessionLocal
    from app.db.models import Picrust2Run

    try:
        run_id_int = int(run_id)
    except (ValueError, TypeError):
        return None, None

    db = SessionLocal()
    try:
        run = db.query(Picrust2Run).filter(Picrust2Run.id == run_id_int).first()
        if not run or not run.biom_path:
            return None, None
        return get_dataset_metadata_df(run.biom_path)
    finally:
        db.close()


def _try_find_meta_by_samples(picrust_dir: str, pred_type: str):
    """Load a PICRUSt2 table and search DB for metadata matching its sample names."""
    try:
        counts_df, _ = load_picrust2_table(picrust_dir, pred_type)
        sample_ids = counts_df.columns.tolist()
        if not sample_ids:
            return None, None, None
        return find_metadata_for_samples(sample_ids)
    except Exception:
        return None, None, None


def get_layout():
    run_opts = [{"label": "-- None --", "value": ""}] + get_picrust2_run_options()
    return dbc.Container([
        html.H3("KEGG Pathway Map", className="mb-2"),
        html.P(
            "Inspect individual KEGG pathways: see which KOs/ECs are detected in your "
            "PICRUSt2 data, view them on the KEGG pathway map, and compare activity between groups.",
            className="text-muted mb-4",
        ),

        # Stores
        dcc.Store(id="km-picrust-dir"),
        dcc.Store(id="km-uploaded-file"),
        dcc.Store(id="km-meta-store"),
        dcc.Store(id="km-sample-id-col"),
        dcc.Store(id="km-pathway-kos"),       # list of pathway member IDs
        dcc.Store(id="km-pathway-id-store"),   # normalized pathway ID
        dcc.Store(id="km-pathway-name"),
        dcc.Store(id="km-ko-csv"),             # coverage table CSV for download

        dbc.Row([
            # ── Left panel ────────────────────────────────────────────
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Input Data"),
                    dbc.CardBody([
                        dbc.Label("PICRUSt2 Results", className="fw-bold"),
                        dbc.Select(
                            id="km-select-run",
                            options=run_opts,
                            placeholder="Select a completed PICRUSt2 run...",
                            className="mb-2",
                        ),
                        html.Div("-- or upload prediction file --",
                                 className="text-center text-muted small mb-2"),
                        dcc.Upload(
                            id="km-upload-file",
                            children=html.Div([
                                "Drag & drop or ", html.A("select prediction file (.tsv / .tsv.gz)"),
                            ]),
                            style={
                                "borderWidth": "2px", "borderStyle": "dashed",
                                "borderRadius": "5px", "borderColor": "#555",
                                "textAlign": "center", "padding": "8px",
                            },
                            multiple=False,
                        ),
                        html.Div(id="km-input-status", className="mt-1 small mb-3"),

                        dbc.Label("Metadata (CSV/TSV)", className="fw-bold"),
                        html.Div(
                            "Auto-loaded from matching pipeline datasets. Upload to override.",
                            className="text-muted small mb-1",
                        ),
                        dcc.Upload(
                            id="km-upload-meta",
                            children=html.Div([
                                "Drag & drop or ", html.A("select metadata file"),
                            ]),
                            style={
                                "borderWidth": "2px", "borderStyle": "dashed",
                                "borderRadius": "5px", "borderColor": "#555",
                                "textAlign": "center", "padding": "8px",
                            },
                            multiple=False,
                        ),
                        html.Div(id="km-meta-status", className="mt-1 small mb-3"),
                    ]),
                ], className="mb-3"),

                dbc.Card([
                    dbc.CardHeader("Pathway Lookup"),
                    dbc.CardBody([
                        dbc.Label("Feature Type", className="fw-bold"),
                        dbc.RadioItems(
                            id="km-feature-type",
                            options=[
                                {"label": "KEGG Orthologs (KO)", "value": "ko"},
                                {"label": "Enzyme Commission (EC)", "value": "ec"},
                            ],
                            value="ko",
                            className="mb-3",
                        ),

                        dbc.Label("Common Pathways", className="fw-bold"),
                        dbc.Select(
                            id="km-common-pathways",
                            options=COMMON_PATHWAYS,
                            value="",
                            className="mb-2",
                        ),

                        dbc.Label("KEGG Pathway ID", className="fw-bold"),
                        dbc.Input(
                            id="km-pathway-id",
                            type="text",
                            placeholder="e.g. 00680",
                            className="mb-2",
                        ),

                        dbc.Button(
                            "Lookup Pathway", id="km-btn-lookup",
                            color="info", className="w-100",
                        ),
                        html.Div(id="km-lookup-status", className="mt-2 small"),
                    ]),
                ], className="mb-3"),

                dbc.Card([
                    dbc.CardHeader("Group Comparison"),
                    dbc.CardBody([
                        dbc.Label("Group Column", className="fw-bold"),
                        dbc.Select(
                            id="km-group-col",
                            placeholder="Select or upload metadata...",
                            className="mb-2",
                        ),
                        dbc.Label("Reference Group", className="fw-bold"),
                        dbc.Select(
                            id="km-ref-group",
                            placeholder="Select group column first...",
                            className="mb-2",
                        ),
                        dbc.Label("Test Group", className="fw-bold"),
                        dbc.Select(
                            id="km-test-group",
                            placeholder="Select group column first...",
                            className="mb-3",
                        ),
                        dbc.Button(
                            "Run Analysis", id="km-btn-run",
                            color="primary", className="w-100", disabled=True,
                        ),
                    ]),
                ], className="mb-3"),
            ], md=4),

            # ── Right panel ───────────────────────────────────────────
            dbc.Col([
                html.Div(id="km-pathway-info", className="mb-2"),
                html.Div(id="km-kegg-link", className="mb-3"),
                dcc.Loading(
                    type="dot",
                    color="#3498db",
                    children=[
                        html.Div(id="km-error", className="mb-2"),
                        dcc.Graph(id="km-boxplot", style=HIDDEN),
                        html.Div(id="km-stats-table", className="mb-3"),
                        html.Div(id="km-ko-table"),
                        html.Div(id="km-download-area"),
                    ],
                ),
            ], md=8),
        ]),
    ], fluid=True)


# ── Callback 1: Input data change ────────────────────────────────────────────


@dash_app.callback(
    Output("km-picrust-dir", "data"),
    Output("km-input-status", "children"),
    Output("km-meta-store", "data"),
    Output("km-sample-id-col", "data"),
    Output("km-meta-status", "children"),
    Output("km-group-col", "options"),
    Output("km-btn-run", "disabled"),
    Output("km-uploaded-file", "data"),
    Input("km-select-run", "value"),
    Input("km-upload-file", "contents"),
    Input("km-upload-meta", "contents"),
    State("km-upload-file", "filename"),
    State("km-upload-meta", "filename"),
    State("km-picrust-dir", "data"),
    State("km-meta-store", "data"),
    State("km-sample-id-col", "data"),
    State("km-uploaded-file", "data"),
    prevent_initial_call=True,
)
def on_input_change(run_id, file_contents, meta_contents,
                    file_filename, meta_filename,
                    prev_picrust_dir, prev_meta_json, prev_sid_col,
                    prev_uploaded_file):
    trigger = callback_context.triggered_id

    picrust_dir = prev_picrust_dir
    input_status = no_update
    meta_json = prev_meta_json
    sid_col = prev_sid_col
    meta_status = no_update
    group_opts = no_update
    btn_disabled = no_update
    uploaded_file = prev_uploaded_file

    if trigger == "km-select-run":
        uploaded_file = None
        if not run_id:
            return None, "", None, None, "", [], True, None

        picrust_dir = run_id
        input_status = dbc.Alert(f"PICRUSt2 run #{run_id} selected", color="success")

        db_meta, db_sid = _try_auto_load_meta_for_run(run_id)
        ds_name = None
        if db_meta is None:
            for pt in ["ko", "ec"]:
                db_meta, db_sid, ds_name = _try_find_meta_by_samples(run_id, pt)
                if db_meta is not None:
                    break

        if db_meta is not None:
            meta_json = db_meta.to_json(date_format="iso", orient="split")
            sid_col = db_sid
            gcols = get_group_columns(db_meta, db_sid)
            group_opts = [{"label": c, "value": c} for c in gcols]
            source = f" (from dataset '{ds_name}')" if ds_name else ""
            meta_status = dbc.Alert(
                f"Metadata auto-loaded{source}: {len(db_meta)} samples, "
                f"{len(gcols)} group columns",
                color="success",
            )
            btn_disabled = len(gcols) == 0
        else:
            meta_json = None
            sid_col = None
            group_opts = []
            meta_status = dbc.Alert(
                "No linked metadata found. Upload a metadata file.", color="info",
            )
            btn_disabled = True

    elif trigger == "km-upload-file" and file_contents:
        picrust_dir = None  # clear run selection

        temp_path, error = parse_uploaded_prediction_file(
            file_contents, file_filename or ""
        )
        if error:
            uploaded_file = None
            input_status = dbc.Alert(error, color="danger")
            btn_disabled = True
        else:
            uploaded_file = temp_path
            input_status = dbc.Alert(
                f"Prediction file loaded: {file_filename}", color="success"
            )

            # Auto-load metadata by matching sample names
            try:
                counts_df, _ = load_prediction_file(temp_path)
                sample_ids = counts_df.columns.tolist()
                db_meta, db_sid, ds_name = (
                    find_metadata_for_samples(sample_ids) if sample_ids
                    else (None, None, None)
                )
            except Exception:
                db_meta, db_sid, ds_name = None, None, None

            if db_meta is not None:
                meta_json = db_meta.to_json(date_format="iso", orient="split")
                sid_col = db_sid
                gcols = get_group_columns(db_meta, db_sid)
                group_opts = [{"label": c, "value": c} for c in gcols]
                source = f" (from dataset '{ds_name}')" if ds_name else ""
                meta_status = dbc.Alert(
                    f"Metadata auto-loaded{source}: {len(db_meta)} samples, "
                    f"{len(gcols)} group columns",
                    color="success",
                )
                btn_disabled = len(gcols) == 0
            else:
                meta_json = None
                sid_col = None
                group_opts = []
                meta_status = dbc.Alert(
                    "No matching metadata found. Upload a metadata file.",
                    color="info",
                )
                btn_disabled = True

    elif trigger == "km-upload-meta" and meta_contents:
        df, s_col, error = parse_uploaded_metadata(meta_contents, meta_filename or "")
        if error:
            meta_json = None
            sid_col = None
            meta_status = dbc.Alert(error, color="danger")
            group_opts = []
            btn_disabled = True
        else:
            meta_json = df.to_json(date_format="iso", orient="split")
            sid_col = s_col
            gcols = get_group_columns(df, s_col)
            group_opts = [{"label": c, "value": c} for c in gcols]
            meta_status = dbc.Alert(f"Metadata loaded: {len(df)} samples", color="success")
            has_data = picrust_dir is not None or uploaded_file is not None
            btn_disabled = (not has_data) or len(gcols) == 0

    return (picrust_dir, input_status, meta_json, sid_col, meta_status, group_opts,
            btn_disabled, uploaded_file)


# ── Callback 2: Group column change ──────────────────────────────────────────


@dash_app.callback(
    Output("km-ref-group", "options"),
    Output("km-test-group", "options"),
    Input("km-group-col", "value"),
    State("km-meta-store", "data"),
    prevent_initial_call=True,
)
def on_group_col_change(group_col, meta_json):
    if not group_col or not meta_json:
        return [], []
    meta_df = pd.read_json(io.StringIO(meta_json), orient="split")
    unique_vals = sorted(meta_df[group_col].dropna().unique().tolist())
    options = [{"label": str(v), "value": str(v)} for v in unique_vals]
    return options, options


# ── Callback 3: Common pathway preset selection ──────────────────────────────


@dash_app.callback(
    Output("km-pathway-id", "value"),
    Input("km-common-pathways", "value"),
    prevent_initial_call=True,
)
def on_common_pathway_select(preset_val):
    if not preset_val:
        return no_update
    return preset_val


# ── Callback 4: Pathway lookup ───────────────────────────────────────────────


@dash_app.callback(
    Output("km-pathway-kos", "data"),
    Output("km-pathway-id-store", "data"),
    Output("km-pathway-name", "data"),
    Output("km-pathway-info", "children"),
    Output("km-kegg-link", "children"),
    Output("km-lookup-status", "children"),
    Output("km-btn-run", "disabled", allow_duplicate=True),
    Input("km-btn-lookup", "n_clicks"),
    State("km-pathway-id", "value"),
    State("km-feature-type", "value"),
    State("km-picrust-dir", "data"),
    State("km-uploaded-file", "data"),
    State("km-meta-store", "data"),
    prevent_initial_call=True,
)
def on_pathway_lookup(n_clicks, pathway_id_raw, feature_type, picrust_dir,
                      uploaded_file, meta_json):
    no7 = (no_update,) * 7
    if not n_clicks or not pathway_id_raw:
        return no7

    try:
        pathway_id = normalize_pathway_id(pathway_id_raw)
    except ValueError as e:
        return (None, None, None, "", "",
                dbc.Alert(str(e), color="danger"), no_update)

    try:
        name = fetch_pathway_name(pathway_id)
    except Exception as e:
        return (None, None, None, "", "",
                dbc.Alert(f"KEGG API error: {e}", color="danger"), no_update)

    try:
        if feature_type == "ec":
            member_ids = fetch_pathway_ecs(pathway_id)
        else:
            member_ids = fetch_pathway_kos(pathway_id)
    except Exception as e:
        return (None, None, None, "", "",
                dbc.Alert(f"KEGG API error fetching {feature_type.upper()}s: {e}",
                          color="danger"), no_update)

    if not member_ids:
        ft_label = "KOs" if feature_type == "ko" else "ECs"
        return ([], pathway_id, name,
                dbc.Alert(f"{name} ({pathway_id}): No {ft_label} found", color="warning"),
                "", dbc.Alert("Pathway has no mapped features.", color="warning"),
                no_update)

    # Cross-reference with PICRUSt2 data if available
    info_children = []
    kegg_link = ""
    detected_ids = []
    btn_disabled = no_update
    has_data = picrust_dir or uploaded_file

    if has_data:
        try:
            if uploaded_file:
                counts_df, desc_df = load_prediction_file(uploaded_file)
            else:
                counts_df, desc_df = load_picrust2_table(picrust_dir, feature_type)
            xref = cross_reference_ids(member_ids, counts_df, feature_type)
            detected_ids = xref["detected"]
            n_det = xref["n_detected"]
            n_total = xref["total"]
            pct = xref["coverage_pct"]

            ft_label = "KOs" if feature_type == "ko" else "ECs"
            info_children = dbc.Alert([
                html.Strong(f"{name}"),
                html.Br(),
                f"{n_det}/{n_total} {ft_label} detected ({pct:.1f}% coverage)",
            ], color="info")

            # Build KEGG color URL — yellow preview (pre-DA)
            if detected_ids:
                color_map = {d: "#ffff00" for d in detected_ids}
                color_url = build_kegg_color_url(pathway_id, color_map)
                kegg_link = html.A(
                    "View colored pathway on KEGG (yellow = detected, pre-DA) ->",
                    href=color_url, target="_blank",
                    className="btn btn-outline-info btn-sm",
                )

            # Enable run button if metadata is available
            if meta_json:
                btn_disabled = False

        except Exception as e:
            info_children = dbc.Alert([
                html.Strong(f"{name}"),
                html.Br(),
                f"{len(member_ids)} {feature_type.upper()}s in pathway (could not load prediction data to cross-reference)",
            ], color="secondary")
    else:
        ft_label = "KOs" if feature_type == "ko" else "ECs"
        info_children = dbc.Alert([
            html.Strong(f"{name}"),
            html.Br(),
            f"{len(member_ids)} {ft_label} in pathway. Load prediction data to check coverage.",
        ], color="secondary")

    status = dbc.Alert(f"Pathway {pathway_id} loaded: {name}", color="success")

    return (member_ids, pathway_id, name, info_children, kegg_link, status, btn_disabled)


# ── Callback 5: Run analysis ─────────────────────────────────────────────────


@dash_app.callback(
    Output("km-error", "children"),
    Output("km-boxplot", "figure"),
    Output("km-boxplot", "style"),
    Output("km-stats-table", "children"),
    Output("km-ko-table", "children"),
    Output("km-download-area", "children"),
    Output("km-ko-csv", "data"),
    Output("km-kegg-link", "children", allow_duplicate=True),
    Output("km-pathway-info", "children", allow_duplicate=True),
    Input("km-btn-run", "n_clicks"),
    State("km-picrust-dir", "data"),
    State("km-uploaded-file", "data"),
    State("km-meta-store", "data"),
    State("km-sample-id-col", "data"),
    State("km-feature-type", "value"),
    State("km-pathway-kos", "data"),
    State("km-pathway-id-store", "data"),
    State("km-pathway-name", "data"),
    State("km-group-col", "value"),
    State("km-ref-group", "value"),
    State("km-test-group", "value"),
    prevent_initial_call=True,
)
def on_run(n_clicks, picrust_dir, uploaded_file, meta_json, sid_col, feature_type,
           pathway_member_ids, pathway_id, pathway_name,
           group_col, ref_group, test_group):
    no9 = (no_update,) * 9

    if not n_clicks:
        return no9

    if not all([meta_json, pathway_member_ids, group_col, ref_group, test_group]):
        return (dbc.Alert("Please fill all inputs and look up a pathway first.", color="warning"),
                *no9[1:])

    if not picrust_dir and not uploaded_file:
        return (dbc.Alert("Please select a PICRUSt2 run or upload a prediction file.",
                          color="warning"), *no9[1:])

    if ref_group == test_group:
        return (dbc.Alert("Reference and test groups must differ.", color="warning"),
                *no9[1:])

    try:
        meta_df = pd.read_json(io.StringIO(meta_json), orient="split")

        if uploaded_file:
            counts_df, desc_df = load_prediction_file(uploaded_file)
        else:
            counts_df, desc_df = load_picrust2_table(picrust_dir, feature_type)

        # Cross-reference
        xref = cross_reference_ids(pathway_member_ids, counts_df, feature_type)
        detected = xref["detected"]
        id_map = xref["id_map"]

        if not detected:
            return (dbc.Alert("No pathway features detected in your data.", color="warning"),
                    *no9[1:])

        # Validate metadata vs samples
        biom_sample_ids = counts_df.columns.tolist()
        match_info = validate_metadata_vs_biom(meta_df, sid_col, biom_sample_ids)
        matched = match_info["matched"]
        if len(matched) < 2:
            return (dbc.Alert("Need at least 2 matched samples.", color="danger"),
                    *no9[1:])

        # ── Subset counts to detected pathway features ──────────────
        original_detected = [id_map[d] for d in detected if d in id_map]
        original_detected = [o for o in original_detected if o in counts_df.index]
        subset_counts = counts_df.loc[original_detected]

        # ALDEx2 needs >= 2 features with non-zero counts
        nonzero_features = (subset_counts.sum(axis=1) > 0).sum()
        if nonzero_features < 2:
            return (dbc.Alert(
                f"ALDEx2 requires at least 2 features with non-zero counts, "
                f"but only {nonzero_features} found in this pathway.",
                color="warning",
            ), *no9[1:])

        # ── Run ALDEx2 via run_pathway_da ───────────────────────────
        da_results = run_pathway_da(
            subset_counts, meta_df, sid_col, matched,
            group_col, ref_group, test_group,
        )

        # ── Build color_map from ALDEx2 results ────────────────────
        # Reverse id_map: original index -> normalized ID
        rev_map = {v: k for k, v in id_map.items()}

        color_map = {}
        for _, row in da_results.iterrows():
            feat = row["feature"]
            norm_id = rev_map.get(feat, feat)
            q = row.get("qvalue", row.get("wi.eBH", 1.0))
            lfc = row.get("log2fc", row.get("diff.btw", 0.0))
            if pd.notna(q) and q < 0.05:
                if lfc < 0:
                    color_map[norm_id] = "#ff0000"   # higher in ref
                else:
                    color_map[norm_id] = "#00ff00"   # higher in test
            else:
                color_map[norm_id] = "#ffff00"       # not significant

        # Also color detected features not in DA results as yellow
        for d in detected:
            if d not in color_map:
                color_map[d] = "#ffff00"

        color_url = build_kegg_color_url(pathway_id, color_map)
        kegg_link = html.A(
            "View DA-colored pathway on KEGG ->",
            href=color_url, target="_blank",
            className="btn btn-outline-info btn-sm",
        )

        # ── Boxplot (summed pathway activity — kept for overview) ───
        activity = compute_pathway_activity(detected, counts_df, id_map)
        meta_df[sid_col] = meta_df[sid_col].astype(str)
        ref_ids = meta_df.loc[
            meta_df[group_col].astype(str) == str(ref_group), sid_col
        ].tolist()
        test_ids = meta_df.loc[
            meta_df[group_col].astype(str) == str(test_group), sid_col
        ].tolist()

        ref_vals = activity.reindex(ref_ids).dropna()
        test_vals = activity.reindex(test_ids).dropna()

        ft_label = "KOs" if feature_type == "ko" else "ECs"
        fig = go.Figure()
        fig.add_trace(go.Box(
            y=ref_vals.values, name=str(ref_group),
            marker_color="#636EFA", boxpoints="all",
        ))
        fig.add_trace(go.Box(
            y=test_vals.values, name=str(test_group),
            marker_color="#EF553B", boxpoints="all",
        ))
        fig.update_layout(
            title=f"Summed Pathway Activity: {pathway_name} ({pathway_id})",
            yaxis_title=f"Summed {ft_label} Abundance",
            template="plotly_dark",
            showlegend=False,
            height=400,
        )

        # ── ALDEx2 summary stats table ──────────────────────────────
        n_tested = len(da_results)
        n_sig = int(da_results.apply(
            lambda r: (r.get("qvalue", r.get("wi.eBH", 1.0)) or 1.0) < 0.05, axis=1
        ).sum())
        n_up = int(da_results.apply(
            lambda r: ((r.get("qvalue", r.get("wi.eBH", 1.0)) or 1.0) < 0.05
                       and (r.get("log2fc", r.get("diff.btw", 0.0)) or 0.0) > 0),
            axis=1,
        ).sum())
        n_down = n_sig - n_up

        stats_table = dbc.Table(
            [
                html.Thead(html.Tr([
                    html.Th("Statistic"), html.Th("Value"),
                ])),
                html.Tbody([
                    html.Tr([html.Td("Method"), html.Td("ALDEx2")]),
                    html.Tr([html.Td("Features tested"), html.Td(str(n_tested))]),
                    html.Tr([html.Td("Significant (q < 0.05)"), html.Td(str(n_sig))]),
                    html.Tr([html.Td(f"Higher in {test_group}"), html.Td(str(n_up))]),
                    html.Tr([html.Td(f"Higher in {ref_group}"), html.Td(str(n_down))]),
                    html.Tr([html.Td("Color legend"), html.Td([
                        html.Span("\u2588 ", style={"color": "#00ff00"}),
                        f"higher in {test_group}  ",
                        html.Span("\u2588 ", style={"color": "#ff0000"}),
                        f"higher in {ref_group}  ",
                        html.Span("\u2588 ", style={"color": "#ffff00"}),
                        "not significant",
                    ])]),
                ]),
            ],
            bordered=True, hover=True, color="dark", size="sm",
            className="w-auto",
        )

        # ── Coverage table with ALDEx2 columns merged ───────────────
        coverage_df = ko_coverage_by_group(
            pathway_member_ids, counts_df, id_map, meta_df, sid_col,
            group_col, ref_group, test_group, desc_df,
        )

        # Merge DA results into coverage table
        if not da_results.empty:
            da_merge = da_results[["feature"]].copy()
            da_merge["log2fc"] = da_results.get(
                "log2fc", da_results.get("diff.btw", pd.Series(dtype=float))
            )
            da_merge["qvalue"] = da_results.get(
                "qvalue", da_results.get("wi.eBH", pd.Series(dtype=float))
            )
            # Map feature (original index) back to normalized ID
            da_merge["feature_id"] = da_merge["feature"].map(
                lambda f: rev_map.get(f, f)
            )
            coverage_df = coverage_df.merge(
                da_merge[["feature_id", "log2fc", "qvalue"]],
                on="feature_id", how="left",
            )

            # Add color indicator column
            def _color_label(row):
                q = row.get("qvalue")
                lfc = row.get("log2fc")
                if pd.isna(q) or pd.isna(lfc):
                    return "-"
                if q < 0.05:
                    return f"UP ({test_group})" if lfc > 0 else f"UP ({ref_group})"
                return "n.s."
            coverage_df["DA_result"] = coverage_df.apply(_color_label, axis=1)

        # Format for display
        display_df = coverage_df.head(100).copy()
        display_df["detected"] = display_df["detected"].map(
            {True: "Yes", False: "No"}
        )
        for col in display_df.columns:
            if col.startswith("mean_"):
                display_df[col] = display_df[col].map(
                    lambda x: f"{x:.2f}" if pd.notna(x) else "-"
                )
        if "log2fc" in display_df.columns:
            display_df["log2fc"] = display_df["log2fc"].map(
                lambda x: f"{x:.3f}" if pd.notna(x) else "-"
            )
        if "qvalue" in display_df.columns:
            display_df["qvalue"] = display_df["qvalue"].map(
                lambda x: f"{x:.4g}" if pd.notna(x) else "-"
            )

        n_det = int(coverage_df["detected"].sum())
        n_total = len(coverage_df)
        ko_table = html.Div([
            html.H5(
                f"{ft_label} Coverage ({n_det}/{n_total} detected)"
                + (f" - showing top 100" if n_total > 100 else ""),
                className="mt-3",
            ),
            dbc.Table.from_dataframe(
                display_df, striped=True, bordered=True, hover=True,
                color="dark", size="sm",
            ),
        ])

        # ── Download ────────────────────────────────────────────────
        csv_data = coverage_df.to_csv(sep="\t", index=False)
        dl = html.Div([
            dcc.Download(id="km-download"),
            dbc.Button(
                "Download Coverage Table (TSV)", id="km-btn-download",
                color="secondary", size="sm", className="mt-2",
            ),
        ])

        # ── Updated pathway info banner ─────────────────────────────
        pathway_info = dbc.Alert([
            html.Strong(f"{pathway_name}"),
            html.Br(),
            f"{len(detected)}/{len(pathway_member_ids)} {ft_label} detected — "
            f"ALDEx2: {n_sig} significant ({n_up} up in {test_group}, "
            f"{n_down} up in {ref_group})",
        ], color="success")

        return ("", fig, {"display": "block"}, stats_table, ko_table, dl, csv_data,
                kegg_link, pathway_info)

    except Exception as e:
        return (
            dbc.Alert(f"Error: {e}\n{traceback.format_exc()}", color="danger"),
            *no9[1:],
        )


# ── Callback 6: Download ─────────────────────────────────────────────────────


@dash_app.callback(
    Output("km-download", "data"),
    Input("km-btn-download", "n_clicks"),
    State("km-ko-csv", "data"),
    prevent_initial_call=True,
)
def on_download(n_clicks, csv_data):
    if not n_clicks or not csv_data:
        return no_update
    return dcc.send_string(csv_data, "kegg_pathway_coverage.tsv")
