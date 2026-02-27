"""
MicrobiomeDash — Server-side file/directory browser modal component.

Click-to-navigate explorer: click folders to enter them, click files to select.
"""
from pathlib import Path

import dash_bootstrap_components as dbc
from dash import ALL, Input, Output, State, ctx, dcc, html, no_update

from app.dashboard.app import app as dash_app


def _get_start_path() -> str:
    """Return a sensible starting path."""
    for p in [Path.home(), Path("/mnt/c"), Path("/tmp")]:
        if p.exists():
            return str(p)
    return "/"


def _list_dir(directory: str, mode: str) -> list:
    """Return ListGroup items for a directory listing.

    mode: "directory" — only show directories (for FASTQ dir selection)
          "file"      — show directories + files (for metadata selection)
    """
    d = Path(directory)
    if not d.is_dir():
        return [dbc.ListGroupItem("Directory not accessible", color="danger")]

    items = []

    # Parent directory
    if d.parent != d:
        items.append(
            dbc.ListGroupItem(
                ".. (up)",
                id={"type": "browser-entry", "path": str(d.parent)},
                action=True,
                className="py-2 text-warning fw-bold",
                style={"cursor": "pointer"},
            )
        )

    try:
        entries = sorted(d.iterdir(), key=lambda e: (not e.is_dir(), e.name.lower()))
    except PermissionError:
        items.append(dbc.ListGroupItem("Permission denied", color="danger"))
        return items

    for entry in entries:
        if entry.name.startswith("."):
            continue
        if entry.is_dir():
            items.append(
                dbc.ListGroupItem(
                    [html.Span("📁 ", style={"marginRight": "6px"}), entry.name],
                    id={"type": "browser-entry", "path": str(entry)},
                    action=True,
                    className="py-2",
                    style={"cursor": "pointer"},
                )
            )
        elif entry.is_file():
            size = entry.stat().st_size
            if size > 1024 * 1024:
                size_str = f"{size / (1024*1024):.1f} MB"
            elif size > 1024:
                size_str = f"{size / 1024:.0f} KB"
            else:
                size_str = f"{size} B"
            if mode == "file":
                # Clickable file entry for file-selection mode
                items.append(
                    dbc.ListGroupItem(
                        [
                            html.Span("📄 ", style={"marginRight": "6px"}),
                            entry.name,
                            html.Small(f"  ({size_str})", className="text-muted"),
                        ],
                        id={"type": "browser-entry", "path": str(entry)},
                        action=True,
                        className="py-2",
                        style={"cursor": "pointer"},
                    )
                )
            else:
                # Non-clickable file preview for directory-selection mode
                items.append(
                    dbc.ListGroupItem(
                        [
                            html.Span("📄 ", style={"marginRight": "6px"}),
                            html.Span(entry.name, className="text-muted"),
                            html.Small(f"  ({size_str})", className="text-muted"),
                        ],
                        className="py-1",
                        style={"opacity": "0.6", "fontSize": "0.85rem"},
                        disabled=True,
                    )
                )

    if not items:
        items.append(dbc.ListGroupItem("(empty directory)", disabled=True))

    return items


def create_browser_modal():
    """Return the browser modal. Include once in the page layout."""
    return html.Div(
        [
            dcc.Store(id="browser-mode", data="directory"),
            dcc.Store(id="browser-current-path", data=""),
            dcc.Store(id="browser-selected-path", data=""),
            dbc.Modal(
                [
                    dbc.ModalHeader(dbc.ModalTitle(id="browser-title")),
                    dbc.ModalBody(
                        [
                            # Current path breadcrumb
                            html.Div(
                                id="browser-path-display",
                                className="mb-2 p-2 rounded",
                                style={
                                    "backgroundColor": "#1a1a2e",
                                    "fontFamily": "monospace",
                                    "fontSize": "0.9rem",
                                    "wordBreak": "break-all",
                                },
                            ),
                            # Directory listing
                            html.Div(
                                id="browser-listing",
                                style={
                                    "maxHeight": "400px",
                                    "overflowY": "auto",
                                },
                            ),
                            # Hint text
                            html.Small(
                                id="browser-hint",
                                className="text-muted mt-2 d-block",
                            ),
                        ]
                    ),
                    dbc.ModalFooter(
                        [
                            dbc.Button(
                                "Select This Directory",
                                id="browser-btn-select",
                                color="primary",
                            ),
                            dbc.Button(
                                "Cancel",
                                id="browser-btn-cancel",
                                color="secondary",
                                className="ms-2",
                            ),
                        ]
                    ),
                ],
                id="browser-modal",
                size="lg",
                is_open=False,
            ),
        ]
    )


# ── Open / close modal ───────────────────────────────────────────────────────


@dash_app.callback(
    Output("browser-modal", "is_open"),
    Output("browser-mode", "data"),
    Output("browser-current-path", "data"),
    Output("browser-title", "children"),
    Input("btn-browse-fastq", "n_clicks"),
    Input("browser-btn-cancel", "n_clicks"),
    Input("browser-btn-select", "n_clicks"),
    State("browser-modal", "is_open"),
    prevent_initial_call=True,
)
def toggle_modal(n_fastq, n_cancel, n_select, is_open):
    triggered = ctx.triggered_id
    if triggered == "btn-browse-fastq":
        return True, "directory", _get_start_path(), "Select FASTQ Directory"
    # Cancel or Select closes
    return False, no_update, no_update, no_update


# ── Populate listing when path changes ────────────────────────────────────────


@dash_app.callback(
    Output("browser-listing", "children"),
    Output("browser-path-display", "children"),
    Output("browser-hint", "children"),
    Output("browser-btn-select", "children"),
    Input("browser-current-path", "data"),
    State("browser-mode", "data"),
    prevent_initial_call=True,
)
def refresh_listing(current_path, mode):
    if not current_path:
        return [], "", "", "Select"

    items = _list_dir(current_path, mode)
    listing = dbc.ListGroup(items, flush=True)

    hint = (
        "Click a folder to navigate into it. Click 'Select This Directory' to choose the current folder."
        if mode == "directory"
        else "Click a folder to navigate. Click a file to select it, then click 'Select This File'."
    )
    btn_text = "Select This Directory" if mode == "directory" else "Select This File"

    return listing, current_path, hint, btn_text


# ── Handle clicking an entry (navigate into dirs or select files) ─────────────


@dash_app.callback(
    Output("browser-current-path", "data", allow_duplicate=True),
    Output("browser-selected-path", "data"),
    Input({"type": "browser-entry", "path": ALL}, "n_clicks"),
    State("browser-mode", "data"),
    prevent_initial_call=True,
)
def on_entry_click(n_clicks_list, mode):
    if not any(n_clicks_list):
        return no_update, no_update

    triggered = ctx.triggered_id
    if not triggered:
        return no_update, no_update

    selected = triggered["path"]
    p = Path(selected)

    if p.is_dir():
        # Navigate into directory
        return str(p), str(p)

    return no_update, no_update


# ── "Select" button → push result to the input fields ────────────────────────


@dash_app.callback(
    Output("input-fastq-dir", "value"),
    Input("browser-btn-select", "n_clicks"),
    State("browser-mode", "data"),
    State("browser-selected-path", "data"),
    State("browser-current-path", "data"),
    prevent_initial_call=True,
)
def on_select(n_clicks, mode, selected_path, current_path):
    if not n_clicks:
        return no_update

    if mode == "directory":
        path = selected_path or current_path
        if path and Path(path).is_dir():
            return path

    return no_update
