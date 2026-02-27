"""
MicrobiomeDash — Dash application initialization.
"""
import dash
import dash_bootstrap_components as dbc
import dash_uploader as du

from app.config import UPLOAD_DIR

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.DARKLY],
    suppress_callback_exceptions=True,
    title="MicrobiomeDash",
    update_title="Loading...",
)

server = app.server

# Allow large file uploads (BIOM files can be several MB; default 500 KB is too low)
server.config["MAX_CONTENT_LENGTH"] = 200 * 1024 * 1024  # 200 MB
server.config["MAX_FORM_MEMORY_SIZE"] = 200 * 1024 * 1024

# Configure dash-uploader for chunked FASTQ uploads
du.configure_upload(app, str(UPLOAD_DIR), use_upload_id=True)
