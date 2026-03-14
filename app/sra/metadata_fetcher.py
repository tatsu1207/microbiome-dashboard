"""
SRA Metadata Fetcher — Fetch BioSample metadata from NCBI for downloaded SRA data.
"""
import json
import logging
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from app.db.database import SessionLocal
from app.db.models import UploadMetadata

logger = logging.getLogger(__name__)


def fetch_biosample_metadata(srr_ids: list[str]) -> dict[str, dict[str, str]]:
    """
    Fetch BioSample metadata for a list of SRR accessions via NCBI E-utilities.
    Returns {sample_name: {attribute_name: attribute_value, ...}}.

    The sample_name is derived from the SRR ID (same as fasterq-dump output naming).
    """
    result = {}

    for srr_id in srr_ids:
        try:
            attrs = _fetch_single_biosample(srr_id)
            if attrs:
                result[srr_id] = attrs
            # Respect NCBI rate limit (3 requests/sec without API key)
            time.sleep(0.4)
        except Exception as e:
            logger.warning(f"Failed to fetch BioSample metadata for {srr_id}: {e}")

    return result


def _fetch_single_biosample(srr_id: str) -> dict[str, str]:
    """Fetch BioSample attributes for a single SRR accession."""
    # Step 1: Get the linked BioSample ID from SRA
    params = urllib.parse.urlencode({
        "db": "sra",
        "term": srr_id,
        "retmode": "json",
    })
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}"

    with urllib.request.urlopen(url, timeout=15) as resp:
        data = json.loads(resp.read())

    uid_list = data.get("esearchresult", {}).get("idlist", [])
    if not uid_list:
        return {}

    # Step 2: Get linked BioSample from SRA record
    sra_uid = uid_list[0]
    link_params = urllib.parse.urlencode({
        "dbfrom": "sra",
        "db": "biosample",
        "id": sra_uid,
        "retmode": "json",
    })
    link_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?{link_params}"

    with urllib.request.urlopen(link_url, timeout=15) as resp:
        link_data = json.loads(resp.read())

    # Extract BioSample UID
    biosample_uid = None
    for linkset in link_data.get("linksets", []):
        for linksetdb in linkset.get("linksetdbs", []):
            if linksetdb.get("dbto") == "biosample":
                links = linksetdb.get("links", [])
                if links:
                    biosample_uid = links[0]
                    break

    if not biosample_uid:
        return {}

    time.sleep(0.4)

    # Step 3: Fetch BioSample XML
    fetch_params = urllib.parse.urlencode({
        "db": "biosample",
        "id": biosample_uid,
        "retmode": "xml",
    })
    fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{fetch_params}"

    with urllib.request.urlopen(fetch_url, timeout=15) as resp:
        xml_text = resp.read().decode("utf-8")

    # Step 4: Parse BioSample attributes
    return _parse_biosample_xml(xml_text)


def _parse_biosample_xml(xml_text: str) -> dict[str, str]:
    """Parse BioSample XML to extract attribute name-value pairs."""
    attrs = {}
    try:
        root = ET.fromstring(xml_text)
        for biosample in root.iter("BioSample"):
            # Get accession
            accession = biosample.get("accession", "")
            if accession:
                attrs["biosample_accession"] = accession

            # Get attributes
            for attr in biosample.iter("Attribute"):
                name = attr.get("attribute_name") or attr.get("harmonized_name", "")
                value = attr.text or ""
                if name and value.strip():
                    attrs[name] = value.strip()
    except ET.ParseError:
        logger.warning("Failed to parse BioSample XML")

    return attrs


def store_biosample_metadata(upload_id: int, metadata: dict[str, dict[str, str]]):
    """
    Store fetched BioSample metadata in the UploadMetadata table.
    metadata = {sample_name: {attr_name: attr_value, ...}}
    """
    if not metadata:
        return

    db = SessionLocal()
    try:
        for sample_name, attrs in metadata.items():
            for key, value in attrs.items():
                # Skip if already exists
                existing = (
                    db.query(UploadMetadata)
                    .filter(
                        UploadMetadata.upload_id == upload_id,
                        UploadMetadata.sample_name == sample_name,
                        UploadMetadata.key == key,
                    )
                    .first()
                )
                if not existing:
                    db.add(
                        UploadMetadata(
                            upload_id=upload_id,
                            sample_name=sample_name,
                            key=key,
                            value=value,
                        )
                    )
        db.commit()
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()
