"""
Methods Text Generator — Auto-generate a Materials & Methods paragraph
from pipeline parameters stored in the Dataset model.
"""
from app.db.database import SessionLocal
from app.db.models import Dataset, Sample


def generate_methods_text(dataset_id: int) -> str:
    """
    Generate a Materials & Methods paragraph for a completed dataset.
    Returns a formatted text string ready for publication.
    """
    db = SessionLocal()
    try:
        dataset = db.query(Dataset).filter(Dataset.id == dataset_id).first()
        if not dataset:
            return "Dataset not found."

        samples = (
            db.query(Sample)
            .filter(Sample.dataset_id == dataset_id)
            .all()
        )

        return _build_text(dataset, samples)
    finally:
        db.close()


def _build_text(dataset, samples: list) -> str:
    """Build the methods paragraph from dataset and sample data."""
    parts = []

    # Sequencing description
    region = dataset.variable_region or "16S rRNA"
    platform = _platform_name(dataset.platform)
    seq_type = dataset.sequencing_type or "paired-end"

    parts.append(
        f"Raw 16S rRNA gene amplicon sequences targeting the {region} region "
        f"were generated using {platform} {seq_type} sequencing."
    )

    # Primer trimming
    if dataset.custom_fwd_primer or dataset.custom_rev_primer:
        fwd = dataset.custom_fwd_primer or "N/A"
        rev = dataset.custom_rev_primer or "N/A"
        parts.append(
            f"Primer sequences (forward: {fwd}, reverse: {rev}) were removed "
            f"using Cutadapt (Martin, 2011)."
        )
    else:
        parts.append(
            "Region-specific primer sequences were removed using Cutadapt (Martin, 2011)."
        )

    # DADA2 processing
    dada2_text = "Quality filtering and amplicon sequence variant (ASV) inference were performed using DADA2 (Callahan et al., 2016)"

    if dataset.platform in ("pacbio", "nanopore"):
        error_model = "PacBioErrfun" if dataset.platform == "pacbio" else "loessErrfun"
        dada2_text += f" with the {error_model} error model for long-read data"
    elif seq_type == "paired-end":
        trunc_f = dataset.trunc_len_f
        trunc_r = dataset.trunc_len_r
        if trunc_f and trunc_r:
            dada2_text += (
                f". Forward reads were truncated at {trunc_f} bp "
                f"and reverse reads at {trunc_r} bp"
            )
            if dataset.min_overlap:
                dada2_text += f" with a minimum overlap of {dataset.min_overlap} bp for read merging"
    elif seq_type == "single-end":
        trunc_f = dataset.trunc_len_f
        if trunc_f:
            dada2_text += f". Reads were truncated at {trunc_f} bp"

    dada2_text += ". Chimeric sequences were removed using the consensus method."
    parts.append(dada2_text)

    # Taxonomy
    parts.append(
        "Taxonomy was assigned using the naive Bayesian classifier (Wang et al., 2007) "
        "against the SILVA NR99 v138.1 reference database (Quast et al., 2013), "
        "with species-level assignment via exact matching."
    )

    # Phylogenetic tree
    parts.append(
        "Representative ASV sequences were aligned using MAFFT (Katoh and Standley, 2013), "
        "and a phylogenetic tree was constructed using FastTree 2 (Price et al., 2010) "
        "under the GTR+CAT model."
    )

    # Summary statistics
    if samples:
        nonchim = [s.read_count_nonchimeric for s in samples if s.read_count_nonchimeric]
        if nonchim:
            mean_reads = int(sum(nonchim) / len(nonchim))
            min_reads = min(nonchim)
            max_reads = max(nonchim)
            parts.append(
                f"A total of {dataset.asv_count or '?'} ASVs were identified "
                f"across {dataset.sample_count or len(samples)} samples, "
                f"with a mean of {mean_reads:,} non-chimeric reads per sample "
                f"(range: {min_reads:,}–{max_reads:,})."
            )
        else:
            parts.append(
                f"A total of {dataset.asv_count or '?'} ASVs were identified "
                f"across {dataset.sample_count or len(samples)} samples."
            )

    # PICRUSt2 (if run)
    if dataset.picrust_dir_path:
        parts.append(
            "Functional potential was predicted using PICRUSt2 (Douglas et al., 2020) "
            "based on phylogenetic placement of ASV sequences against a reference tree "
            "of sequenced genomes."
        )

    # Pipeline credit
    parts.append(
        "All analyses were performed using 16S-Pipeline "
        "(https://github.com/tatsu1207/16S-Pipeline), "
        "an open-source web-based platform for end-to-end 16S rRNA amplicon analysis."
    )

    return " ".join(parts)


def _platform_name(platform: str | None) -> str:
    """Convert internal platform name to publication-friendly name."""
    names = {
        "illumina": "Illumina",
        "pacbio": "PacBio HiFi",
        "nanopore": "Oxford Nanopore",
    }
    return names.get(platform or "illumina", "Illumina")
