"""
Load GSE275334 (immune exhaustion in ME/CFS and Long COVID) into tidy tables.

Source: Gene Expression Omnibus, accession GSE275334.
Assay:  NanoString nCounter panel on PBMCs, 635 genes x 47 samples.
Groups: 18 healthy control, 15 Long COVID, 14 ME/CFS.

Produces, in data/:
  samples.csv     one row per sample: id, group, age, sex, QC metrics
  expression.csv  genes x samples, normalised counts
  expr_log2.csv   log2(normalised + 1), the matrix used for PCA

Run:  python src/load_data.py
"""

from pathlib import Path
import gzip
import re
import subprocess
import urllib.request

import numpy as np
import pandas as pd

GSE = "GSE275334"
DATA = Path(__file__).resolve().parent.parent / "data"

XLSX_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    f"?acc={GSE}&format=file&file={GSE}%5FFile%5F1%5FNormalised%2Exlsx"
)
MATRIX_URL = (
    f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GSE[:-3]}nnn/{GSE}"
    f"/matrix/{GSE}_series_matrix.txt.gz"
)

# The expression sheet labels samples HC01/LC01/ME01; the GEO series matrix
# labels the same samples "Healthy Control #1" etc. Map one to the other.
GROUP_PREFIX = {"Healthy Control": "HC", "Long COVID": "LC", "ME/CFS": "ME"}


def download(url: str, dest: Path) -> Path:
    """Fetch url to dest if not already present."""
    if dest.exists():
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.suffix == ".xlsx":
        # The FTP path 403s for this file; the GEO download endpoint serves it,
        # but redirects, so follow with curl -L.
        subprocess.run(["curl", "-sSL", "-o", str(dest), url], check=True)
    else:
        urllib.request.urlretrieve(url, dest)
    return dest


def parse_series_matrix(path: Path) -> pd.DataFrame:
    """Pull per-sample metadata out of a GEO series matrix file.

    Each !Sample_characteristics_ch1 line is one field across all samples,
    formatted "field: value". Returns one row per sample.
    """
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            key, *values = line.rstrip("\n").split("\t")
            values = [v.strip('"') for v in values]
            if key == "!Sample_characteristics_ch1":
                name = values[0].split(":")[0].strip()
                fields[name] = [
                    v.split(": ", 1)[1] if ": " in v else v for v in values
                ]
            elif key in ("!Sample_title", "!Sample_geo_accession"):
                fields[key.replace("!Sample_", "")] = values
    return pd.DataFrame(fields)


def sample_id_from_title(title: str) -> str:
    """'Long COVID #7' -> 'LC07'."""
    label, number = title.rsplit("#", 1)
    return f"{GROUP_PREFIX[label.strip()]}{int(number):02d}"


def main() -> None:
    xlsx = download(XLSX_URL, DATA / f"{GSE}_normalised.xlsx")
    matrix = download(MATRIX_URL, DATA / f"{GSE}_series_matrix.txt.gz")

    sheets = pd.ExcelFile(xlsx)
    norm = sheets.parse("Normalised counts")
    qc = sheets.parse("Sample QC")

    gene_cols = ["Symbol", "GeneID", "Alias", "Description"]
    expr = norm.set_index("Symbol").drop(columns=[c for c in gene_cols if c != "Symbol"])

    # --- sample table -----------------------------------------------------
    meta = parse_series_matrix(matrix)
    meta["sample_id"] = meta["title"].map(sample_id_from_title)
    meta["group"] = meta["sample_id"].str[:2].map(
        {"HC": "Healthy control", "LC": "Long COVID", "ME": "ME/CFS"}
    )
    meta["age"] = pd.to_numeric(meta["age (years)"], errors="coerce")
    meta["sex"] = meta["Sex"]

    qc_keep = {
        "Sample name": "sample_id",
        "Imaging quality (%FOV captured)": "imaging_quality_pct",
        "Binding density": "binding_density",
        "housekeepingGeneScalingFactor": "housekeeping_scaling_factor",
        "QC Issues": "qc_issues",
    }
    qc_small = qc[list(qc_keep)].rename(columns=qc_keep)

    samples = (
        meta[["sample_id", "geo_accession", "group", "age", "sex"]]
        .merge(qc_small, on="sample_id", how="left", validate="one_to_one")
        .sort_values(["group", "sample_id"])
        .reset_index(drop=True)
    )

    # --- checks -----------------------------------------------------------
    # The expression columns are not in sorted order (LC05 precedes LC01), so
    # every downstream step must align on sample_id rather than on position.
    assert set(expr.columns) == set(samples.sample_id), (
        "expression columns and sample table disagree: "
        f"{set(expr.columns) ^ set(samples.sample_id)}"
    )
    assert samples.sample_id.is_unique
    assert not samples.group.isna().any()
    assert expr.notna().all().all(), "unexpected missing values in expression matrix"
    assert len(samples) == 47 and expr.shape[0] == 635

    expr = expr[samples.sample_id]  # enforce the shared column order
    expr_log2 = np.log2(expr + 1)

    samples.to_csv(DATA / "samples.csv", index=False)
    expr.to_csv(DATA / "expression.csv")
    expr_log2.to_csv(DATA / "expr_log2.csv")

    print(f"genes x samples: {expr.shape}")
    print(samples.group.value_counts().to_string())
    print("\nage by group:")
    print(samples.groupby("group")["age"].agg(["median", "min", "max"]).to_string())
    print("\nsex by group:")
    print(pd.crosstab(samples.group, samples.sex).to_string())
    print("\nsamples flagged by the vendor QC sheet:")
    flagged = samples.qc_issues.notna() & (samples.qc_issues.astype(str) != "nan")
    print(f"  {int(flagged.sum())} of {len(samples)}")


if __name__ == "__main__":
    main()
