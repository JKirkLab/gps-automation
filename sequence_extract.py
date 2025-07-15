import pandas as pd
import re

def parse_modifications(row):
    """Extracts phospho site info from one row."""
    mod_str = row["Modifications in Master Proteins"]
    desc_str = row["Master Protein Descriptions"]

    if not isinstance(mod_str, str) or not isinstance(desc_str, str):
        return []

    match = re.match(r"^(?P<accession>\w+)\s\d+xPhospho\s+\[(?P<sites>[^\]]+)\]", mod_str)
    if not match:
        return []

    accession = match.group("accession")
    sites_str = match.group("sites")

    gene_match = re.search(r"GN=([\w\-\.]+)", desc_str)
    gene_name = gene_match.group(1) if gene_match else "gene"

    parsed = []
    for site in sites_str.split(";"):
        site = site.strip()
        site_match = re.match(r"(?P<residue>[A-Z])(?P<position>\d+)\((?P<conf>[\d.]+)\)", site)
        if site_match:
            parsed.append({
                "accession": accession,
                "residue": site_match.group("residue"),
                "position": int(site_match.group("position")),
                "confidence": float(site_match.group("conf")),
                "gene_name": gene_name
            })

    return parsed


def generate_cleaned_df(all_entries):
    flat_entries = [item for sublist in all_entries for item in sublist]
    parsed_df = pd.DataFrame(flat_entries)
    return parsed_df
    
