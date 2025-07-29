import pandas as pd
import re



def parse_modifications(row):
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
            confidence = float(site_match.group("conf"))
        else:
            site_match = re.match(r"(?P<residue>[A-Z])(?P<position>\d+)", site)
            if site_match:
                confidence = None
            else:
                continue 

        parsed.append({
            "accession": accession,
            "residue": site_match.group("residue"),
            "position": int(site_match.group("position")),
            "confidence": confidence,
            "gene_name": gene_name
        })

    return parsed

# def parse_modifications(row):
#     """
#     Extract phospho-site information from a single DataFrame row.

#     This function parses the "Modifications in Master Proteins" and "Master Protein Descriptions"
#     columns of a proteomics dataset to extract phosphorylation site details, including residue, 
#     position, localization confidence, and gene name.

#     Parameters:
#     row : pd.Series
#         A single row from a DataFrame containing at least the columns:
#         - "Modifications in Master Proteins"
#         - "Master Protein Descriptions"

#     Returns:
#     list of dict
#         A list of dictionaries, each containing:
#         - 'accession' : str – the protein accession ID
#         - 'residue'   : str – the amino acid residue (e.g., 'S', 'T', or 'Y')
#         - 'position'  : int – the site position on the protein
#         - 'confidence': float – localization confidence (e.g., 0.99)
#         - 'gene_name' : str – gene name (if available), otherwise "gene"
    
#     Notes:
#     - Returns an empty list if the modification or description fields are missing or malformed.
#     - Assumes phosphorylation entries are in the format:
#       "ACCESSION #xPhospho [S123(0.95); T456(0.87)]"
#     """
#     mod_str = row["Modifications in Master Proteins"]
#     desc_str = row["Master Protein Descriptions"]

#     if not isinstance(mod_str, str) or not isinstance(desc_str, str):
#         return []

#     match = re.match(r"^(?P<accession>\w+)\s\d+xPhospho\s+\[(?P<sites>[^\]]+)\]", mod_str)
#     if not match:
#         return []

#     accession = match.group("accession")
#     sites_str = match.group("sites")

#     gene_match = re.search(r"GN=([\w\-\.]+)", desc_str)
#     gene_name = gene_match.group(1) if gene_match else "gene"

#     parsed = []
#     for site in sites_str.split(";"):
#         site = site.strip()
#         site_match = re.match(r"(?P<residue>[A-Z])(?P<position>\d+)\((?P<conf>[\d.]+)\)", site)
#         if site_match:
#             parsed.append({
#                 "accession": accession,
#                 "residue": site_match.group("residue"),
#                 "position": int(site_match.group("position")),
#                 "confidence": float(site_match.group("conf")),
#                 "gene_name": gene_name
#             })

#     return parsed


def generate_cleaned_df(all_entries):
    """
    Flattens a nested list of parsed phospho-site entries into a clean DataFrame.

    This utility function takes the parsed output from multiple rows (e.g., from
    `parse_modifications`) and consolidates them into a single pandas DataFrame.

    Parameters:
    - all_entries : list of list of dict
        A list where each element is a list of parsed modification dicts 
        (as returned by `parse_modifications`).

    Returns:
    - pd.DataFrame
        A flat DataFrame containing one row per phospho-site entry, with columns:
        - 'accession'
        - 'residue'
        - 'position'
        - 'confidence'
        - 'gene_name'
    """
    flat_entries = [item for sublist in all_entries for item in sublist]
    parsed_df = pd.DataFrame(flat_entries)
    return parsed_df
    
