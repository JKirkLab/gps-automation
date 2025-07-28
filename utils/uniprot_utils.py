import pandas as pd
import requests
import re
def chunk_list(lst, size):
    """
    Helper function to chunk arrays into workable sizes for uniprot.
    """
    return [lst[i:i+size] for i in range(0, len(lst), size)]

def extract_accession(header):
    """
    Helper function to extract protein accession from a header.
    """
    match = re.match(r'(?:\w+\|)?(\w+)\|?', header)
    return match.group(1) if match else None

def parse_fasta_entry(entry):
    """
    Helper function to parse a FASTA entry into its accession and sequence.

    Parameters:
    - entry (str): FASTA-formatted string starting with a header line

    Returns:
    - accession (str or None): Extracted UniProt accession ID
    - sequence (str or None): Full amino acid sequence (no line breaks)
    """
    lines = entry.splitlines()
    if not lines:
        return None, None
    header = lines[0]
    sequence = "".join(lines[1:])
    accession = extract_accession(header)
    return accession, sequence


def query_full_seq(df):
    """
    Queries UniProt’s REST API to fetch full amino acid sequences for accessions
    found in the input DataFrame. Requests are made in chunks of 100 accessions.

    Parameters:
    - df (pd.DataFrame): DataFrame containing a column named 'accession'

    Returns:
    - all_fasta_data (List[str]): List of FASTA-formatted response strings
    """
    accession_list = df['accession'].to_list()
    unique_accessions = list(set(accession_list))
    
    query_string = ""
    chunks = chunk_list(unique_accessions, 100)

    all_fasta_data = []

    for chunk in chunks: 
        query_string = " OR ".join(f"accession:{acc}" for acc in chunk)
        url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {
        "format": "fasta",
        "query": query_string 
        }
        response = requests.get(url, params=params)
        if response.status_code == 200:
            fasta_data = response.text
            all_fasta_data.append(fasta_data)
        else:
            print("request failed with status code:", response.status_code)

    return all_fasta_data

def process_fasta_data(fasta_data):
    """
    Parses a list of FASTA-formatted strings into a dictionary of accessions to sequences,
    and also returns the raw FASTA entries for downstream analysis.

    Parameters:
    - fasta_data (List[str]): List of FASTA responses (strings) from UniProt

    Returns:
    - fasta_dict (dict): Mapping from accession → sequence
    - entries (List[str]): List of raw FASTA entry strings (split on '>')
    """
    fasta_text = "\n".join(fasta_data)
    entries = fasta_text.strip().split(">")

    fasta_dict = {}

    for entry in entries:
        if not entry:
            continue
        accession, sequence = parse_fasta_entry(entry)
        if accession and sequence:
            fasta_dict[accession] = sequence
    
    return fasta_dict, entries


def req_obsolete_accessions(unique_accessions, entries):
    """
    Identifies accessions missing from initial FASTA results and performs fallback requests
    to fetch them individually (e.g., for obsolete or redirected UniProt entries).

    Parameters:
    - unique_accessions (List[str]): List of originally requested accessions
    - entries (List[str]): FASTA entries returned from the initial query

    Returns:
    - missing_fasta_dict (dict): Mapping of missing accession → (new_accession, sequence) tuples
    """
    returned_accessions = []

    for entry in entries:
        if not entry:
            continue
        accession, _sequence = parse_fasta_entry(entry)
        returned_accessions.append(accession)
    missing_in_fasta = set(unique_accessions) - set(returned_accessions)
    missing_fasta_data = []
    missing_fasta_dict = {}
    for accession in list(missing_in_fasta):
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession}.fasta", allow_redirects=True)
        if response.status_code == 200:
            entry = response.text.strip().lstrip(">")
            new_accession, sequence = parse_fasta_entry(entry)
            if new_accession and sequence:
                missing_fasta_dict[accession] = (new_accession, sequence)
        else:
            print("request failed with status code:", response.status_code)
    return missing_fasta_dict

def fetch_all_sequences(original_df):
    """
    Fetches full amino acid sequences for all UniProt accessions in the input DataFrame.

    This function performs the following:
    - Queries UniProt in batches (max 100) for accession sequences
    - Parses FASTA-formatted responses into accession-sequence mappings
    - Identifies and resolves missing/obsolete accessions
    - Appends the retrieved sequences as a new column in the original DataFrame

    Parameters:
    - original_df (pd.DataFrame): DataFrame with an 'accession' column containing UniProt IDs

    Returns:
    - updated_df (pd.DataFrame): Input DataFrame with an added 'sequence' column
    - missing_fasta_dict (dict): Dictionary of accessions requiring fallback queries,
                                 mapping to (new_accession, sequence)
    - fasta_dict (dict): All successfully retrieved accession to sequence mappings
    """
    all_fasta_data = query_full_seq(original_df)
    fasta_dict, entries  =  process_fasta_data(all_fasta_data)
    missing_fasta_dict = {}

    fasta_keys = set(fasta_dict.keys())
    original_keys = set(original_df['accession'].unique())
    missing = original_keys - fasta_keys
    if missing:
        missing_fasta_dict = req_obsolete_accessions(original_df['accession'].unique(), entries)
        # for k, v in missing_fasta_dict.items():
        #     fasta_dict[k] = v[1]
    original_df['sequence'] = original_df['accession'].map(fasta_dict)
    return original_df, missing_fasta_dict, fasta_dict

    



