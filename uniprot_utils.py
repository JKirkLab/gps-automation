import pandas as pd
import requests
import re
def chunk_list(lst, size):
    return [lst[i:i+size] for i in range(0, len(lst), size)]

def extract_accession(header):
    match = re.match(r'(?:\w+\|)?(\w+)\|?', header)
    return match.group(1) if match else None

def parse_fasta_entry(entry):
    lines = entry.splitlines()
    if not lines:
        return None, None
    header = lines[0]
    sequence = "".join(lines[1:])
    accession = extract_accession(header)
    return accession, sequence


def query_full_seq(df):
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
    returned_accessions = []

    for entry in entries:
        if not entry:
            continue
        accession, _sequence = parse_fasta_entry(entry)
        returned_accessions.append(accession)
    missing_in_fasta = set(unique_accessions) - set(returned_accessions)
    print("Missing obsolete accesions:", missing_in_fasta)

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
    all_fasta_data = query_full_seq(original_df)
    fasta_dict, entries  =  process_fasta_data(all_fasta_data)
    missing_fasta_dict = {}

    fasta_keys = set(fasta_dict.keys())
    original_keys = set(original_df['accession'].unique())
    missing = original_keys - fasta_keys
    if missing:
        print("different number of accessions in fasta and original df")
        missing_fasta_dict = req_obsolete_accessions(original_df['accession'].unique(), entries)
        # for k, v in missing_fasta_dict.items():
        #     fasta_dict[k] = v[1]
    original_df['sequence'] = original_df['accession'].map(fasta_dict)
    return original_df, missing_fasta_dict, fasta_dict

    



