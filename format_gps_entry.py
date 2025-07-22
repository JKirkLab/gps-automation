import streamlit as st
import pandas as pd
from io import StringIO, BytesIO
import csv
def generate_gps_input(df):
    """
    Generates GPS input format from mass spectrometry data.

    The output follows the format required by GPS (Group-based Prediction System):
        >gene|Center = #
        [extracted_sequence]

    Where:
    - `gene` is the gene or protein identifier
    - `Center` indicates the zero-based index of the modification site relative to the extracted sequence
    - `[extracted_sequence]` is the 21-mer amino acid sequence window surrounding the modification site

    Returns:
    - A string or file content suitable for input to GPS prediction tools.
    """
    fasta_buffer = StringIO()
    for _, row in df.iterrows():
        fasta_buffer.write(f">{row['gene_name']}|Center = {row['center_index']}\n")
        fasta_buffer.write(f"{row['extracted_sequence']}\n")
    return fasta_buffer.getvalue()

def prepare_excel_download(df, sheet_name="Sheet1"):
    """
    Helper function to allow dataframes to be downloaded in streamlit.
    """
    output = BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
    output.seek(0)
    return output