import streamlit as st
import pandas as pd
from io import StringIO, BytesIO
import csv
def generate_gps_input(df):
    fasta_buffer = StringIO()
    for _, row in df.iterrows():
        fasta_buffer.write(f">{row['gene_name']}\n")
        fasta_buffer.write(f"{row['sequence']}\n")
    return fasta_buffer.getvalue()

def prepare_excel_download(df, sheet_name="Sheet1"):
    output = BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False, sheet_name=sheet_name)
    output.seek(0)
    return output