import streamlit as st
import pandas as pd
from uniprot_utils import fetch_all_sequences
import sequence_extract
import align_sequence
import format_gps_entry
import traceback

def main():
    st.title("Protein Sequence Extractor")
    st.write("Upload an Excel file with UniProt accessions and modification positions.")

    uploaded_file = st.file_uploader("Upload Excel File", type=["xlsx"])

    if uploaded_file:
        try:
            df = pd.read_excel(uploaded_file)
            required_columns = {"Master Protein Descriptions", "Modifications in Master Proteins", "Annotated Sequence"}
            if not required_columns.issubset(df.columns):
                st.error(f"Input file must contain the following columns: {required_columns}")
                return

            df = df.dropna(subset=["Modifications in Master Proteins"])

            all_entries = df.apply(sequence_extract.parse_modifications, axis=1)
            flat_entries = [item for sublist in all_entries for item in sublist]
            parsed_df = pd.DataFrame(flat_entries)
            with st.spinner("Fetching sequences from UniProt..."):
                complete_df, missing_fasta_dict, fasta_dict = fetch_all_sequences(parsed_df)
            
            complete_df = complete_df.dropna(subset=["sequence"])
            st.success("Sequences fetched successfully!")
            if missing_fasta_dict:
                st.warning("Some accessions were not found in the UniProt database. Please check the following:")
                for accession, (new_accession, sequence) in missing_fasta_dict.items():
                    st.write(f"Accession: {accession}, New Accession: {new_accession}")
                
            st.info("Aligning peptide sequences...")
            aligned_df = align_sequence.align_peptide_sequence(complete_df)
            st.success("Peptide sequences aligned successfully!")

            st.dataframe(complete_df)

            st.info("Generating GPS input format and csv file...")

            
            gps_input = format_gps_entry.generate_gps_input(aligned_df)
            st.download_button(
                label="Download GPS Input",
                data=gps_input,
                file_name="gps_input.txt",
                mime="text/plain"
            )

            excel_buffer = format_gps_entry.prepare_excel_download(complete_df)

            st.download_button(
                label="Download Full Excel File",
                data=excel_buffer,
                file_name="full_data.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )

            st.success("GPS input format generated successfully!")
        except:
            st.error("An error occurred while processing the file. Please ensure it is formatted correctly.")
            st.text(traceback.format_exc())
if __name__ == "__main__":
    main()
