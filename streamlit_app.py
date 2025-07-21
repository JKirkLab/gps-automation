import traceback
import streamlit as st
import pandas as pd
import sequence_extract
import align_sequence
import format_gps_entry
import process_output
import plot_utils
from uniprot_utils import fetch_all_sequences



def main():
    st.title("Downstream Output Grapher & Bounded Amino-Acid Region for Kinase prediction (DOGBARK)")
    docs_tab, results_tab = st.tabs(["Documentation", "Results"])
    with docs_tab:
        st.header("Documentation")
        st.write("There are several assumptions that the input files require in order for the app to function:")
        st.markdown("""
            - The following columns must be present in the excel file in exact case and spelling:
                - Modifications in Master Proteins
                - Master Protein Descriptions
            - Additionally, the file must be for singular amino acid modifications only (does not support multiple different AA modifications in the same file)
        """)
        st.markdown("""
            - Each entry in `Modifications in Master Proteins` is expected to follow this standardized format:

                ```
                Accession  #xPhospho [AApos1(confidence); AApos2(confidence); ...]
                ```

            - For example:

                ```
                Q62261  3xPhospho [S2315(97.6); S2318(100); S2322(97.6)]
                ```
            """)

        st.header("Usage")
        st.subheader("Preparing the Input")
        st.markdown("""
            1. Begin by uploading mass spectrometry input data.
            2. The results will appear in the results section with several intermediate log comments.
            3. A table will be displayed with the following columns:
                - `accession`: Accession of the protein.
                - `residue`: The modification that was present in the mass spec data.
                - `confidence`: The confidence in the detected modification.
                - `gene_name`: The gene name for the corresponding accession.
                - `sequence`: The full extracted sequence from UniProt for the accession.
                - `extracted_sequence`: The 21 AA sequence centered around the modification.
                - `center_index`: The position of the modification within the 21 AA sequence. This is primarily used for edge cases.
            4. A downloadable TXT file button will appear.
                - The TXT file contains output in the following format:
                    ```text
                    >gene
                    [extracted_sequence]
                    ```
            5. A downloadable Excel file button will appear.
                - This Excel file is the table described in step 3.
            """)
        st.subheader("Processing the Output")
        st.markdown("""
            1. Begin by uploading multiple csv files from the output of GPS.
                - This will aggregate all the csv files that are uploaded into one large table.
            2. The results will appear in the results section with several intermediate log comments. 
            3. A pie chart will be displayed for the top level of kinase groupings.
            4. A dropdown will appear that allows the user to select one of the primary kinase groupings
            5. A second pie chart will display the distribution within the selected primary kinase grouping.
            6. A table will be displayed with the following columns:
                - `Position`: Position of the modification within the extracted 21 AA sequence (usually position 11).
                - `Code`: The one letter representation of the AA.
                - `Kinase`: The full information on the predicted kinase.
                - `Peptide`: The 21 AA peptide sequence centered around the modification.
                - `Score`: The score assigned by GPS for this prediction.
                - `Cutoff`: I don't actually know what this represents...
                - `Gene`: The gene name for the protein.
                - `Kinase_Group`: The top level (primary) kinase prediction.
                - `Kinase_Subgroup`: The secondary level kinase prediction.
            7. A download button will appear to download the table described in 6. 
        """)
    st.markdown("---")

    with results_tab:

        st.header("Results")
        st.write("The results will appear here once a file is uploaded:")

        with st.sidebar:
            uploaded_file = st.file_uploader("Upload Mass Spec Excel File", type=["xlsx"])
            output_files = st.file_uploader("Upload one or multiple GPS Output File(s)", type=["csv"], accept_multiple_files=True)
        if uploaded_file:
            with st.expander("Mass Spec Input File Processing", expanded=True):
                try:
                    df = pd.read_excel(uploaded_file)
                    required_columns = {"Master Protein Descriptions", "Modifications in Master Proteins", "Annotated Sequence"}
                    if not required_columns.issubset(df.columns):
                        st.error(f"Input file must contain the following columns: {required_columns}")
                        return

                    df = df.dropna(subset=["Modifications in Master Proteins"]).copy()

                    all_entries = df.apply(sequence_extract.parse_modifications, axis=1)
                    flat_entries = [item for sublist in all_entries for item in sublist]
                    parsed_df = pd.DataFrame(flat_entries)
                    with st.spinner("Fetching sequences from UniProt..."):
                        complete_df, missing_fasta_dict, fasta_dict = fetch_all_sequences(parsed_df)

                    complete_df = complete_df.dropna(subset=["sequence"]).copy()
                    st.success("Sequences fetched successfully!")
                    if missing_fasta_dict:
                        st.warning("Some accessions were not found in the UniProt database. Please check the following:")
                        for accession, (new_accession, sequence) in missing_fasta_dict.items():
                            st.write(f"Obsolete Accession: {accession}, New Accession: {new_accession}")

                    st.info("Aligning peptide sequences...")
                    aligned_df = align_sequence.align_peptide_sequence(complete_df)
                    st.success("Peptide sequences aligned successfully!")

                    st.dataframe(aligned_df)
                    st.info("Generating GPS input format and csv file...")

                    gps_input = format_gps_entry.generate_gps_input(aligned_df)
                    st.download_button(
                        label="Download GPS Input",
                        data=gps_input,
                        file_name="gps_input.txt",
                        mime="text/plain",
                        icon=":material/download:"
                    )

                    excel_buffer = format_gps_entry.prepare_excel_download(aligned_df)

                    st.download_button(
                        label="Download Full Excel File",
                        data=excel_buffer,
                        file_name="full_data.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        icon=":material/download:"
                    )

                    st.success("GPS input format generated successfully!")
                except:
                    st.error("An error occurred while processing the file. Please ensure it is formatted correctly.")
                    st.text(traceback.format_exc())
        
        
        if output_files:
            with st.expander("GPS Output File Processing", expanded=True):
                aggregate_df = pd.DataFrame()
                st.info("Processing Output file(s)")
                for file in output_files:
                    try:
                        output = pd.read_csv(file)
                        processed_df = process_output.process_custom_csv(output)
                        processed_df = plot_utils.split_kinase_hierarchy(processed_df)
                        aggregate_df = pd.concat([aggregate_df, processed_df])
                        aggregate_df = aggregate_df.reset_index(drop=True)
                    except:
                        st.error("An error occured while processing the output. Please ensure it is formatted correctly.")
                        st.text(traceback.format_exc())

                st.success("Successfully Processed Output File!")
                unique_groups = aggregate_df["Kinase_Group"].dropna().unique()
                st.info("Plotting Kinase Distribution")
                plot_utils.plot_kinase_pie_chart(aggregate_df, group_col="Kinase_Group")

                selected_group = st.selectbox("Select Kinase Group to explore subfamilies:", sorted(unique_groups))
    
                filtered_df = aggregate_df[aggregate_df["Kinase_Group"] == selected_group]
                st.info(f"Subfamily distribution within {selected_group}")
                plot_utils.plot_kinase_pie_chart(filtered_df, group_col="Kinase_Subgroup")

                output_cleaned = format_gps_entry.prepare_excel_download(aggregate_df)

                st.dataframe(aggregate_df)

                st.download_button(
                    label="Download Processed Output Excel",
                    data=output_cleaned,
                    file_name="processed_output.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    icon=":material/download:"
                )

        
if __name__ == "__main__":
    main()
