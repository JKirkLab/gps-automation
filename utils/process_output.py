import pandas as pd
import streamlit as st

def process_custom_csv(df):
    """
    Processes the output of GPS(Group Based Prediction System) for kinase prediction

    Parameters:
    - df (pd.DataFrame): DataFrame containing Kinase prediction data

    Returns:
    - processed_df (pd.DataFrame): A cleaned dataframe where the genes have been exploded into its own column and filtered for modifications that match the relative center position.
    """
    current_gene = None
    output_rows = []

    for _, row in df.iterrows():
        first_col = str(row.iloc[0])
        if first_col.startswith(">"):
            header = first_col.lstrip(">")
            parts = header.split("|")
            current_gene = parts[0]
            center_str = parts[1].split('=')[1].strip()
            center = int(center_str)
        else:
            try:
                position = int(row['Position'])
                if position - 1 == center:
                    new_row = row.copy()
                    new_row["Gene"] = current_gene
                    output_rows.append(new_row)
            except (KeyError, ValueError, TypeError):
                continue 
    processed_df = pd.DataFrame(output_rows)
    processed_df.sort_values(by="Gene", inplace=True)
    processed_df.reset_index(drop=True, inplace=True)
    return processed_df


def filter_output(df, absolute_cutoff, relative_cutoff):
    if 'Score' in df.columns and 'Cutoff' in df.columns:
        df['abs_diff'] = df['Score'] - df['Cutoff']
        df['rel_diff'] = (df['Score'] - df['Cutoff']) / (1 - df['Cutoff'])
        return df[(df['abs_diff'] > absolute_cutoff) & (df['rel_diff'] > relative_cutoff)]
    else:
        st.error("No score and cutoff columns in final output!")
        return df
    

def filter_top_kinase_mod(df, k):
    df_sorted = df.sort_values("Score", ascending=False)
    kth_scores = df_sorted.groupby("Peptide")["Score"].transform(
        lambda x: x.nlargest(k).min()
    )
    return df_sorted[df_sorted["Score"] >= kth_scores]
