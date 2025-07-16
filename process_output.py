import pandas as pd

def process_custom_csv(df):
    current_gene = None
    output_rows = []

    for _, row in df.iterrows():
        first_col = str(row.iloc[0])
        if first_col.startswith(">"):
            current_gene = first_col.lstrip(">")
        else:
            try:
                position = int(row['Position'])
                if position == 11:
                    new_row = row.copy()
                    new_row["Gene"] = current_gene
                    output_rows.append(new_row)
            except (KeyError, ValueError, TypeError):
                continue 
    processed_df = pd.DataFrame(output_rows)
    processed_df.sort_values(by="Gene", inplace=True)
    processed_df.reset_index(drop=True, inplace=True)
    return processed_df