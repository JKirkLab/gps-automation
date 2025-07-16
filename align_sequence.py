import pandas as pd

def extract_surrounding_sequence(sequence, position):
    position -= 1
    l , r = max(position - 10, 0), min(position + 11, len(sequence))
    return sequence[l:r]

def align_peptide_sequence(df):
    df['extracted_sequence'] = df.apply(lambda row: extract_surrounding_sequence(row['sequence'], row['position']), axis=1)
    return df
