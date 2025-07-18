import pandas as pd

def extract_surrounding_sequence(sequence, position):
    """
    Helper function that extracts the 21 AA chain surrounding the specified post translational modification.
    
    Parameters:
    - sequence (str): Sequence of the full protein.
    - position (int): One-based index of the position of the post translational modification.
    Returns:
    - extracted_sequence (str): String containing the 21 AA chain centered around the post translational modification.
    """
    position -= 1
    l , r = max(position - 10, 0), min(position + 11, len(sequence))
    return sequence[l:r]

def align_peptide_sequence(df):
    """
    Aligns the peptide sequences around the position of the modified residue.
    
    Parameters:
    - df (pd.DataFrame): Dataframe containing residue modification information.

    Returns:
    - df (pd.DataFrame): Dataframe containing the new extracted sequence as a column.
    """
    df['extracted_sequence'] = df.apply(lambda row: extract_surrounding_sequence(row['sequence'], row['position']), axis=1)
    return df
