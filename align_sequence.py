import pandas as pd

def extract_surrounding_sequence(sequence, position):
    """
    Helper function that extracts the 21 AA chain surrounding the specified post translational modification.
    
    Parameters:
    - sequence (str): Sequence of the full protein.
    - position (int): One-based index of the position of the post translational modification.
    Returns:
    - extracted_sequence (str): String containing the 21 AA chain centered around the post translational modification.
    - relative_pos (int): The relative position of the modification within the 21-mer (or shorter if an edge case).
    """
    position -= 1
    l , r = max(position - 10, 0), min(position + 11, len(sequence))
    extracted = sequence[l:r]
    relative_pos = position - l
    return extracted, relative_pos

def align_peptide_sequence(df):
    """
    Aligns the peptide sequences around the position of the modified residue.
    
    Parameters:
    - df (pd.DataFrame): Dataframe containing residue modification information.

    Returns:
    - df (pd.DataFrame): Dataframe containing the new extracted sequence as a column.
    """

    def apply_fn(row):
        extracted, rel_pos = extract_surrounding_sequence(row['sequence'], row['position'])
        return pd.Series({'extracted_sequence': extracted, 'center_index': rel_pos})

    df_out = df.copy()
    df_out[['extracted_sequence', 'center_index']] = df_out.apply(apply_fn, axis=1)
    return df_out
