import matplotlib.pyplot as plt
import streamlit as st

def plot_kinase_pie_chart(df, group_col, kinase_column="Kinase"):
    """
    Plots a pie chart of kinase classifications from a DataFrame and displays it in Streamlit.

    Parameters:
    - df (pd.DataFrame): DataFrame containing Kinase prediction data
    - group_col (str): Column name to use for grouping (e.g., 'Kinase_Group', 'Kinase_Subgroup')
    - kinase_column (str): (Optional) Original column containing raw kinase annotations (default: "Kinase")

    Returns:
    - None, displays chart in streamlit
    """
    if group_col not in df.columns:
        st.warning(f"Column '{group_col}' not found in the data.")
        return

    df_filtered = df.dropna(subset=[group_col]).copy()
    group_counts = df_filtered[group_col].value_counts()
    n = len(df_filtered[group_col])
    if group_counts.empty:
        st.warning("No data available for the selected kinase group.")
        return

    fig, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie(
        group_counts,
        labels=None,
        autopct='%1.1f%%',
        startangle=100
    )
    plt.text(0.95, 0.9, f'n={n}', transform=ax.transAxes,
         fontsize=12, bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))
    ax.axis('equal')
    ax.legend(wedges, group_counts.index, title=group_col, loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    
    st.subheader(f"{group_col} Distribution")
    st.pyplot(fig)

def split_kinase_hierarchy(df):
    """
    Splits the primary group and subgroup for the predicted kinase. 
    
    Parameters:
    - df (pd.DataFrame): Dataframe containing predicted kinase data

    Returns:
    - df (pd.DataFrame): Dataframe containing primary group and subgroup kinase data
    """
    df = df.copy()
    df['Kinase_Group'] = df['Kinase'].str.split('/').str[0]
    df['Kinase_Subgroup'] = df['Kinase'].str.split("/").str[1]
    return df