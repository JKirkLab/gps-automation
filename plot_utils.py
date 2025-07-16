import matplotlib.pyplot as plt
import streamlit as st

def plot_kinase_pie_chart(df, kinase_column="Kinase"):
    """
    Plots a pie chart of top-level kinase groups from a DataFrame and displays it in Streamlit.

    Parameters:
    - df (pd.DataFrame): DataFrame containing a 'Kinase' column with kinase groupings (e.g., "CAMK/RAD53")
    - kinase_column (str): Column name to extract kinase groups from

    Returns:
    - None (displays chart directly in Streamlit)
    """
    df_filtered = df.dropna(subset=[kinase_column]).copy()
    df_filtered['Kinase Group'] = df_filtered[kinase_column].str.split("/").str[0]

    group_counts = df_filtered['Kinase Group'].value_counts()

    if group_counts.empty:
        st.warning("No kinase data found for plotting.")
        return

    fig, ax = plt.subplots()
    wedges, texts, autotexts = ax.pie(
        group_counts,
        labels=None,  
        autopct='%1.1f%%',
        startangle = 100 
    )
    
    ax.axis('equal')  # Equal aspect ratio to ensure it's a circle
    ax.legend(wedges, group_counts.index, title="Kinase Groups", loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    st.subheader("Kinase Group Distribution")
    st.pyplot(fig)
