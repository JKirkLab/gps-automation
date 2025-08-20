import matplotlib.pyplot as plt
import streamlit as st
import numpy as np
import plotly.express as px
def plot_kinase_pie_chart(df, group_col, kinase_column="Kinase", pct=False, legend=False):
    """
    Plots a pie chart of kinase classifications from a DataFrame and displays it in Streamlit.
    """

    if group_col not in df.columns:
        st.warning(f"Column '{group_col}' not found in the data.")
        return

    df_filtered = df.dropna(subset=[group_col]).copy()
    group_counts = df_filtered[group_col].value_counts().reset_index()
    group_counts.columns = [group_col, "count"]
    n = len(df_filtered[group_col])

    if group_counts.empty:
        st.warning("No data available for the selected kinase group.")
        return

    fig = px.pie(
        group_counts,
        values="count",
        names=group_col,
        title=f"{group_col} Distribution (n={n})",
    )

    fig.update_traces(textinfo="none")
    fig.update_layout(showlegend=legend, legend_title_text=group_col)

    hover = f"{group_col}=%{{label}}<br>count=%{{value:,.0f}}"
    if pct:
        hover += "<br>percent=%{percent:.2%}"
    hover += "<extra></extra>"
    fig.update_traces(hovertemplate=hover)

    st.subheader(f"{group_col} Distribution")
    st.plotly_chart(fig)
    # fig, ax = plt.subplots()
    # wedges, texts, autotexts = ax.pie(
    #     group_counts,
    #     labels=None,
    #     autopct='%1.1f%%',
    #     startangle=100
    # )
    # plt.text(0.95, 0.9, f'n={n}', transform=ax.transAxes,
    #      fontsize=12, bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))
    # ax.axis('equal')
    # # ax.legend(wedges, group_counts.index, title=group_col, loc="center left", bbox_to_anchor=(1, 0.5))
    # plt.tight_layout()
    
    # st.subheader(f"{group_col} Distribution")
    # st.pyplot(fig)

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


def percent_contour(df, abs_col='abs_diff', rel_col='rel_diff',
                    score_col='Score', cutoff_col='Cutoff',
                    n_abs_bins=200, n_rel_bins=200, clip_rel=None,
                    levels=(25, 50, 75)):
    """
    Make a contour plot where lines correspond to given % of rows surviving.
    Default levels are 25%, 50%, 75%.
    """
    # --- get abs_diff / rel_diff arrays ---
    if abs_col in df.columns and rel_col in df.columns:
        x = df[abs_col].to_numpy()
        y = df[rel_col].to_numpy()
    elif score_col in df.columns and cutoff_col in df.columns:
        den = (1 - df[cutoff_col]).replace(0, np.nan)
        x = (df[score_col] - df[cutoff_col]).to_numpy()
        y = ((df[score_col] - df[cutoff_col]) / den).to_numpy()
        y = np.nan_to_num(y, nan=-np.inf)
    else:
        raise ValueError("Need (abs_diff & rel_diff) or (Score & Cutoff) columns.")

    if clip_rel is not None:
        y = np.clip(y, clip_rel[0], clip_rel[1])

    # --- 2D histogram ---
    x_edges = np.linspace(x.min(), x.max(), n_abs_bins+1)
    y_edges = np.linspace(y.min(), y.max(), n_rel_bins+1)
    H, xe, ye = np.histogram2d(x, y, bins=[x_edges, y_edges])

    # --- suffix cumulative counts ---
    H_rev = H[::-1, ::-1]
    C_rev = H_rev.cumsum(0).cumsum(1)
    C = C_rev[::-1, ::-1]

    # --- % surface ---
    Z = 100.0 * C / max(1, len(x))
    X, Y = np.meshgrid(xe[:-1], ye[:-1], indexing='ij')

    # --- plot contours for given % levels ---
    fig, ax = plt.subplots()
    cs = ax.contour(X, Y, Z, levels=levels)
    ax.clabel(cs, inline=True, fmt="%.0f%%")
    ax.set_xlabel("absolute_cutoff")
    ax.set_ylabel("relative_cutoff")
    ax.set_title("Rows surviving filter (% contours)")
    return fig