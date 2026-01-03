# mcpath/interactive_plots.py
import os
import glob
from typing import Dict, List, Optional, Tuple

import pandas as pd


# ------------------------- Plotly setup -------------------------
def _ensure_plotly_renderer():
    """
    Sets a renderer that works well in Colab / notebooks.
    Safe to call multiple times.
    """
    try:
        import plotly.io as pio
        # Prefer Colab renderer if available
        if "colab" in pio.renderers:
            pio.renderers.default = "colab"
        elif "notebook_connected" in pio.renderers:
            pio.renderers.default = "notebook_connected"
        elif "notebook" in pio.renderers:
            pio.renderers.default = "notebook"
    except Exception:
        # If plotly isn't installed or renderer can't be set, we fail later with a clearer message.
        pass


def _display_fig(fig):
    """
    Display plotly figure in a way that behaves well inside ipywidgets.Output.
    """
    try:
        from IPython.display import display
        display(fig)
        return
    except Exception:
        pass
    try:
        fig.show()
    except Exception:
        pass


# ------------------------- File discovery -------------------------
def find_centrality_files(work_dir: str) -> Dict[str, List[str]]:
    """
    Find candidate closeness/betweenness output files in work_dir.
    Returns a dict: {"closeness":[...], "betweenness":[...]}
    """
    work_dir = os.path.abspath(work_dir)
    if not os.path.isdir(work_dir):
        raise FileNotFoundError(f"work_dir not found: {work_dir}")

    all_files = []
    for ext in ("*.txt", "*.dat", "*.out", "*.csv", "*"):
        all_files += glob.glob(os.path.join(work_dir, ext))

    # Keep only regular files
    all_files = [f for f in all_files if os.path.isfile(f)]

    def _score(name: str) -> Tuple[int, int]:
        """
        Higher score -> more likely.
        Primary: keyword hits
        Secondary: shorter names (less noisy)
        """
        n = name.lower()
        hit = 0
        hit += 3 if "closeness" in n else 0
        hit += 2 if "close" in n else 0
        hit += 3 if "betweenness" in n else 0
        hit += 2 if "between" in n else 0
        hit += 1 if "central" in n else 0
        return (hit, -len(n))

    close_cands = []
    betw_cands = []
    for f in all_files:
        base = os.path.basename(f).lower()
        if any(k in base for k in ["closeness", "close"]):
            close_cands.append(f)
        if any(k in base for k in ["betweenness", "between", "betw"]):
            betw_cands.append(f)

    # Sort by heuristic score, then newest
    def _sort(files: List[str]) -> List[str]:
        files = sorted(files, key=lambda f: (_score(os.path.basename(f)), os.path.getmtime(f)), reverse=True)
        # de-dup
        out = []
        seen = set()
        for f in files:
            if f not in seen:
                seen.add(f)
                out.append(f)
        return out

    return {
        "closeness": _sort(close_cands),
        "betweenness": _sort(betw_cands),
    }


# ------------------------- Parsing -------------------------
def load_metric_table(path: str) -> pd.DataFrame:
    """
    Robustly load a metric file and return a dataframe with:
      residue_index (int), value (float)
    Accepts whitespace-delimited or CSV-ish inputs.
    Drops non-numeric header rows automatically.
    """
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

    # Try whitespace first
    df = None
    read_errors = []

    for kwargs in [
        dict(sep=r"\s+", engine="python", header=None, comment="#"),
        dict(sep=",", engine="python", header=None, comment="#"),
        dict(sep=None, engine="python", header=None, comment="#"),  # auto-sep
    ]:
        try:
            df = pd.read_csv(path, **kwargs)
            if df is not None and df.shape[1] >= 2:
                break
        except Exception as e:
            read_errors.append(str(e))
            df = None

    if df is None or df.shape[1] < 2:
        raise ValueError(
            f"Could not parse at least 2 columns from {os.path.basename(path)}.\n"
            f"Tried multiple parsers. Last errors: {read_errors[-2:]}"
        )

    # Take first two columns, coerce numeric
    out = df.iloc[:, :2].copy()
    out.columns = ["residue_index", "value"]
    out["residue_index"] = pd.to_numeric(out["residue_index"], errors="coerce")
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    out = out.dropna(subset=["residue_index", "value"]).copy()

    # Make residue_index int-like if possible
    out["residue_index"] = out["residue_index"].round().astype("int64")

    return out


# ------------------------- Plotting -------------------------
def plot_metric(df: pd.DataFrame, title: str, file_label: Optional[str] = None):
    """
    Interactive plotly line+markers with hover.
    """
    _ensure_plotly_renderer()
    try:
        import plotly.express as px
    except Exception as e:
        raise ImportError(
            "Plotly is not available. In Colab run: !pip -q install plotly"
        ) from e

    hover = ["residue_index", "value"]
    if file_label:
        df = df.copy()
        df["file"] = file_label
        hover.append("file")

    fig = px.scatter(
        df,
        x="residue_index",
        y="value",
        hover_data=hover,
        title=title
    )
    fig.update_traces(mode="lines+markers")
    fig.update_layout(xaxis_title="Residue index", yaxis_title="Centrality value")
    _display_fig(fig)
    return fig


def plot_centrality_results(
    work_dir: str,
    closeness_file: Optional[str] = None,
    betweenness_file: Optional[str] = None,
    take_top: int = 1,
) -> Dict[str, List]:
    """
    Plot closeness and betweenness interactively.
    - If explicit files are provided, uses them.
    - Otherwise, discovers candidates in work_dir and plots up to `take_top` for each metric.

    Returns: {"closeness":[figs...], "betweenness":[figs...]}
    """
    figs = {"closeness": [], "betweenness": []}

    if closeness_file or betweenness_file:
        if closeness_file:
            dfc = load_metric_table(closeness_file)
            figs["closeness"].append(
                plot_metric(dfc, "Closeness (interactive)", os.path.basename(closeness_file))
            )
        if betweenness_file:
            dfb = load_metric_table(betweenness_file)
            figs["betweenness"].append(
                plot_metric(dfb, "Betweenness (interactive)", os.path.basename(betweenness_file))
            )
        return figs

    found = find_centrality_files(work_dir)
    for metric in ["closeness", "betweenness"]:
        files = found.get(metric, [])[: max(1, int(take_top))]
        for f in files:
            try:
                df = load_metric_table(f)
                figs[metric].append(
                    plot_metric(df, f"{metric.title()} (interactive)", os.path.basename(f))
                )
            except Exception:
                # Skip files that match name patterns but aren't actually the metric tables
                continue

    if not figs["closeness"] and not figs["betweenness"]:
        # Give a helpful error showing what is in the directory
        names = sorted(os.listdir(work_dir))[:200]
        raise FileNotFoundError(
            "No plottable closeness/betweenness tables found.\n"
            f"work_dir={work_dir}\n"
            f"Files (first 200): {names}\n\n"
            "If your output filenames don't include 'close/closeness' or 'betw/betweenness', "
            "either rename them or pass explicit closeness_file / betweenness_file."
        )

    return figs
