"""Common matplotlib style for PMG algorithm documentation figures."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# palette
C_FINE   = "#3B6FB6"   # fine level / rank blue
C_COARSE = "#D9822B"   # coarse level orange
C_GC     = "#B03A48"   # global-coarse red
C_GHOST  = "#9AA5B1"   # ghost gray
C_SEND   = "#2E8B57"   # send green
C_RECV   = "#8E44AD"   # recv purple
C_BG     = "#F7F8FA"

RANK_COLORS = ["#3B6FB6", "#D9822B", "#2E8B57", "#8E44AD"]
RANK_LIGHT  = ["#C6D5EC", "#F2D5B4", "#BFE0CE", "#DCC8E8"]

def setup(figsize=(10, 6), dpi=150):
    plt.rcParams.update({
        "font.size": 11,
        "font.family": "DejaVu Sans",
        "axes.linewidth": 0.8,
        "savefig.bbox": "tight",
        "savefig.dpi": dpi,
        "savefig.facecolor": "white",
    })
    return plt.subplots(figsize=figsize)

def save(fig, name):
    import os
    out = os.path.join(os.path.dirname(__file__), "..", "figures", name)
    fig.savefig(out)
    print("saved:", os.path.abspath(out))
