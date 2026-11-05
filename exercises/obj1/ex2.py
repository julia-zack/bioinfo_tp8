"""Objetivo 1, ejercicio 2: volcano plot de todas las sondas de DE_table."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from degs import LOGFC_MIN, P_COL, P_MAX, load_de_table

FIG_DIR = Path(__file__).parents[2] / "figures"

# Mismas categorías que la Fig. 1 de la consigna: qué cortes pasa cada sonda.
CATEGORIES = [
    # (nombre, color, pasa el corte de p, pasa el corte de logFC)
    ("NS", "lightgray", False, False),
    (f"solo |logFC| > {LOGFC_MIN}", "tab:green", False, True),
    (f"solo {P_COL} < {P_MAX}", "tab:blue", True, False),
    ("ambos (DEGs)", "tab:red", True, True),
]


def volcano(ax, de, masks, y_col, y_cut):
    """Dibuja el volcano con −log10(y_col) en el eje Y y la línea de corte en y_cut."""
    y = -np.log10(de[y_col])
    for (name, color, _, _), mask in zip(CATEGORIES, masks):
        ax.scatter(de.loc[mask, "logFC"], y[mask], s=4, color=color,
                   label=f"{name} ({mask.sum()})", rasterized=True)
    ax.axhline(-np.log10(y_cut), color="black", linestyle="--", linewidth=0.8)
    for x in (-LOGFC_MIN, LOGFC_MIN):
        ax.axvline(x, color="black", linestyle="--", linewidth=0.8)
    ax.set_xlabel("logFC (WT − PTC)\n←  sube en el tumor   |   baja en el tumor  →")
    ax.set_ylabel(f"−log10({y_col})")
    ax.legend(loc="upper left", fontsize=8, markerscale=3)
    print(f"  {y_col}: línea en {y_col} = {y_cut:.2e} (−log10 = {-np.log10(y_cut):.2f}), "
          f"eje Y de {y.min():.2f} a {y.max():.2f}")


def run():
    de = load_de_table()
    pass_p = de[P_COL] < P_MAX
    pass_fc = de["logFC"].abs() > LOGFC_MIN
    masks = [(pass_p == in_p) & (pass_fc == in_fc) for _, _, in_p, in_fc in CATEGORIES]

    print(f"{'categoría':28} {'sondas':>7}")
    for (name, *_), mask in zip(CATEGORIES, masks):
        print(f"{name:28} {mask.sum():>7}")
    print(f"Rango de logFC: {de['logFC'].min():.2f} a {de['logFC'].max():.2f}\n")

    # Panel izquierdo: p crudo, como pide la consigna. El corte es por el ajustado,
    # pero como el ajuste por FDR (Benjamini-Hochberg) no cambia el orden de los
    # p-valores, equivale a un corte en P.Value: el mayor p crudo que todavía pasa.
    # Panel derecho: p ajustado, donde la línea queda directamente en P_MAX.
    p_raw_cut = de.loc[pass_p, "P.Value"].max()
    fig, (ax_raw, ax_adj) = plt.subplots(1, 2, figsize=(15, 6), sharex=True)
    volcano(ax_raw, de, masks, "P.Value", p_raw_cut)
    volcano(ax_adj, de, masks, "adj.P.Val", P_MAX)
    ax_raw.set_title("p-valor crudo")
    ax_adj.set_title("p-valor ajustado (FDR)")
    fig.suptitle("GSE3467: volcano plot (54.675 sondas)")

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "volcano.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nFigura guardada en {out}")
