"""Objetivo 2, PEA, ejercicio 3: gráficos del PEA (dotplot y barplot de gseapy)."""

from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt

from degs import filter_degs, load_de_table
from exercises.obj2.pea1 import GENE_SETS
from exercises.obj2.pea2 import ADJ_P_MAX, run_enrichr

FIG_DIR = Path(__file__).parents[2] / "figures"
TOP_TERM = 7  # términos por base, como en la consigna


def deg_direction():
    """Para cada símbolo DEG, si sube o baja en el tumor (logFC es WT − PTC)."""
    degs = filter_degs(load_de_table()).dropna(subset=["SYMBOL"])
    degs = degs.assign(SYMBOL=degs["SYMBOL"].str.upper())
    # Si un gen tiene varias sondas DEG, usamos la de mayor |logFC|.
    best = degs.loc[degs["logFC"].abs().groupby(degs["SYMBOL"]).idxmax()]
    return best.set_index("SYMBOL")["logFC"].map(lambda x: "sube" if x < 0 else "baja")


def run():
    results = run_enrichr()
    significant = results[results["Adjusted P-value"] < ADJ_P_MAX]
    FIG_DIR.mkdir(exist_ok=True)

    gp.dotplot(significant, column="Adjusted P-value", x="Gene_set", top_term=TOP_TERM,
               size=5, figsize=(6, 16), xticklabels_rot=45, show_ring=True,
               title="PEA: términos más significativos por base",
               ofname=str(FIG_DIR / "pea_dotplot.png"))
    gp.barplot(significant, column="Adjusted P-value", group="Gene_set",
               top_term=TOP_TERM, size=10, figsize=(4, 14),
               color=dict(zip(GENE_SETS, plt.get_cmap("tab10").colors)),
               title="PEA: términos más significativos por base",
               ofname=str(FIG_DIR / "pea_barplot.png"))
    plt.close("all")
    print(f"Figuras guardadas en {FIG_DIR}/pea_dotplot.png y pea_barplot.png\n")

    # El PEA no distingue dirección: contamos cuántos genes de cada término
    # suben y cuántos bajan en el tumor.
    direction = deg_direction()
    print(f"Dirección de los genes en los {TOP_TERM} términos principales de cada base:")
    for db, group in significant.groupby("Gene_set", sort=False):
        print(f"\n== {db} ==")
        for _, row in group.nsmallest(TOP_TERM, "Adjusted P-value").iterrows():
            genes = row["Genes"].split(";")
            up = sum(direction.get(g) == "sube" for g in genes)
            print(f"  {row['Term'][:65]:65} {up:>3} suben / {len(genes) - up:>3} bajan")
