"""Objetivo 1, ejercicio 4: heatmap con clustering de las sondas DEG y las muestras."""

from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import fcluster

from degs import load_de_table
from exercises.obj1.ex3 import COLORS, LABELS, load_deg_matrix

FIG_DIR = Path(__file__).parents[2] / "figures"


def run():
    deg_matrix = load_deg_matrix()
    conditions = deg_matrix.columns.str.split("_").str[0]
    print(f"Matriz de DEGs: {deg_matrix.shape[0]} sondas x {deg_matrix.shape[1]} muestras\n")

    # Igual que la consigna: z_score=0 estandariza cada sonda (fila) por separado, así
    # el color indica si la sonda está más o menos expresada en esa muestra que en el
    # promedio de las 18, sin importar su nivel de expresión absoluto.
    # Como en el PCA, la condición va en una barra de color y rotulamos solo el número.
    grid = sns.clustermap(
        deg_matrix, method="complete", cmap="RdBu_r", z_score=0, figsize=(8, 8),
        col_colors=pd.Series(conditions.map(COLORS), index=deg_matrix.columns, name=""),
        yticklabels=False, cbar_kws={"label": "z-score"}, tree_kws={"linewidths": 0.4},
    )
    order = [deg_matrix.columns[i] for i in grid.dendrogram_col.reordered_ind]
    grid.ax_heatmap.set_xticklabels([name.split("_")[1].lstrip("0") for name in order],
                                    rotation=0)
    grid.ax_heatmap.set_xlabel("muestra (número de paciente)")
    grid.ax_col_colors.set_yticks([])
    grid.ax_heatmap.set_ylabel(f"{deg_matrix.shape[0]} sondas DEG")
    grid.ax_col_dendrogram.legend(
        handles=[Patch(color=COLORS[c], label=LABELS[c]) for c in COLORS],
        loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8, frameon=False)
    grid.figure.suptitle("GSE3467: heatmap de los DEGs", y=1.02)

    print("Orden de las muestras en el heatmap (según el clustering):")
    print("  " + ", ".join(order) + "\n")

    # Los dos grupos que se forman al cortar el dendrograma de muestras en la raíz.
    col_groups = fcluster(grid.dendrogram_col.linkage, t=2, criterion="maxclust")
    for g in (1, 2):
        members = deg_matrix.columns[col_groups == g]
        print(f"Grupo {g} de muestras: {', '.join(members)}")

    # Lo mismo para las sondas: cuántas hay en cada grupo y hacia dónde cambian.
    logfc_all = load_de_table().loc[deg_matrix.index, "logFC"].to_numpy()
    row_groups = fcluster(grid.dendrogram_row.linkage, t=2, criterion="maxclust")
    print()
    for g in (1, 2):
        logfc = logfc_all[row_groups == g]
        print(f"Grupo {g} de sondas: {len(logfc)} sondas, "
              f"{(logfc < 0).sum()} suben en el tumor y {(logfc > 0).sum()} bajan")

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "heatmap.png"
    grid.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nFigura guardada en {out}")

