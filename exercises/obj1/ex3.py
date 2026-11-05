"""Objetivo 1, ejercicio 3: PCA de las muestras usando solo los DEGs."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from degs import DATA_DIR, filter_degs, load_de_table

FIG_DIR = Path(__file__).parents[2] / "figures"
COLORS = {"WT": "tab:blue", "PTC": "tab:red"}
LABELS = {"WT": "WT (tejido sano)", "PTC": "PTC (tumor)"}


def load_deg_matrix():
    """Expresión de las sondas DEG: una fila por sonda, una columna por muestra."""
    matrix = pd.read_csv(DATA_DIR / "full_matrix.csv", index_col=0)
    return matrix.loc[filter_degs(load_de_table()).index]


def run():
    deg_matrix = load_deg_matrix()
    print(f"Matriz de DEGs: {deg_matrix.shape[0]} sondas x {deg_matrix.shape[1]} muestras\n")

    # Transponemos para que cada fila sea una muestra: así el PCA busca las
    # combinaciones de genes que más diferencian a las muestras, y cada muestra
    # queda como un punto. StandardScaler estandariza cada sonda (media 0, desvío 1).
    samples = deg_matrix.T
    pca = PCA(n_components=2)
    coords = pca.fit_transform(StandardScaler().fit_transform(samples))
    var = pca.explained_variance_ratio_ * 100
    coords = pd.DataFrame(coords, index=samples.index, columns=["PC1", "PC2"])
    coords["condición"] = coords.index.str.split("_").str[0]

    print(f"Varianza explicada: PC1 {var[0]:.1f} %, PC2 {var[1]:.1f} %\n")
    print(coords.sort_values("PC1").round(1).to_string(), "\n")

    # Para comparar: el código de la consigna, sin transponer (cada sonda es un punto
    # del PCA y las muestras aparecen recién en components_).
    pca_consigna = PCA(n_components=2).fit(StandardScaler().fit_transform(deg_matrix))
    var_c = pca_consigna.explained_variance_ratio_ * 100
    pc1_c = pca_consigna.components_[0]
    print(f"Código de la consigna (sin transponer): PC1 {var_c[0]:.1f} %, PC2 {var_c[1]:.1f} %; "
          f"PC1 de las muestras entre {pc1_c.min():.2f} y {pc1_c.max():.2f}")

    fig, ax = plt.subplots(figsize=(8, 6))
    for cond, color in COLORS.items():
        sub = coords[coords["condición"] == cond]
        ax.scatter(sub["PC1"], sub["PC2"], color=color, s=40, label=LABELS[cond])
    # La condición ya la indica el color, así que rotulamos solo el número de muestra
    # (que es el del paciente: WT_003 y PTC_003 son del mismo paciente).
    for name, row in coords.iterrows():
        number = name.split("_")[1].lstrip("0")
        ax.annotate(number, (row["PC1"], row["PC2"]), textcoords="offset points",
                    xytext=(5, -3), fontsize=8)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)
    ax.set_xlabel(f"PC1 ({var[0]:.1f} %)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f} %)")
    ax.set_title(f"GSE3467: PCA de las 18 muestras ({deg_matrix.shape[0]} sondas DEG)")
    ax.legend(title="condición")

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "pca.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"\nFigura guardada en {out}")
