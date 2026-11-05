"""Objetivo 2, GSEA, ejercicio 2: gráfico de GSEA para un término (síntesis de hormona tiroidea)."""

from pathlib import Path

import matplotlib.pyplot as plt

from exercises.obj2.gsea1 import run_prerank

FIG_DIR = Path(__file__).parents[2] / "figures"
TERM = "KEGG_2021_Human__Thyroid hormone synthesis"


def run():
    result = run_prerank()
    row = result.res2d.set_index("Term").loc[TERM]
    print(f"{TERM}: ES {float(row['ES']):.3f}, NES {float(row['NES']):.3f}, "
          f"FDR q-val {float(row['FDR q-val']):.3f}, "
          f"genes del término en el ranking: {row['Tag %']}")

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "gsea_thyroid.png"
    result.plot(terms=TERM, show_ranking=True, ofname=str(out))
    plt.close("all")
    print(f"Figura guardada en {out}")
