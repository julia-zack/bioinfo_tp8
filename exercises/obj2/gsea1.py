"""Objetivo 2, GSEA, ejercicio 1: ranking de todos los genes y GSEA (prerank)."""

import pickle
from pathlib import Path

import gseapy as gp
import numpy as np

from degs import load_de_table
from exercises.obj2.pea1 import GENE_SETS

RESULTS_DIR = Path(__file__).parents[2] / "results"
CACHE = RESULTS_DIR / "gsea_prerank.pkl"
FDR_MAX = 0.05
SEED = 42  # semilla de las permutaciones, para que los resultados sean reproducibles


def ranking(collapse="aveexpr"):
    """Un puntaje por gen, ordenado de mayor a menor (positivo = sube en el tumor).

    collapse: cómo elegir la sonda de los genes que tienen varias. "aveexpr" (la que
    usamos) toma la de mayor expresión media; "extremo" toma la de mayor |rank|.
    """
    de = load_de_table()
    # La consigna usa −log10(P.Value) × logFC. Nuestro logFC es WT − PTC, así que lo
    # invertimos (PTC − WT) para que un puntaje positivo signifique "sube en el tumor".
    de["rank"] = -np.log10(de["P.Value"]) * -de["logFC"]
    de = de.dropna(subset=["SYMBOL"])
    de["SYMBOL"] = de["SYMBOL"].str.upper()
    # Si un gen tiene varias sondas, nos quedamos con la de mayor expresión media
    # (AveExpr): la que mejor mide el gen, elegida sin mirar el resultado.
    key = de["AveExpr"] if collapse == "aveexpr" else de["rank"].abs()
    best = de.loc[key.groupby(de["SYMBOL"]).idxmax()]
    return (best[["SYMBOL", "rank"]].dropna()
            .sort_values("rank", ascending=False).reset_index(drop=True))


def run_prerank(force=False):
    """Corre el GSEA (tarda varios minutos) o lo lee del cache en results/."""
    if CACHE.exists() and not force:
        with open(CACHE, "rb") as f:
            return pickle.load(f)
    result = gp.prerank(rnk=ranking(), gene_sets=GENE_SETS, permutation_num=1000,
                        organism="human", seed=SEED, threads=4, outdir=None)
    RESULTS_DIR.mkdir(exist_ok=True)
    with open(CACHE, "wb") as f:
        pickle.dump(result, f)
    return result


def run():
    rnk = ranking()
    print(f"Ranking: {len(rnk)} genes únicos "
          f"({(rnk['rank'] > 0).sum()} con puntaje positivo, {(rnk['rank'] < 0).sum()} negativo)")
    print("Principio de la fila (suben en el tumor):")
    print(rnk.head(5).round(2).to_string(index=False))
    print("Final de la fila (suben en el tejido sano):")
    print(rnk.tail(5).round(2).to_string(index=False), "\n")

    res = run_prerank().res2d
    res["FDR q-val"] = res["FDR q-val"].astype(float)
    res["NES"] = res["NES"].astype(float)
    sig = res[res["FDR q-val"] < FDR_MAX]
    print(f"Términos evaluados: {len(res)}; con FDR q-val < {FDR_MAX}: {len(sig)} "
          f"({(sig['NES'] > 0).sum()} más activos en el tumor, "
          f"{(sig['NES'] < 0).sum()} más activos en el tejido sano)\n")

    cols = ["Term", "NES", "FDR q-val", "Lead_genes"]
    up, down = sig[sig["NES"] > 0], sig[sig["NES"] < 0]
    for label, part in [("Más activos en el tumor (NES > 0)", up.nlargest(10, "NES")),
                        ("Más activos en el tejido sano (NES < 0)", down.nsmallest(10, "NES"))]:
        print(f"== {label} ==")
        part = part[cols].assign(Lead_genes=part["Lead_genes"].str.split(";").str.len())
        print(part.rename(columns={"Lead_genes": "n_lead"}).round(4).to_string(index=False), "\n")

    thyroid = res[res["Term"].str.contains("Thyroid hormone synthesis")]
    print("Thyroid hormone synthesis:")
    print(thyroid[["Term", "NES", "FDR q-val"]].to_string(index=False))
