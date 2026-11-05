"""Objetivo 2, PEA, ejercicio 2: PEA (análisis de enriquecimiento de vías) con Enrichr."""

import gseapy as gp
import pandas as pd

from degs import filter_degs, load_de_table
from exercises.obj2.pea1 import GENE_SETS

ADJ_P_MAX = 0.05


def deg_symbols():
    """Símbolos únicos de los DEGs, en mayúscula y sin sondas sin gen asignado."""
    degs = filter_degs(load_de_table())
    return sorted(degs["SYMBOL"].dropna().str.upper().unique())


def run_enrichr():
    """Corre el PEA contra Enrichr y devuelve la tabla de resultados (todas las vías)."""
    res = gp.enrichr(gene_list=deg_symbols(), gene_sets=GENE_SETS,
                     organism="human", outdir=None)
    return res.results


def run():
    symbols = deg_symbols()
    print(f"Lista de entrada: {len(symbols)} símbolos de gen únicos\n")

    results = run_enrichr()
    significant = results[results["Adjusted P-value"] < ADJ_P_MAX]
    print(f"Términos evaluados: {len(results)}; con Adjusted P-value < {ADJ_P_MAX}: "
          f"{len(significant)}\n")
    print(significant.groupby("Gene_set").size().reindex(GENE_SETS, fill_value=0)
          .to_string(), "\n")

    cols = ["Term", "Overlap", "Adjusted P-value"]
    with pd.option_context("display.max_colwidth", 70, "display.width", 200):
        for db in GENE_SETS:
            top = significant[significant["Gene_set"] == db].nsmallest(8, "Adjusted P-value")
            print(f"== {db} ==")
            print(top[cols].to_string(index=False) if len(top) else "  (ninguno)", "\n")
