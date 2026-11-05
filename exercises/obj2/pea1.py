"""Objetivo 2, PEA, ejercicio 1: bases de datos de conjuntos de genes disponibles en Enrichr."""

import re

import gseapy as gp

# Bases elegidas para el PEA y el GSEA: las que nombra la consigna, en su versión más
# nueva, salvo KEGG: usamos la 2021 porque el ejercicio de GSEA de la consigna grafica
# el término "KEGG_2021_Human__Thyroid hormone synthesis".
GENE_SETS = [
    "GO_Biological_Process_2026",
    "GO_Molecular_Function_2026",
    "GO_Cellular_Component_2026",
    "KEGG_2021_Human",
    "Reactome_Pathways_2024",
    "MSigDB_Hallmark_2020",
]

FAMILIES = ["GO_", "KEGG", "Reactome", "MSigDB", "WikiPathways"]


def run():
    libraries = gp.get_library_name(organism="Human")
    print(f"Enrichr tiene {len(libraries)} bases de datos para humanos. "
          f"Las de las familias más usadas:\n")
    for family in FAMILIES:
        matches = [lib for lib in libraries if re.match(family, lib)]
        print(f"  {family}: {', '.join(matches)}")

    missing = [db for db in GENE_SETS if db not in libraries]
    if missing:
        raise SystemExit(f"\nNo están disponibles: {', '.join(missing)}")
    print("\nBases elegidas (todas disponibles):")
    for db in GENE_SETS:
        print(f"  {db}")
