"""Carga de las tablas de GSE3467 y filtrado de DEGs."""

from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).parent / "data" / "gse3467"


def load_de_table():
    """DE_table indexada por ID de sonda (la primera columna, sin nombre)."""
    return pd.read_csv(DATA_DIR / "DE_table.csv", index_col=0)


def passes(de, p_col, p_max, logfc_min):
    """Máscara de las sondas con p < p_max y |logFC| > logfc_min."""
    return (de[p_col] < p_max) & (de["logFC"].abs() > logfc_min)


# Criterio elegido en el ejercicio 1.1: adj.P.Val < 0.05 y |logFC| > 1 (2 veces).
P_COL = "adj.P.Val"
P_MAX = 0.05
LOGFC_MIN = 1


def filter_degs(de):
    """Sondas que pasan el criterio de corte elegido."""
    return de[passes(de, P_COL, P_MAX, LOGFC_MIN)]
