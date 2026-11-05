"""Objetivo 1, ejercicio 1: cuántos DEGs quedan con distintos criterios de corte."""

from degs import LOGFC_MIN, P_COL, P_MAX, filter_degs, load_de_table, passes

P_CUTS = [("P.Value", 0.05), ("adj.P.Val", 0.05), ("adj.P.Val", 0.01)]
LOGFC_CUTS = [0, 0.58, 1, 1.5]  # 0 = sin corte por logFC; 0.58, 1, 1.5 = 1.5x, 2x, 3x


def run():
    de = load_de_table()
    print(f"DE_table: {len(de)} sondas, {de['SYMBOL'].nunique()} símbolos únicos, "
          f"{de['SYMBOL'].isna().sum()} sondas sin SYMBOL\n")

    print(f"  {'p-valor':18} {'|logFC| >':>9} {'sondas':>8} {'símbolos':>9} {'logFC>0':>8} {'logFC<0':>8}")
    for p_col, p_max in P_CUTS:
        for logfc_min in LOGFC_CUTS:
            hits = de[passes(de, p_col, p_max, logfc_min)]
            up = (hits["logFC"] > 0).sum()
            print(f"  {p_col + ' < ' + str(p_max):18} {logfc_min:>9} {len(hits):>8} "
                  f"{hits['SYMBOL'].nunique():>9} {up:>8} {len(hits) - up:>8}")
        print()

    degs = filter_degs(de)
    print(f"Criterio elegido: {P_COL} < {P_MAX} y |logFC| > {LOGFC_MIN} -> "
          f"{len(degs)} sondas, {degs['SYMBOL'].nunique()} genes únicos, "
          f"{degs['SYMBOL'].isna().sum()} sondas sin SYMBOL")
