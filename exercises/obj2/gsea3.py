"""Objetivo 2, GSEA, ejercicio 3: gráfico de GSEA de 4 términos significativos."""

from pathlib import Path

import matplotlib.pyplot as plt

from exercises.obj2.gsea1 import FDR_MAX, run_prerank

FIG_DIR = Path(__file__).parents[2] / "figures"

# 2 términos más activos en el tumor (NES > 0), uno por cada gran tema (respuesta
# inmune; invasión y matriz extracelular), y 2 más activos en el tejido sano (NES < 0).
# Del lado sano solo hay 4 términos significativos: elegimos la fosforilación oxidativa,
# el único con interpretación biológica clara, y el cilio móvil, el que tiene más genes
# líderes de los otros tres.
MIS_TERMS = [
    "MSigDB_Hallmark_2020__Interferon Gamma Response",
    "MSigDB_Hallmark_2020__Epithelial Mesenchymal Transition",
    "MSigDB_Hallmark_2020__Oxidative Phosphorylation",
    "GO_Cellular_Component_2026__9+2 Motile Cilium (GO:0097729)",
]


def run():
    result = run_prerank()
    res = result.res2d.set_index("Term")
    sig = res[res["FDR q-val"].astype(float) < FDR_MAX]
    print(f"Términos con FDR q-val < {FDR_MAX}: {len(sig)} "
          f"({(sig['NES'].astype(float) > 0).sum()} con NES > 0, "
          f"{(sig['NES'].astype(float) < 0).sum()} con NES < 0)\n")

    for term in MIS_TERMS:
        row = sig.loc[term]  # falla si el término no es significativo
        print(f"{term}\n  ES {float(row['ES']):.3f}, NES {float(row['NES']):.3f}, "
              f"FDR q-val {float(row['FDR q-val']):.4f}, genes líderes / genes del término: "
              f"{row['Tag %']}\n  primeros genes líderes: "
              f"{', '.join(row['Lead_genes'].split(';')[:10])}")

    FIG_DIR.mkdir(exist_ok=True)
    out = FIG_DIR / "gsea_mis_terms.png"
    result.plot(terms=MIS_TERMS, show_ranking=True, figsize=(6, 6), ofname=str(out))
    plt.close("all")
    print(f"\nFigura guardada en {out}")
