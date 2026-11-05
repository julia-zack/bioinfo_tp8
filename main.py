"""Entry point for the TP. Runs one exercise at a time, by name:

    venv/bin/python main.py 1.1
    venv/bin/python main.py          # lists the available exercises

The names are <objetivo>.<ejercicio>, matching the numbering in the TP PDF.
"""

import sys

from exercises.obj1 import ex1 as o1_ex1
from exercises.obj1 import ex2 as o1_ex2
from exercises.obj1 import ex3 as o1_ex3
from exercises.obj1 import ex4 as o1_ex4

EXERCISES = {
    "1.1":   o1_ex1.run,   # cuántos DEGs quedan con cada criterio de corte
    "1.2":   o1_ex2.run,   # volcano plot
    "1.3":   o1_ex3.run,   # PCA de las muestras con los DEGs
    "1.4":   o1_ex4.run,   # heatmap con clustering
}


def usage():
    print("usage: python3 main.py <exercise>")
    print("available exercises:")
    for name in EXERCISES:
        print(f"  {name}")


def main():
    if len(sys.argv) != 2:
        usage()
        return
    name = sys.argv[1]
    if name not in EXERCISES:
        print(f"'{name}' is not a known exercise.\n")
        usage()
        return
    EXERCISES[name]()


if __name__ == "__main__":
    main()
