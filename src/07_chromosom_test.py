#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 15:23:27 2026

@author: Heiko Han
"""

# %% 07_chromosom_test.py - Test: ein ganzes Chromosom mit Spark (nur Ausgabe, nichts wird geschrieben)
# Kontrolle: Im OAS-Fenster muessen dieselben Gruppenzahlen herauskommen wie in Skript 04.

import sys
import time
from pathlib import Path

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
except NameError:
    sys.path.insert(0, str(Path.cwd()))

import importlib
cfg = importlib.import_module("00_config")
import genom_funktionen as gf
import genom_spark as gs
import genom_chrom as gc
for m in (gf, gs, gc):
    importlib.reload(m)

from pyspark.sql import SparkSession, functions as F

CHROM = 12

# %% Spark: begrenzt auf 6 Kerne und 5 GB (Laptop hat 12 Kerne, 15 GiB RAM)
spark = (SparkSession.builder.master("local[4]")
         .appName("genom_07_chromosom_test")
         .config("spark.driver.memory", "4g")
         .config("spark.sql.shuffle.partitions", "16")
         .getOrCreate())

# %% Pipeline fuer das ganze Chromosom; eine einzige Aktion am Ende
t0 = time.time()
df = gc.baue_chromosomtabelle(spark, CHROM)

oas = cfg.REGIONEN["OAS"]
ergebnis = (df.groupBy("gruppe")
              .agg(F.count("*").alias("n_chromosom"),
                   F.sum(F.when(F.col("pos").between(oas["start"], oas["ende"]), 1).otherwise(0)).alias("n_im_OAS_fenster"))
              .orderBy(F.desc("n_chromosom"))
              .collect())

print(f"\nChromosom {CHROM}: Laufzeit {time.time() - t0:.0f} Sekunden")
print(f"{'gruppe':<42}{'ganzes Chromosom':>18}{'OAS-Fenster':>14}")
for r in ergebnis:
    print(f"{r['gruppe']:<42}{r['n_chromosom']:>18}{r['n_im_OAS_fenster']:>14}")
print("Summe:", sum(r["n_chromosom"] for r in ergebnis))
print("Erwartet im OAS-Fenster: 745 / 51 / 50 / 60 / 1065 / 18")

spark.stop()