#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 08:34:26 2026

@author: Heiko Han
"""

# %% 11 - Export der Kandidaten-Tabellen als CSV fuer die Datenbanken
# Liest Ausgabe/chromosomen.parquet (nur Kandidaten) und Ausgabe/kandidaten_pop.parquet.
# Schreibt: Ausgabe/export/kandidaten.csv und Ausgabe/export/haeufigkeiten.csv
# (neu angelegter Ordner; ein Neulauf ueberschreibt nur diese beiden Dateien).

import importlib
from pathlib import Path
from pyspark.sql import SparkSession, functions as F

cfg = importlib.import_module("00_config")
importlib.reload(cfg)

# %% Einstellungen
KANDIDAT = "geteilt_Neandertaler_nicht_afrikanisch"
AUSGABE = Path(cfg.AUSGABE)
CHROM_PQ = AUSGABE / "chromosomen.parquet"
POP_PQ = AUSGABE / "kandidaten_pop.parquet"
EXPORT = AUSGABE / "export"
EXPORT.mkdir(exist_ok=True)

spark = (SparkSession.builder
         .master("local[4]")
         .appName("11_export_fuer_datenbanken")
         .config("spark.driver.memory", "4g")
         .config("spark.sql.shuffle.partitions", "16")
         .getOrCreate())

def fertige_ordner(basis):
    return sorted(p for p in basis.glob("chrom=chr*") if (p / "_SUCCESS").exists())

# %% Kandidaten (alle Spalten aus Skript 08 plus Chromosom)
ordner = fertige_ordner(CHROM_PQ)
print("Chromosomen mit fertigem Parquet:", len(ordner))
alle = spark.read.option("basePath", str(CHROM_PQ)).parquet(*[str(p) for p in ordner])

kand = (alle.filter(F.col("gruppe") == KANDIDAT)
            .withColumn("chrom_nr", F.regexp_extract("chrom", r"chr(\d+)", 1).cast("int"))
            .orderBy("chrom_nr", "pos"))
print("Spalten:", kand.columns)
n = kand.count()
print("Kandidaten:", n, "(erwartet 88632)")
kand.toPandas().to_csv(EXPORT / "kandidaten.csv", index=False)

# %% Haeufigkeiten je Population
ordner_pop = fertige_ordner(POP_PQ)
print("Chromosomen mit Populationsdaten:", len(ordner_pop))
pop = (spark.read.option("basePath", str(POP_PQ)).parquet(*[str(p) for p in ordner_pop])
            .withColumn("chrom_nr", F.regexp_extract("chrom", r"chr(\d+)", 1).cast("int"))
            .select("chrom", "chrom_nr", "pos", "ref", "alt", "pop", "superpop",
                    "n_abg", "n_all", "daf")
            .orderBy("chrom_nr", "pos", "pop"))
m = pop.count()
print("Haeufigkeiten:", m, "(erwartet 2304432)")
pop.toPandas().to_csv(EXPORT / "haeufigkeiten.csv", index=False)

print("Geschrieben nach:", EXPORT)