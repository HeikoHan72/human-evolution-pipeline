#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 15:35:04 2026

@author: Heiko Han
"""

# %% 08_alle_chromosomen.py - Spark ueber alle vollstaendig geladenen Chromosomen, Ausgabe als Parquet
# SCHREIBT nach Ausgabe/chromosomen.parquet/chrom=chrN (neu angelegt). Vorhandene Chromosomen werden uebersprungen.

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

# %% Vollstaendigkeit: erwartete Dateigroesse laut Server (Bytes). Chromosomen ohne Eintrag werden nicht verarbeitet.
ERWARTET = {1: 2390515911, 2: 2506332497, 3: 2097106049, 4: 2086491071, 5: 1899165551,
            6: 1878333223, 7: 1751799644, 8: 1617138777, 9: 1320051150, 10: 1493904721,
            11: 1444152921, 12: 1418304697, 13: 1073594120, 14: 966241761, 15: 885687136,
            16: 961378545, 17: 861943967, 18: 834387264, 19: 717737562, 20: 687320087,
            21: 427136458, 22: 445701977}
REIHENFOLGE = list(range(1, 23))   # alle 22 Chromosomen; fertige werden uebersprungen
ZIEL = cfg.AUSGABE / "chromosomen.parquet"

# %% Spark (sparsam wegen Akademie-Browser)
spark = (SparkSession.builder.master("local[4]")
         .appName("genom_08_alle_chromosomen")
         .config("spark.driver.memory", "4g")
         .config("spark.sql.shuffle.partitions", "16")
         .getOrCreate())

# %% Schleife ueber Chromosomen
cfg.AUSGABE.mkdir(exist_ok=True)
ZIEL.mkdir(exist_ok=True)
for chrom in REIHENFOLGE:
    vcf = cfg.tausend_genome_vcf(chrom)
    ordner = ZIEL / f"chrom=chr{chrom}"
    if not vcf.exists() or vcf.stat().st_size != ERWARTET[chrom]:
        print(f"chr{chrom}: VCF fehlt oder unvollstaendig -> uebersprungen")
        continue
    if (ordner / "_SUCCESS").exists():
        print(f"chr{chrom}: Ausgabe existiert schon -> uebersprungen")
        continue
    t0 = time.time()
    df = gc.baue_chromosomtabelle(spark, chrom).drop("chrom")      # chrom steckt im Ordnernamen
    df.write.mode("errorifexists").parquet(str(ordner))
    n = spark.read.parquet(str(ordner)).count()
    print(f"chr{chrom}: {n} Zeilen geschrieben in {time.time() - t0:.0f} Sekunden")

# %% Gesamtuebersicht aus den Parquet-Dateien
alle = spark.read.parquet(str(ZIEL))
print("\nKandidaten und Gruppen je Chromosom:")
(alle.groupBy("chrom", "gruppe").count()
     .orderBy("chrom", "gruppe").show(60, truncate=False))
print("Zeilen gesamt:", alle.count())

spark.stop()