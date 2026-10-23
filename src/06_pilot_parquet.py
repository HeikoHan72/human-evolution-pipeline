#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 15:12:14 2026

@author: Heiko Han
"""

# %% 06_pilot_parquet.py - drei Pilotregionen in einer Tabelle, als Parquet speichern (partitioniert nach Region)
# SCHREIBT: legt den Ordner "Ausgabe" an und schreibt dort "pilot_varianten.parquet".
# Existiert das Ziel schon, bricht Spark ab (nichts wird ueberschrieben).

import sys
from pathlib import Path
from functools import reduce

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
except NameError:
    sys.path.insert(0, str(Path.cwd()))

import importlib
cfg = importlib.import_module("00_config")
import genom_funktionen as gf
import genom_spark as gs
importlib.reload(gf)
importlib.reload(gs)

from pyspark.sql import SparkSession, DataFrame, functions as F

# %% Spark starten
spark = (SparkSession.builder.master("local[*]")
         .appName("genom_06_pilot_parquet")
         .config("spark.driver.memory", "4g")
         .getOrCreate())

# %% Regionen berechnen und zusammenfuegen
teile = []
for name in ["OAS", "STAT2", "TLR"]:
    t = gs.baue_regionstabelle(spark, name)
    print(f"{name}: {t.count()} Zeilen")
    teile.append(t)
alle = reduce(DataFrame.unionByName, teile)

# %% Parquet schreiben
cfg.AUSGABE.mkdir(exist_ok=True)
ziel = cfg.AUSGABE / "pilot_varianten.parquet"
alle.write.mode("errorifexists").partitionBy("region").parquet(str(ziel))
print("geschrieben nach:", ziel)

# %% Kontrolle: wieder einlesen
gelesen = spark.read.parquet(str(ziel))
print("Zeilen gesamt:", gelesen.count())
gelesen.groupBy("region", "gruppe").count().orderBy("region", "gruppe").show(40, truncate=False)

spark.stop()