#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 12:00:29 2026

@author: Heiko Han
"""

# %% 03_oas_haeufigkeiten.py - OAS-Cluster: Neandertaler + Ahnenzustand + 1000-Genomes-Haeufigkeiten (Spark-Join)
# Es wird nichts geschrieben, nur ausgegeben.

import sys
from pathlib import Path

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
except NameError:
    sys.path.insert(0, str(Path.cwd()))

import importlib
cfg = importlib.import_module("00_config")
import genom_funktionen as gf
importlib.reload(gf)        # damit die neuen Funktionen sicher geladen sind

# %% Daten holen
R = cfg.REGIONEN["OAS"]
nean = gf.lese_neandertaler(R["chrom"], R["start"], R["ende"])
seq = gf.lese_ahnenzustand(R["chrom"], R["start"], R["ende"])
nean["anc"] = [seq[p - R["start"]] for p in nean["pos"]]

var = nean[nean["alt"] != "."].copy()
var["status"] = [gf.ahnenzustand_status(r, a, n) for r, a, n in zip(var["ref"], var["alt"], var["anc"])]
print("Neandertaler-Variantenpositionen:", len(var))

kg = gf.lese_tausend_genome_info(R["chrom"], R["start"], R["ende"])
print("1000-Genomes-SNVs (biallelisch) im Fenster:", len(kg))
print(kg.head(3).to_string())

# %% Spark-Join ueber Position + Referenz + Alt-Base
from pyspark.sql import SparkSession

spark = (SparkSession.builder.master("local[*]")
         .appName("genom_03_oas_haeufigkeiten")
         .config("spark.driver.memory", "4g")
         .getOrCreate())

s_nean = spark.createDataFrame(
    var[["pos", "ref", "alt", "anc", "status", "altai", "vindija", "denisova", "chagyrskaya"]],
    schema="pos long, ref string, alt string, anc string, status string, "
           "altai string, vindija string, denisova string, chagyrskaya string")
s_kg = spark.createDataFrame(
    kg, schema="pos long, ref string, alt string, id string, af_alle double, "
               "af_afr double, af_amr double, af_eas double, af_eur double, af_sas double")

joined = s_nean.join(s_kg.drop("id"), on=["pos", "ref", "alt"], how="left")
print("Zeilen nach Join:", joined.count(),
      "| davon mit 1000-Genomes-Treffer:", joined.filter("af_alle is not null").count())

# %% Kandidaten: Alt-Allel in mindestens einem der drei Neandertaler
kand = joined.filter("altai like '%1%' or vindija like '%1%' or chagyrskaya like '%1%'").orderBy("pos")
kand.select("pos", "ref", "alt", "anc", "status", "altai", "vindija", "denisova", "chagyrskaya",
            "af_afr", "af_eur", "af_eas", "af_sas", "af_amr").show(40, truncate=False)

spark.stop()