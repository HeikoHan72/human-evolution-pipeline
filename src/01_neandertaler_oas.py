#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 10:50:42 2026

@author: Heiko Han
"""

# %% 01_neandertaler_oas.py - Neandertaler-Genotypen des OAS-Clusters lesen (pysam) und an Spark uebergeben
# Es wird nichts geschrieben, nur ausgegeben.

import sys
import importlib
from pathlib import Path

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
except NameError:                      # Zellenmodus: kein __file__
    sys.path.insert(0, str(Path.cwd()))

cfg = importlib.import_module("00_config")

import pandas as pd
from pysam import VariantFile

# %% Einstellungen
REGION = cfg.REGIONEN["OAS"]
KURZ = {"AltaiNeandertal": "altai", "Vindija33.19": "vindija",
        "Denisova": "denisova", "Chagyrskaya-Phalanx": "chagyrskaya"}

def gt_text(gt):
    """Genotyp als Text. Fehlende Werte werden './.' (nicht '0/0')."""
    if gt is None or any(a is None for a in gt):
        return "./."
    return "/".join(str(a) for a in gt)

# %% BCF lesen
bcf_pfad = cfg.neandertaler_bcf(REGION["chrom"])
zeilen = []
with VariantFile(str(bcf_pfad)) as vf:
    # fetch: 0-basiert, Ende exklusiv
    for rec in vf.fetch(f"chr{REGION['chrom']}", REGION["start"] - 1, REGION["ende"]):
        zeile = {"chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref,
                 "alt": ",".join(rec.alts) if rec.alts else "."}
        for probe, kurz in KURZ.items():
            zeile[kurz] = gt_text(rec.samples[probe]["GT"])
        zeilen.append(zeile)

df = pd.DataFrame(zeilen)
print("Eintraege im OAS-Fenster:", len(df))
print(df.head(5).to_string())

# %% Kontrolle: fehlende Genotypen pro Probe, Variantenpositionen
for kurz in KURZ.values():
    print(f"{kurz:<12} fehlend: {(df[kurz] == './.').sum():>5}   "
          f"mit Alt-Allel: {df[kurz].str.contains('1').sum():>4}")
print("Positionen mit ALT != '.':", (df["alt"] != ".").sum())

# %% An Spark uebergeben (lokaler Modus wie im Kurs)
from pyspark.sql import SparkSession

spark = (SparkSession.builder.master("local[*]")
         .appName("genom_01_neandertaler_oas")
         .config("spark.driver.memory", "4g")
         .getOrCreate())

schema = "chrom string, pos long, ref string, alt string, altai string, vindija string, denisova string, chagyrskaya string"
sdf = spark.createDataFrame(df, schema=schema)
sdf.printSchema()
print("Zeilen in Spark:", sdf.count())
sdf.filter("alt != '.'").orderBy("pos").show(10, truncate=False)

spark.stop()