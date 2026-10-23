#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 14:49:05 2026

@author: Heiko Han
"""

# %% 04_oas_gruppen.py - OAS-Cluster: Varianten nach Herkunft in Gruppen einteilen (Spark)
# Startpunkt sind die 1000-Genomes-Varianten (heutiger Mensch); daran haengen wir Ahnenzustand und
# archaische Genotypen. Es wird nichts geschrieben, nur ausgegeben.

import sys
from pathlib import Path

try:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
except NameError:
    sys.path.insert(0, str(Path.cwd()))

import importlib
cfg = importlib.import_module("00_config")
import genom_funktionen as gf
importlib.reload(gf)

# %% Parameter (eigene Festlegung, in der Dokumentation nennen)
AFR_MAX = 0.03          # Haeufigkeit der abgeleiteten Base in Afrika hoechstens
NICHTAFR_MIN = 0.05     # Mittel der anderen vier Grossgruppen mindestens
REGION_NAME = "OAS"     # "OAS", "TLR" oder "STAT2"

# %% Daten holen
R = cfg.REGIONEN[REGION_NAME]
nean = gf.lese_neandertaler(R["chrom"], R["start"], R["ende"])
seq = gf.lese_ahnenzustand(R["chrom"], R["start"], R["ende"])
kg = gf.lese_tausend_genome_info(R["chrom"], R["start"], R["ende"])
kg["anc"] = [seq[p - R["start"]] for p in kg["pos"]]
print(f"Region {REGION_NAME}: 1000-Genomes-SNVs {len(kg)}, Neandertaler-Eintraege {len(nean)}")

# %% Spark: DataFrames anlegen
from pyspark.sql import SparkSession, functions as F

spark = (SparkSession.builder.master("local[*]")
         .appName("genom_04_gruppen")
         .config("spark.driver.memory", "4g")
         .getOrCreate())

s_kg = spark.createDataFrame(
    kg, schema="pos long, ref string, alt string, id string, af_alle double, af_afr double, "
               "af_amr double, af_eas double, af_eur double, af_sas double, anc string")
arch = nean[["pos", "ref", "alt", "altai", "vindija", "denisova", "chagyrskaya"]].rename(
    columns={"pos": "pos_n", "ref": "ref_n", "alt": "alt_n"})
s_arch = spark.createDataFrame(
    arch, schema="pos_n long, ref_n string, alt_n string, altai string, vindija string, "
                 "denisova string, chagyrskaya string")

# %% Join: pro 1000-Genomes-Position die archaischen Genotypen (Position und Referenz muessen passen)
df = s_kg.join(s_arch, (s_kg.pos == s_arch.pos_n) & (s_kg.ref == s_arch.ref_n), "left")

# %% Abgeleitete Base und Haeufigkeit der abgeleiteten Base (DAF)
df = df.withColumn("abgeleitet_ist",
                   F.when(F.upper("anc") == F.col("ref"), "alt")
                    .when(F.upper("anc") == F.col("alt"), "ref"))
GRUPPEN = ["afr", "amr", "eas", "eur", "sas"]
for g in GRUPPEN:
    df = df.withColumn(f"daf_{g}",
                       F.when(F.col("abgeleitet_ist") == "alt", F.col(f"af_{g}"))
                        .when(F.col("abgeleitet_ist") == "ref", 1 - F.col(f"af_{g}")))
df = df.withColumn("daf_nichtafr",
                   (F.col("daf_amr") + F.col("daf_eas") + F.col("daf_eur") + F.col("daf_sas")) / 4)

# %% Traegt eine archaische Probe die abgeleitete Base?  (None = fehlt/nicht abgedeckt)
def hat_info(gt):
    return F.col(gt).isNotNull() & (F.col(gt) != "./.")

def traegt_abgeleitet(gt):
    alt_gleich = F.col("alt_n") == F.col("alt")       # archaische Alt-Base = 1000-Genomes-Alt-Base
    tr_alt = F.when(alt_gleich, F.col(gt).contains("1")).otherwise(F.lit(False))
    tr_ref = F.col(gt).contains("0")
    return (F.when(~hat_info(gt), F.lit(None).cast("boolean"))
             .when(F.col("abgeleitet_ist") == "alt", tr_alt)
             .when(F.col("abgeleitet_ist") == "ref", tr_ref))

PROBEN = ["altai", "vindija", "denisova", "chagyrskaya"]
for p in PROBEN:
    df = df.withColumn(f"d_{p}", traegt_abgeleitet(p))

df = df.withColumn("n_abgedeckt", sum(F.when(hat_info(p), 1).otherwise(0) for p in PROBEN))
df = df.withColumn("nean_traegt", F.coalesce("d_altai", F.lit(False)) | F.coalesce("d_vindija", F.lit(False))
                                  | F.coalesce("d_chagyrskaya", F.lit(False)))
df = df.withColumn("deni_traegt", F.coalesce("d_denisova", F.lit(False)))

# %% Gruppe zuweisen
df = df.withColumn(
    "gruppe",
    F.when(F.col("abgeleitet_ist").isNull(), "unklar_Ahnenzustand")
     .when(F.col("n_abgedeckt") == 0, "archaisch_nicht_abgedeckt")
     .when(F.col("nean_traegt") & (F.col("daf_afr") <= AFR_MAX) & (F.col("daf_nichtafr") >= NICHTAFR_MIN),
           "geteilt_Neandertaler_nicht_afrikanisch")
     .when(F.col("nean_traegt"), "geteilt_Neandertaler_sonst")
     .when(F.col("deni_traegt"), "nur_Denisova_geteilt")
     .otherwise("nur_Mensch"))

# %% Ergebnis
print("\nAnzahl 1000-Genomes-Varianten je Gruppe:")
df.groupBy("gruppe").count().orderBy(F.desc("count")).show(truncate=False)

print("Kandidaten: geteilt_Neandertaler_nicht_afrikanisch")
(df.filter("gruppe = 'geteilt_Neandertaler_nicht_afrikanisch'").orderBy("pos")
   .select("pos", "ref", "alt", "anc", "abgeleitet_ist", "altai", "vindija", "denisova", "chagyrskaya",
           F.round("daf_afr", 3).alias("afr"), F.round("daf_eur", 3).alias("eur"),
           F.round("daf_eas", 3).alias("eas"), F.round("daf_sas", 3).alias("sas"),
           F.round("daf_amr", 3).alias("amr"))
   .show(60, truncate=False))

spark.stop()