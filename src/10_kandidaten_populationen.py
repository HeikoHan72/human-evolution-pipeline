#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 17:00:16 2026

@author: Heiko Han
"""

# %% 10 - Haeufigkeit der Kandidaten-Allele je Population (26)
# Fuer jede Kandidaten-Position (Gruppe geteilt_Neandertaler_nicht_afrikanisch)
# wird die Haeufigkeit der abgeleiteten Base je Population aus den Genotypen
# der 1000-Genomes-VCF berechnet.
# Schreibt: Ausgabe/kandidaten_pop.parquet/chrom=chrN  (je Chromosom, errorifexists)
#           Ausgabe/kandidaten_pop_mittel.csv           (wird bei Neulauf ueberschrieben)

import importlib
from pathlib import Path
import pandas as pd
import pysam
from pyspark.sql import SparkSession, functions as F

cfg = importlib.import_module("00_config")
importlib.reload(cfg)

# %% Einstellungen (alles an einer Stelle)
CHROMOSOMEN = None          # Test mit chr12; spaeter None = alle fertigen Chromosomen
NUR_GRUENDER = True         # True: nur Personen ohne Eltern in der Tabelle (FatherID und MotherID = 0)
KANDIDAT = "geteilt_Neandertaler_nicht_afrikanisch"
POP_TABELLE = None          # None = automatisch suchen; sonst Pfad eintragen

AUSGABE = Path(cfg.AUSGABE)
EINGABE = AUSGABE / "chromosomen.parquet"
ZIEL = AUSGABE / "kandidaten_pop.parquet"
ZIEL_CSV = AUSGABE / "kandidaten_pop_mittel.csv"

# %% Spark starten (gleiche sparsame Einstellung wie 08 und 09)
spark = (SparkSession.builder
         .master("local[4]")
         .appName("10_kandidaten_populationen")
         .config("spark.driver.memory", "4g")
         .config("spark.sql.shuffle.partitions", "16")
         .getOrCreate())

# %% Welche Chromosomen? Nur mit fertigem Kandidaten-Parquet (_SUCCESS)
def fertige_chromosomen():
    nummern = []
    for p in EINGABE.glob("chrom=chr*"):
        if (p / "_SUCCESS").exists():
            nummern.append(int(p.name.replace("chrom=chr", "")))
    return sorted(nummern)

chroms = CHROMOSOMEN if CHROMOSOMEN is not None else fertige_chromosomen()
print("Chromosomen:", chroms)

# %% Populationstabelle finden und einlesen
if POP_TABELLE is None:
    suche = list(Path(cfg.TAUSEND_GENOME).glob("*ped_population*.txt")) + \
            list(Path(cfg.TAUSEND_GENOME).parent.glob("*ped_population*.txt"))
    if not suche:
        raise FileNotFoundError("Populationstabelle nicht gefunden: bitte POP_TABELLE setzen")
    POP_TABELLE = suche[0]
print("Populationstabelle:", POP_TABELLE)

ped = pd.read_csv(POP_TABELLE, sep=r"\s+", dtype=str)
print(ped.columns.tolist())
print(ped.head(3))
noetig = {"SampleID", "FatherID", "MotherID", "Population", "Superpopulation"}
if not noetig.issubset(ped.columns):
    raise ValueError(f"Spalten fehlen, gefunden: {ped.columns.tolist()}")

ped["gruender"] = (ped["FatherID"] == "0") & (ped["MotherID"] == "0")
print("Personen gesamt:", len(ped), "| Gruender:", int(ped["gruender"].sum()))

# %% Reihenfolge der Personen aus dem VCF-Kopf, Zuordnung Spaltenindex -> Population
vcf_erste = Path(str(cfg.tausend_genome_vcf(chroms[0])))
proben = list(pysam.VariantFile(str(vcf_erste)).header.samples)
print("Personen im VCF:", len(proben))

nachschlag = ped.set_index("SampleID")
zeilen = []
ohne_treffer = 0
for i, s in enumerate(proben):
    if s not in nachschlag.index:
        ohne_treffer += 1
        continue
    r = nachschlag.loc[s]
    if NUR_GRUENDER and not r["gruender"]:
        continue
    zeilen.append((i, r["Population"], r["Superpopulation"]))
print("Ohne Treffer in der Tabelle:", ohne_treffer, "| verwendete Personen:", len(zeilen))
print(pd.DataFrame(zeilen, columns=["idx", "pop", "superpop"]).groupby("superpop").size())

pop_df = spark.createDataFrame(zeilen, "idx int, pop string, superpop string")

# %% Funktion: ein Chromosom verarbeiten
def verarbeite(chrom):
    ordner_aus = ZIEL / f"chrom=chr{chrom}"
    if (ordner_aus / "_SUCCESS").exists():
        print(f"chr{chrom}: schon fertig, uebersprungen")
        return

    # Kandidaten mit Ahnenbase; abgeleitet = die Base, die nicht dem Ahnenzustand entspricht
    kand = (spark.read.parquet(str(EINGABE / f"chrom=chr{chrom}"))
            .filter(F.col("gruppe") == KANDIDAT)
            .select(F.col("pos").cast("int").alias("pos"), "ref", "alt",
                    F.upper(F.col("anc")).alias("anc"))
            .withColumn("abg_alt",
                        F.when(F.col("anc") == F.col("ref"), F.lit(True))
                         .when(F.col("anc") == F.col("alt"), F.lit(False))))
    n_kand = kand.count()
    kand = kand.filter(F.col("abg_alt").isNotNull()) \
               .withColumn("abg_symbol", F.when(F.col("abg_alt"), "1").otherwise("0"))
    print(f"chr{chrom}: Kandidaten {n_kand}, davon mit eindeutiger Ahnenbase {kand.count()}")

    # VCF-Zeilen lesen, erst nach Position/Ref/Alt filtern, dann Genotypen aufteilen
    roh = (spark.read.text(str(cfg.tausend_genome_vcf(chrom)))
           .filter(~F.col("value").startswith("#")))
    muster = r"^[^\t]*\t([^\t]*)\t[^\t]*\t([^\t]*)\t([^\t]*)\t"
    zeilen_df = roh.select(
        "value",
        F.regexp_extract("value", muster, 1).cast("int").alias("pos"),
        F.regexp_extract("value", muster, 2).alias("ref"),
        F.regexp_extract("value", muster, 3).alias("alt"))
    treffer = zeilen_df.join(F.broadcast(kand), ["pos", "ref", "alt"], "inner")

    gt = (treffer.repartition(16)
          .withColumn("felder", F.split("value", "\t"))
          .withColumn("gts", F.expr("slice(felder, 10, size(felder) - 9)"))
          .select("pos", "ref", "alt", "abg_symbol",
                  F.posexplode("gts").alias("idx", "gt"))
          .withColumn("a1", F.split("gt", r"[|/]").getItem(0))
          .withColumn("a2", F.split("gt", r"[|/]").getItem(1))
          .join(F.broadcast(pop_df), "idx"))

    n_abg = (F.coalesce((F.col("a1") == F.col("abg_symbol")).cast("int"), F.lit(0)) +
             F.coalesce((F.col("a2") == F.col("abg_symbol")).cast("int"), F.lit(0)))
    n_all = (F.coalesce(F.col("a1").isin("0", "1").cast("int"), F.lit(0)) +
             F.coalesce(F.col("a2").isin("0", "1").cast("int"), F.lit(0)))

    res = (gt.groupBy("pos", "ref", "alt", "pop", "superpop")
             .agg(F.sum(n_abg).alias("n_abg"), F.sum(n_all).alias("n_all"))
             .withColumn("daf", F.round(F.col("n_abg") / F.col("n_all"), 4)))

    res.write.mode("errorifexists").parquet(str(ordner_aus))
    print(f"chr{chrom}: {spark.read.parquet(str(ordner_aus)).count()} Zeilen geschrieben")

# %% Lauf (ein Fehler bei einem Chromosom stoppt die anderen nicht)
import time
fehler = []
for c in chroms:
    t = time.time()
    try:
        verarbeite(c)
    except Exception as e:
        fehler.append(c)
        print(f"chr{c}: FEHLER {type(e).__name__}: {str(e)[:300]}")
    print(f"chr{c}: {time.time() - t:.0f} Sekunden")
print("Fehlgeschlagen:", fehler)

# %% Auswertung: mittlere Haeufigkeit der Kandidaten-Allele je Population
fertig = sorted(p for p in ZIEL.glob("chrom=chr*") if (p / "_SUCCESS").exists())
alle = spark.read.option("basePath", str(ZIEL)).parquet(*[str(p) for p in fertig])

mittel = (alle.groupBy("superpop", "pop")
              .agg(F.countDistinct("chrom", "pos").alias("kandidaten"),
                   F.round(F.avg("daf"), 4).alias("mittel_daf"))
              .orderBy("superpop", F.col("mittel_daf").desc()))
mittel.show(30, truncate=False)
mittel.toPandas().to_csv(ZIEL_CSV, index=False)
print("Gespeichert:", ZIEL_CSV)

# %% Kontrolle am OAS-Locus (Vergleich mit Skript 05: AFR ca. 1 %, EUR ca. 31-37 %)
if 12 in [int(p.name.replace("chrom=chr", "")) for p in fertig]:
    o = cfg.REGIONEN["OAS"]
    (alle.filter((F.col("chrom") == "chr12") & F.col("pos").between(o["start"], o["ende"]))
         .groupBy("superpop")
         .agg(F.countDistinct("pos").alias("kandidaten"),
              F.round(F.avg("daf"), 3).alias("mittel_daf"))
         .orderBy("superpop").show())