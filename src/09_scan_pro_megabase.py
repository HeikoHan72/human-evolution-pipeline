#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 16:51:37 2026

@author: Heiko Han
"""

# %% 09 - Genomweiter Scan: Kandidaten pro Megabase
# Liest die fertigen Chromosomen-Ordner aus Ausgabe/chromosomen.parquet,
# teilt jedes Chromosom in 1-Mb-Fenster und zaehlt die Gruppen je Fenster.
# Schreibt: Ausgabe/scan_1mb.parquet und Ausgabe/scan_1mb.csv
# (ein erneuter Lauf ueberschreibt diese beiden eigenen Ausgaben).

import importlib
from pathlib import Path
from pyspark.sql import SparkSession, Window, functions as F

cfg = importlib.import_module("00_config")
importlib.reload(cfg)

# %% Einstellungen (alles an einer Stelle)
FENSTER = 1_000_000                       # Fenstergroesse in Basen (1 Mb)
KANDIDAT = "geteilt_Neandertaler_nicht_afrikanisch"

AUSGABE = Path(cfg.AUSGABE)
EINGABE = AUSGABE / "chromosomen.parquet"
ZIEL_PARQUET = AUSGABE / "scan_1mb.parquet"
ZIEL_CSV = AUSGABE / "scan_1mb.csv"

# %% Spark starten (gleiche sparsame Einstellung wie Skript 08)
spark = (SparkSession.builder
         .master("local[4]")
         .appName("09_scan_pro_megabase")
         .config("spark.driver.memory", "4g")
         .config("spark.sql.shuffle.partitions", "16")
         .getOrCreate())

# %% Nur vollstaendige Chromosomen-Ordner lesen
fertig = sorted(p for p in EINGABE.glob("chrom=chr*") if (p / "_SUCCESS").exists())
print("Gelesene Chromosomen:", [p.name.replace("chrom=", "") for p in fertig])

df = (spark.read.option("basePath", str(EINGABE))
      .parquet(*[str(p) for p in fertig])
      .select("chrom", "pos", "gruppe")
      .withColumn("chrom_nr", F.regexp_extract("chrom", r"chr(\d+)", 1).cast("int"))
      .withColumn("fenster", F.floor((F.col("pos") - 1) / FENSTER).cast("int")))

df.createOrReplaceTempView("varianten")

# %% Aggregation mit Spark SQL: Zaehlungen je Chromosom und Fenster
scan = spark.sql(f"""
    SELECT chrom, chrom_nr, fenster,
           fenster * {FENSTER}              AS start_pos,
           (fenster + 1) * {FENSTER}        AS ende_pos,
           COUNT(*)                         AS n_snv,
           SUM(CASE WHEN gruppe = '{KANDIDAT}' THEN 1 ELSE 0 END)
                                            AS n_kandidaten,
           SUM(CASE WHEN gruppe = 'geteilt_Neandertaler_sonst' THEN 1 ELSE 0 END)
                                            AS n_neandertaler_sonst,
           SUM(CASE WHEN gruppe = 'nur_Denisova_geteilt' THEN 1 ELSE 0 END)
                                            AS n_denisova,
           SUM(CASE WHEN gruppe = 'archaisch_nicht_abgedeckt' THEN 1 ELSE 0 END)
                                            AS n_nicht_abgedeckt
    FROM varianten
    GROUP BY chrom, chrom_nr, fenster
""")

# Anteil Kandidaten je 1000 SNVs (macht Fenster mit unterschiedlicher SNV-Dichte vergleichbar)
scan = scan.withColumn("kandidaten_je_1000_snv",
                       F.round(F.col("n_kandidaten") / F.col("n_snv") * 1000, 2))

# Window Function: Rang des Fensters innerhalb seines Chromosoms
w = Window.partitionBy("chrom").orderBy(F.col("n_kandidaten").desc())
scan = scan.withColumn("rang_im_chrom", F.rank().over(w)).cache()

# %% Kontrolle: Summe der SNVs muss zur Ausgabe von Skript 08 passen
summe = scan.agg(F.sum("n_snv").alias("snv"), F.sum("n_kandidaten").alias("kand")).first()
print("Fenster gesamt:", scan.count())
print("SNVs gesamt:", summe["snv"], "| Kandidaten gesamt:", summe["kand"])

# %% Zusammenfassung je Chromosom
(scan.groupBy("chrom_nr")
     .agg(F.count("*").alias("fenster"),
          F.sum("n_kandidaten").alias("kandidaten"),
          F.round(F.avg("n_kandidaten"), 1).alias("mittel_je_mb"),
          F.max("n_kandidaten").alias("max_je_mb"))
     .orderBy("chrom_nr").show(30, truncate=False))

# %% Die 15 Fenster mit den meisten Kandidaten
(scan.orderBy(F.col("n_kandidaten").desc())
     .select("chrom", "start_pos", "n_snv", "n_kandidaten", "kandidaten_je_1000_snv")
     .show(15, truncate=False))

# %% Ergebnis speichern (klein: einige tausend Zeilen)
scan.orderBy("chrom_nr", "fenster").write.mode("overwrite").parquet(str(ZIEL_PARQUET))
scan.orderBy("chrom_nr", "fenster").toPandas().to_csv(ZIEL_CSV, index=False)
print("Gespeichert:", ZIEL_PARQUET, "und", ZIEL_CSV)