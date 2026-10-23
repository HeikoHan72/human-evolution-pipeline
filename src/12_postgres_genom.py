#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 08:42:13 2026

@author: Heiko Han
"""

# %% 12 - PostgreSQL: Datenbank genom aufbauen und abfragen
# Legt (falls noetig) die Datenbank genom im Container bde-postgres an,
# erstellt Tabellen, laedt die CSV aus Ausgabe/export und scan_1mb.csv und fuehrt Abfragen aus.
# Zugangsdaten: Umgebungsvariablen POSTGRES_* (Vorlage: .env.example).
# Schreibt nur in die Datenbank genom (Schema public). Dateien werden nur gelesen.

import os
import io
import importlib
from pathlib import Path
import pandas as pd
import psycopg2
from sqlalchemy import create_engine, text

cfg = importlib.import_module("00_config")
importlib.reload(cfg)

# %% Einstellungen (alles an einer Stelle)
DB_HOST = cfg.env("POSTGRES_HOST", "localhost")
DB_PORT = int(cfg.env("POSTGRES_PORT", "5432"))
DB_USER = cfg.env("POSTGRES_USER", "postgres")
DB_PASS = cfg.env("POSTGRES_PASSWORD")          # Pflicht, kommt aus .env oder der Umgebung
DB_NAME = "genom"
NEU_AUFBAUEN = os.environ.get("GENOM_NEU_AUFBAUEN") == "1"        # True: Tabellen in genom loeschen und neu aufbauen
                            # (eine leere/halb gebaute Datenbank wird automatisch neu aufgebaut)
EXPORT = Path(cfg.AUSGABE) / "export"
SCAN_CSV = Path(cfg.AUSGABE) / "scan_1mb.csv"

# %% Stammdaten
POPULATIONEN = {
    "AFR": ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"],
    "AMR": ["CLM", "MXL", "PEL", "PUR"],
    "EAS": ["CDX", "CHB", "CHS", "JPT", "KHV"],
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "SAS": ["BEB", "GIH", "ITU", "PJL", "STU"],
}
# Gene (GRCh38, aus der NCBI-Annotation GCF_000001405.40): gen, region, chrom, start, ende, strang
GENE = [
    ("TLR10", "TLR", "chr4", 38772238, 38782990, "-"),
    ("TLR1", "TLR", "chr4", 38787569, 38805644, "-"),
    ("TLR6", "TLR", "chr4", 38822897, 38868390, "-"),
    ("STAT2", "STAT2", "chr12", 56341597, 56360107, "-"),
    ("OAS1", "OAS", "chr12", 112906962, 112933219, "+"),
    ("OAS3", "OAS", "chr12", 112938474, 112973251, "+"),
    ("OAS2", "OAS", "chr12", 112978519, 113011723, "+"),
]

# %% Datenbank anlegen, falls sie fehlt
admin = psycopg2.connect(host=DB_HOST, port=DB_PORT, user=DB_USER, password=DB_PASS, dbname="postgres")
admin.autocommit = True
with admin.cursor() as c:
    c.execute("SELECT 1 FROM pg_database WHERE datname = %s", (DB_NAME,))
    if not c.fetchone():
        c.execute(f"CREATE DATABASE {DB_NAME}")
        print("Datenbank angelegt:", DB_NAME)
    else:
        print("Datenbank vorhanden:", DB_NAME)
admin.close()

conn = psycopg2.connect(host=DB_HOST, port=DB_PORT, user=DB_USER, password=DB_PASS, dbname=DB_NAME)
cur = conn.cursor()
cur.execute("SELECT to_regclass('public.kandidaten')")
if cur.fetchone()[0] is not None:
    cur.execute("SELECT count(*) FROM kandidaten")
    n_vorhanden = cur.fetchone()[0]
    if n_vorhanden > 0 and not NEU_AUFBAUEN:
        raise SystemExit(f"Tabellen existieren schon ({n_vorhanden} Kandidaten). "
                         "Zum Neuaufbau NEU_AUFBAUEN = True setzen.")
    cur.execute("DROP SCHEMA public CASCADE; CREATE SCHEMA public;")
    conn.commit()
    print("Tabellen in genom geloescht, Neuaufbau")

# %% Tabellen anlegen
cur.execute("""
CREATE TABLE regionen (
    region TEXT PRIMARY KEY, chrom TEXT NOT NULL, start_pos INTEGER NOT NULL, ende_pos INTEGER NOT NULL);
CREATE TABLE gene (
    gen TEXT PRIMARY KEY, region TEXT REFERENCES regionen(region), chrom TEXT NOT NULL,
    start_pos INTEGER NOT NULL, ende_pos INTEGER NOT NULL, strang CHAR(1));
CREATE TABLE populationen (pop TEXT PRIMARY KEY, superpop TEXT NOT NULL);
CREATE TABLE kandidaten (
    chrom TEXT NOT NULL, chrom_nr SMALLINT NOT NULL, pos INTEGER NOT NULL, ref TEXT NOT NULL, alt TEXT NOT NULL,
    id TEXT, anc TEXT, abgeleitet_ist TEXT,
    af_alle DOUBLE PRECISION, af_afr DOUBLE PRECISION, af_amr DOUBLE PRECISION, af_eas DOUBLE PRECISION,
    af_eur DOUBLE PRECISION, af_sas DOUBLE PRECISION,
    daf_afr DOUBLE PRECISION, daf_amr DOUBLE PRECISION, daf_eas DOUBLE PRECISION, daf_eur DOUBLE PRECISION,
    daf_sas DOUBLE PRECISION, daf_nichtafr DOUBLE PRECISION,
    altai TEXT, vindija TEXT, denisova TEXT, chagyrskaya TEXT, n_abgedeckt SMALLINT,
    PRIMARY KEY (chrom, pos, ref, alt));
CREATE TABLE haeufigkeiten (
    chrom TEXT NOT NULL, chrom_nr SMALLINT NOT NULL, pos INTEGER NOT NULL, ref TEXT NOT NULL, alt TEXT NOT NULL,
    pop TEXT NOT NULL, superpop TEXT NOT NULL, n_abg INTEGER, n_all INTEGER, daf DOUBLE PRECISION);
CREATE TABLE scan_1mb (
    chrom TEXT NOT NULL, chrom_nr SMALLINT NOT NULL, fenster INTEGER NOT NULL, start_pos INTEGER, ende_pos INTEGER,
    n_snv INTEGER, n_kandidaten INTEGER, n_neandertaler_sonst INTEGER, n_denisova INTEGER, n_nicht_abgedeckt INTEGER,
    kandidaten_je_1000_snv DOUBLE PRECISION, rang_im_chrom INTEGER, PRIMARY KEY (chrom, fenster));
""")
conn.commit()

# %% Stammdaten laden
for name, r in cfg.REGIONEN.items():
    chrom = "chr" + str(r["chrom"]).replace("chr", "")
    cur.execute("INSERT INTO regionen VALUES (%s, %s, %s, %s)", (name, chrom, int(r["start"]), int(r["ende"])))
for g in GENE:
    cur.execute("INSERT INTO gene VALUES (%s, %s, %s, %s, %s, %s)", g)
for sup, pops in POPULATIONEN.items():
    for p in pops:
        cur.execute("INSERT INTO populationen VALUES (%s, %s)", (p, sup))
conn.commit()
print("Stammdaten geladen: regionen, gene, populationen")

# %% Hilfsfunktion: DataFrame per COPY in eine Tabelle laden
def copy_df(df, tabelle, spalten):
    puffer = io.StringIO()
    df[spalten].to_csv(puffer, index=False, header=False)
    puffer.seek(0)
    cur.copy_expert(f"COPY {tabelle} ({', '.join(spalten)}) FROM STDIN WITH (FORMAT csv, NULL '')", puffer)

# %% Kandidaten laden
spalten_k = ["chrom", "chrom_nr", "pos", "ref", "alt", "id", "anc", "abgeleitet_ist",
             "af_alle", "af_afr", "af_amr", "af_eas", "af_eur", "af_sas",
             "daf_afr", "daf_amr", "daf_eas", "daf_eur", "daf_sas", "daf_nichtafr",
             "altai", "vindija", "denisova", "chagyrskaya", "n_abgedeckt"]
text_k = {c: str for c in ["chrom", "ref", "alt", "id", "anc", "abgeleitet_ist",
                           "altai", "vindija", "denisova", "chagyrskaya"]}
k = pd.read_csv(EXPORT / "kandidaten.csv", dtype=text_k, keep_default_na=False, na_values=[""])
k["n_abgedeckt"] = k["n_abgedeckt"].astype("Int64")
copy_df(k, "kandidaten", spalten_k)
conn.commit()
print("kandidaten geladen:", len(k))
del k

# %% Haeufigkeiten laden (in Stuecken, 2,3 Mio Zeilen)
spalten_h = ["chrom", "chrom_nr", "pos", "ref", "alt", "pop", "superpop", "n_abg", "n_all", "daf"]
text_h = {c: str for c in ["chrom", "ref", "alt", "pop", "superpop"]}
n_h = 0
for stueck in pd.read_csv(EXPORT / "haeufigkeiten.csv", dtype=text_h, keep_default_na=False,
                          na_values=[""], chunksize=500_000):
    stueck["n_abg"] = stueck["n_abg"].astype("Int64")
    stueck["n_all"] = stueck["n_all"].astype("Int64")
    copy_df(stueck, "haeufigkeiten", spalten_h)
    n_h += len(stueck)
    print("  haeufigkeiten:", n_h)
conn.commit()

# %% Scan laden
s = pd.read_csv(SCAN_CSV, dtype={"chrom": str})
spalten_s = list(s.columns)
copy_df(s, "scan_1mb", spalten_s)
conn.commit()
print("scan_1mb geladen:", len(s))

# %% Schluessel und Indizes nach dem Laden (schneller als vorher)
cur.execute("""
ALTER TABLE haeufigkeiten ADD PRIMARY KEY (chrom, pos, ref, alt, pop);
ALTER TABLE haeufigkeiten ADD FOREIGN KEY (chrom, pos, ref, alt) REFERENCES kandidaten (chrom, pos, ref, alt);
ALTER TABLE haeufigkeiten ADD FOREIGN KEY (pop) REFERENCES populationen (pop);
CREATE INDEX idx_haeufigkeiten_pop ON haeufigkeiten (pop);
CREATE INDEX idx_kandidaten_chrom_pos ON kandidaten (chrom, pos);
CREATE INDEX idx_gene_chrom ON gene (chrom, start_pos, ende_pos);
""")
conn.commit()
cur.execute("ANALYZE")
conn.commit()
print("Schluessel und Indizes angelegt")

# %% View: mittlere Haeufigkeit je Kandidat, Afrika gegen Nicht-Afrika
# Mittelwert ueber die Populationen (jede Population zaehlt gleich, nicht nach Personenzahl gewichtet).
cur.execute("""
CREATE VIEW v_kandidat_mittel AS
SELECT chrom, pos, ref, alt,
       AVG(daf) FILTER (WHERE superpop = 'AFR')  AS mittel_afr,
       AVG(daf) FILTER (WHERE superpop <> 'AFR') AS mittel_nichtafr
FROM haeufigkeiten
GROUP BY chrom, pos, ref, alt;
""")
conn.commit()

# %% Kontrolle der Zeilenzahlen
for t, erwartet in [("kandidaten", 88632), ("haeufigkeiten", 2304432), ("scan_1mb", 2777),
                    ("regionen", None), ("gene", 7), ("populationen", 26)]:
    cur.execute(f"SELECT count(*) FROM {t}")
    print(f"{t:15s}", cur.fetchone()[0], "" if erwartet is None else f"(erwartet {erwartet})")
cur.close()
conn.close()

# %% Abfragen
engine = create_engine(f"postgresql+psycopg2://{DB_USER}:{DB_PASS}@{DB_HOST}:{DB_PORT}/{DB_NAME}")

def abfrage(titel, sql):
    print("\n" + "=" * 70 + "\n" + titel + "\n" + "=" * 70)
    with engine.connect() as con:
        print(pd.read_sql(text(sql), con).to_string(index=False))

abfrage("1) Kandidaten je Gen (JOIN ueber Position; STAT2 zeigt den negativen Befund)", """
SELECT g.gen, g.region, g.chrom, COUNT(k.pos) AS kandidaten
FROM gene g
LEFT JOIN kandidaten k ON k.chrom = g.chrom AND k.pos BETWEEN g.start_pos AND g.ende_pos
GROUP BY g.gen, g.region, g.chrom
ORDER BY g.region, g.gen""")

abfrage("2) Mittlere Haeufigkeit der abgeleiteten Base je Region und Superpopulation (Mittel ueber Populationen)", """
SELECT r.region, h.superpop,
       COUNT(DISTINCT h.pos) AS kandidaten,
       ROUND(AVG(h.daf)::numeric, 3) AS mittel_daf
FROM haeufigkeiten h
JOIN regionen r ON r.chrom = h.chrom AND h.pos BETWEEN r.start_pos AND r.ende_pos
GROUP BY r.region, h.superpop
ORDER BY r.region, h.superpop""")

abfrage("3) Kandidaten je Chromosom mit Rang (Window Function)", """
SELECT chrom, SUM(n_kandidaten) AS kandidaten, SUM(n_snv) AS snv,
       ROUND(1000.0 * SUM(n_kandidaten) / SUM(n_snv), 2) AS je_1000_snv,
       RANK() OVER (ORDER BY SUM(n_kandidaten) DESC) AS rang
FROM scan_1mb
GROUP BY chrom, chrom_nr
ORDER BY chrom_nr""")

abfrage("4) Top 10 Megabase-Fenster nach Kandidatenzahl", """
SELECT chrom, start_pos, ende_pos, n_snv, n_kandidaten, kandidaten_je_1000_snv
FROM scan_1mb
ORDER BY n_kandidaten DESC
LIMIT 10""")

abfrage("5) Top 10 Kandidaten nach mittlerer Haeufigkeit bei Nicht-Afrikanern (View, mit Gen falls vorhanden)", """
SELECT v.chrom, v.pos, v.ref, v.alt, g.gen,
       ROUND(v.mittel_afr::numeric, 3) AS mittel_afr,
       ROUND(v.mittel_nichtafr::numeric, 3) AS mittel_nichtafr,
       k.altai, k.vindija, k.denisova, k.chagyrskaya
FROM v_kandidat_mittel v
JOIN kandidaten k USING (chrom, pos, ref, alt)
LEFT JOIN gene g ON g.chrom = v.chrom AND v.pos BETWEEN g.start_pos AND g.ende_pos
ORDER BY v.mittel_nichtafr DESC
LIMIT 10""")
