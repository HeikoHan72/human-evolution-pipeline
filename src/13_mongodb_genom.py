#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Heiko Han
"""

# %% 13 - MongoDB: Kandidaten als Dokumente mit eingebetteten Populationshaeufigkeiten
# Legt in der Datenbank genom (Container/Dienst auf localhost:27017) die Collection kandidaten an:
# ein Dokument je Kandidat, darin die Haeufigkeiten der 26 Populationen als Liste.
# Liest: Ausgabe/export/kandidaten.csv und haeufigkeiten.csv (nur lesen).
# Schreibt: nur die Collection genom.kandidaten in MongoDB. Bestehende andere Datenbanken bleiben unberuehrt.

import os
import importlib
from pathlib import Path
import pandas as pd
from pymongo import MongoClient, ASCENDING

cfg = importlib.import_module("00_config")
importlib.reload(cfg)

# %% Einstellungen
MONGO_URI = cfg.env("MONGO_URI", "mongodb://localhost:27017")
DB_NAME = "genom"
COLL = "kandidaten"
NEU_AUFBAUEN = os.environ.get("GENOM_NEU_AUFBAUEN") == "1"        # True: Collection genom.kandidaten loeschen und neu laden
EXPORT = Path(cfg.AUSGABE) / "export"
BATCH = 5000

# Gene wie in Skript 12 (spaeter nach 00_config verschieben)
GENE = [
    ("TLR10", "TLR", "chr4", 38772238, 38782990, "-"),
    ("TLR1", "TLR", "chr4", 38787569, 38805644, "-"),
    ("TLR6", "TLR", "chr4", 38822897, 38868390, "-"),
    ("STAT2", "STAT2", "chr12", 56341597, 56360107, "-"),
    ("OAS1", "OAS", "chr12", 112906962, 112933219, "+"),
    ("OAS3", "OAS", "chr12", 112938474, 112973251, "+"),
    ("OAS2", "OAS", "chr12", 112978519, 113011723, "+"),
]

# %% Verbindung
client = MongoClient(MONGO_URI, serverSelectionTimeoutMS=5000)
client.admin.command("ping")
coll = client[DB_NAME][COLL]
vorhanden = coll.estimated_document_count()
if vorhanden > 0 and not NEU_AUFBAUEN:
    raise SystemExit(f"Collection {DB_NAME}.{COLL} enthaelt schon {vorhanden} Dokumente. "
                     "Zum Neuaufbau NEU_AUFBAUEN = True setzen.")
if vorhanden > 0:
    coll.drop()
    print("Collection geloescht, Neuaufbau")

# %% Hilfsfunktionen
def gen_zu(chrom, pos):
    """Gen (nur innerhalb der Genkoordinaten) und Pilotregion (Fenster aus 00_config, wie in PostgreSQL)."""
    gen = next((g[0] for g in GENE if g[2] == chrom and g[3] <= pos <= g[4]), None)
    region = next((name for name, r in cfg.REGIONEN.items()
                   if "chr" + str(r["chrom"]).replace("chr", "") == chrom
                   and int(r["start"]) <= pos <= int(r["ende"])), None)
    return gen, region

def n_abgeleitet(gt, abgeleitet_ist):
    """Anzahl abgeleiteter Allele (0, 1, 2) aus einem Genotyp wie '1/1'.
    abgeleitet_ist = 'alt': das ALT-Allel (1) ist abgeleitet; 'ref': das REF-Allel (0) ist abgeleitet.
    Fehlend ('./.') oder unklar: None."""
    if not isinstance(gt, str) or "." in gt or abgeleitet_ist not in ("alt", "ref"):
        return None
    ziel = "1" if abgeleitet_ist == "alt" else "0"
    return sum(1 for a in gt.replace("|", "/").split("/") if a == ziel)

# %% Daten lesen
textspalten = {c: str for c in ["chrom", "ref", "alt", "id", "anc", "abgeleitet_ist",
                                "altai", "vindija", "denisova", "chagyrskaya"]}
k = pd.read_csv(EXPORT / "kandidaten.csv", dtype=textspalten, keep_default_na=False, na_values=[""])
h = pd.read_csv(EXPORT / "haeufigkeiten.csv", dtype={c: str for c in ["chrom", "ref", "alt", "pop", "superpop"]},
                keep_default_na=False, na_values=[""])
print("Gelesen: Kandidaten", len(k), "| Haeufigkeiten", len(h))

def nn(v):
    """pandas-Fehlwerte (NaN/NA) zu None, numpy-Zahlen zu Python-Zahlen."""
    if pd.isna(v):
        return None
    return v.item() if hasattr(v, "item") else v

# %% Haeufigkeiten je Kandidat einsammeln
freq = {}
for chrom, pos, ref, alt, pop, sup, n_abg, n_all, daf in zip(
        h["chrom"], h["pos"], h["ref"], h["alt"], h["pop"], h["superpop"], h["n_abg"], h["n_all"], h["daf"]):
    freq.setdefault((chrom, int(pos), ref, alt), []).append(
        {"pop": pop, "superpop": sup, "n_abg": nn(n_abg), "n_all": nn(n_all), "daf": nn(daf)})
del h

# %% Dokumente bauen und einfuegen
def dokument(z):
    chrom, pos, ref, alt = z.chrom, int(z.pos), z.ref, z.alt
    gen, region = gen_zu(chrom, pos)
    ist = nn(z.abgeleitet_ist)
    archaisch = {n: nn(getattr(z, n)) for n in ("altai", "vindija", "denisova", "chagyrskaya")}
    return {
        "_id": f"{chrom}:{pos}:{ref}:{alt}",
        "chrom": chrom, "chrom_nr": int(z.chrom_nr), "pos": pos, "ref": ref, "alt": alt,
        "rsid_oder_id": nn(z.id), "ahnenzustand": nn(z.anc), "abgeleitet_ist": ist,
        "gen": gen, "region": region,
        "af": {s: nn(getattr(z, "af_" + s)) for s in ("alle", "afr", "amr", "eas", "eur", "sas")},
        "daf": {s: nn(getattr(z, "daf_" + s)) for s in ("afr", "amr", "eas", "eur", "sas", "nichtafr")},
        "archaisch_genotyp": archaisch,
        "archaisch_abgeleitet": {n: n_abgeleitet(g, ist) for n, g in archaisch.items()},
        "n_abgedeckt": nn(z.n_abgedeckt),
        "haeufigkeiten": sorted(freq.get((chrom, pos, ref, alt), []), key=lambda d: d["pop"]),
    }

stapel, n = [], 0
for z in k.itertuples(index=False):
    stapel.append(dokument(z))
    if len(stapel) == BATCH:
        coll.insert_many(stapel, ordered=False)
        n += len(stapel)
        stapel = []
        print("  eingefuegt:", n)
if stapel:
    coll.insert_many(stapel, ordered=False)
    n += len(stapel)
print("Dokumente eingefuegt:", n, "(erwartet 88632)")

# %% Indizes
coll.create_index([("chrom", ASCENDING), ("pos", ASCENDING)])
coll.create_index([("gen", ASCENDING)])
coll.create_index([("region", ASCENDING)])
print("Indizes:", [i["name"] for i in coll.list_indexes()])

# %% Beispieldokument
bsp = coll.find_one({"gen": "OAS1"}, {"haeufigkeiten": {"$slice": 2}})
print("\nBeispieldokument (Haeufigkeitsliste auf 2 Eintraege gekuerzt):")
for kk, v in bsp.items():
    print(f"  {kk}: {v}")

# %% Abfragen (Aggregation Pipeline)
def zeige(titel, pipeline):
    print("\n" + "=" * 70 + "\n" + titel + "\n" + "=" * 70)
    print(pd.DataFrame(list(coll.aggregate(pipeline))).to_string(index=False))

zeige("1) Mittlere Haeufigkeit je Region und Superpopulation (Gegenprobe zu PostgreSQL Abfrage 2)", [
    {"$match": {"region": {"$ne": None}}},
    {"$unwind": "$haeufigkeiten"},
    {"$group": {"_id": {"region": "$region", "superpop": "$haeufigkeiten.superpop"},
                "kandidaten": {"$addToSet": "$_id"},
                "mittel_daf": {"$avg": "$haeufigkeiten.daf"}}},
    {"$project": {"_id": 0, "region": "$_id.region", "superpop": "$_id.superpop",
                  "kandidaten": {"$size": "$kandidaten"}, "mittel_daf": {"$round": ["$mittel_daf", 3]}}},
    {"$sort": {"region": 1, "superpop": 1}},
])

zeige("2) Kandidaten je Chromosom", [
    {"$group": {"_id": {"chrom": "$chrom", "nr": "$chrom_nr"}, "kandidaten": {"$sum": 1}}},
    {"$sort": {"_id.nr": 1}},
    {"$project": {"_id": 0, "chrom": "$_id.chrom", "kandidaten": 1}},
])

zeige("3) Top 10 Kandidaten nach mittlerer Haeufigkeit bei Nicht-Afrikanern (Mittel ueber Populationen)", [
    {"$addFields": {
        "mittel_nichtafr": {"$avg": {"$map": {
            "input": {"$filter": {"input": "$haeufigkeiten", "as": "x", "cond": {"$ne": ["$$x.superpop", "AFR"]}}},
            "as": "y", "in": "$$y.daf"}}},
        "mittel_afr": {"$avg": {"$map": {
            "input": {"$filter": {"input": "$haeufigkeiten", "as": "x", "cond": {"$eq": ["$$x.superpop", "AFR"]}}},
            "as": "y", "in": "$$y.daf"}}}}},
    {"$sort": {"mittel_nichtafr": -1}},
    {"$limit": 10},
    {"$project": {"_id": 1, "gen": 1, "mittel_afr": {"$round": ["$mittel_afr", 3]},
                  "mittel_nichtafr": {"$round": ["$mittel_nichtafr", 3]},
                  "archaisch_abgeleitet": 1}},
])

zeige("4) Wer traegt die abgeleitete Base? Neandertaler-Individuen (von 3) und Denisova (Genotyp auf abgeleitet umgerechnet)", [
    {"$addFields": {
        "nean_traeger": {"$size": {"$filter": {
            "input": [{"$ifNull": ["$archaisch_abgeleitet.altai", 0]},
                      {"$ifNull": ["$archaisch_abgeleitet.vindija", 0]},
                      {"$ifNull": ["$archaisch_abgeleitet.chagyrskaya", 0]}],
            "as": "g", "cond": {"$gt": ["$$g", 0]}}}},
        "denisova_traegt": {"$gt": [{"$ifNull": ["$archaisch_abgeleitet.denisova", 0]}, 0]}}},
    {"$group": {"_id": {"nean_traeger": "$nean_traeger", "denisova_traegt": "$denisova_traegt"},
                "kandidaten": {"$sum": 1}}},
    {"$sort": {"_id.nean_traeger": -1, "_id.denisova_traegt": -1}},
    {"$project": {"_id": 0, "nean_traeger": "$_id.nean_traeger", "denisova_traegt": "$_id.denisova_traegt",
                  "kandidaten": 1}},
])

client.close()
