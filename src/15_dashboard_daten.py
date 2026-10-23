#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Heiko Han
"""

# %% 15 - Dashboard-Daten: aus den Ergebnissen eine kleine JSON-Datei bauen
# Liest (nur lesen): Ausgabe/chromosomen.parquet (Gruppenzahlen), Ausgabe/scan_1mb.csv,
#                    Ausgabe/export/kandidaten.csv, Ausgabe/export/haeufigkeiten.csv
# Schreibt: Ausgabe/dashboard/daten.json (neuer Ordner; ein Neulauf ueberschreibt nur diese Datei)
# Kein Spark noetig: von den 61,6 Mio SNVs werden nur die Spalten pos und gruppe je Chromosom gelesen.

import importlib
import json
from datetime import datetime
from pathlib import Path
import pandas as pd

cfg = importlib.import_module("00_config")
importlib.reload(cfg)

# %% Einstellungen
AUSGABE = Path(cfg.AUSGABE)
CHROM_PQ = AUSGABE / "chromosomen.parquet"
EXPORT = AUSGABE / "export"
ZIEL = AUSGABE / "dashboard"
KANDIDAT = "geteilt_Neandertaler_nicht_afrikanisch"
AFR_MAX, NICHTAFR_MIN = 0.03, 0.05      # eigene Schwellenwerte der Kandidatengruppe (siehe Dokumentation)
TOP_TABELLE = 300                       # Kandidaten mit der hoechsten mittleren Haeufigkeit bei Nicht-Afrikanern
ARCHAISCH = ["altai", "vindija", "chagyrskaya", "denisova"]

# %% Hilfsfunktionen
def n_abgeleitet(gt, abgeleitet_ist):
    """Anzahl abgeleiteter Allele (0, 1, 2) aus einem Genotyp; fehlend/unklar: None."""
    if not isinstance(gt, str) or "." in gt or abgeleitet_ist not in ("alt", "ref"):
        return None
    ziel = "1" if abgeleitet_ist == "alt" else "0"
    return sum(1 for a in gt.replace("|", "/").split("/") if a == ziel)

def region_von(chrom, pos):
    for name, r in cfg.REGIONEN.items():
        if "chr" + str(r["chrom"]).replace("chr", "") == chrom and int(r["start"]) <= pos <= int(r["ende"]):
            return name
    return None

def gen_von(chrom, pos):
    for name, g in cfg.GENE.items():
        if "chr" + str(g["chrom"]).replace("chr", "") == chrom and int(g["start"]) <= pos <= int(g["ende"]):
            return name
    return None

# %% 1) Gruppenzahlen je Chromosom und in den Pilotregionen (aus dem Parquet, nur zwei Spalten)
gruppen_chrom, gruppen_pilot = {}, {name: {} for name in cfg.REGIONEN}
for nr in range(1, 23):
    ordner = CHROM_PQ / f"chrom=chr{nr}"
    if not (ordner / "_SUCCESS").exists():
        raise SystemExit(f"chr{nr}: Parquet fehlt oder ist nicht fertig ({ordner})")
    df = pd.read_parquet(ordner, columns=["pos", "gruppe"])
    gruppen_chrom[nr] = {str(k): int(v) for k, v in df["gruppe"].value_counts().items()}
    for name, r in cfg.REGIONEN.items():
        if int(str(r["chrom"]).replace("chr", "")) == nr:
            m = df["pos"].between(int(r["start"]), int(r["ende"]))
            gruppen_pilot[name] = {str(k): int(v) for k, v in df.loc[m, "gruppe"].value_counts().items()}
    print(f"chr{nr}: {len(df):>9} SNVs")
    del df
snv_gesamt = sum(sum(g.values()) for g in gruppen_chrom.values())
print("SNVs gesamt:", snv_gesamt, "(erwartet 61599150)")

# %% 2) Scan (1-Mb-Fenster)
scan = pd.read_csv(AUSGABE / "scan_1mb.csv", dtype={"chrom": str})
scan_zeilen = [[int(r.chrom_nr), int(r.start_pos // 1_000_000), int(r.n_snv), int(r.n_kandidaten)]
               for r in scan.itertuples(index=False)]
print("Scan-Fenster:", len(scan_zeilen))

# %% 3) Kandidaten und Haeufigkeiten
textspalten = {c: str for c in ["chrom", "ref", "alt", "id", "anc", "abgeleitet_ist",
                                "altai", "vindija", "denisova", "chagyrskaya"]}
k = pd.read_csv(EXPORT / "kandidaten.csv", dtype=textspalten, keep_default_na=False, na_values=[""])
h = pd.read_csv(EXPORT / "haeufigkeiten.csv",
                dtype={c: str for c in ["chrom", "ref", "alt", "pop", "superpop"]},
                keep_default_na=False, na_values=[""])
print("Kandidaten:", len(k), "| Haeufigkeitszeilen:", len(h))

schluessel = ["chrom", "pos", "ref", "alt"]
je_super = h.groupby(schluessel + ["superpop"])["daf"].mean().unstack("superpop")
nichtafr = h[h["superpop"] != "AFR"].groupby(schluessel)["daf"].mean().rename("nichtafr")
k = k.merge(je_super.reset_index(), on=schluessel, how="left").merge(nichtafr.reset_index(), on=schluessel, how="left")

# Archaische Individuen: Anzahl abgeleiteter Allele (Genotyp umgerechnet auf abgeleitet)
for name in ARCHAISCH:
    k["a_" + name] = [n_abgeleitet(g, ist) for g, ist in zip(k[name], k["abgeleitet_ist"])]
k["nean"] = sum((k["a_" + n].fillna(0) > 0).astype(int) for n in ["altai", "vindija", "chagyrskaya"])
k["den"] = (k["a_denisova"].fillna(0) > 0).astype(int)
k["region"] = [region_von(c, int(p)) for c, p in zip(k["chrom"], k["pos"])]
k["gen"] = [gen_von(c, int(p)) for c, p in zip(k["chrom"], k["pos"])]

# Staerke des Hinweises: wie viele der 3 Neandertaler-Individuen tragen die abgeleitete Base, Denisova ja/nein
staerke = (k.groupby(["nean", "den"]).size().reset_index(name="n")
             .sort_values(["nean", "den"], ascending=False))
staerke_liste = [[int(r.nean), int(r.den), int(r.n)] for r in staerke.itertuples(index=False)]

# Histogramm der mittleren Haeufigkeit bei Nicht-Afrikanern (Klassen von 5 Prozentpunkten)
klassen = pd.cut(k["nichtafr"], bins=[i / 20 for i in range(0, 21)], include_lowest=True)
hist = [int(x) for x in klassen.value_counts(sort=False).values]

# Pilotregionen: Kandidaten je Gen
pilot = {}
for name in cfg.REGIONEN:
    sub = k[k["region"] == name]
    gene = {g: int((sub["gen"] == g).sum()) for g, v in cfg.GENE.items()
            if "chr" + str(v["chrom"]).replace("chr", "") == "chr" + str(cfg.REGIONEN[name]["chrom"]).replace("chr", "")
            and v["start"] >= cfg.REGIONEN[name]["start"] and v["ende"] <= cfg.REGIONEN[name]["ende"]}
    pilot[name] = {"gruppen": gruppen_pilot[name], "gene": gene,
                   "rand": int(sub["gen"].isna().sum()), "kandidaten": int(len(sub))}

# %% 4) Populationen: mittlere Haeufigkeit der abgeleiteten Base (Genom, OAS, TLR)
h = h.merge(k[schluessel + ["region"]], on=schluessel, how="left")
pops = h.groupby(["superpop", "pop"])["daf"].mean().rename("genom").reset_index()
for name in ("OAS", "TLR"):
    m = h[h["region"] == name].groupby("pop")["daf"].mean().rename(name)
    pops = pops.merge(m.reset_index(), on="pop", how="left")
pop_liste = [{"pop": r.pop, "superpop": r.superpop, "genom": round(float(r.genom), 4),
              "OAS": round(float(r.OAS), 4), "TLR": round(float(r.TLR), 4)} for r in pops.itertuples(index=False)]

# %% 5) Tabelle: alle Kandidaten der Pilotregionen plus die Top-Kandidaten
top = k.sort_values("nichtafr", ascending=False).head(TOP_TABELLE)
tab = pd.concat([k[k["region"].notna()], top]).drop_duplicates(subset=schluessel).sort_values("nichtafr", ascending=False)
def r3(x):
    return None if pd.isna(x) else round(float(x), 3)
tabelle = [{"id": f"{r.chrom}:{int(r.pos)}:{r.ref}>{r.alt}", "chrom": r.chrom, "pos": int(r.pos),
            "gen": r.gen if isinstance(r.gen, str) else None, "region": r.region if isinstance(r.region, str) else None,
            "afr": r3(r.AFR), "amr": r3(r.AMR), "eas": r3(r.EAS), "eur": r3(r.EUR), "sas": r3(r.SAS),
            "nichtafr": r3(r.nichtafr), "nean": int(r.nean), "den": int(r.den)}
           for r in tab.itertuples(index=False)]
print("Tabellenzeilen:", len(tabelle))

# Megabase-Fenster der Pilotregionen (fuer die Markierung im Genomscan)
region_fenster = {name: {"chrom": int(str(r["chrom"]).replace("chr", "")),
                         "mb": list(range(int(r["start"]) // 1_000_000, int(r["ende"]) // 1_000_000 + 1))}
                  for name, r in cfg.REGIONEN.items()}

# %% 6) JSON schreiben
daten = {
    "meta": {"erstellt": datetime.now().strftime("%d.%m.%Y %H:%M"), "kandidatengruppe": KANDIDAT,
             "afr_max": AFR_MAX, "nichtafr_min": NICHTAFR_MIN, "top_tabelle": TOP_TABELLE},
    "kpi": {"snv": int(snv_gesamt), "kandidaten": int(len(k)), "fenster": len(scan_zeilen),
            "populationen": int(h["pop"].nunique()), "chromosomen": 22,
            "pilot_kandidaten": int(sum(p["kandidaten"] for p in pilot.values()))},
    "gruppen_chrom": gruppen_chrom,
    "pilot": pilot,
    "region_fenster": region_fenster,
    "scan": scan_zeilen,
    "pop": pop_liste,
    "staerke": staerke_liste,
    "hist_nichtafr": hist,
    "tabelle": tabelle,
}
ZIEL.mkdir(exist_ok=True)
ausgabe = ZIEL / "daten.json"
ausgabe.write_text(json.dumps(daten, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
print(f"Geschrieben: {ausgabe} ({ausgabe.stat().st_size / 1024:.0f} KB)")
