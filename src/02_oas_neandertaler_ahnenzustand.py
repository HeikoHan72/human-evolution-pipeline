#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 11:05:11 2026

@author: Heiko Han
"""

# %% 02_oas_neandertaler_ahnenzustand.py - OAS-Cluster: Neandertaler-Genotypen plus Ahnenzustand
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

# %% Daten holen
R = cfg.REGIONEN["OAS"]
df = gf.lese_neandertaler(R["chrom"], R["start"], R["ende"])
seq = gf.lese_ahnenzustand(R["chrom"], R["start"], R["ende"])
print("Neandertaler-Eintraege:", len(df), "| Laenge Ahnenzustand:", len(seq),
      "| erwartet:", R["ende"] - R["start"] + 1)

# %% Ahnenzustand an jede Position haengen
df["anc"] = [seq[p - R["start"]] for p in df["pos"]]
print("\nZeichen des Ahnenzustands an unseren Positionen:")
print(df["anc"].value_counts().to_string())

# %% Kontrolle: stimmt die Koordinate?  An Nicht-Variantenpositionen sollte ref == Ahnenzustand sein
ohne = df[(df["alt"] == ".") & df["anc"].isin(list("ACGT"))]
anteil = (ohne["ref"] == ohne["anc"]).mean()
print(f"\nKontrolle: ref == Ahnenzustand an {len(ohne)} Nicht-Variantenpositionen: {anteil:.1%}")
print("(Wenn das nur ca. 25 % waere, waeren die Koordinaten verschoben.)")

# %% Variantenpositionen: welche Base ist ancestral?
var = df[df["alt"] != "."].copy()

def status(z):
    a = z["anc"].upper()
    if a not in "ACGT" or a == "":
        return "kein_Ahnenzustand"
    if a == z["ref"]:
        return "Referenz=ancestral"
    if a in z["alt"].split(","):
        return "Alt=ancestral"
    return "weder_noch"

var["status"] = var.apply(status, axis=1)
print("\nVariantenpositionen nach Ahnenzustand:")
print(var["status"].value_counts().to_string())

# %% Die 36 Positionen mit Alt-Allel in einem Neandertaler
nean = var[var.altai.str.contains("1") | var.vindija.str.contains("1")
           | var.chagyrskaya.str.contains("1")]
print("\nPositionen mit Alt-Allel in Altai/Vindija/Chagyrskaya:", len(nean))
print(nean[["pos", "ref", "alt", "anc", "status", "altai", "vindija", "denisova",
            "chagyrskaya"]].to_string(index=False))