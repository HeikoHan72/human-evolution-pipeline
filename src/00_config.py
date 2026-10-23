#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 10:48:27 2026

@author: Heiko Han
"""

# %% 00_config.py - zentrale Einstellungen fuer das Human-Genom-Projekt
# Alle Pfade, Regionen und Probennamen stehen NUR hier. Alle anderen Skripte importieren diese Datei.

import os
from pathlib import Path

# %% Pfade
# CODE: Ordner dieser Datei (src/). Bei zellenweiser Ausfuehrung ohne __file__ (z. B. Spyder) der aktuelle Ordner.
try:
    CODE = Path(__file__).resolve().parent
except NameError:
    CODE = Path.cwd()
REPO = CODE.parent


def lade_dotenv(datei=REPO / ".env"):
    """Liest KEY=WERT-Zeilen aus .env (falls vorhanden) in os.environ. Bereits gesetzte Variablen gewinnen."""
    if not datei.exists():
        return
    for zeile in datei.read_text(encoding="utf-8").splitlines():
        zeile = zeile.strip()
        if zeile and not zeile.startswith("#") and "=" in zeile:
            schluessel, wert = zeile.split("=", 1)
            os.environ.setdefault(schluessel.strip(), wert.strip().strip('"').strip("'"))


lade_dotenv()


def env(name, default=None):
    """Umgebungsvariable lesen. Ohne Standardwert ist sie Pflicht (klare Fehlermeldung statt Passwort im Code)."""
    wert = os.environ.get(name, default)
    if wert is None:
        raise SystemExit(f"Umgebungsvariable {name} fehlt. Vorlage: .env.example nach .env kopieren und ausfuellen.")
    return wert


# GENOM_BASIS: Ordner mit Datensaetze/, Ausgabe/ und Zwischen/ (Standard: <Repository>/data)
BASIS = Path(os.environ.get("GENOM_BASIS", REPO / "data"))
DATEN = BASIS / "Datensätze"
AUSGABE = BASIS / "Ausgabe"          # hier landen Parquet- und JSON-Ergebnisse

NEANDERTALER = DATEN / "Neandertaler_hg38"
UST_ISHIM = DATEN / "UstIshim_hg19"
AHNENZUSTAND = DATEN / "Ahnenzustand_GRCh38"
TAUSEND_GENOME = DATEN / "1000Genomes_GRCh38"

# %% Dateinamen pro Chromosom (chrom = Zahl 1 bis 22)
def neandertaler_bcf(chrom):
    return NEANDERTALER / f"individuals_highcov.{chrom}.bcf"

def ust_ishim_vcf(chrom):
    return UST_ISHIM / f"chr{chrom}_mq25_mapab100.vcf.gz"

def tausend_genome_vcf(chrom):
    return TAUSEND_GENOME / (f"1kGP_high_coverage_Illumina.chr{chrom}."
                             "filtered.SNV_INDEL_SV_phased_panel.vcf.gz")

# %% Proben der archaischen Genome (Reihenfolge wie in der BCF-Datei)
ARCHAISCHE_PROBEN = ["AltaiNeandertal", "Vindija33.19", "Denisova", "Chagyrskaya-Phalanx"]

# %% Gene (GRCh38, 1-basiert, aus NCBI GCF_000001405.40, GTF)
GENE = {
    "TLR10": {"chrom": 4,  "start": 38772238,  "ende": 38782990,  "strang": "-"},
    "TLR1":  {"chrom": 4,  "start": 38787569,  "ende": 38805644,  "strang": "-"},
    "TLR6":  {"chrom": 4,  "start": 38822897,  "ende": 38868390,  "strang": "-"},
    "STAT2": {"chrom": 12, "start": 56341597,  "ende": 56360107,  "strang": "-"},
    "OAS1":  {"chrom": 12, "start": 112906962, "ende": 112933219, "strang": "+"},
    "OAS3":  {"chrom": 12, "start": 112938474, "ende": 112973251, "strang": "+"},
    "OAS2":  {"chrom": 12, "start": 112978519, "ende": 113011723, "strang": "+"},
}

# %% Pilot-Regionen (Cluster-Fenster ueber die Gene, ohne Flanken)
REGIONEN = {
    "TLR":   {"chrom": 4,  "start": 38772238,  "ende": 38868390},
    "STAT2": {"chrom": 12, "start": 56341597,  "ende": 56360107},
    "OAS":   {"chrom": 12, "start": 112906962, "ende": 113011723},
}

# %% Selbsttest: Dateien vorhanden?  (laeuft nur, wenn die Datei direkt gestartet wird)
if __name__ == "__main__":
    print("BASIS vorhanden:", BASIS.exists())
    for chrom in (4, 12):
        for name, pfad in [("Neandertaler", neandertaler_bcf(chrom)),
                           ("Ust-Ishim", ust_ishim_vcf(chrom)),
                           ("1000 Genomes", tausend_genome_vcf(chrom))]:
            if pfad.exists():
                print(f"chr{chrom:>2} {name:<13} OK  {pfad.stat().st_size / 1e6:10.1f} MB")
            else:
                print(f"chr{chrom:>2} {name:<13} FEHLT  {pfad}")