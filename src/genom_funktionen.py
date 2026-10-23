#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 11:01:09 2026

@author: Heiko Han
"""

# genom_funktionen.py - gemeinsame Lesefunktionen fuer das Human-Genom-Projekt
# Wird von den nummerierten Skripten importiert.

import importlib
import pandas as pd
from pysam import VariantFile

cfg = importlib.import_module("00_config")

# Kurznamen der archaischen Proben (Spark-tauglich: keine Punkte/Bindestriche)
KURZ = {"AltaiNeandertal": "altai", "Vindija33.19": "vindija",
        "Denisova": "denisova", "Chagyrskaya-Phalanx": "chagyrskaya"}


def gt_text(gt):
    """Genotyp als Text. Fehlende Werte werden './.' und nie '0/0'."""
    if gt is None or any(a is None for a in gt):
        return "./."
    return "/".join(str(a) for a in gt)


def lese_neandertaler(chrom, start, ende):
    """Archaische Genotypen aus dem Zenodo-BCF (hg38). start/ende 1-basiert, inklusive."""
    zeilen = []
    with VariantFile(str(cfg.neandertaler_bcf(chrom))) as vf:
        for rec in vf.fetch(f"chr{chrom}", start - 1, ende):
            zeile = {"chrom": rec.chrom, "pos": rec.pos, "ref": rec.ref,
                     "alt": ",".join(rec.alts) if rec.alts else "."}
            for probe, kurz in KURZ.items():
                zeile[kurz] = gt_text(rec.samples[probe]["GT"])
            zeilen.append(zeile)
    return pd.DataFrame(zeilen)


def lese_ahnenzustand(chrom, start, ende):
    """Ahnenzustand (Ensembl EPO) als String, ein Zeichen pro Position von start bis ende (1-basiert).
    Zeichen: ACGT sicher, acgt unsicher, N unbekannt, '-' Einfuegung, '.' nicht abgedeckt."""
    pfad = (cfg.AHNENZUSTAND / "entpackt" / "homo_sapiens_ancestor_GRCh38"
            / f"homo_sapiens_ancestor_{chrom}.fa")
    with open(pfad, "rb") as f:
        kopf = f.readline()                        # Kopfzeile ">ANCESTOR_for_chromosome..."
        zeilenlaenge = len(f.readline().rstrip(b"\n"))   # Zeichen pro Zeile (hier 100)

        def versatz(pos):                          # Byte-Position von pos (1-basiert)
            k = pos - 1
            return len(kopf) + k + k // zeilenlaenge

        f.seek(versatz(start))
        roh = f.read(versatz(ende) - versatz(start) + 1)
    return roh.replace(b"\n", b"").decode()




# ---------------------------------------------------------------------------
# Ergaenzung (Skript 03): Ahnenzustand-Status und 1000-Genomes-Haeufigkeiten
# ---------------------------------------------------------------------------
GROSSGRUPPEN = ["AFR", "AMR", "EAS", "EUR", "SAS"]


def ahnenzustand_status(ref, alt, anc):
    """Welche Base ist laut Ahnenzustand ancestral? (anc = Zeichen aus der Ensembl-Datei)"""
    a = anc.upper()
    if a not in ("A", "C", "G", "T"):
        return "kein_Ahnenzustand"
    if a == ref:
        return "Referenz=ancestral"
    if a in alt.split(","):
        return "Alt=ancestral"
    return "weder_noch"


def _erster_wert(wert):
    """INFO-Felder mit Number=A kommen als Tuple; wir nehmen den ersten Wert."""
    if wert is None:
        return None
    if isinstance(wert, (tuple, list)):
        return wert[0] if len(wert) > 0 else None
    return wert


def lese_tausend_genome_info(chrom, start, ende):
    """Allelhaeufigkeiten je Grossgruppe (AF_xxx_unrel = 2504 unverwandte Personen).
    Nur biallelische SNVs; Indels, Strukturvarianten und Multiallele werden uebersprungen."""
    zeilen = []
    with VariantFile(str(cfg.tausend_genome_vcf(chrom))) as vf:
        for rec in vf.fetch(f"chr{chrom}", start - 1, ende):
            if rec.alts is None or len(rec.alts) != 1:
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue
            zeile = {"pos": rec.pos, "ref": rec.ref, "alt": rec.alts[0], "id": rec.id}
            zeile["af_alle"] = _erster_wert(rec.info.get("AF"))
            for g in GROSSGRUPPEN:
                zeile[f"af_{g.lower()}"] = _erster_wert(rec.info.get(f"AF_{g}_unrel"))
            zeilen.append(zeile)
    return pd.DataFrame(zeilen)




# ---------------------------------------------------------------------------
# Ergaenzung (Skript 05): Ust'-Ishim (hg19) lesen und nach hg38 umrechnen
# ---------------------------------------------------------------------------
def lese_ust_ishim(chrom, start38, ende38, min_dp=10, min_gq=30):
    """Ust'-Ishim-Genotypen fuer eine hg38-Region. Positionen werden von hg19 nach hg38 umgerechnet.
    Genotypen mit zu geringer Tiefe/Qualitaet werden './.' (nicht abgedeckt). Nur Strang '+'."""
    from pyliftover import LiftOver
    kette = cfg.DATEN / "Liftover_Chain"
    nach19 = LiftOver(str(kette / "hg38ToHg19.over.chain.gz"))
    nach38 = LiftOver(str(kette / "hg19ToHg38.over.chain.gz"))

    a = nach19.convert_coordinate(f"chr{chrom}", start38 - 1)
    b = nach19.convert_coordinate(f"chr{chrom}", ende38 - 1)
    if not a or not b:
        raise ValueError("Regionsrand laesst sich nicht nach hg19 umrechnen")
    start19, ende19 = a[0][1] + 1, b[0][1] + 1          # wieder 1-basiert

    zeilen = []
    with VariantFile(str(cfg.ust_ishim_vcf(chrom))) as vf:
        for rec in vf.fetch(str(chrom), start19 - 1, ende19):     # Kontigname hier ohne 'chr'
            res = nach38.convert_coordinate(f"chr{chrom}", rec.pos - 1)
            if not res or res[0][2] != "+":
                continue
            s = rec.samples["Ust_Ishim"]
            dp = s["DP"] if "DP" in s else None
            gq = s["GQ"] if "GQ" in s else None
            if dp is None or gq is None or dp < min_dp or gq < min_gq:
                gt = "./."
            else:
                gt = gt_text(s["GT"])
            zeilen.append({"pos_u": res[0][1] + 1, "ref_u": rec.ref,
                           "alt_u": ",".join(rec.alts) if rec.alts else ".",
                           "ust": gt, "dp_u": float(dp) if dp is not None else None,
                           "gq_u": float(gq) if gq is not None else None})
    return pd.DataFrame(zeilen)