#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Heiko Han
"""

# %% 14 - Neo4j: kleiner Graph der Pilotregionen (Gen - Variante - archaisches Individuum - Population)
# Nur die Kandidaten in den Pilotregionen (OAS, TLR; STAT2 ohne Kandidaten) kommen in den Graphen.
# Grund: Pfadabfragen sind hier anschaulich; die Massendaten liegen in PostgreSQL und MongoDB.
# Liest: Ausgabe/export/kandidaten.csv und haeufigkeiten.csv (nur lesen).
# Schreibt: nur Knoten mit dem zusaetzlichen Label :Genom in Neo4j. Andere Daten in Neo4j (z.B. der
# Kurs-Graph) werden nie angefasst, auch nicht bei NEU_AUFBAUEN = True.

import os
import importlib
from pathlib import Path
import pandas as pd
from neo4j import GraphDatabase

cfg = importlib.import_module("00_config")
importlib.reload(cfg)

# %% Einstellungen
NEO_URI = cfg.env("NEO4J_URI", "bolt://localhost:7687")
NEO_USER = cfg.env("NEO4J_USER", "neo4j")
NEO_PASS = cfg.env("NEO4J_PASSWORD")             # Pflicht, kommt aus .env oder der Umgebung
NEU_AUFBAUEN = os.environ.get("GENOM_NEU_AUFBAUEN") == "1"        # True: nur die Knoten mit Label :Genom loeschen und neu aufbauen
EXPORT = Path(cfg.AUSGABE) / "export"

POPULATIONEN = {
    "AFR": ["ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI"],
    "AMR": ["CLM", "MXL", "PEL", "PUR"],
    "EAS": ["CDX", "CHB", "CHS", "JPT", "KHV"],
    "EUR": ["CEU", "FIN", "GBR", "IBS", "TSI"],
    "SAS": ["BEB", "GIH", "ITU", "PJL", "STU"],
}
# Gene wie in Skript 12 und 13 (spaeter nach 00_config verschieben)
GENE = [
    ("TLR10", "TLR", "chr4", 38772238, 38782990, "-"),
    ("TLR1", "TLR", "chr4", 38787569, 38805644, "-"),
    ("TLR6", "TLR", "chr4", 38822897, 38868390, "-"),
    ("STAT2", "STAT2", "chr12", 56341597, 56360107, "-"),
    ("OAS1", "OAS", "chr12", 112906962, 112933219, "+"),
    ("OAS3", "OAS", "chr12", 112938474, 112973251, "+"),
    ("OAS2", "OAS", "chr12", 112978519, 113011723, "+"),
]
ARCHAISCH = [("altai", "Neandertaler"), ("vindija", "Neandertaler"),
             ("chagyrskaya", "Neandertaler"), ("denisova", "Denisovaner")]

# %% Hilfsfunktionen
def gen_zu(chrom, pos):
    gen = next((g[0] for g in GENE if g[2] == chrom and g[3] <= pos <= g[4]), None)
    region = next((name for name, r in cfg.REGIONEN.items()
                   if "chr" + str(r["chrom"]).replace("chr", "") == chrom
                   and int(r["start"]) <= pos <= int(r["ende"])), None)
    return gen, region

def n_abgeleitet(gt, abgeleitet_ist):
    """Anzahl abgeleiteter Allele (0, 1, 2) aus einem Genotyp; fehlend/unklar: None."""
    if not isinstance(gt, str) or "." in gt or abgeleitet_ist not in ("alt", "ref"):
        return None
    ziel = "1" if abgeleitet_ist == "alt" else "0"
    return sum(1 for a in gt.replace("|", "/").split("/") if a == ziel)

# %% Kandidaten der Pilotregionen aus der CSV waehlen
textspalten = {c: str for c in ["chrom", "ref", "alt", "id", "anc", "abgeleitet_ist",
                                "altai", "vindija", "denisova", "chagyrskaya"]}
k = pd.read_csv(EXPORT / "kandidaten.csv", dtype=textspalten, keep_default_na=False, na_values=[""])
varianten = []
for z in k.itertuples(index=False):
    gen, region = gen_zu(z.chrom, int(z.pos))
    if region is None:
        continue
    vid = f"{z.chrom}:{int(z.pos)}:{z.ref}:{z.alt}"
    ist = None if pd.isna(z.abgeleitet_ist) else z.abgeleitet_ist
    traeger = []
    for name, _ in ARCHAISCH:
        kopien = n_abgeleitet(getattr(z, name), ist)
        if kopien:                               # nur Individuen mit mindestens einer abgeleiteten Base
            traeger.append({"individuum": name, "kopien": kopien})
    varianten.append({"id": vid, "chrom": z.chrom, "pos": int(z.pos), "ref": z.ref, "alt": z.alt,
                      "abgeleitet_ist": ist, "daf_nichtafr": None if pd.isna(z.daf_nichtafr) else float(z.daf_nichtafr),
                      "gen": gen, "region": region, "traeger": traeger})
ids = {v["id"] for v in varianten}
print("Kandidaten in den Pilotregionen:", len(varianten), "(erwartet 81: OAS 51 + TLR 30)")

# Haeufigkeiten nur fuer diese Kandidaten (Datei in Stuecken lesen)
teile = []
for stueck in pd.read_csv(EXPORT / "haeufigkeiten.csv", dtype={c: str for c in ["chrom", "ref", "alt", "pop", "superpop"]},
                          keep_default_na=False, na_values=[""], chunksize=500_000):
    stueck["id"] = stueck["chrom"] + ":" + stueck["pos"].astype(str) + ":" + stueck["ref"] + ":" + stueck["alt"]
    teile.append(stueck[stueck["id"].isin(ids)])
h = pd.concat(teile)
haeufig = [{"id": r.id, "pop": r.pop, "daf": None if pd.isna(r.daf) else float(r.daf),
            "n_all": None if pd.isna(r.n_all) else int(r.n_all)} for r in h.itertuples(index=False)]
print("Haeufigkeiten fuer diese Kandidaten:", len(haeufig), "(erwartet 81 * 26 = 2106)")

# %% Verbindung und Schutz
driver = GraphDatabase.driver(NEO_URI, auth=(NEO_USER, NEO_PASS))
driver.verify_connectivity()

def lies(cypher, **p):
    with driver.session() as s:
        return s.run(cypher, **p).data()

def schreibe(cypher, **p):
    with driver.session() as s:
        s.execute_write(lambda tx: tx.run(cypher, **p).consume())

vorhanden = lies("MATCH (n:Genom) RETURN count(n) AS n")[0]["n"]
if vorhanden > 0 and not NEU_AUFBAUEN:
    raise SystemExit(f"Es gibt schon {vorhanden} Knoten mit Label :Genom. Zum Neuaufbau NEU_AUFBAUEN = True setzen.")
if vorhanden > 0:
    schreibe("MATCH (n:Genom) DETACH DELETE n")
    print("Knoten mit Label :Genom geloescht, Neuaufbau")

# %% Constraints (Eindeutigkeit, dazu automatisch ein Index)
for label, prop in [("Variante", "id"), ("Gen", "name"), ("Region", "name"),
                    ("Population", "name"), ("Superpopulation", "name"), ("Archaisch", "name")]:
    schreibe(f"CREATE CONSTRAINT IF NOT EXISTS FOR (n:{label}) REQUIRE n.{prop} IS UNIQUE")

# %% Stammdaten: Regionen, Gene, Populationen, archaische Individuen
regionen = [{"name": n, "chrom": "chr" + str(r["chrom"]).replace("chr", ""),
             "start": int(r["start"]), "ende": int(r["ende"])} for n, r in cfg.REGIONEN.items()]
schreibe("UNWIND $r AS x CREATE (:Region:Genom {name: x.name, chrom: x.chrom, start: x.start, ende: x.ende})", r=regionen)

gene = [{"name": g[0], "region": g[1], "chrom": g[2], "start": g[3], "ende": g[4], "strang": g[5]} for g in GENE]
schreibe("""UNWIND $g AS x
            MATCH (r:Region {name: x.region})
            CREATE (:Gen:Genom {name: x.name, chrom: x.chrom, start: x.start, ende: x.ende, strang: x.strang})
                   -[:GEHOERT_ZU]->(r)""", g=gene)
# (CREATE mit Pfad legt Knoten und Beziehung in einem Schritt an)

pops = [{"pop": p, "sup": s} for s, ps in POPULATIONEN.items() for p in ps]
schreibe("""UNWIND $p AS x
            MERGE (s:Superpopulation:Genom {name: x.sup})
            CREATE (:Population:Genom {name: x.pop, superpop: x.sup})-[:TEIL_VON]->(s)""", p=pops)
# Beziehung Population -> Superpopulation: der Pfad im CREATE verwendet den gefundenen Knoten s.

schreibe("UNWIND $a AS x CREATE (:Archaisch:Genom {name: x.name, art: x.art})",
         a=[{"name": n, "art": a} for n, a in ARCHAISCH])

# %% Varianten und ihre Beziehungen
schreibe("""UNWIND $v AS x
            MATCH (r:Region {name: x.region})
            CREATE (v:Variante:Genom {id: x.id, chrom: x.chrom, pos: x.pos, ref: x.ref, alt: x.alt,
                                      abgeleitet_ist: x.abgeleitet_ist, daf_nichtafr: x.daf_nichtafr})
            CREATE (v)-[:LIEGT_IN_REGION]->(r)""", v=varianten)
schreibe("""UNWIND $v AS x
            WITH x WHERE x.gen IS NOT NULL
            MATCH (v:Variante {id: x.id}), (g:Gen {name: x.gen})
            CREATE (v)-[:LIEGT_IN]->(g)""", v=varianten)
schreibe("""UNWIND $v AS x
            UNWIND x.traeger AS t
            MATCH (v:Variante {id: x.id}), (a:Archaisch {name: t.individuum})
            CREATE (v)-[:TRAEGT_ABGELEITET {kopien: t.kopien}]->(a)""", v=varianten)
schreibe("""UNWIND $h AS x
            MATCH (p:Population {name: x.pop}), (v:Variante {id: x.id})
            CREATE (p)-[:HAT_ALLEL {daf: x.daf, n_all: x.n_all}]->(v)""", h=haeufig)

# %% Kontrolle
print("\nKnoten je Label:")
for z in lies("MATCH (n:Genom) UNWIND labels(n) AS l WITH l WHERE l <> 'Genom' RETURN l AS label, count(*) AS anzahl ORDER BY label"):
    print(f"  {z['label']:16s}{z['anzahl']}")
print("Beziehungen je Typ:")
for z in lies("MATCH (:Genom)-[r]->(:Genom) RETURN type(r) AS typ, count(r) AS anzahl ORDER BY typ"):
    print(f"  {z['typ']:16s}{z['anzahl']}")

# %% Abfragen (Cypher)
def zeige(titel, cypher):
    print("\n" + "=" * 70 + "\n" + titel + "\n" + "=" * 70)
    print(pd.DataFrame(lies(cypher)).to_string(index=False))

zeige("1) Kandidaten je Gen (STAT2 = negativer Befund, OPTIONAL MATCH)", """
MATCH (g:Gen)
OPTIONAL MATCH (v:Variante)-[:LIEGT_IN]->(g)
RETURN g.name AS gen, count(v) AS kandidaten ORDER BY gen""")

zeige("2) Welche archaischen Individuen tragen die abgeleitete Base? (Pfad Gen - Variante - Individuum)", """
MATCH (g:Gen)<-[:LIEGT_IN]-(v:Variante)-[t:TRAEGT_ABGELEITET]->(a:Archaisch)
RETURN g.name AS gen, a.name AS individuum, count(v) AS varianten ORDER BY gen, individuum""")

zeige("3) Mittlere Haeufigkeit je Region und Superpopulation (Gegenprobe zu PostgreSQL und MongoDB)", """
MATCH (r:Region)<-[:LIEGT_IN_REGION]-(v:Variante)<-[h:HAT_ALLEL]-(p:Population)-[:TEIL_VON]->(s:Superpopulation)
RETURN r.name AS region, s.name AS superpop, count(DISTINCT v) AS kandidaten, round(avg(h.daf), 3) AS mittel_daf
ORDER BY region, superpop""")

zeige("4) Varianten in Genen, die alle drei Neandertaler tragen: mittlere Haeufigkeit in Europa", """
MATCH (v:Variante)-[:LIEGT_IN]->(g:Gen)
MATCH (v)-[:TRAEGT_ABGELEITET]->(a:Archaisch {art: 'Neandertaler'})
WITH v, g, count(a) AS neandertaler WHERE neandertaler = 3
MATCH (:Superpopulation {name: 'EUR'})<-[:TEIL_VON]-(p:Population)-[h:HAT_ALLEL]->(v)
RETURN g.name AS gen, v.id AS variante, round(avg(h.daf), 3) AS mittel_daf_eur
ORDER BY mittel_daf_eur DESC LIMIT 10""")

driver.close()
