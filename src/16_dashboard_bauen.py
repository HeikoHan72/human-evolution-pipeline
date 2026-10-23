#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Heiko Han
"""

# %% 16 - Dashboard bauen: eine einzige HTML-Datei, die ohne Internet funktioniert
# Liest (nur lesen): src/dashboard_vorlage.html und Ausgabe/dashboard/daten.json (aus Skript 15)
# Schreibt: Ausgabe/dashboard/dashboard.html (ein Neulauf ueberschreibt nur diese Datei)
# Die Vorlage enthaelt Gestaltung, Programmcode (JavaScript, ohne externe Bibliotheken) und die Weltkarte
# (Natural Earth, gemeinfrei). Die Daten werden an der Stelle __DATEN__ eingesetzt.

import importlib
import json
from pathlib import Path

cfg = importlib.import_module("00_config")
importlib.reload(cfg)

VORLAGE = Path(cfg.CODE) / "dashboard_vorlage.html"
DATEN = Path(cfg.AUSGABE) / "dashboard" / "daten.json"
ZIEL = Path(cfg.AUSGABE) / "dashboard" / "dashboard.html"

vorlage = VORLAGE.read_text(encoding="utf-8")
if vorlage.count("__DATEN__") != 1:
    raise SystemExit("Platzhalter __DATEN__ in der Vorlage nicht genau einmal gefunden")
daten = json.loads(DATEN.read_text(encoding="utf-8"))      # prueft, dass die JSON-Datei lesbar ist
html = vorlage.replace("__DATEN__", json.dumps(daten, ensure_ascii=False, separators=(",", ":")))
ZIEL.write_text(html, encoding="utf-8")
print(f"Geschrieben: {ZIEL} ({ZIEL.stat().st_size / 1024:.0f} KB), Kandidaten: {daten['kpi']['kandidaten']}")
