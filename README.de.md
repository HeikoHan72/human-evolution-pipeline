# Project Human Evolution (Kurzfassung)

Abschlussprojekt des Lehrgangs AI Engineer, Modul 5 Big Data Engineering. Autor: Heiko Han.
Englische Hauptdatei: [README.md](README.md). Vollständige Dokumentation: [docs/Projekt_Evolution_Dokumentation.pdf](docs/Projekt_Evolution_Dokumentation.pdf).

**Frage:** Welche DNA-Varianten heutiger Menschen sehen nach Neandertaler-Erbe aus? Eine Spark-Pipeline verknüpft
61,6 Mio. SNVs der 1000 Genomes (3.202 Personen, 26 Populationen) mit dem Ahnenzustand (Ensembl) und vier
archaischen Genomen. Getestet werden die Immun-Genregionen OAS1/2/3 (Chr 12, Positivkontrolle), TLR1/6/10 (Chr 4)
und STAT2 (Chr 12), dazu ein Scan über alle 22 Autosomen.

**Ergebnis:** 88.632 Kandidaten-Allele (0,14 %), davon 81 in den Pilotregionen (OAS 51, TLR 30, STAT2 0).
Die Kandidaten sind bei Nicht-Afrikanern im Mittel mit 8,6 bis 11,8 % häufig, bei afrikanischen Populationen mit
0,3 bis 2,6 %. Das sind Hinweise auf Kandidaten-Allele in Regionen, kein genomweiter Neandertaler-Anteil.

**Aufbau:** Skripte 08 bis 14 (Spark, PostgreSQL, MongoDB, Neo4j), Airflow-DAG `src/genom_pipeline.py`,
Offline-Dashboard `docs/index.html` (Skripte 15 und 16). Die Rohdaten (ca. 80 GB) sind nicht enthalten, Quellen
stehen in der Dokumentation. Zugangsdaten kommen aus einer `.env`-Datei (Vorlage `.env.example`).

**Grenzen:** zwei Personenmengen (2.504 Unverwandte für die Regel, 2.590 Gründer für die Häufigkeiten), Literatur
nur für OAS geprüft, Spark-Tasks nicht über den DAG getestet, Schwellen 3 % und 5 % sind eigene Entscheidungen.
Die vollständige Liste steht im Abschnitt "Limitations" der englischen README.
