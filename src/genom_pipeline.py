# -*- coding: utf-8 -*-
"""
DAG genom_pipeline: Projekt Evolution (Human-Genom-Projekt)

Ruft die nummerierten Skripte 08 bis 14 nacheinander auf. Die Skripte selbst laufen im Environment
bde-kurs (Spark, pandas, Datenbanktreiber); Airflow laeuft getrennt in bde-kurs-airflow
(SQLAlchemy-Konflikt, siehe Tag 13). Deshalb startet jeder Task das Skript als eigenen Prozess
mit dem Python aus bde-kurs, statt es im Airflow-Prozess zu importieren.

Wiederholbar:
- 08 und 10 ueberspringen Chromosomen, deren Parquet-Ordner fertig ist (_SUCCESS).
- 09 und 11 ueberschreiben nur ihre eigenen Ausgabedateien.
- 12, 13, 14 bauen ihre Datenbestaende neu auf (GENOM_NEU_AUFBAUEN=1). Geloescht wird dabei nur,
  was die Skripte selbst angelegt haben: Schema public in der Datenbank genom (PostgreSQL),
  die Collection genom.kandidaten (MongoDB), Knoten mit Label :Genom (Neo4j).

Ein erster Lauf ueber alle 22 Chromosomen dauert viele Stunden (08 und 10 sind Spark-Laeufe).
Fuer Tests einzelne Tasks verwenden: airflow tasks test genom_pipeline neo4j_14 2026-10-21
"""

import os
from datetime import datetime, timedelta
from pathlib import Path

from airflow import DAG
from airflow.operators.bash import BashOperator

# %% Einstellungen (nur hier anpassen)
CODE = os.environ.get("GENOM_CODE", str(Path(__file__).resolve().parent))   # Ordner mit den Skripten (src/)

# Python-Umgebung der Skripte. Standard: Schwester-Environment bde-kurs neben dem aktiven Airflow-Environment.
# Ueberschreibbar mit der Umgebungsvariable GENOM_ENV (Pfad zum Environment-Ordner).
ENV_SETUP = (
    'ENV="${GENOM_ENV:-$CONDA_PREFIX/../bde-kurs}"; '
    'export PATH="$ENV/bin:$PATH" JAVA_HOME="$ENV"; '
)

SKRIPTE = {
    "chromosomen_08": ("08_alle_chromosomen.py", timedelta(hours=14)),
    "scan_09": ("09_scan_pro_megabase.py", timedelta(hours=1)),
    "kandidaten_pop_10": ("10_kandidaten_populationen.py", timedelta(hours=14)),
    "export_11": ("11_export_fuer_datenbanken.py", timedelta(hours=1)),
    "postgres_12": ("12_postgres_genom.py", timedelta(minutes=30)),
    "mongodb_13": ("13_mongodb_genom.py", timedelta(minutes=30)),
    "neo4j_14": ("14_neo4j_genom.py", timedelta(minutes=30)),
}
DB_TASKS = {"postgres_12", "mongodb_13", "neo4j_14"}

# %% DAG
with DAG(
    dag_id="genom_pipeline",
    description="Human-Genom-Projekt: Spark-Auswertung, Export, PostgreSQL, MongoDB, Neo4j",
    start_date=datetime(2026, 10, 1),
    schedule=None,          # nur manuell starten
    catchup=False,
    default_args={"retries": 0},
    tags=["genom", "projekt"],
    doc_md=__doc__,
) as dag:

    tasks = {}
    for task_id, (skript, limit) in SKRIPTE.items():
        tasks[task_id] = BashOperator(
            task_id=task_id,
            bash_command=ENV_SETUP + f'cd "{CODE}" && python "{skript}"',
            env={"GENOM_NEU_AUFBAUEN": "1"} if task_id in DB_TASKS else {},
            append_env=True,    # bestehende Umgebung (PATH, CONDA_PREFIX) behalten, nur ergaenzen
            execution_timeout=limit,
        )

    # Reihenfolge: bewusst nacheinander (15 GiB Arbeitsspeicher, kein paralleler Spark- oder Datenbank-Betrieb)
    reihenfolge = list(SKRIPTE)
    for a, b in zip(reihenfolge, reihenfolge[1:]):
        tasks[a] >> tasks[b]
