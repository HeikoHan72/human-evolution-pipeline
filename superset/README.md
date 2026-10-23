# Superset dashboard

Apache Superset 6.0.0 on the PostgreSQL database `genom` (built by `src/12_postgres_genom.py`).
It is an addition to the self-contained offline dashboard in `docs/index.html`, which needs no database.

![Dashboard](dashboard.jpg)

## Charts

| Chart | Source | Query idea |
|---|---|---|
| Mean candidate-allele frequency per population (%) | table `haeufigkeiten` | `AVG(daf)` per `pop` (26 populations) |
| Candidate alleles per chromosome | table `scan_1mb` | `SUM(n_kandidaten)` per `chrom_nr` |
| Candidate alleles per pilot gene | tables `gene`, `kandidaten` | join on chromosome and position within the gene |

The per-gene counts only cover positions inside the gene bodies. The region totals (OAS 51, TLR 30) are
higher because they cover the whole region, so the gene bars do not add up to them.

The SQL of each chart is stored in `dashboard_export_20261023T120515/datasets/genom/*.yaml`.

## Run it

Requirements: Docker, the database `genom` filled by script 12, and about 5 to 8 GB of disk space
(rough estimate).

```bash
git clone https://github.com/apache/superset superset-app
cd superset-app
git checkout tags/6.0.0
docker compose -f docker-compose-image-tag.yml up -d
```

This is the quickstart from the [Superset documentation](https://superset.apache.org/user-docs/quickstart).
Open <http://localhost:8088> (default login of that quickstart: `admin` / `admin`, for local use only).

Make the PostgreSQL container reachable from Superset. If your database runs in a container, attach it to the
Superset network (the container name `bde-postgres` is from my setup, replace it with yours):

```bash
docker network connect superset-app_default bde-postgres
```

Then add the database in Superset: **Settings, Database Connections, + Database, PostgreSQL**, host
`bde-postgres`, port `5432`, database `genom`, user and password of your database.

## Import the exported dashboard

The folder `dashboard_export_20261023T120515/` is a Superset export (database passwords are masked).
Zip it and import it under **Dashboards, Import**:

```bash
cd superset
zip -r dashboard_export.zip dashboard_export_20261023T120515
```

Superset asks for the database password during the import. The export refers to the database connection
`postgresql+psycopg2://postgres:...@bde-postgres:5432/genom`; if your host name differs, adjust the connection
afterwards under Database Connections.

I did not test the import into a fresh Superset instance.

## Stop and clean up

```bash
cd superset-app
docker compose -f docker-compose-image-tag.yml down                       # stop, keep data
docker compose -f docker-compose-image-tag.yml down --volumes --rmi all   # remove everything
```

## Limits

- Three charts only. Candidate strength (how many Neandertals carry the allele) and the genome scan per 1 Mb
  window from the offline dashboard are not rebuilt here.
- The chart bars are sorted alphabetically by population, not by value.
- The Superset quickstart images are not meant for production use.
