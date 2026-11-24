"""SQLite persistence: campaigns, batches, measurements, model runs.

A campaign is the unit of work — it survives restarts, and every measurement is
attributable to the cycle and batch that requested it.
"""
from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from pathlib import Path

SCHEMA = """
CREATE TABLE IF NOT EXISTS campaign (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  name TEXT NOT NULL UNIQUE,
  target TEXT,
  profile TEXT NOT NULL,
  settings TEXT NOT NULL,
  created_at TEXT NOT NULL DEFAULT (datetime('now'))
);
CREATE TABLE IF NOT EXISTS batch (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  campaign_id INTEGER NOT NULL REFERENCES campaign(id),
  cycle INTEGER NOT NULL,
  acquisition TEXT NOT NULL,
  plan TEXT NOT NULL,
  diversity TEXT NOT NULL,
  created_at TEXT NOT NULL DEFAULT (datetime('now'))
);
CREATE TABLE IF NOT EXISTS measurement (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  campaign_id INTEGER NOT NULL REFERENCES campaign(id),
  batch_id INTEGER REFERENCES batch(id),
  compound_id TEXT NOT NULL,
  smiles TEXT NOT NULL,
  endpoint TEXT NOT NULL DEFAULT 'pIC50',
  value REAL NOT NULL,
  predicted REAL,
  uncertainty REAL,
  source TEXT NOT NULL,
  created_at TEXT NOT NULL DEFAULT (datetime('now')),
  UNIQUE (campaign_id, compound_id, endpoint)
);
CREATE TABLE IF NOT EXISTS model_run (
  id INTEGER PRIMARY KEY AUTOINCREMENT,
  campaign_id INTEGER NOT NULL REFERENCES campaign(id),
  cycle INTEGER NOT NULL,
  report TEXT NOT NULL,
  prospective TEXT,
  created_at TEXT NOT NULL DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_meas_campaign ON measurement(campaign_id);
"""


class Store:
    def __init__(self, path: str | Path = "campaigns.db"):
        self.path = str(path)
        with self.conn() as cx:
            cx.executescript(SCHEMA)

    @contextmanager
    def conn(self):
        cx = sqlite3.connect(self.path)
        cx.row_factory = sqlite3.Row
        try:
            yield cx
            cx.commit()
        finally:
            cx.close()

    # ------------------------------------------------------- campaigns
    def create_campaign(self, name: str, target: str, profile: list[dict], settings: dict) -> int:
        with self.conn() as cx:
            cur = cx.execute(
                "INSERT INTO campaign (name, target, profile, settings) VALUES (?,?,?,?)",
                (name, target, json.dumps(profile), json.dumps(settings)))
            return int(cur.lastrowid)

    def get_campaign(self, name: str) -> dict | None:
        with self.conn() as cx:
            row = cx.execute("SELECT * FROM campaign WHERE name=?", (name,)).fetchone()
        if not row:
            return None
        d = dict(row)
        d["profile"] = json.loads(d["profile"])
        d["settings"] = json.loads(d["settings"])
        return d

    def list_campaigns(self) -> list[dict]:
        with self.conn() as cx:
            rows = cx.execute(
                "SELECT c.id, c.name, c.target, c.created_at,"
                " (SELECT COUNT(*) FROM measurement m WHERE m.campaign_id=c.id) AS measured,"
                " (SELECT COUNT(*) FROM batch b WHERE b.campaign_id=c.id) AS cycles"
                " FROM campaign c ORDER BY c.id DESC").fetchall()
        return [dict(r) for r in rows]

    def delete_campaign(self, campaign_id: int) -> None:
        with self.conn() as cx:
            cx.execute("DELETE FROM measurement WHERE campaign_id=?", (campaign_id,))
            cx.execute("DELETE FROM model_run WHERE campaign_id=?", (campaign_id,))
            cx.execute("DELETE FROM batch WHERE campaign_id=?", (campaign_id,))
            cx.execute("DELETE FROM campaign WHERE id=?", (campaign_id,))

    # ---------------------------------------------------------- batches
    def add_batch(self, campaign_id: int, cycle: int, acquisition: str,
                  plan: dict, diversity: dict) -> int:
        with self.conn() as cx:
            cur = cx.execute(
                "INSERT INTO batch (campaign_id, cycle, acquisition, plan, diversity) VALUES (?,?,?,?,?)",
                (campaign_id, cycle, acquisition, json.dumps(plan), json.dumps(diversity)))
            return int(cur.lastrowid)

    def batches(self, campaign_id: int) -> list[dict]:
        with self.conn() as cx:
            rows = cx.execute("SELECT * FROM batch WHERE campaign_id=? ORDER BY cycle", (campaign_id,)).fetchall()
        out = []
        for r in rows:
            d = dict(r)
            d["plan"] = json.loads(d["plan"])
            d["diversity"] = json.loads(d["diversity"])
            out.append(d)
        return out

    def next_cycle(self, campaign_id: int) -> int:
        with self.conn() as cx:
            row = cx.execute("SELECT MAX(cycle) AS c FROM batch WHERE campaign_id=?", (campaign_id,)).fetchone()
        return 0 if row["c"] is None else int(row["c"]) + 1

    # ----------------------------------------------------- measurements
    def add_measurements(self, campaign_id: int, batch_id: int | None, rows: list[dict], source: str) -> int:
        with self.conn() as cx:
            cx.executemany(
                "INSERT OR REPLACE INTO measurement"
                " (campaign_id, batch_id, compound_id, smiles, endpoint, value, predicted, uncertainty, source)"
                " VALUES (?,?,?,?,?,?,?,?,?)",
                [(campaign_id, batch_id, r["compound_id"], r["smiles"], r.get("endpoint", "pIC50"),
                  r["value"], r.get("predicted"), r.get("uncertainty"), source) for r in rows])
        return len(rows)

    def measurements(self, campaign_id: int) -> list[dict]:
        with self.conn() as cx:
            rows = cx.execute(
                "SELECT m.*, b.cycle FROM measurement m LEFT JOIN batch b ON m.batch_id=b.id"
                " WHERE m.campaign_id=? ORDER BY m.id", (campaign_id,)).fetchall()
        return [dict(r) for r in rows]

    # ------------------------------------------------------- model runs
    def add_model_run(self, campaign_id: int, cycle: int, report: dict, prospective: dict | None) -> None:
        with self.conn() as cx:
            cx.execute("INSERT INTO model_run (campaign_id, cycle, report, prospective) VALUES (?,?,?,?)",
                       (campaign_id, cycle, json.dumps(report), json.dumps(prospective) if prospective else None))

    def model_runs(self, campaign_id: int) -> list[dict]:
        with self.conn() as cx:
            rows = cx.execute("SELECT * FROM model_run WHERE campaign_id=? ORDER BY cycle", (campaign_id,)).fetchall()
        out = []
        for r in rows:
            d = dict(r)
            d["report"] = json.loads(d["report"])
            d["prospective"] = json.loads(d["prospective"]) if d["prospective"] else None
            out.append(d)
        return out
