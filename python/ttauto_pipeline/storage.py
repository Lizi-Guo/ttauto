# SPDX-License-Identifier: GPL-3.0-or-later
"""Atomic, content-bound checkpoints and transactional candidate records."""
from hashlib import sha256
import json
import os
from pathlib import Path
import sqlite3


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def digest(value):
    return sha256(canonical(value).encode("utf-8")).hexdigest()


def file_digest(path):
    result = sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            result.update(block)
    return result.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name+".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def source_digest():
    root = Path(__file__).resolve().parent.parent
    return digest({p.relative_to(root).as_posix(): file_digest(p)
                   for package in ("ttauto_pipeline", "ttauto_cert")
                   for p in sorted((root/package).rglob("*.py"))})


class Journal:
    def __init__(self, path, identity, resume=False):
        path = Path(path)
        if path.exists() != resume:
            raise ValueError("existing journal requires --resume; --resume requires an existing journal")
        self.connection = sqlite3.connect(path)
        self.connection.execute("PRAGMA synchronous=FULL")
        if not resume:
            self.connection.executescript("""
                CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL, digest TEXT NOT NULL);
                CREATE TABLE candidates(id INTEGER PRIMARY KEY, payload TEXT NOT NULL, digest TEXT NOT NULL);
            """)
            self.put("identity", identity)
            self.connection.commit()
        else:
            try:
                if self.get("identity") != identity:
                    raise ValueError("resume refused: configuration, source code, automaton or backend changed")
            except BaseException:
                self.connection.close()
                raise

    def get(self, key):
        row = self.connection.execute("SELECT value,digest FROM meta WHERE key=?", (key,)).fetchone()
        if row is None:
            return None
        value = json.loads(row[0])
        if digest(value) != row[1]:
            raise ValueError(f"corrupt journal metadata: {key}")
        return value

    def put(self, key, value):
        self.connection.execute("INSERT OR REPLACE INTO meta VALUES(?,?,?)", (key, canonical(value), digest(value)))

    def add_candidate(self, number, record):
        self.connection.execute("INSERT INTO candidates VALUES(?,?,?)", (number, canonical(record), digest(record)))

    def checkpoint(self, value):
        self.put("checkpoint", value)
        self.connection.commit()

    def records(self):
        for index, payload, checksum in self.connection.execute("SELECT id,payload,digest FROM candidates ORDER BY id"):
            record = json.loads(payload)
            if digest(record) != checksum:
                raise ValueError(f"corrupt candidate record: {index}")
            yield index, record

    def close(self):
        self.connection.close()
