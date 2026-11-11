"""Durable execution records; these are checkpoints, not mathematical proofs.

SQLite transactions prevent a partial write from marking a state complete.
Every reusable item is bound to its input and checked for accidental corruption.
The search still verifies mathematical results when reading cached backend data.
"""
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import tempfile


def canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=True, allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical_json(value).encode('ascii')).hexdigest()


def atomic_json(path, value):
    """Replace a report only after its complete contents have been flushed."""
    path = Path(path)
    encoded = json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n'
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='\n',
                                         dir=path.parent, prefix=path.name+'.',
                                         suffix='.tmp', delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


class SearchJournal:
    """Single-writer resume log, bound to an exact run identity.

    Identity must contain surface, bound, source hashes and executable hashes.
    Resource allowances may change on resume; mathematical inputs may not.
    A changed implementation requires a new run directory, not a stale resume.
    """
    def __init__(self, path, identity, *, resume=False):
        path = Path(path)
        exists = path.exists()
        if exists and not resume:
            raise ValueError('search journal already exists; use --resume or a new output directory')
        if resume and not exists:
            raise ValueError('cannot resume: no search journal in this output directory')
        self.connection = sqlite3.connect(path)
        self.connection.execute('PRAGMA synchronous=FULL')
        self.connection.execute('CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)')
        self.connection.execute('CREATE TABLE IF NOT EXISTS items '
            '(name TEXT PRIMARY KEY, input_hash TEXT NOT NULL, output TEXT NOT NULL, output_hash TEXT NOT NULL)')
        expected = canonical_json(identity)
        row = self.connection.execute("SELECT value FROM metadata WHERE key='identity'").fetchone()
        if row is None:
            if exists:
                self.close()
                raise ValueError('search journal is missing its committed identity')
            with self.connection:
                self.connection.execute("INSERT INTO metadata VALUES ('identity', ?)", (expected,))
        elif row[0] != expected:
            self.close()
            raise ValueError('resume identity mismatch: surface, bound, sources or executables changed')

    def get(self, name, inputs):
        row = self.connection.execute(
            'SELECT input_hash, output, output_hash FROM items WHERE name=?', (name,)).fetchone()
        if row is None:
            return None
        if row[0] != digest(inputs):
            raise ValueError('checkpoint input mismatch: '+name)
        value = json.loads(row[1])
        if digest(value) != row[2]:
            raise ValueError('checkpoint checksum mismatch: '+name)
        return value

    def put(self, name, inputs, output, *, replace=False):
        previous = self.get(name, inputs)
        if previous is not None:
            if canonical_json(previous) == canonical_json(output):
                return
            if not replace:
                raise ValueError('deterministic checkpoint changed: '+name)
        with self.connection:
            self.connection.execute('INSERT OR REPLACE INTO items VALUES (?, ?, ?, ?)',
                (name, digest(inputs), canonical_json(output), digest(output)))

    def names(self, prefix):
        # Length then lexical order sorts nonnegative decimal ordinals without
        # relying on SQLite's fixed-width integer conversion.
        rows = self.connection.execute(
            'SELECT name FROM items WHERE substr(name,1,?)=? ORDER BY length(name),name',
            (len(prefix), prefix))
        for row in rows:
            yield row[0]

    def close(self):
        self.connection.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()
