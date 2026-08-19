import contextlib
import csv
import os
import sys


def ensure(path):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    return path


def read_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def write_rows(path, header, rows):
    ensure(path)
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def write_dicts(path, rows, fieldnames=None):
    ensure(path)
    fields = list(fieldnames or rows[0])
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_text(path, text):
    ensure(path)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def save(fig, path, **kw):
    ensure(path)
    fig.savefig(path, **kw)
    return path


class _Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
        return len(data)

    def flush(self):
        for s in self.streams:
            s.flush()

    def isatty(self):
        return False


@contextlib.contextmanager
def captured(path):
    ensure(path)
    with open(path, "w") as fh:
        keep = sys.stdout
        sys.stdout = _Tee(keep, fh)
        try:
            yield path
        finally:
            sys.stdout = keep


def note(*args):
    print(*args, file=sys.stderr)
