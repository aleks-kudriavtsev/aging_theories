"""Verified numeric-column FWF intake without executing or interpreting SAS.

Only declared INPUT ranges are supported. SAS FORMAT statements are output
formatting and never rescale the actual numeric literals in the source file.
The calling cohort must supply its locked expected_rows for real releases.
"""
from __future__ import annotations
from collections import Counter
import hashlib
from pathlib import Path
import re
import numpy as np
import pandas as pd

NUMBER = re.compile(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z')
NAME = re.compile(r'[A-Za-z_][A-Za-z0-9_]*\Z')


def _metadata(fields, nominal_width, expected_rows):
    if type(nominal_width) is not int or nominal_width <= 0:
        raise ValueError('Positive integer nominal width required')
    if expected_rows is not None and (type(expected_rows) is not int or expected_rows <= 0):
        raise ValueError('Locked expected row count must be a positive integer or None')
    if not isinstance(fields, dict) or not fields or 'SEQN' not in fields:
        raise ValueError('Nonempty field mapping including SEQN required')
    normalized = {}
    occupied = set()
    for name, span in fields.items():
        if not isinstance(name, str) or NAME.fullmatch(name) is None:
            raise ValueError('Invalid field name')
        if (not isinstance(span, (list, tuple)) or len(span) != 2
                or any(type(v) is not int for v in span)):
            raise ValueError('Two integer INPUT endpoints required')
        start, end = span
        if not 1 <= start <= end <= nominal_width:
            raise ValueError('INPUT range outside declared source width')
        columns = set(range(start, end + 1))
        if occupied & columns:
            raise ValueError('Overlapping selected INPUT fields')
        occupied.update(columns)
        normalized[name] = (start, end)
    return normalized


def _without_comments(sas):
    # Do not execute SAS. Drop its two ordinary comment forms before locating
    # statements; quoted strings are handled separately for INFILE keywords.
    sas = re.sub(r'/\*.*?\*/', '', sas, flags=re.S)
    return re.sub(r'(?ms)^[ \t]*\*[^;]*;', '', sas)


def _without_quoted_strings(value):
    return re.sub(r'''"(?:[^"]|"")*"|'(?:[^']|'')*' ''', ' ', value, flags=re.X)


def _input_contract(sas, fields):
    clean = _without_comments(sas.replace('\r\n', '\n').replace('\r', '\n'))
    blocks = re.findall(r'(?ims)^[ \t]*INPUT\b(.*?);', clean)
    if len(blocks) != 1:
        raise ValueError('Exactly one SAS INPUT block required')
    for name, expected in fields.items():
        matches = re.findall(r'(?im)^[ \t]*' + re.escape(name)
            + r'[ \t]+(\d+)(?:[ \t]*-[ \t]*(\d+))?[ \t]*$', blocks[0])
        positions = [(int(start), int(end or start)) for start, end in matches]
        if positions != [expected]:
            raise ValueError('Exact INPUT positions unresolved: ' + name)
    in_files = re.findall(r'(?ims)^[ \t]*INFILE\b([^;]*);', clean)
    if len(in_files) != 1:
        raise ValueError('Exactly one SAS INFILE statement required')
    return bool(re.search(r'\bMISSOVER\b', _without_quoted_strings(in_files[0]), re.I))


def _numeric(token, name):
    token = token.strip(' ')
    if token in ('', '.'):
        return np.nan
    if NUMBER.fullmatch(token) is None:
        raise ValueError('Non-numeric selected source field: ' + name)
    value = float(token)
    if not np.isfinite(value):
        raise ValueError('Nonfinite selected source field: ' + name)
    return value


def read_checked(data_path, sas_path, fields, nominal_width, expected_rows=None):
    """Return a dataframe and aggregate audit after two exact numeric reads.

    The manual pass streams source lines and retains only selected numeric
    fields. Unselected trailing characters are never padded or manufactured.
    Missing numeric codes (e.g. 8888) remain numbers for the cohort definition;
    this parser treats only literal dot/blank as missing. No source rows or IDs
    are printed. The source SHA covers exact original ASCII bytes/newlines.
    """
    fields = _metadata(fields, nominal_width, expected_rows)
    sas_bytes = Path(sas_path).read_bytes()
    missover = _input_contract(sas_bytes.decode('ascii'), fields)
    minimum = max(end for start, end in fields.values())
    widths = Counter(); rows = []; digest = hashlib.sha256(); source_bytes = 0
    with Path(data_path).open('r', encoding='ascii', newline='') as stream:
        for source_line in stream:
            encoded = source_line.encode('ascii')
            digest.update(encoded); source_bytes += len(encoded)
            if source_line.endswith('\r\n'):
                line = source_line[:-2]
            elif source_line.endswith(('\r', '\n')):
                line = source_line[:-1]
            else:
                line = source_line
            width = len(line); widths[width] += 1
            supported = minimum <= width <= nominal_width if missover else width == nominal_width
            if not supported:
                raise ValueError('Selected field truncation or unexpected source record width')
            rows.append({name: _numeric(line[start-1:end], name)
                         for name, (start, end) in fields.items()})
    if not rows:
        raise ValueError('Empty source file')
    if expected_rows is not None and len(rows) != expected_rows:
        raise ValueError('Source row count differs from locked expected count')
    manual = pd.DataFrame(rows, columns=list(fields))
    ids = manual.SEQN.to_numpy()
    if (not np.isfinite(ids).all() or np.any(ids <= 0)
            or np.any(ids != np.floor(ids)) or manual.SEQN.duplicated().any()):
        raise ValueError('Invalid or duplicate source SEQN')
    independent = pd.read_fwf(data_path, colspecs=[(start-1, end) for start, end in fields.values()],
        names=list(fields), header=None, dtype=float, keep_default_na=False, na_values=['', '.'])
    try:
        pd.testing.assert_frame_equal(manual, independent, check_exact=True, check_dtype=False)
    except AssertionError:
        raise ValueError('Independent source readers disagree; no source values exported') from None
    audit = dict(rows=len(manual), fields=len(fields), cell_comparisons=int(manual.size),
        exact_match=True, source_record_width_counts=dict(sorted(widths.items())),
        SAS_MISSOVER=missover, selected_fields_end=minimum, nominal_content_width=nominal_width,
        expected_rows=expected_rows,
        rows_with_absent_trailing_unselected_characters=sum(count for width, count in widths.items()
                                                          if width < nominal_width),
        source_bytes=source_bytes, source_sha256=digest.hexdigest(),
        sas_sha256=hashlib.sha256(sas_bytes).hexdigest(),
        FORMAT_scaling_applied=False)
    return manual, audit
