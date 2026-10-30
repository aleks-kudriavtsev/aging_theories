"""Fetch age/sex-specific mortality rates for USA, Germany and Russia from UN WPP 2024
and digitize them into data/mortality/ as plain CSV (mx and derived qx by single age).

Source: United Nations, DESA, Population Division, World Population Prospects 2024
(https://population.un.org/wpp), redistributed as the `wpp2024` R package
(https://github.com/PPgp/wpp2024, CC BY 3.0 IGO). The official mortality.org (HMD) site
requires a login and is unreachable from this environment's network whitelist (only
GitHub/PyPI are reachable, see docs/07); UN WPP is the closest openly-licensed substitute
that already covers all three countries with a consistent single-year-of-age (1x1) series
1950-2023 (historical estimates) plus UN median projections 2024-2100.

`mx1dt.rda` is the full global dataset (~85 MB, all countries/years/ages) -- this script
downloads it once, filters to the three countries requested, and discards the rest; the
raw file is not committed (see data/mortality/README.md).

Requires `pyreadr` (only for this one-off fetch script, not a runtime dependency of the
model/web app -- see analysis/requirements.txt):

    pip install pyreadr
    python tools/fetch_wpp_mortality.py
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import tempfile
import urllib.request

import pyreadr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(ROOT, "data", "mortality")
OUT_CSV = os.path.join(OUT_DIR, "mx_qx_usa_deu_rus.csv")

SOURCE_URL = "https://raw.githubusercontent.com/PPgp/wpp2024/master/data/mx1dt.rda"

# UN M49 numeric country codes (standard, used as-is by wpp2024's country_code column).
COUNTRIES = {
    840: "USA",
    276: "DEU",  # Germany -- WPP treats it as one continuous historical series
    643: "RUS",  # Russian Federation -- WPP reconstructs a pre-1991 series too
}

LAST_HISTORICAL_YEAR = 2023  # 2024+ in this package are UN median projections, not observed


def mx_to_qx(mx: float) -> float:
    """Probability of death within one year from the central death rate, under the
    standard actuarial assumption of a constant force of mortality within the single-year
    age interval: qx = 1 - exp(-mx). Exact under that assumption; a simplification at the
    extremes (age 0, and the oldest ages where mx is large) where the true qx is usually
    computed with an age-specific a_x (average fraction of the interval lived by those who
    die) instead of the implicit a_x=0.5 this formula assumes -- flagged here, not silently
    hidden, per this repo's convention of marking approximations explicitly."""
    return 1.0 - math.exp(-mx)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--with-projections", action="store_true",
                    help="also include UN median projections (2024-2100), not just "
                         "historical estimates (<=2023). Off by default so the committed "
                         "file only contains observed/estimated mortality, not projections.")
    ap.add_argument("--both-sexes", action="store_true",
                    help="also include the combined-sex (B) series. Off by default -- the "
                         "request this was built for wants sex-specific (M/F) rates.")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    with tempfile.NamedTemporaryFile(suffix=".rda", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        print(f"Downloading {SOURCE_URL} ...")
        urllib.request.urlretrieve(SOURCE_URL, tmp_path)
        print("Parsing .rda ...")
        result = pyreadr.read_r(tmp_path)
        df = result["mx1dt"]
        df = df[df["country_code"].isin(COUNTRIES)]
        if not args.with_projections:
            df = df[df["year"] <= LAST_HISTORICAL_YEAR]

        sexes = [("M", "mxM"), ("F", "mxF")] + ([("B", "mxB")] if args.both_sexes else [])
        rows = []
        for _, r in df.iterrows():
            code = int(r["country_code"])
            year = int(r["year"])
            age = int(r["age"])
            for sex, col in sexes:
                mx = float(r[col])
                rows.append({
                    "country": COUNTRIES[code],
                    "sex": sex,
                    "year": year,
                    "age": age,
                    "mx": round(mx, 6),
                    "qx": round(mx_to_qx(mx), 6),
                    "series": "estimate" if year <= LAST_HISTORICAL_YEAR else "projection",
                })
        rows.sort(key=lambda r: (r["country"], r["sex"], r["year"], r["age"]))

        with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=["country", "sex", "year", "age", "mx", "qx", "series"])
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {OUT_CSV} -- {len(rows)} rows "
              f"({len(COUNTRIES)} countries x {len(sexes)} sexes x "
              f"{df['year'].nunique()} years x {df['age'].nunique()} ages)")
    finally:
        os.unlink(tmp_path)


if __name__ == "__main__":
    main()
