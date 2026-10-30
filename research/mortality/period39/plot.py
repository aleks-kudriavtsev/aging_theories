"""Optional static chart of a frozen aggregate reference; no risk calculations."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


def series(payload):
    if (payload.get("schema") != "period39-results/1" or payload.get("reference_rows") != 164
            or payload.get("clinical_validation") is not False or payload.get("current_calendar_forecast") is not False
            or type(payload.get("source_year")) is not int or payload["source_year"] != 2023
            or payload.get("period_assumption") != "hold_2023_attained_age_schedule"
            or payload.get("scope") != "provisional_local_derived_period_reference_under_declared_units"):
        raise ValueError("Only the named provisional aggregate period reference can be plotted")
    points = {}
    for row in payload["rows"]:
        r = row["request"]
        fixed = {"schema": "period39-request/1", "model": "RUS_WPP2023_period_reference",
                 "country": "RUS", "source_year": 2023, "age_origin": "exact_birthday",
                 "period_assumption": "hold_2023_attained_age_schedule", "purpose": "population_period_reference_only"}
        if (type(r) is not dict or set(r) != set(fixed) | {"source_sex", "age_exact_years", "horizon_years"}
                or any(r[k] != v for k,v in fixed.items())
                or type(r["source_year"]) is not int or type(r["age_exact_years"]) is not int
                or type(r["horizon_years"]) is not int or type(r["source_sex"]) is not str):
            raise ValueError("Plot labels require the exact Russian2023 birthday-period input contract")
        key = (r["source_sex"], r["age_exact_years"], r["horizon_years"])
        value = row["methods"]["constant_force"]["probability"]
        if key in points or type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("Unique finite aggregate probabilities required")
        points[key] = value
    expected = {(s, a, h) for s in ("F", "M") for a in range(40, 81) for h in (5, 10)}
    if set(points) != expected:
        raise ValueError("Exact 164-cell age/sex/horizon grid required")
    return points


def render(results, svg, png_preview=None):
    from .analyze import verify_freeze
    verify_freeze()
    results, svg = Path(results), Path(svg)
    preview = Path(png_preview) if png_preview else None
    if svg.exists() or (preview is not None and preview.exists()):
        raise ValueError("New image paths required; no overwrite")
    raw = results.read_bytes()
    source_sha = hashlib.sha256(raw).hexdigest()
    points = series(json.loads(raw))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    with plt.rc_context({"svg.hashsalt": "period39", "font.family": "DejaVu Sans"}):
        fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.5), sharey=True)
        ages = list(range(40, 81))
        for axis, horizon in zip(axes, (5, 10)):
            for sex, label, color in (("F", "F — женская серия", "#8d3b72"), ("M", "M — мужская серия", "#176b87")):
                axis.plot(ages, [points[sex, age, horizon] for age in ages], color=color,
                          label=label, linewidth=2, marker=".", markersize=3)
            axis.set_title(f"Горизонт {horizon} лет")
            axis.set_xlabel("Точный возраст на исходный день рождения, лет")
            axis.set_xlim(40, 80)
            axis.set_ylim(0, 1)
            axis.yaxis.set_major_formatter(PercentFormatter(xmax=1))
            axis.grid(alpha=.2)
            axis.legend(loc="upper left", frameon=False, fontsize=9)
        axes[0].set_ylabel("Условная периодная вероятность смерти")
        fig.suptitle("Россия: локальный периодный ориентир по ставкам 2023 года", fontsize=13)
        fig.text(.5, .08, "Constant-force approximation; единицы и источник — по сохранённому derived CSV.", ha="center", fontsize=9)
        fig.text(.5, .045, "Не индивидуальный прогноз и не disease-free baseline. Линии соединяют целочисленные возрасты.", ha="center", fontsize=9)
        fig.tight_layout(rect=(0, .12, 1, .94))
        svg.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(svg, format="svg", metadata={"Date": None, "Creator": "period39 local research",
            "Title": "Provisional Russian population period reference",
            "Description": "source_results_sha256="+source_sha+"; not individual prognosis or clinical validation"})
        if preview:
            preview.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(preview, dpi=150, metadata={"Description": "period39 provisional population reference; source_sha256="+source_sha})
        plt.close(fig)
    return {"schema": "period39-chart-provenance/1", "source_results_sha256": source_sha,
            "svg_sha256": hashlib.sha256(svg.read_bytes()).hexdigest(),
            "png_sha256": hashlib.sha256(preview.read_bytes()).hexdigest() if preview else None,
            "matplotlib_version": matplotlib.__version__, "clinical_validation": False}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--svg", type=Path, required=True)
    p.add_argument("--png-preview", type=Path)
    x = p.parse_args()
    print(json.dumps(render(x.results, x.svg, x.png_preview)))
