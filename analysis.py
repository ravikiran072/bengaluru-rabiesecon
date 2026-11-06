"""Bengaluru analysis and portable results exports."""
from __future__ import annotations

import io
import json
from pathlib import Path
import zipfile

import pandas as pd

from model.program import normalize, run_analysis

ROOT = Path(__file__).resolve().parent


def default_inputs():
    return json.loads((ROOT / "bengaluru_inputs.json").read_text(encoding="utf-8"))


def calculate(inputs):
    """Simulate disease burden, vaccination and economic outcomes."""
    return run_analysis(inputs)


def is_report_scenario(inputs):
    return normalize(inputs) == normalize(default_inputs())


def summary_figure(result, formats=("svg", "pdf")):
    """Four-panel figure with zero-based outcome axes and editable vector text."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.ticker import MaxNLocator, StrMethodFormatter

    a = pd.DataFrame(result["annual"])
    m = pd.DataFrame(result["monthly_immunity"])
    styles = {"font.size": 9, "font.family": "DejaVu Sans", "svg.fonttype": "none",
              "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False}
    with matplotlib.rc_context(styles):
        fig = Figure(figsize=(9, 6.8), layout="constrained")
        FigureCanvasAgg(fig)
        axes = fig.subplots(2, 2)
        panels = [
            (axes[0, 0], "a  Human rabies deaths", a, "year", "human_deaths", "Deaths per year"),
            (axes[0, 1], "b  Dog rabies deaths", a, "year", "dog_rabies_deaths", "Deaths per year"),
            (axes[1, 0], "c  Dog population immunity", m, "month", "immune_percent", "Immune dogs (%)"),
            (axes[1, 1], "d  PEP demand", a, "year", "pep_courses", "Courses per year"),
        ]
        for ax, title, data, x, y, ylabel in panels:
            marker = None if x == "month" else "o"
            ax.plot(data[x], data[y], color="#21618c", label="Program", marker=marker,
                    markersize=3, linewidth=1.7)
            ax.plot(data[x], data["comparator_" + y], color="#666666", label="Comparator",
                    linestyle="--", marker=marker, markersize=3, linewidth=1.4)
            ax.set_title(title, loc="left", fontsize=10)
            ax.set_xlabel("Month after campaign start" if x == "month" else "Campaign year")
            ax.set_ylabel(ylabel)
            ax.set_ylim(bottom=0)
            if x == "month":
                ax.set_xlim(0, int(m.month.max()))
            else:
                ax.set_xticks(a.year)
            ax.yaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
            ax.grid(axis="y", color="#dddddd", linewidth=0.5)
            if y == "pep_courses":
                ax.plot(data[x], data.pep_true_exposures, color="#a05a2c", label="True-exposure PEP",
                        marker="s", markersize=3, linewidth=1.3)
            ax.legend(frameon=False, fontsize=8)
        files = {}
        for fmt in formats:
            buffer = io.BytesIO()
            fig.savefig(buffer, format=fmt, dpi=160)
            files[f"campaign_outcomes.{fmt}"] = buffer.getvalue()
        return files


def burden_comparison_data(result):
    """Annual outcomes and unrounded cumulative totals over campaign years only."""
    fields = ["clinical_dog_cases", "comparator_clinical_dog_cases",
              "human_deaths", "comparator_human_deaths"]
    data = pd.DataFrame(result["annual"])[["year", *fields]].copy()
    for field in fields:
        data["cumulative_" + field] = data[field].cumsum()
    return data


def burden_comparison_figure(result, formats=("svg", "pdf")):
    """Annual and cumulative disease burden, using the selected comparator."""
    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.ticker import MaxNLocator, StrMethodFormatter

    data = burden_comparison_data(result)
    p = result["inputs"]
    unchanged_baseline = (
        p["history_mode"] == p["comparator_mode"]
        and p["history_weeks"] == p["comparator_weeks"]
        and all(value == p["history_value"] for value in p["comparator_values"])
    )
    comparator_label = "Continued baseline vaccination" if unchanged_baseline else "Comparator vaccination schedule"
    styles = {"font.size": 9, "font.family": "DejaVu Sans", "svg.fonttype": "none",
              "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False}
    with matplotlib.rc_context(styles):
        fig = Figure(figsize=(9, 6.8), layout="constrained")
        FigureCanvasAgg(fig)
        axes = fig.subplots(2, 2)
        panels = [
            (axes[0, 0], "a  Clinical dog cases (annual)", "clinical_dog_cases", "Cases per year", False),
            (axes[0, 1], "b  Clinical dog cases (cumulative)", "clinical_dog_cases", "Cumulative cases", True),
            (axes[1, 0], "c  Human rabies deaths (annual)", "human_deaths", "Deaths per year", False),
            (axes[1, 1], "d  Human rabies deaths (cumulative)", "human_deaths", "Cumulative deaths", True),
        ]
        for ax, title, field, ylabel, cumulative in panels:
            prefix = "cumulative_" if cumulative else ""
            ax.plot(data.year, data[prefix + field], color="#21618c", label="Vaccination program",
                    marker="o", markersize=3, linewidth=1.7)
            ax.plot(data.year, data[prefix + "comparator_" + field], color="#666666",
                    label=comparator_label, linestyle="--", marker="o", markersize=3, linewidth=1.4)
            ax.set_title(title, loc="left", fontsize=10)
            ax.set_xlabel("Campaign year")
            ax.set_ylabel(ylabel)
            ax.set_xticks(data.year)
            ax.set_ylim(bottom=0)
            ax.yaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.0f}"))
            ax.grid(axis="y", color="#dddddd", linewidth=0.5)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False, fontsize=9)
        files = {}
        for fmt in formats:
            buffer = io.BytesIO()
            fig.savefig(buffer, format=fmt, dpi=160)
            files[f"burden_comparison.{fmt}"] = buffer.getvalue()
        return files


def export_files(result):
    """Unrounded numerical outputs, exact inputs and vector figures."""
    def encoded(obj):
        return (json.dumps(obj, indent=2, allow_nan=False) + "\n").encode("utf-8")

    files = {
        "inputs.json": encoded(result["inputs"]),
        "engine_config.json": encoded(result["engine_config"]),
        "result.json": encoded(result),
        "summary.json": encoded({"summary": result["summary"], "economics": result["economics"]}),
    }
    for name in ("annual", "full_timeline", "monthly_immunity", "weekly_immunity", "benefit_cost"):
        files[name + ".csv"] = pd.DataFrame(result[name]).to_csv(index=False).encode("utf-8")
    files["burden_comparison.csv"] = burden_comparison_data(result).to_csv(index=False).encode("utf-8")
    files.update(summary_figure(result))
    files.update(burden_comparison_figure(result))
    return files


def results_zip(result):
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
        for name, content in export_files(result).items():
            archive.writestr(name, content)
    return buffer.getvalue()


def write_results(result, directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for name, content in export_files(result).items():
        (directory / name).write_bytes(content)
