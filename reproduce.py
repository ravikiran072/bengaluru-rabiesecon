"""Recalculate the Bengaluru scenario and export its results and figures."""
import argparse
import json
from pathlib import Path

from analysis import ROOT, calculate, default_inputs, is_report_scenario, write_results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, help="Alternative complete or partial input JSON")
    parser.add_argument("--out", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    inputs = default_inputs()
    if args.inputs:
        inputs.update(json.loads(args.inputs.read_text(encoding="utf-8")))
    result = calculate(inputs)
    write_results(result, args.out)
    scenario = "Report scenario" if is_report_scenario(inputs) else "Modified scenario"
    print(f"{scenario}; run {result['run_id']}")
    for label, value in result["summary"].items():
        if label in ("human_deaths", "pep_courses", "incremental_deaths_averted", "total_cost"):
            print(f"{label}: {value:,.6f}")
    print(f"Discounted net additional cost: {result['economics']['net_additional_cost_pv']:,.6f}")
    print(f"Results: {args.out.resolve()}")


if __name__ == "__main__":
    main()
