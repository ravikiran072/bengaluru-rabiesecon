# Rabies vaccination in Bangalore

Supplementary analysis software for the Bangalore dog vaccination scenario. The
application recalculates disease burden and economic outcomes from the supplied
GDREP V2.1.25 model and opens with the parameters used in the Bangalore burden
report updated 4 November 2026.

The five-year program delivers 567,353 dog vaccinations after a 50-year lead-in
with 10,000 vaccinations per year. Both arms use 90% PEP uptake. The comparator
continues 10,000 dog vaccinations annually.

| Five-year outcome | Report scenario |
|---|---:|
| Human rabies deaths | 645.4 |
| Deaths averted relative to comparator | 2,165.7 |
| PEP courses | 752,937 |
| Total program cost, undiscounted | USD 43,226,291 |
| Net additional cost, discounted at 3% | USD 1,422,572 |
| Cost per model DALY averted | USD 26.48 |

These are deterministic scenario projections. Human deaths denote fatal clinical
onsets in the source model. Costs use the supplied unit prices. Values in this
table are rounded; exported results retain full numerical precision.

## Installation and launch

Use Python 3.12 or newer. The supplied environment was tested with Python 3.14.3.
From this directory:

```bash
python -m venv .venv
```

Activate the environment on Windows:

```powershell
.venv\Scripts\Activate.ps1
```

Or on macOS and Linux:

```bash
source .venv/bin/activate
```

Install the dependencies and start the application:

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

On Windows, `start_dashboard.bat` creates the environment and installs the
dependencies on first use. Installation requires internet access; the installed
analysis runs locally without an external data service or Excel.

## Use

The interface has two tabs. **Rabies burden and vaccination** shows annual and
cumulative dog cases and human deaths, vaccination coverage, current immunity,
PEP use, dog rabies deaths and expected detected human cases. The cumulative
curves sum campaign years only. Annual tables, lead-in trajectories and methods
are available in supporting panels.

**Costs and economics** shows annual cost components, additional spending and
direct savings, cost-effectiveness results, unit costs and economic assumptions.
Annual costs and discounted economic results are labeled separately.

Use **Edit inputs** in the sidebar to open the scenario form. Submit edits together
with *Apply inputs and run*. *Restore report defaults* reloads the accepted inputs.
Modified scenarios are labeled separately from the report scenario.

The results download contains exact inputs, the evaluated engine configuration,
unrounded CSV tables and two four-panel figures in PDF and SVG.
`burden_comparison` contains annual and cumulative clinical dog cases and
human deaths. `campaign_outcomes` contains annual human and dog deaths, immunity
and PEP demand. Figures use campaign years 1–5 and months 0–60 for current dog
immunity. The full 55-year trajectory is available separately. The comparator
continues baseline vaccination; it is distinct from the zero-intervention null arm.

## Reproduce the report results

```bash
python reproduce.py
```

This command runs the report scenario and writes its outputs to `results/`.
Outputs are generated on demand and are not stored in the repository.

For an alternative scenario, use a separate JSON file. Unspecified fields inherit
the accepted Bangalore inputs:

```bash
python reproduce.py --inputs alternative_inputs.json --out alternative_results
```

Changing the number of campaign years requires schedules of equal length for
program vaccination, comparator vaccination, program PEP and comparator PEP.

## Repository contents

| Path | Contents |
|---|---|
| `app.py` | Streamlit interface |
| `bangalore_inputs.json` | Complete accepted Bangalore inputs |
| `analysis.py` | Calculation and export functions |
| `reproduce.py` | Command-line reproduction |
| `model/program.py` | Bangalore scenario adapter |
| `model/base_config.json` | Supplied engine configuration |
| `model/GDREP_V2_1_25years/` | Source engine modules and model notes |
| `METHODS.md` | Scenario specification and interpretation |
| `requirements.txt`, `start_dashboard.bat`, `.streamlit/` | Dependencies, Windows launcher and application configuration |

The source engine equations and base configuration are unchanged. The package
contains the model and inputs needed to run independently of the parent
RabiesEcon workspace and the original ZIP and Word files.

## Source and citation

The source package is `GDREP_Bangalore_50_Year_Full_Python (1).zip`. Scenario
parameters and reported outcomes come from
`Bangalore_GDREP_V2_1_25_Burden_Report (3).docx`, updated 4 November 2026.

The associated manuscript citation, repository DOI and reuse license should be
added by the authors when the supplementary repository is released.
