"""Streamlit interface for the Bengaluru rabies analysis."""
import copy
import json

import altair as alt
import pandas as pd
import streamlit as st

from analysis import ROOT, burden_comparison_data, calculate, default_inputs, is_report_scenario, results_zip

st.set_page_config(page_title="RabiesEcon | Bengaluru", layout="wide")


@st.cache_data(show_spinner=False, max_entries=8)
def cached_analysis(serialized):
    """Calculate the selected scenario from the packaged model and inputs."""
    return calculate(json.loads(serialized))


@st.cache_data(show_spinner=False, max_entries=8)
def cached_download(serialized):
    """Generate numerical results and figures for the selected scenario."""
    return results_zip(cached_analysis(serialized))


def restore_defaults():
    st.session_state.inputs = default_inputs()
    st.session_state.revision += 1
    st.session_state.results_tab = "Rabies burden and vaccination"
    st.session_state.editing_inputs = False


def close_inputs():
    st.session_state.editing_inputs = False


def number(value, decimals=0):
    return "Not defined" if value is None else f"{value:,.{decimals}f}"


def money(value, decimals=0):
    return "Not defined" if value is None else "$" + number(value, decimals)


def show_table(data, columns, decimals=None, total=False, percent=()):
    """Format a copy for display; downloads retain unrounded numerical values."""
    frame = data[list(columns)].copy()
    if total:
        row = {key: ("Total" if key == "year" else frame[key].sum()) for key in frame}
        frame["year"] = frame.year.astype(str)
        frame = pd.concat([frame, pd.DataFrame([row])], ignore_index=True)
    for key in frame:
        if key == "year":
            continue
        digits = (decimals or {}).get(key, 0)
        if key in percent:
            frame[key] = frame[key].map(lambda value: f"{value:.{digits}%}")
        else:
            frame[key] = frame[key].map(lambda value: number(value, digits))
    st.dataframe(frame.rename(columns=columns), hide_index=True, width="stretch",
                 alt="Annual model results with labeled outcome columns")


def comparison_chart(data, x, y, xlabel, ylabel, title, *, comparator=None, percent=False, key=None):
    comparator = comparator or "comparator_" + y
    frame = data[[x, y, comparator]].rename(
        columns={y: "Program", comparator: "Comparator"})
    frame = frame.melt(id_vars=x, var_name="Arm", value_name="Value")
    chart = alt.Chart(frame).mark_line(point=x == "year", strokeWidth=2).encode(
        x=alt.X(f"{x}:Q", title=xlabel, axis=alt.Axis(tickMinStep=1),
                scale=alt.Scale(zero=x == "month")),
        y=alt.Y("Value:Q", title=ylabel, scale=alt.Scale(zero=True),
                axis=alt.Axis(format=".0%" if percent else ",.0f")),
        color=alt.Color("Arm:N", scale=alt.Scale(domain=["Program", "Comparator"],
                        range=["#21618c", "#777777"]), legend=alt.Legend(title=None, orient="bottom")),
        strokeDash=alt.StrokeDash("Arm:N", scale=alt.Scale(domain=["Program", "Comparator"],
                                  range=[[1, 0], [6, 4]]), legend=None),
        tooltip=[alt.Tooltip(f"{x}:Q", title=xlabel, format=".0f"), "Arm:N",
                 alt.Tooltip("Value:Q", title=ylabel, format=".1%" if percent else ",.2f")],
    ).properties(height=250)
    st.markdown(f"**{title}**")
    st.altair_chart(chart, width="stretch", key=key, alt=f"{title}: program and comparator over time")


@st.dialog("Scenario inputs", width="large", on_dismiss=close_inputs)
def edit_inputs():
    st.markdown("The report inputs load automatically. Edit a scenario and select **Apply inputs and run**. "
                "Restore report defaults returns to the accepted analysis.")
    draft = copy.deepcopy(st.session_state.inputs)
    rev = st.session_state.revision

    def fields(specs):
        columns = st.columns(2)
        for idx, (key, label) in enumerate(specs):
            with columns[idx % 2]:
                val = draft[key]
                fmt = "%.7f" if key == "z" else "%.9f" if key in ("test_probability", "waning_per_week") else "%g"
                draft[key] = st.number_input(label, min_value=0 if type(val) is int else 0.0,
                    value=val, step=1 if type(val) is int else None,
                    format="%d" if type(val) is int else fmt, key=f"{rev}:{key}")

    with st.form("scenario_inputs"):
        with st.expander("Population and history", expanded=True):
            fields([("human_population", "Humans at campaign start"), ("dog_population", "Dogs at campaign start"),
                    ("area_km2", "Transmission area (km²)"), ("lead_in_years", "Lead-in years"),
                    ("history_value", "Historical vaccinations per year"), ("history_weeks", "Historical delivery weeks per year")])
        with st.expander("Annual vaccination schedules", expanded=True):
            program_text = st.text_input("Program vaccinations, one count per year", ", ".join(f"{x:g}" for x in draft["program_values"]), key=f"{rev}:program_values")
            comparator_text = st.text_input("Comparator vaccinations, one count per year", ", ".join(f"{x:g}" for x in draft["comparator_values"]), key=f"{rev}:comparator_values")
            fields([("program_weeks", "Program delivery weeks per year"), ("comparator_weeks", "Comparator delivery weeks per year"), ("start_week", "First vaccination week")])
        with st.expander("Transmission and vaccine protection"):
            fields([("r0", "Requested R₀"), ("r0_cap", "R₀ cap"), ("z", "Z (km² per infectious dog per week)"),
                ("incubation_days", "Dog incubation (days)"), ("dog_recruitment", "Dog recruitment per dog-year"),
                ("dog_mortality", "Dog natural mortality per dog-year"), ("human_birth_rate", "Human births per 1,000 per year"),
                ("life_expectancy", "Human life expectancy (years)"), ("dog_vaccine_efficacy", "Dog vaccine efficacy (fraction)"),
                ("exposed_protection", "Protection in vaccinated exposed dogs (fraction)"), ("protected_weeks", "Protected period (weeks)"),
                ("revaccination_weeks", "Revaccination interval (weeks)"), ("waning_per_week", "Immunity waning per week"),
                ("dog_clinical_probability", "Dog clinical probability"), ("dog_exit_rate", "Dog infectious exit rate per week")])
        with st.expander("PEP and case detection"):
            program_pep_text = st.text_input("Program PEP uptake, one fraction per year", ", ".join(f"{x:g}" for x in draft["program_pep"]), key=f"{rev}:program_pep")
            comparator_pep_text = st.text_input("Comparator PEP uptake, one fraction per year", ", ".join(f"{x:g}" for x in draft["comparator_pep"]), key=f"{rev}:comparator_pep")
            fields([("history_pep", "Historical PEP uptake (fraction)"), ("pep_efficacy", "PEP efficacy (fraction)"),
                    ("pep_multiplier", "PEP demand multiplier"), ("cdr_before", "Detection before campaign (fraction)"),
                    ("cdr_after", "Detection during campaign (fraction)")])
        with st.expander("Unit costs and economic assumptions"):
            st.caption("All monetary inputs are USD. Investigations are a fixed annual workload in both arms.")
            fields([("dog_cost", "Cost per dog vaccination"), ("pep_material", "PEP materials per injection"),
                    ("pep_overhead", "PEP overhead per visit"), ("pep_vaccine", "PEP vaccine per dose"),
                    ("pep_doses_per_visit", "PEP doses per visit"), ("pep_visits", "PEP visits per course"),
                    ("rig_cost", "RIG cost per recipient"), ("rig_fraction", "Fraction receiving RIG"),
                    ("quarantine_probability", "Quarantine probability per suspected exposure"), ("quarantine_cost", "Cost per quarantine"),
                    ("test_probability", "Testing probability per suspected exposure"), ("test_cost", "Cost per laboratory test"),
                    ("investigation_cost", "Cost per investigation"), ("cost_discount", "Annual cost discount rate"),
                    ("outcome_discount", "Annual health discount rate"), ("average_age_at_death", "Average age at death (years)")])
            current_unit = st.session_state.inputs["investigation_cost"]
            workload = st.session_state.inputs["other_program_cost"] / current_unit if current_unit else 200.0
            investigations = st.number_input("Investigations per year in each arm", min_value=0.0,
                value=workload, key=f"{rev}:investigations")
        submitted = st.form_submit_button("Apply inputs and run", type="primary")
    if submitted:
        try:
            for key, text in [("program_values", program_text), ("comparator_values", comparator_text),
                              ("program_pep", program_pep_text), ("comparator_pep", comparator_pep_text)]:
                draft[key] = [float(value.strip()) for value in text.split(",")]
            draft["other_program_cost"] = draft["other_comparator_cost"] = investigations * draft["investigation_cost"]
            draft["investigation_probability"] = 0.0
            cached_analysis(json.dumps(draft, sort_keys=True))
            st.session_state.inputs = draft
            st.session_state.revision += 1
            st.session_state.editing_inputs = False
            st.rerun()
        except (ValueError, TypeError, KeyError, ArithmeticError, AssertionError) as exc:
            st.error(f"Inputs were not applied: {exc}")
    st.caption("The inputs download contains the complete parameter set used for the displayed result.")


if "inputs" not in st.session_state:
    st.session_state.inputs = default_inputs()
    st.session_state.revision = 0

with st.sidebar:
    st.title("RabiesEcon")
    st.caption("Bengaluru analysis")
    if st.button("Edit inputs", width="stretch", icon=":material/tune:"):
        st.session_state.editing_inputs = True
    st.button("Restore report defaults", on_click=restore_defaults, width="stretch")

st.title("Bengaluru rabies burden")
st.caption("Dog vaccination, human health and economic outcomes | GDREP V2.1.25")

serialized = json.dumps(st.session_state.inputs, sort_keys=True)
try:
    with st.spinner("Calculating the vaccination history and campaign outcomes…"):
        result = cached_analysis(serialized)
except (ValueError, TypeError, KeyError, ArithmeticError, AssertionError) as exc:
    st.error(f"The analysis could not be completed: {exc}")
    st.stop()

p, s, e = result["inputs"], result["summary"], result["economics"]
annual = pd.DataFrame(result["annual"])
monthly = pd.DataFrame(result["monthly_immunity"])
benefit = pd.DataFrame(result["benefit_cost"])
burden = burden_comparison_data(result)
years = s["years"]
annual["comparator_pep_true_exposures"] = annual.comparator_human_exposures * pd.Series(p["comparator_pep"])
annual["dog_deaths_averted"] = annual.comparator_dog_rabies_deaths - annual.dog_rabies_deaths

st.caption(f"{number(p['human_population'])} people · {number(p['dog_population'])} dogs at campaign start · "
           f"{p['area_km2']:g} km² · {p['lead_in_years']}-year lead-in · {years}-year campaign")
if is_report_scenario(p):
    st.caption("Report scenario. Calculated using the Bengaluru report inputs.")
else:
    st.info("Modified scenario. Results use the edited inputs and are not the report estimates.")

comparator_values = p["comparator_values"]
if len(set(comparator_values)) == 1:
    value = comparator_values[0]
    comparator_schedule = (f"{number(value)} dog vaccinations per year" if p["comparator_mode"] == "count"
                           else f"{value:.1%} annual dose coverage")
else:
    comparator_schedule = "the year-specific vaccination schedule in the inputs"
comparator_pep = (f"{p['comparator_pep'][0]:.0%} PEP uptake" if len(set(p["comparator_pep"])) == 1
                  else "year-specific PEP uptake")
program_pep = (f"{p['program_pep'][0]:.0%}" if len(set(p["program_pep"])) == 1 else "year-specific")

with st.sidebar:
    st.markdown("**Scenario**")
    st.caption(f"R₀ used: {s['actual_r0']:g}\n\nZ: {p['z']:g}\n\nProgram PEP uptake: {program_pep}")
    st.markdown("**Comparator**")
    st.caption(f"{comparator_schedule.capitalize()}; {comparator_pep}.")
    st.download_button("Download results and figures", cached_download(serialized),
        f"bengaluru_{result['run_id']}.zip", "application/zip", width="stretch",
        help="Exact inputs, unrounded CSV tables, and figures in PDF and SVG.")
    st.download_button("Download inputs", json.dumps(p, indent=2),
        "bengaluru_inputs.json", "application/json", width="stretch")
    st.caption(f"Run {result['run_id']}")

burden_tab, economics_tab = st.tabs(
    ["Rabies burden and vaccination", "Costs and economics"], key="results_tab", on_change="rerun")

if burden_tab.open:
    with burden_tab:
        st.markdown(f"**Campaign years 1–{years}.** The comparator uses {comparator_schedule} and {comparator_pep}. "
                    "Both arms share the same lead-in.")
        top = st.columns(3)
        top[0].metric("Human rabies deaths", number(s["human_deaths"], 1), border=True)
        top[1].metric("Human deaths averted", number(s["incremental_deaths_averted"], 1), border=True)
        top[2].metric("Clinical dog cases", number(s["clinical_dog_cases"]), border=True)
        second = st.columns(3)
        second[0].metric("Dog vaccinations", number(s["dog_doses"]), border=True)
        second[1].metric("PEP courses", number(s["pep_courses"]), border=True)
        second[2].metric("Final dog immunity", f"{s['year_end_immunity_percent']:.1f}%", border=True)
        if s["comparator_human_deaths"] > 0:
            reduction = s["incremental_deaths_averted"] / s["comparator_human_deaths"]
            effect = f"{abs(reduction):.1%} {'fewer' if reduction >= 0 else 'more'} human deaths"
            st.caption(f"{number(s['comparator_human_deaths'], 1)} comparator deaths; {effect} under the program. "
                       "Event counts are campaign totals. Immunity is the proportion at the end of the final year.")

        st.subheader("Annual and cumulative disease burden")
        left, right = st.columns(2)
        with left:
            comparison_chart(burden, "year", "clinical_dog_cases", "Campaign year", "Cases per year",
                             "Clinical dog cases — annual", key="annual_dog_cases")
        with right:
            comparison_chart(burden, "year", "cumulative_clinical_dog_cases", "Campaign year", "Cumulative cases",
                             "Clinical dog cases — cumulative", comparator="cumulative_comparator_clinical_dog_cases",
                             key="cumulative_dog_cases")
        left, right = st.columns(2)
        with left:
            comparison_chart(burden, "year", "human_deaths", "Campaign year", "Deaths per year",
                             "Human rabies deaths — annual", key="annual_human_deaths")
        with right:
            comparison_chart(burden, "year", "cumulative_human_deaths", "Campaign year", "Cumulative deaths",
                             "Human rabies deaths — cumulative", comparator="cumulative_comparator_human_deaths",
                             key="cumulative_human_deaths")
        st.caption(f"Cumulative curves sum unrounded outcomes from campaign year 1; they exclude the lead-in. "
                   f"Clinical dog cases averted: {number(s['program_dog_cases_averted'])}. "
                   "Human deaths denote fatal clinical onsets. Counts are deterministic expected values.")
        with st.expander("Annual disease burden table"):
            show_table(annual, {"year": "Year", "clinical_dog_cases": "Dog cases", "comparator_clinical_dog_cases": "Comparator dog cases",
                "human_deaths": "Human deaths", "comparator_human_deaths": "Comparator human deaths",
                "incremental_deaths_averted": "Human deaths averted"},
                {"human_deaths": 1, "comparator_human_deaths": 1, "incremental_deaths_averted": 1}, total=True)

        st.subheader("Vaccination coverage and population immunity")
        left, right = st.columns(2)
        with left:
            comparison_chart(annual, "year", "coverage", "Campaign year", "Annual dose coverage",
                             "Annual vaccination coverage", percent=True, key="vaccination_coverage")
        with right:
            comparison_chart(monthly, "month", "immune_percent", "Month after campaign start", "Immune dogs (%)",
                             "Current dog population immunity", key="population_immunity")
        st.caption(f"Coverage uses each year's opening population. Dog population is {number(annual.opening_dogs.iloc[0])} "
                   f"at campaign start and {number(annual.closing_dogs.iloc[-1])} at the final year end. "
                   "Current immunity includes protected and still-immune waning dogs; it is not vaccination coverage.")
        with st.expander("Vaccination and immunity table"):
            show_table(annual, {"year": "Year", "opening_dogs": "Opening dogs", "dog_doses": "Vaccinations",
                "coverage": "Coverage", "mean_dog_immunity": "Mean immunity", "end_dog_immunity": "Year-end immunity",
                "closing_dogs": "Year-end dogs"}, {"coverage": 2, "mean_dog_immunity": 1, "end_dog_immunity": 1},
                percent=("coverage", "mean_dog_immunity", "end_dog_immunity"))

        st.subheader("PEP use and true exposure")
        left, right = st.columns(2)
        with left:
            comparison_chart(annual, "year", "pep_courses", "Campaign year", "Courses per year",
                             "Total PEP demand", key="total_pep")
        with right:
            comparison_chart(annual, "year", "pep_true_exposures", "Campaign year", "Courses per year",
                             "PEP following true exposure", key="true_exposure_pep")
        st.caption(f"True-exposure PEP totals {number(s['pep_true_exposures'])} courses; other PEP totals "
                   f"{number(s['pep_non_exposures'])}. The panels have separate vertical scales. "
                   f"Total PEP demand = {p['pep_multiplier']:g} × max(actual exposures, D00 null exposures) × PEP uptake.")
        show_table(annual, {"year": "Year", "pep_courses": "Total PEP", "pep_true_exposures": "True-exposure PEP",
            "pep_non_exposures": "Other PEP", "comparator_pep_courses": "Comparator PEP",
            "expected_detected_cases": "Expected detected cases"}, {"expected_detected_cases": 1}, total=True)
        st.caption("Other PEP is residual modeled demand, not a judgment about clinical appropriateness. "
                   "Courses are not injections or unique recipients. Totals use unrounded values.")

        st.subheader("Dog deaths and human case detection")
        left, right = st.columns(2)
        with left:
            comparison_chart(annual, "year", "dog_rabies_deaths", "Campaign year", "Deaths per year",
                             "Dog rabies deaths", key="dog_deaths")
            st.caption(f"Program: {number(annual.dog_rabies_deaths.sum())}; comparator: "
                       f"{number(annual.comparator_dog_rabies_deaths.sum())}. These are rabies death exits, distinct from clinical dog cases.")
        with right:
            comparison_chart(annual, "year", "expected_detected_cases", "Campaign year", "Expected cases per year",
                             "Expected detected human cases", key="detected_cases")
            st.caption(f"Campaign total: {number(s['expected_detected_cases'], 1)}. Detection is "
                       f"{p['cdr_before']:.0%} before the campaign and {p['cdr_after']:.0%} afterward; it changes reporting only.")

        history_panel = st.expander("Lead-in trajectories", key="lead_in", on_change="rerun")
        if history_panel.open:
            with history_panel:
                timeline = pd.DataFrame(result["full_timeline"])
                historical = timeline[timeline.model_year <= p["lead_in_years"]]
                for field, label in [("clinical_dog_cases", "Clinical dog cases"), ("human_deaths", "Human rabies deaths")]:
                    path = timeline[["model_year", field]].rename(columns={field: "Value"}).assign(Arm="Program pathway")
                    comp = annual[["year", "comparator_" + field]].rename(columns={"year": "model_year", "comparator_" + field: "Value"})
                    comp["model_year"] += p["lead_in_years"]
                    if len(historical):
                        comp = pd.concat([historical[["model_year", field]].tail(1).rename(columns={field: "Value"}), comp])
                    comp["Arm"] = "Comparator after lead-in"
                    data = pd.concat([path, comp], ignore_index=True)
                    lines = alt.Chart(data).mark_line().encode(
                        x=alt.X("model_year:Q", title="Model year"),
                        y=alt.Y("Value:Q", title=label + " per year", scale=alt.Scale(zero=True)),
                        color=alt.Color("Arm:N", scale=alt.Scale(domain=["Program pathway", "Comparator after lead-in"], range=["#21618c", "#777777"]), legend=alt.Legend(title=None, orient="bottom")),
                        strokeDash=alt.StrokeDash("Arm:N", legend=None),
                        tooltip=["model_year", "Arm", alt.Tooltip("Value:Q", format=",.2f")])
                    boundary = alt.Chart(pd.DataFrame({"model_year": [p["lead_in_years"]]})).mark_rule(color="#999999", strokeDash=[3, 3]).encode(x="model_year:Q")
                    st.altair_chart((lines + boundary).properties(height=230), width="stretch",
                                   alt=f"{label} through the lead-in and campaign")
                st.caption(f"The dashed vertical line marks campaign start after model year {p['lead_in_years']}. "
                           "Annual totals are plotted at year end. The lead-in is numerical initialization; year zero is not a zero-burden observation.")

        with st.expander("Methods and interpretation"):
            st.markdown(f"The supplied Bengaluru adapter runs the GDREP V2.1.25 dog, human and economic equations. "
                f"The {p['lead_in_years']}-year vaccination lead-in follows 10,000 unvaccinated weekly initialization steps. "
                "Populations are anchored at campaign start. Vaccination counts are converted to coverage using each year's opening dog population.")
            st.markdown("Inputs define a deterministic scenario; they are not fitted to observed Bengaluru annual rabies or PEP counts. "
                "No uncertainty interval is implied. D00 has no dog vaccination and no PEP and is distinct from the comparator shown here. "
                "Monthly immunity interpolates immune and total counts before division. Annual mean immunity uses immune dog-years divided by all dog-years.")
            for warning in result["warnings"]:
                st.markdown(f"- {warning}")
            stability = result["stability"]
            if "case_relative_range" in stability:
                st.caption(f"Five-year continuation of baseline vaccination: clinical dog case range "
                    f"{stability['case_relative_range']:.2%}; opening dog population range {stability['population_relative_range']:.3%}. "
                    "The stability criterion is 1% for each.")
            st.download_button("Download model notes", (ROOT / "model/GDREP_V2_1_25years/MODEL_NOTES.md").read_text(encoding="utf-8"),
                               "MODEL_NOTES.md", "text/markdown")

if economics_tab.open:
    with economics_tab:
        st.markdown(f"**{years}-year economic comparison.** Costs are in {p['currency']} at supplied unit prices, without escalation. "
                    f"Costs are discounted at {p['cost_discount']:.1%} and health outcomes at {p['outcome_discount']:.1%}, "
                    "with year-end timing from campaign start.")
        top = st.columns(4)
        top[0].metric("Total program cost", money(s["total_cost"]), border=True)
        top[1].metric("Discounted net additional cost", money(e["net_additional_cost_pv"]), border=True)
        top[2].metric("Cost per death averted", money(e["cost_per_death_averted"]), border=True)
        top[3].metric("Cost per model DALY averted", money(e["cost_per_model_daly_averted"], 2), border=True)
        st.caption("Total program cost is undiscounted. Incremental cost and cost-effectiveness ratios use discounted outcomes relative to the comparator.")

        left, right = st.columns(2)
        with left:
            st.markdown("**Annual program cost by component**")
            cost_labels = {"dog_vaccination_cost": "Dog vaccination", "pep_cost": "PEP",
                           "suspect_exposure_cost": "Quarantine and tests", "other_cost": "Investigations"}
            data = annual[["year", *cost_labels]].rename(columns=cost_labels).melt(id_vars="year", var_name="Component", value_name="Cost")
            chart = alt.Chart(data).mark_bar().encode(
                x=alt.X("year:O", title="Campaign year"),
                y=alt.Y("Cost:Q", title="Undiscounted cost (USD)", axis=alt.Axis(format="$,.0f")),
                color=alt.Color("Component:N", scale=alt.Scale(domain=list(cost_labels.values()),
                    range=["#21618c", "#97b4c7", "#a05a2c", "#777777"]), legend=alt.Legend(title=None, orient="bottom", columns=2)),
                tooltip=["year", "Component", alt.Tooltip("Cost:Q", format="$,.0f")]).properties(height=270)
            st.altair_chart(chart, width="stretch", key="annual_cost_components", alt="Annual program costs stacked by component")
            pep_share = s["pep_cost"] / s["total_cost"] if s["total_cost"] else 0
            st.caption(f"PEP accounts for {pep_share:.1%} of the total program cost.")
        with right:
            st.markdown("**Additional spending and direct savings**")
            labels = {"additional_program_cost": "Additional spending", "direct_savings": "Direct savings", "net_additional_cost": "Net additional cost"}
            data = benefit[["year", *labels]].rename(columns=labels).melt(id_vars="year", var_name="Measure", value_name="Cost")
            chart = alt.Chart(data).mark_bar().encode(
                x=alt.X("year:O", title="Campaign year"), xOffset="Measure:N",
                y=alt.Y("Cost:Q", title="Undiscounted cost (USD)", axis=alt.Axis(format="$,.0f")),
                color=alt.Color("Measure:N", scale=alt.Scale(domain=list(labels.values()), range=["#21618c", "#a05a2c", "#777777"]),
                                legend=alt.Legend(title=None, orient="bottom", columns=1)),
                tooltip=["year", "Measure", alt.Tooltip("Cost:Q", format="$,.0f")]).properties(height=270)
            st.altair_chart(chart, width="stretch", key="incremental_costs", alt="Annual additional spending, direct savings and net additional cost")
            st.caption("Net additional cost equals additional spending minus direct savings.")

        st.subheader("Annual program and comparator costs")
        show_table(annual, {"year": "Year", "dog_vaccination_cost": "Dog vaccination (USD)", "pep_cost": "PEP (USD)",
            "suspect_exposure_cost": "Quarantine and tests (USD)", "other_cost": "Investigations (USD)",
            "total_cost": "Program total (USD)", "comparator_total_cost": "Comparator total (USD)"}, total=True)
        st.subheader("Incremental costs and health gains")
        show_table(benefit, {"year": "Year", "additional_dog_cost": "Additional vaccination (USD)",
            "direct_savings": "Direct savings (USD)", "net_additional_cost": "Net additional cost (USD)",
            "deaths_averted": "Deaths averted"}, {"deaths_averted": 1}, total=True)
        st.caption("Annual tables are undiscounted. Totals use unrounded values, so displayed components may not sum exactly.")

        st.subheader("Discounted economic outcomes")
        rows = [("Additional program spending (USD)", number(e["additional_program_cost_pv"])),
                ("Avoided PEP expenditure (USD)", number(e["pep_savings_pv"])),
                ("Avoided response cost (USD)", number(e["response_savings_pv"])),
                ("Total direct savings (USD)", number(e["direct_savings_pv"])),
                ("Net additional cost (USD)", number(e["net_additional_cost_pv"])),
                ("Net financial benefit (USD)", number(e["net_financial_benefit_pv"])),
                ("Direct financial benefit–cost ratio", number(e["direct_financial_bcr"], 3)),
                ("Deaths averted", number(e["discounted_deaths_averted"], 1)),
                ("Model DALYs averted", number(e["discounted_dalys_averted"])),
                ("Mortality valuation breakeven (USD per death averted)", number(e["mortality_breakeven"], 2))]
        st.dataframe(pd.DataFrame(rows, columns=["Measure", "Discounted result"]), hide_index=True, alt="Discounted economic outcomes and benefit cost ratio")
        st.markdown("Direct financial benefits comprise avoided PEP, quarantine and testing costs. Equal investigation costs cancel. "
                    "The financial benefit–cost ratio divides direct savings by additional program spending. "
                    "Mortality benefits remain unpriced; the breakeven value is a threshold, not an estimate of willingness to pay.")

        st.subheader("Supporting unit costs")
        left, right = st.columns(2)
        with left:
            st.markdown("**PEP course components**")
            rows = [("Materials per injection", money(p["pep_material"], 2)), ("Overhead per visit", money(p["pep_overhead"], 2)),
                    ("Vaccine per dose", money(p["pep_vaccine"], 2)), ("Doses per visit", number(p["pep_doses_per_visit"], 1)),
                    ("Visits per course", number(p["pep_visits"], 1)), ("RIG per recipient", money(p["rig_cost"], 2)),
                    ("Recipients receiving RIG", f"{p['rig_fraction']:.0%}"), ("Calculated cost per course", money(result["derived"]["pep_cost_per_course"], 3))]
            st.dataframe(pd.DataFrame(rows, columns=["Component", "Value"]), hide_index=True, alt="PEP unit costs and calculated course cost")
        with right:
            st.markdown("**Dog vaccination and response**")
            workload = p["other_program_cost"] / p["investigation_cost"] if p["investigation_cost"] else None
            rows = [("Dog vaccination per dose", money(p["dog_cost"], 2)), ("Quarantine per event", money(p["quarantine_cost"], 2)),
                    ("Quarantine probability", f"{p['quarantine_probability']:.7g}"), ("Laboratory test per event", money(p["test_cost"], 2)),
                    ("Testing probability", f"{p['test_probability']:.9f}"), ("Investigations per year", number(workload)),
                    ("Cost per investigation", money(p["investigation_cost"], 2)), ("Annual investigations cost per arm", money(p["other_program_cost"]))]
            st.dataframe(pd.DataFrame(rows, columns=["Item", "Value"]), hide_index=True, alt="Dog vaccination, quarantine, testing and fixed investigation costs")
        with st.expander("Economic assumptions and interpretation"):
            st.markdown("Quarantine and test probabilities apply to suspected exposures. Investigations are a fixed annual workload in both arms; "
                "the exposure-based investigation probability is zero. Costs use the supplied price basis, with no price escalation. "
                "Model DALYs retain the source model's age-weighted years-of-life-lost calculation. "
                "Travel, household time, productivity, nonfatal bite morbidity and animal welfare are not monetized. "
                "No benefits beyond the selected campaign horizon are credited.")
            st.caption(f"Average age at death: {p['average_age_at_death']:g} years; life expectancy: {p['life_expectancy']:g} years. "
                       "Shared historical costs cancel from the incremental comparison.")

if st.session_state.get("editing_inputs", False):
    edit_inputs()
