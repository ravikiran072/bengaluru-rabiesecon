# Bangalore model: equations and use

This application adapts the RabiesEcon framework described in the publications
below to the Bangalore vaccination scenario. See [METHODS.md](../../METHODS.md)
for scenario assumptions and outcome definitions.

## Run

After following the installation steps in [README.md](../../README.md), run these
commands from the `bangalore_rabiesecon` directory:

```bash
# Open the dashboard
python -m streamlit run app.py

# Recalculate and export the report scenario
python reproduce.py
```

On Windows, `start_dashboard.bat` installs the dependencies and opens the dashboard.

## Inputs

The app and reproduction script load `bangalore_inputs.json`: a 50-year lead-in
followed by five campaign years, with 90% PEP uptake in both arms.

Use **Edit inputs** in the app, or supply a JSON file:

```bash
python reproduce.py --inputs alternative_inputs.json --out alternative_results
```

Omitted fields inherit the report inputs. When editing a scenario:

- Keep `program_values`, `comparator_values`, `program_pep` and
  `comparator_pep` the same length: one value per campaign year.
- Vaccination values are annual dose counts when the corresponding mode is
  `count`, or fractions when it is `coverage`. PEP uptake and vaccine efficacy
  are fractions between 0 and 1.
- Vaccination counts must satisfy population and dose-eligibility constraints.
  Invalid schedules raise an error.
- Human and dog populations refer to campaign start, after the lead-in.
  Area is in km², incubation in days and campaign duration in model weeks.
  The model uses 52 weekly steps per year and continuous expected counts.

## Core equations

These are discrete weekly updates. Dog states are normalized to carrying capacity
`K`; human states are counts. Annual outcomes sum the weekly flows.

### Dog transmission and vaccination

`S`, `E`, `I`, `Vp` and `Vw` denote susceptible, exposed, infectious, protected
immune and waning immune dogs. Weekly rates and the force of infection are:

```text
n     = S + E + I + Vp + Vw
b     = dog_births_per_year / 52
d     = dog_natural_mortality_per_year / 52
sigma = days_per_model_week / incubation_days
r     = dog_clinical_probability
mu    = dog_infectious_exit_per_week
omega = waning_per_week
m     = d + (b - d)*n
betaK = R0*(sigma + b)*(mu + b)/(sigma*r)
lambda = betaK*I

R0 = betaK*sigma*r / ((sigma + b)*(mu + b))
```

Vaccination occurs before each weekly update. The annual dose budget is the
entered count, or coverage multiplied by opening dog population. Weekly allocation
is that budget divided by campaign weeks, and zero outside the scheduled campaign.
Eligible dogs were never vaccinated or received their last dose at least
`revaccination_weeks` earlier. For normalized eligible population `N_eligible`,
the selected amount in compartment `j` and transfers when `protected_weeks > 0` are:

```text
q_j = min(weekly_doses/K, N_eligible)*X_j,eligible/N_eligible

S  -> Vp: efficacy*q_S
E  -> Vp: exposed_protection*q_E
Vw -> Vp: efficacy*q_Vw
```

Selection is zero when no dogs are eligible. All recipients reset their
revaccination clock, including vaccine failures. Infectious dogs receive no
therapeutic benefit. With zero protected weeks, new immunity enters `Vw` directly
and waning immune recipients remain in `Vw`.

Using the post-vaccination states, the compartment updates are:

```text
S_next  = S*(1 - m) - lambda*S + sigma*(1 - r)*E + omega*Vw + b*n
E_next  = E*(1 - m - sigma) + lambda*S
I_next  = I*(1 - m - mu) + sigma*r*E
Vp_next = Vp*(1 - m) - G
Vw_next = Vw*(1 - m - omega) + G
```

`G` is the surviving protected cohort completing `protected_weeks` after the
weekly update. It enters `Vw` and first wanes the following week. Births enter
susceptible, vaccine-eligible dogs. In dog counts, weekly clinical cases are
`K*sigma*r*E` and rabies deaths are `K*mu*I`.

### Dog-to-human transmission

`SH`, `EH`, `IH` and `RH` denote susceptible, exposed, clinically infectious and
protected humans. `u` is PEP uptake and `eH` is PEP efficacy. The weekly clinical,
nonclinical and infectious-exit rates are `pi`, `pn` and `muH`.

```text
bH = human_births_per_1000_year / 1000 / 52
dH = 1 / human_life_expectancy_years / 52
infectious_dog_density = K*I / area_km2
exposure = Z*SH*infectious_dog_density
onset = pi*EH*(1 - u*eH)
protected = pi*EH*u*eH
births = bH*(SH + EH + RH)

SH_next = SH + births + pn*EH - dH*SH - exposure
EH_next = EH + exposure - (dH + pi + pn)*EH
IH_next = IH + onset - (dH + muH)*IH
RH_next = RH + protected - dH*RH
```

`Z` has units km² per infectious dog per week. Reported human deaths sum `onset`;
the separate death-exit flow is `muH*IH`. Protection has no waning in this model.

### PEP demand and intervention attribution

For annual exposure totals `X_y` and `X00_y`, suspected exposures and PEP use are:

```text
suspected_exposures_y = multiplier*max(X_y, X00_y)
PEP_courses_y = u_y*suspected_exposures_y
true_exposure_PEP_y = u_y*X_y
```

`X00_y` is exposure without dog vaccination or PEP throughout the modeled history.
The maximum is applied to annual totals. Thus PEP demand can remain high as actual
exposure falls.

Let `Dij` denote annual human deaths with dog vaccination `i` and PEP `j`, where
1 means the entered history and program and 0 means absent throughout. Matched
counterfactuals share historical starting populations:

```text
total_averted = D00 - D11
PEP_attribution = 0.5*((D00 - D01) + (D10 - D11))
vaccination_attribution = 0.5*((D00 - D10) + (D01 - D11))
```

The economic comparator instead shares the program's lead-in, then follows its
own vaccination and PEP schedules.

### Costs and cost-effectiveness

Unit costs and annual cost in each arm are:

```text
dog_unit_cost = (vaccine + materials + marking)*(1 + wastage_markup)
                + campaign_fixed_cost/reference_campaign_doses
PEP_course_cost = (PEP_materials + PEP_overhead + PEP_vaccine)*doses_per_visit*visits
                  + RIG_cost*RIG_fraction
suspect_unit_cost = quarantine_probability*quarantine_cost
                    + testing_probability*testing_cost
                    + investigation_probability*investigation_cost
annual_cost = planned_dog_doses*dog_unit_cost + PEP_courses*PEP_course_cost
              + suspected_exposures*suspect_unit_cost + other_annual_cost
```

DALYs retain the source model's age-weighted years of life lost (YLL), with no
nonfatal disability component. For age at death `a`, remaining lifespan
`L = life_expectancy - a`, outcome discount rate `rH`, and age-weight parameters
`kA`, `cA`, `betaA`:

```text
YLL = integral from t=0 to L of:
      exp(-rH*t)*[kA*cA*(a + t)*exp(-betaA*(a + t)) + (1 - kA)] dt

deaths_averted_y = comparator_deaths_y - program_deaths_y
DALYs_averted_y = deaths_averted_y*YLL
incremental_cost_y = program_cost_y - comparator_cost_y
PV_cost = sum_y incremental_cost_y/(1 + rC)^y
PV_deaths = sum_y deaths_averted_y/(1 + rH)^y
PV_DALYs = sum_y DALYs_averted_y/(1 + rH)^y
cost_per_death_averted = PV_cost/PV_deaths
cost_per_DALY_averted = PV_cost/PV_DALYs
```

`rC` is the cost discount rate and `y` starts at campaign year 1. Discounting is
at year end; within-death YLL discounting is separate. Ratios are reported only
for positive health gains. The direct financial benefit–cost ratio is discounted
PEP and response-cost savings divided by discounted additional vaccination and
other program spending, when that spending is positive.

## Outputs

`reproduce.py` writes inputs, model configuration, numerical results, CSV tables
and PDF/SVG figures to `results/`, or the directory selected with `--out`.
The same outputs are available through **Download results and figures**.

Human deaths represent fatal clinical onsets. In the report scenario, the comparator
continues baseline vaccination after the shared lead-in. Costs use the supplied USD prices;
discounted outputs are labeled separately. See [METHODS.md](../../METHODS.md)
for the full definitions.

## Model references

1. Borse RH, et al. (2018). Cost-effectiveness of dog rabies vaccination programs
   in East Africa. *PLoS Neglected Tropical Diseases*, **12**(5), e0006490.
   [doi:10.1371/journal.pntd.0006490](https://doi.org/10.1371/journal.pntd.0006490).

2. Jeon S, et al. (2019). Determining the post-elimination level of vaccination
   needed to prevent re-establishment of dog rabies. *PLoS Neglected Tropical
   Diseases*, **13**(12), e0007869.
   [doi:10.1371/journal.pntd.0007869](https://doi.org/10.1371/journal.pntd.0007869).

3. Kunkel A, et al. (2021). The urgency of resuming disrupted dog rabies
   vaccination campaigns: a modeling and cost-effectiveness analysis.
   *Scientific Reports*, **11**, 12476.
   [doi:10.1038/s41598-021-92067-5](https://doi.org/10.1038/s41598-021-92067-5).
