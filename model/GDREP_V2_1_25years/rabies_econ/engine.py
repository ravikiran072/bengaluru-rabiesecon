"""Country-independent RabiesEcon simulation and sensitivity analysis."""
from dataclasses import dataclass, replace
from functools import lru_cache
import io
import json
import math
import zipfile
import numpy as np
import pandas as pd
from .config import Config, ModelError
from .dog import DogModel
from .economics import evaluate, ratio_and_status

VERSION = '2.1.25years'

@dataclass
class Result:
    config: Config
    annual: pd.DataFrame
    summary: dict
    checks: dict
    derived: dict

    def download_zip(self):
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, 'w', zipfile.ZIP_DEFLATED) as archive:
            archive.writestr('inputs.json', self.config.to_json())
            archive.writestr('annual_results.csv', self.annual.to_csv(index=False))
            archive.writestr('summary.json', json.dumps(self.summary, indent=2, allow_nan=False))
            archive.writestr('run_metadata.json', json.dumps({'model_version': VERSION, 'checks': self.checks, 'derived': self.derived}, indent=2, allow_nan=False))
        return buffer.getvalue()


def human_run(config, infectious_density, pep_by_year, initial_population):
    """Human counts; dog input is infectious dogs/km². Return annual clinical onsets and death exits separately."""
    t = config.transmission
    x = np.array([initial_population, 0., 0., 0.])
    mh, b = 1/t.human_life_expectancy_years/52, t.human_births_per_1000_year/1000/52
    pi, pn = t.human_clinical_rate_per_week, t.human_nonclinical_rate_per_week
    rows = []; minimum = 0.; residual = 0.; max_hazard = 0.
    for y, pep in enumerate(pep_by_year):
        opening = x.sum(); totals = np.zeros(6)
        for w in range(52):
            s, e, i, immune = x
            hazard = t.z*infectious_density[y*52+w]
            max_hazard = max(max_hazard, float(hazard))
            if hazard+mh > 1:
                raise ModelError('Human exposure probability exceeds one per week. No result was clipped.')
            exposure = hazard*s
            onset = e*pi*(1-pep*t.pep_efficacy)
            protected = e*pi*pep*t.pep_efficacy
            births = b*(s+e+immune)
            death = t.human_infectious_exit_per_week*i
            natural = mh*x.sum()
            totals += [onset, exposure, x.sum()/52, births, death, natural]
            x += np.array([births+pn*e-mh*s-exposure, exposure-(mh+pi+pn)*e,
                           onset-mh*i-death, protected-mh*immune])
            minimum = min(minimum, float(x.min()))
        residual = max(residual, abs(opening+totals[3]-totals[4]-totals[5]-x.sum()))
        rows.append(dict(zip(['human_deaths', 'human_exposures', 'human_person_years', 'human_births',
                              'human_death_exit_flow', 'human_natural_deaths'], totals)))
    if not np.isfinite(x).all() or minimum < -1e-12*max(1, initial_population) or residual > 1e-9*max(1, initial_population):
        raise ModelError('Human state or population conservation check failed; these inputs exceed the weekly solver domain.')
    return pd.DataFrame(rows), x, {'min_state': minimum, 'max_mass_balance_error': residual, 'max_exposure_hazard': max_hazard}


@lru_cache(maxsize=16)
def _prepare_dogs(setting, transmission, vaccination):
    # Exclude Z and all human/cost inputs from the caller's cache key; no fitted data.
    config = Config(setting=setting, transmission=transmission, vaccination=vaccination)
    m = DogModel(config)
    initial, notes = m.initial_state()
    lead = m.run(initial, [setting.historical_coverage]*setting.lead_in_years, setting.historical_campaign_weeks)
    opening = lead['state'].sum()
    if opening <= 0:
        raise ModelError('No dogs remain at program entry.')
    unit = setting.dog_population/opening
    return m, initial, lead, unit, notes


def prepare_dogs(config):
    # Normalize irrelevant inputs for cache reuse. All dog-specific values remain keys.
    default = Config().transmission
    keep = {'actual_r0', 'dog_births_per_year', 'dog_natural_mortality_per_year', 'incubation_days',
            'days_per_model_week', 'dog_clinical_probability', 'dog_infectious_exit_per_week'}
    t = replace(default, **{k: getattr(config.transmission, k) for k in keep})
    return _prepare_dogs(config.setting, t, config.vaccination)


def effective_config(config):
    """Optional stationary dog-vaccination baseline; preserve the saved pathway inputs."""
    if not config.setting.estimate_at_stabilization:
        return config
    s, v = config.setting, config.vaccination
    return replace(config,
        setting=replace(s, lead_in_years=s.stabilization_years,
                        historical_coverage=s.stabilization_coverage,
                        historical_campaign_weeks=v.campaign_weeks,
                        historical_pep_uptake=v.pep_uptake),
        vaccination=replace(v, coverage=s.stabilization_coverage,
                            comparator_coverage=s.stabilization_coverage,
                            comparator_campaign_weeks=v.campaign_weeks,
                            coverage_schedule=(), comparator_schedule=(),
                            pep_uptake_schedule=(), comparator_pep_uptake_schedule=()))


def stabilization_check(model, lead, unit, coverage, campaign_weeks):
    """Check a five-year continuation at identical annual campaign timing.

    This is a convergence diagnostic, not fitting or flattening of outputs.
    A half-case denominator keeps numerical extinction tails from dominating.
    """
    continuation = model.run(lead['state'], [coverage]*5, campaign_weeks)
    cases = np.array([r['clinical_dog_cases']*unit for r in continuation['annual']])
    populations = np.array([r['opening_dogs']*unit for r in continuation['annual']])
    case_variation = float(np.ptp(cases)/max(float(cases.max()), 0.5))
    population_variation = float(np.ptp(populations)/max(float(populations.max()), 1.0))
    passed = case_variation <= .01 and population_variation <= .01
    return {'passed': bool(passed), 'status': 'Stability check passed' if passed else 'Still stabilizing',
            'tolerance': .01, 'check_years': 5,
            'dog_case_relative_range': case_variation,
            'dog_population_relative_range': population_variation,
            'dog_cases_next_five_years': cases.tolist(),
            'case_denominator_floor': 0.5,
            'interpretation': 'Tests dog transmission and population under unchanged annual vaccination. Human population growth can still change total human exposures. No outputs are flattened.'}


def run(config=None):
    requested = (config or Config()).validate()
    c = effective_config(requested).validate()
    s, t, v = c.setting, c.transmission, c.vaccination
    model, initial, lead, unit, initial_notes = prepare_dogs(c)
    intervention = model.run(lead['state'], v.coverage_schedule or [v.coverage]*s.years, v.campaign_weeks)
    comparator = model.run(lead['state'], v.comparator_schedule or [v.comparator_coverage]*s.years, v.comparator_campaign_weeks)
    null = model.run(initial, [0]*(s.lead_in_years+s.years), v.campaign_weeks)
    density_scale = unit/s.area_km2
    dog_i = np.concatenate([lead['trajectory'], intervention['trajectory']])*density_scale
    dog_c = np.concatenate([lead['trajectory'], comparator['trajectory']])*density_scale
    dog_zero = null['trajectory']*density_scale
    # Anchor BOTH entered populations at program start; historical states are regenerated.
    _, human_lead_unit, _ = human_run(c, lead['trajectory']*density_scale,
                                      [s.historical_pep_uptake]*s.lead_in_years, 1.0)
    if human_lead_unit.sum() <= 0:
        raise ModelError('No humans remain at program entry.')
    starting_humans = s.human_population/human_lead_unit.sum()
    program_pep = list(v.pep_uptake_schedule or [v.pep_uptake]*s.years)
    comparator_pep = list(v.comparator_pep_uptake_schedule or [v.comparator_pep_uptake]*s.years)
    pep = [s.historical_pep_uptake]*s.lead_in_years+program_pep
    cpep = [s.historical_pep_uptake]*s.lead_in_years+comparator_pep
    nop = [0]*len(pep)
    arms = {}; checks = {}
    for name, trajectory, uptake in [('D11', dog_i, pep), ('D10', dog_i, nop),
                                      ('D01', dog_zero, pep), ('D00', dog_zero, nop),
                                      ('comparator', dog_c, cpep)]:
        frame, _, check = human_run(c, trajectory, uptake, starting_humans)
        arms[name] = frame.iloc[s.lead_in_years:].reset_index(drop=True)
        checks[name] = check
    annual = pd.DataFrame(intervention['annual'])*unit
    comp = pd.DataFrame(comparator['annual'])*unit
    for col in arms['D11']:
        annual[col] = arms['D11'][col]
        annual['comparator_'+col] = arms['comparator'][col]
    for col in comp:
        annual['comparator_'+col] = comp[col]
    for arm in ['D10', 'D01', 'D00']:
        annual[arm+'_deaths'] = arms[arm].human_deaths
    annual['total_deaths_averted_vs_no_intervention'] = annual.D00_deaths-annual.human_deaths
    annual['pep_attributed_deaths_averted'] = .5*((annual.D00_deaths-annual.D01_deaths)+(annual.D10_deaths-annual.human_deaths))
    annual['vaccine_attributed_deaths_averted'] = .5*((annual.D00_deaths-annual.D10_deaths)+(annual.D01_deaths-annual.human_deaths))
    for prefix, arm, uptake in [('', 'D11', np.asarray(program_pep)), ('comparator_', 'comparator', np.asarray(comparator_pep))]:
        annual[prefix+'suspected_exposures'] = np.maximum(arms[arm].human_exposures, arms['D00'].human_exposures)*t.suspected_exposure_multiplier
        annual[prefix+'pep_courses'] = annual[prefix+'suspected_exposures']*uptake
        annual[prefix+'mean_dog_immunity'] = annual[prefix+'immune_dog_years']/annual[prefix+'dog_years']
        annual[prefix+'end_dog_immunity'] = annual[prefix+'closing_immune_dogs']/annual[prefix+'closing_dogs']
        annual[prefix+'human_incidence_per_100k'] = annual[prefix+'human_deaths']/annual[prefix+'human_person_years']*1e5
    annual.insert(0, 'year', np.arange(1, s.years+1))
    annual['in_reporting_window'] = annual.year.between(s.report_first_year, s.report_last_year)
    annual, costs = evaluate(annual, c)
    report = annual[annual.in_reporting_window]
    means = {k: float(report[k].mean()) for k in ['human_deaths', 'comparator_human_deaths', 'clinical_dog_cases',
                'dog_rabies_deaths', 'dog_doses', 'planned_dog_doses', 'unused_dog_doses', 'pep_courses',
                'total_cost', 'comparator_total_cost', 'incremental_deaths_averted',
                'total_deaths_averted_vs_no_intervention', 'pep_attributed_deaths_averted', 'vaccine_attributed_deaths_averted']}
    means['mean_dog_immunity'] = float(report.immune_dog_years.sum()/report.dog_years.sum())
    means['human_incidence_per_100k'] = float(report.human_deaths.sum()/report.human_person_years.sum()*1e5)
    totals = {k: float(report[k].sum()) for k in ['discounted_incremental_cost', 'discounted_deaths_averted',
                'discounted_dalys_averted', 'discounted_total_cost', 'discounted_comparator_total_cost']}
    ratio, status = ratio_and_status(totals['discounted_incremental_cost'], totals['discounted_dalys_averted'])
    death_ratio, _ = ratio_and_status(totals['discounted_incremental_cost'], totals['discounted_deaths_averted'])
    summary = {'name': s.name, 'reporting_years': [s.report_first_year, s.report_last_year],
               'annual_means': means, 'reporting_window_totals': totals,
               'cost_per_daly_averted': ratio, 'cost_per_death_averted': death_ratio, 'economic_interpretation': status}
    attribution_error = float(np.abs(annual.total_deaths_averted_vs_no_intervention-
                                   annual.pep_attributed_deaths_averted-annual.vaccine_attributed_deaths_averted).max())
    if attribution_error > 1e-8*max(1, annual.D00_deaths.max()):
        raise ModelError('Matched counterfactual attribution did not reconcile.')
    derived = {**costs, **initial_notes, 'carrying_capacity_dogs': float(unit),
               'carrying_capacity_dogs_per_km2': float(density_scale),
               'beta_density_per_week': float(model.beta_k/density_scale),
               'actual_r0': t.actual_r0, 'human_density_per_km2': s.human_population/s.area_km2,
               'dog_density_per_km2': s.dog_population/s.area_km2,
               'historical_start_humans': float(starting_humans),
               'death_measure': 'Clinical human rabies onsets, retained as the burden proxy; death exits are exported separately.'}
    all_checks = {'dog_lead_in': lead['checks'], 'dogs': intervention['checks'], 'comparator_dogs': comparator['checks'],
                  'no_intervention_dogs': null['checks'], 'humans': checks, 'attribution_error': attribution_error,
                  'dog_population_anchor_error': float(abs(annual.opening_dogs.iloc[0]-s.dog_population)),
                  'human_population_anchor_error': float(abs(human_lead_unit.sum()*starting_humans-s.human_population))}
    if requested.setting.estimate_at_stabilization:
        all_checks['stabilization'] = stabilization_check(model, lead, unit, v.coverage, v.campaign_weeks)
        derived['stabilization'] = {
            'enabled': True, 'numerical_lead_in_years': s.lead_in_years,
            'dog_coverage': v.coverage, 'campaign_weeks': v.campaign_weeks,
            'pep_uptake': v.pep_uptake, 'comparator_pep_uptake': v.comparator_pep_uptake,
            'scope': 'Current dog coverage held constant during initialization and all assessment years; no GDREP progression. PEP uptake is constant within each arm. Prior schedules are retained in saved inputs but not applied in this mode.',
            'lead_in_interpretation': 'Numerical initialization, not an assumed historical vaccination record.'}
    return Result(requested, annual, summary, all_checks, derived)


def set_parameter(config, path, value):
    if config.setting.estimate_at_stabilization and path == 'vaccination.coverage':
        path = 'setting.stabilization_coverage'
    parts = path.split('.')
    if len(parts) != 2 or parts[0] not in config.to_dict():
        raise ModelError('Parameter paths have the form transmission.z or vaccination.coverage.')
    section, name = parts
    if name not in config.to_dict()[section]:
        raise ModelError(f'Unknown parameter: {path}')
    baseline = getattr(getattr(config, section), name)
    if type(baseline) is int and isinstance(value, float) and value.is_integer():
        value = int(value)
    return replace(config, **{section: replace(getattr(config, section), **{name: value})}).validate()


def _sensitivity_row(result):
    return {**result.summary['annual_means'], **result.summary['reporting_window_totals'],
            'cost_per_daly_averted': result.summary['cost_per_daly_averted'],
            'economic_interpretation': result.summary['economic_interpretation']}


def sensitivity_grid(config, z_values=None, campaign_weeks=(4, 10, 26), progress=None):
    """Deterministic joint scenarios. Include the exact entered baseline even if off-grid."""
    z_values = list(z_values if z_values is not None else [config.transmission.z/2, config.transmission.z, config.transmission.z*2])
    z_values = list(dict.fromkeys([config.transmission.z]+z_values))
    durations = list(dict.fromkeys([config.vaccination.campaign_weeks]+list(campaign_weeks)))
    rows = []
    for weeks in durations:
        for z in z_values:
            cfg = set_parameter(set_parameter(config, 'transmission.z', z), 'vaccination.campaign_weeks', weeks)
            result = run(cfg)
            rows.append({'z': z, 'campaign_weeks': weeks,
                         'is_entered_baseline': z == config.transmission.z and weeks == config.vaccination.campaign_weeks,
                         **_sensitivity_row(result)})
            if progress:
                progress(len(rows), len(durations)*len(z_values))
    return pd.DataFrame(rows)


def one_way(config, parameter, values, progress=None):
    """Change only one supplied input; all other entered values remain fixed."""
    if parameter == 'vaccination.coverage' and config.vaccination.coverage_schedule and not config.setting.estimate_at_stabilization:
        raise ModelError('Clear the annual coverage schedule before varying constant coverage.')
    for field, schedule in [('pep_uptake', config.vaccination.pep_uptake_schedule),
                            ('comparator_pep_uptake', config.vaccination.comparator_pep_uptake_schedule)]:
        if parameter == 'vaccination.'+field and schedule and not config.setting.estimate_at_stabilization:
            raise ModelError('Clear the annual PEP schedule before varying constant PEP uptake.')
    parts = parameter.split('.')
    if len(parts) != 2 or parts[0] not in config.to_dict() or parts[1] not in config.to_dict()[parts[0]]:
        raise ModelError(f'Unknown parameter: {parameter}')
    base = getattr(getattr(config, parameter.split('.')[0]), parameter.split('.')[1])
    if config.setting.estimate_at_stabilization and parameter == 'vaccination.coverage':
        base = config.setting.stabilization_coverage
    values = list(dict.fromkeys([base]+list(values)))
    rows = []
    for value in values:
        result = run(set_parameter(config, parameter, value))
        rows.append({'parameter': parameter, 'value': value, 'is_entered_baseline': value == base, **_sensitivity_row(result)})
        if progress:
            progress(len(rows), len(values))
    return pd.DataFrame(rows)


def critical_coverage(config):
    """First passing annual coverage at whole-percentage resolution, fully enumerated.

    Apply coverage during every program year, then maintenance coverage. Require
    expected dog cases < threshold in the last program year and every maintenance year.
    Exhaustive enumeration avoids assuming monotonicity near elimination.
    """
    config.validate()
    if config.setting.estimate_at_stabilization:
        raise ModelError('Switch off stabilization mode before evaluating a changing vaccination program with critical coverage.')
    s, v = config.setting, config.vaccination
    model, _, lead, unit, _ = prepare_dogs(config)
    grid = []; found = None; best_annual = None
    for pct in range(101):
        coverage = pct/100
        schedule = [coverage]*s.years+[v.maintenance_coverage]*v.maintenance_years
        result = model.run(lead['state'], schedule, v.campaign_weeks)
        annual = pd.DataFrame(result['annual'])*unit
        worst = float(annual.clinical_dog_cases.iloc[s.years-1:].max())
        passes = worst < v.elimination_threshold_cases
        grid.append({'coverage': coverage, 'worst_cases_in_evaluation_window': worst, 'passes': passes})
        if passes:
            found = coverage; best_annual = annual
            break
    if best_annual is not None:
        best_annual.insert(0, 'year', np.arange(1, len(best_annual)+1))
    return {'coverage': found, 'resolution': .01, 'threshold_cases': v.elimination_threshold_cases,
            'evaluation_years': [s.years, s.years+v.maintenance_years],
            'grid': pd.DataFrame(grid), 'annual': best_annual,
            'status': 'No feasible coverage through 100%' if found is None else 'Minimum passing whole-percentage coverage'}
