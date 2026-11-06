"""Explicit government/program cost accounting and original YLL formula."""
import math
import numpy as np
import pandas as pd


def unit_costs(e):
    dog = (e.vaccine_cost_per_dose+e.dog_materials_cost_per_dose+e.dog_marking_cost_per_dose)*(1+e.vaccine_wastage_markup)
    dog += e.campaign_fixed_cost/e.reference_campaign_doses
    # Source workbook convention: overhead is also multiplied by doses/visit.
    pep = (e.pep_material_cost_per_injection+e.pep_overhead_cost+e.pep_vaccine_cost_per_dose)*e.pep_doses_per_visit*e.pep_visits
    pep += e.rig_cost*e.rig_fraction
    suspect = e.quarantine_probability*e.quarantine_unit_cost+e.testing_probability*e.testing_unit_cost+e.investigation_probability*e.investigation_unit_cost
    return {'dog_cost_per_budgeted_dose': dog, 'pep_cost_per_course': pep, 'suspect_cost_per_exposure': suspect}


def yll_per_death(e):
    """Integral from 0..L of exp(-r*t) * [K*C*(a+t)*exp(-beta*(a+t)) + (1-K)]."""
    a = e.average_age_at_death; L = e.yll_life_expectancy-a
    r, beta, K, C = e.outcome_discount_rate, e.yll_age_beta, e.yll_age_weight, e.yll_age_constant
    q = r+beta
    if q < 1e-7:
        # Stable convergent expansion of integrals exp(-q*t), t*exp(-q*t).
        int0 = sum((-q)**j*L**(j+1)/(math.factorial(j)*(j+1)) for j in range(9))
        int1 = sum((-q)**j*L**(j+2)/(math.factorial(j)*(j+2)) for j in range(9))
    else:
        int0 = -math.expm1(-q*L)/q
        int1 = (-math.expm1(-q*L)-q*L*math.exp(-q*L))/(q*q)
    unweighted = L if r == 0 else -math.expm1(-r*L)/r
    return K*C*math.exp(-beta*a)*(a*int0+int1)+(1-K)*unweighted


def evaluate(annual, config):
    """Compare intervention with the explicitly configured status-quo program."""
    out = annual.copy(); e = config.economics; u = unit_costs(e); yll = yll_per_death(e)
    for prefix, other in [('', e.intervention_other_annual_cost), ('comparator_', e.comparator_other_annual_cost)]:
        out[prefix+'dog_vaccination_cost'] = out[prefix+'planned_dog_doses']*u['dog_cost_per_budgeted_dose']
        out[prefix+'pep_cost'] = out[prefix+'pep_courses']*u['pep_cost_per_course']
        out[prefix+'suspect_exposure_cost'] = out[prefix+'suspected_exposures']*u['suspect_cost_per_exposure']
        out[prefix+'other_cost'] = other
        out[prefix+'total_cost'] = sum(out[prefix+k] for k in ['dog_vaccination_cost', 'pep_cost', 'suspect_exposure_cost', 'other_cost'])
    out['incremental_cost'] = out.total_cost-out.comparator_total_cost
    out['incremental_deaths_averted'] = out.comparator_human_deaths-out.human_deaths
    out['incremental_dalys_averted'] = out.incremental_deaths_averted*yll
    out['discounted_incremental_cost'] = out.incremental_cost/(1+e.cost_discount_rate)**out.year
    out['discounted_deaths_averted'] = out.incremental_deaths_averted/(1+e.outcome_discount_rate)**out.year
    out['discounted_dalys_averted'] = out.incremental_dalys_averted/(1+e.outcome_discount_rate)**out.year
    for key in ['incremental_cost', 'incremental_deaths_averted', 'incremental_dalys_averted',
                'discounted_incremental_cost', 'discounted_deaths_averted', 'discounted_dalys_averted']:
        out['cumulative_'+key] = out[key].cumsum()
    out['discounted_total_cost'] = out.total_cost/(1+e.cost_discount_rate)**out.year
    out['discounted_comparator_total_cost'] = out.comparator_total_cost/(1+e.cost_discount_rate)**out.year
    return out, {**u, 'yll_per_death': yll}


def ratio_and_status(cost, benefit):
    if benefit > 1e-12:
        return cost/benefit, 'Cost-saving with health gain' if cost < 0 else 'Additional cost and health gain'
    if benefit < -1e-12:
        return None, 'More costly and less effective' if cost >= 0 else 'Lower cost with health loss'
    return None, 'No incremental health benefit'
