"""Serializable inputs; no workbook, country table or fitted state is required."""
from dataclasses import dataclass, field, asdict, fields
import json
import math

class ModelError(ValueError):
    pass

@dataclass(frozen=True)
class Setting:
    name: str = 'Synthetic example — replace inputs'
    human_population: float = 100000.0
    dog_population: float = 10000.0
    area_km2: float = 1000.0
    years: int = 13
    report_first_year: int = 1
    report_last_year: int = 13
    lead_in_years: int = 10
    historical_coverage: float = 0.1
    historical_pep_uptake: float = 0.5
    historical_campaign_weeks: int = 10
    initial_mode: str = 'burn_in'
    burn_in_weeks: int = 10000
    initial_capacity_ratio: float = 1.1
    initial_exposed_fraction: float = 0.0002
    initial_infectious_fraction: float = 0.0001
    initial_immune_fraction: float = 0.0
    estimate_at_stabilization: bool = False
    stabilization_years: int = 25
    stabilization_coverage: float = 0.1

@dataclass(frozen=True)
class Transmission:
    actual_r0: float = 1.2
    dog_births_per_year: float = 0.55
    dog_natural_mortality_per_year: float = 0.33
    incubation_days: float = 40.0
    days_per_model_week: float = 7.17703349282297
    dog_clinical_probability: float = 0.45
    dog_infectious_exit_per_week: float = 0.7
    z: float = 0.00015
    human_births_per_1000_year: float = 17.0
    human_life_expectancy_years: float = 69.4
    human_clinical_rate_per_week: float = 0.025089371358398886
    human_nonclinical_rate_per_week: float = 0.09672050116405931
    human_infectious_exit_per_week: float = 1.0
    pep_efficacy: float = 0.93
    suspected_exposure_multiplier: float = 10.0

@dataclass(frozen=True)
class Vaccination:
    efficacy: float = 0.95
    exposed_protection: float = 0.50
    protected_weeks: int = 26
    revaccination_interval_weeks: int = 52
    waning_per_week: float = 1 / 52
    first_campaign_week: int = 22
    campaign_weeks: int = 10
    comparator_campaign_weeks: int = 10
    coverage: float = 0.7
    comparator_coverage: float = 0.1
    pep_uptake: float = 0.5
    comparator_pep_uptake: float = 0.5
    # Empty schedules use the constant coverage above; otherwise one value/year.
    coverage_schedule: tuple = ()
    comparator_schedule: tuple = ()
    pep_uptake_schedule: tuple = ()
    comparator_pep_uptake_schedule: tuple = ()
    maintenance_years: int = 10
    maintenance_coverage: float = 0.4
    elimination_threshold_cases: float = 0.5

@dataclass(frozen=True)
class Economics:
    currency: str = 'currency units'
    # Synthetic demonstration costs; replace with costs in one currency/price year.
    cost_discount_rate: float = 0.03
    outcome_discount_rate: float = 0.03
    vaccine_cost_per_dose: float = 1.0
    dog_materials_cost_per_dose: float = 0.2
    dog_marking_cost_per_dose: float = 0.05
    vaccine_wastage_markup: float = 0.1
    campaign_fixed_cost: float = 1000.0
    reference_campaign_doses: float = 10000.0
    pep_material_cost_per_injection: float = 1.0
    pep_overhead_cost: float = 5.0
    pep_vaccine_cost_per_dose: float = 3.0
    pep_doses_per_visit: float = 1.0
    pep_visits: float = 4.0
    rig_cost: float = 15.0
    rig_fraction: float = 0.5
    quarantine_probability: float = 0.0
    quarantine_unit_cost: float = 0.0
    testing_probability: float = 0.0
    testing_unit_cost: float = 0.0
    investigation_probability: float = 0.0
    investigation_unit_cost: float = 0.0
    intervention_other_annual_cost: float = 0.0
    comparator_other_annual_cost: float = 0.0
    average_age_at_death: float = 15.0
    yll_life_expectancy: float = 69.4
    yll_age_weight: float = 1.0
    yll_age_beta: float = 0.02
    yll_age_constant: float = 0.0634

@dataclass(frozen=True)
class Config:
    setting: Setting = field(default_factory=Setting)
    transmission: Transmission = field(default_factory=Transmission)
    vaccination: Vaccination = field(default_factory=Vaccination)
    economics: Economics = field(default_factory=Economics)

    def to_dict(self):
        return asdict(self)

    def to_json(self):
        return json.dumps(self.to_dict(), indent=2, allow_nan=False)

    @classmethod
    def from_dict(cls, obj):
        if not isinstance(obj, dict):
            raise ModelError('Inputs must be a JSON object.')
        known = {f.name for f in fields(cls)}
        if set(obj) - known:
            raise ModelError(f'Unknown input sections: {sorted(set(obj) - known)}')
        result = {}
        for name, kind in [('setting', Setting), ('transmission', Transmission),
                           ('vaccination', Vaccination), ('economics', Economics)]:
            section = obj.get(name, {})
            if not isinstance(section, dict):
                raise ModelError(f'{name} must be an object.')
            extras = set(section) - {f.name for f in fields(kind)}
            if extras:
                raise ModelError(f'Unknown {name} inputs: {sorted(extras)}')
            section = dict(section)
            for key in ('coverage_schedule', 'comparator_schedule', 'pep_uptake_schedule', 'comparator_pep_uptake_schedule'):
                if key in section:
                    if not isinstance(section[key], (list, tuple)):
                        raise ModelError(f'{key} must be an array.')
                    section[key] = tuple(section[key])
            result[name] = kind(**section)
        config = cls(**result)
        config.validate()
        return config

    @classmethod
    def from_json(cls, text):
        try:
            return cls.from_dict(json.loads(text))
        except (json.JSONDecodeError, TypeError) as exc:
            raise ModelError(f'Invalid input JSON: {exc}') from exc

    def validate(self):
        for section, values in self.to_dict().items():
            for name, value in values.items():
                default = getattr(getattr(Config(), section), name)
                if isinstance(default, bool):
                    if type(value) is not bool:
                        raise ModelError(f'{section}.{name} must be true or false.')
                elif isinstance(default, str):
                    if not isinstance(value, str) or not value.strip():
                        raise ModelError(f'{section}.{name} must be nonempty text.')
                elif isinstance(default, tuple):
                    if not isinstance(value, (tuple, list)):
                        raise ModelError(f'{name} must be a sequence.')
                    if any(isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x) or not 0 <= x <= 1 for x in value):
                        raise ModelError(f'{name} must contain proportions between 0 and 1.')
                else:
                    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                        raise ModelError(f'{section}.{name} must be finite and nonnegative.')
                    if isinstance(default, int) and type(value) is not int:
                        raise ModelError(f'{section}.{name} must be an integer.')
        s, t, v, e = self.setting, self.transmission, self.vaccination, self.economics
        for value, name in [(s.human_population, 'Human population'), (s.dog_population, 'Dog population'),
                            (s.area_km2, 'Area'), (s.initial_capacity_ratio, 'Initial capacity ratio'),
                            (t.incubation_days, 'Incubation'), (t.days_per_model_week, 'Days/model week'),
                            (t.dog_clinical_probability, 'Dog clinical probability'),
                            (t.dog_infectious_exit_per_week, 'Dog infectious exit rate'),
                            (t.human_life_expectancy_years, 'Human life expectancy'),
                            (e.reference_campaign_doses, 'Reference campaign doses')]:
            if value <= 0:
                raise ModelError(f'{name} must be positive.')
        if not 1 <= s.report_first_year <= s.report_last_year <= s.years <= 100 or s.lead_in_years > 100:
            raise ModelError('Choose 1–100 program years and a reporting window within them; lead-in maximum is 100 years.')
        if not 0 <= s.burn_in_weeks <= 100000:
            raise ModelError('Burn-in must be 0–100,000 weeks.')
        if not 1 <= s.stabilization_years <= 100:
            raise ModelError('Stabilization lead-in must be 1–100 years.')
        if not 0 <= s.stabilization_coverage <= 1:
            raise ModelError('Stabilization coverage must be between 0 and 1.')
        if s.initial_mode not in ('burn_in', 'explicit'):
            raise ModelError('Initial mode must be burn_in or explicit.')
        if s.initial_exposed_fraction + s.initial_infectious_fraction + s.initial_immune_fraction > 1:
            raise ModelError('Initial exposed, infectious and immune fractions must sum to at most 1.')
        if s.initial_mode == 'burn_in' and (s.initial_immune_fraction or s.initial_exposed_fraction + s.initial_infectious_fraction == 0):
            raise ModelError('Burn-in requires an infection seed and zero initial immunity; use explicit mode for known state fractions.')
        if t.dog_births_per_year < t.dog_natural_mortality_per_year:
            raise ModelError('This density-dependent model requires dog births ≥ baseline dog mortality.')
        for value, name in [(s.historical_coverage, 'Historical coverage'), (s.historical_pep_uptake, 'Historical PEP'),
                            (t.dog_clinical_probability, 'Clinical probability'), (t.pep_efficacy, 'PEP efficacy'),
                            (v.efficacy, 'Vaccine efficacy'), (v.exposed_protection, 'Protection in exposed dogs'), (v.coverage, 'Coverage'), (v.comparator_coverage, 'Comparator coverage'),
                            (v.pep_uptake, 'PEP uptake'), (v.comparator_pep_uptake, 'Comparator PEP'),
                            (v.maintenance_coverage, 'Maintenance coverage'), (e.rig_fraction, 'RIG fraction'),
                            (e.quarantine_probability, 'Quarantine probability'), (e.testing_probability, 'Testing probability'),
                            (e.investigation_probability, 'Investigation probability'), (e.yll_age_weight, 'YLL age weight')]:
            if not 0 <= value <= 1:
                raise ModelError(f'{name} must be between 0 and 1.')
        if not 1 <= v.revaccination_interval_weeks <= 260 or not 0 <= v.protected_weeks <= v.revaccination_interval_weeks:
            raise ModelError('Revaccination interval must be 1–260 weeks; protected duration must be between 0 and that interval.')
        if not 1 <= v.first_campaign_week <= 52:
            raise ModelError('Campaign first week must be 1–52.')
        for duration in (v.campaign_weeks, v.comparator_campaign_weeks, s.historical_campaign_weeks):
            if not 1 <= duration <= 53 - v.first_campaign_week:
                raise ModelError('Every campaign must fit within the 52-week year.')
        for schedule in (v.coverage_schedule, v.comparator_schedule, v.pep_uptake_schedule, v.comparator_pep_uptake_schedule):
            if schedule and len(schedule) != s.years:
                raise ModelError('A coverage schedule must contain exactly one value per program year.')
        if v.maintenance_years > 100 or v.elimination_threshold_cases <= 0:
            raise ModelError('Maintenance must be 0–100 years and the elimination threshold positive.')
        if e.average_age_at_death > e.yll_life_expectancy:
            raise ModelError('YLL life expectancy must be at least average age at death.')
        if t.days_per_model_week / t.incubation_days + t.dog_births_per_year / 52 > 1:
            raise ModelError('Incubation progression is too fast for this weekly solver.')
        if t.human_clinical_rate_per_week + t.human_nonclinical_rate_per_week + 1 / t.human_life_expectancy_years / 52 > 1:
            raise ModelError('Human exposed-state exit rates exceed one per week.')
        if max(v.waning_per_week, t.dog_infectious_exit_per_week, t.human_infectious_exit_per_week) > 1:
            raise ModelError('Weekly waning and infectious exit rates must be at most one.')
        return self
