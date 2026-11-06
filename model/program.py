"""GDREP V2.1.25.bangalore: program inputs around the frozen V2.1.25years engine.

Only schedules, initialization scale, observation, and reporting are adapted.
All disease, vaccination, human, PEP-demand, attribution, and cost equations
execute the original source files.
"""
from pathlib import Path
from functools import lru_cache
from datetime import datetime, timezone
import copy, hashlib, json, math, sys, types
import numpy as np
import pandas as pd
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parent
ENGINE=ROOT/'GDREP_V2_1_25years'
sys.path.insert(0,str(ENGINE))
from rabies_econ import Config
from rabies_econ.config import ModelError
from rabies_econ.dog import DogModel
from rabies_econ import engine

VERSION='V2.1.25.bangalore'
DEFAULTS={
 'name':'Bangalore India','scenario_name':'Final Bangalore','area_km2':712.0,
 'human_population':13000000.0,'dog_population':309000.0,
 'dog_input_mode':'count','humans_per_dog':42.0,'lead_in_years':25,
 'history_mode':'count','history_value':10000.0,'history_weeks':52,
 'history_pep':0.75,'r0':2.0,'r0_cap':2.0,'z':0.00002,
 'dog_recruitment':0.55,'dog_mortality':0.33,'incubation_days':40.0,
 'human_birth_rate':15.0,'life_expectancy':69.8,'pep_efficacy':0.93,
 'pep_multiplier':10.0,'program_weeks':52,'comparator_weeks':52,'start_week':1,
 'program_mode':'count','comparator_mode':'count',
 'program_values':[112537.,105581.,104123.,89211.,155901.],
 'comparator_values':[10000.]*5,'program_pep':[0.75]*5,'comparator_pep':[0.75]*5,
 'cdr_before':0.01,'cdr_after':0.10,'dog_vaccine_efficacy':0.95,
 'exposed_protection':0.50,'protected_weeks':26,'revaccination_weeks':52,
 'waning_per_week':1/52,'dog_clinical_probability':0.45,'dog_exit_rate':0.70,
 'dog_cost':3.21,'pep_material':1.06,'pep_overhead':5.0,'pep_vaccine':3.0,
 'pep_doses_per_visit':1.0,'pep_visits':5.0,'rig_cost':15.94,'rig_fraction':0.60,
 'quarantine_probability':0.0008,'quarantine_cost':74.38,
 'test_probability':0.011333333,'test_cost':3.19,
 'investigation_probability':0.466666667,'investigation_cost':15.94,
 'other_program_cost':0.0,'other_comparator_cost':0.0,
 'cost_discount':0.03,'outcome_discount':0.03,'average_age_at_death':15.0,
 'currency':'USD','price_year':'Supplied costs','mortality_value':None,
 'observed_annual_pep':None
}

class ScheduleError(ModelError): pass

def normalize(inputs=None):
    p=copy.deepcopy(DEFAULTS)
    if inputs:
        unknown=set(inputs)-set(p)
        if unknown: raise ModelError('Unknown input fields: '+', '.join(sorted(unknown)))
        p.update(copy.deepcopy(inputs))
    for key,default in DEFAULTS.items():
        val=p[key]
        if isinstance(default,str):
            if not isinstance(val,str) or not val.strip(): raise ModelError(key+' must be nonempty text')
            if len(val)>200: raise ModelError(key+' is too long')
        elif isinstance(default,list):
            if not isinstance(val,list): raise ModelError(key+' must be a list')
            for v in val:
                if isinstance(v,bool) or not isinstance(v,(float,int)) or not math.isfinite(v) or v<0: raise ModelError(key+' must contain finite nonnegative numbers')
            p[key]=[float(v) for v in val]
        elif default is None:
            if val is not None and (isinstance(val,bool) or not isinstance(val,(float,int)) or not math.isfinite(val) or val<0): raise ModelError(key+' must be empty or nonnegative')
            if val is not None: p[key]=float(val)
        else:
            if isinstance(val,bool) or not isinstance(val,(int,float)) or not math.isfinite(val) or val<0: raise ModelError(key+' must be a finite nonnegative number')
            if isinstance(default,int):
                if int(val)!=val: raise ModelError(key+' must be a whole number')
                p[key]=int(val)
            else: p[key]=float(val)
    n=len(p['program_values'])
    if not 1<=n<=30: raise ModelError('Choose 1–30 program years in the app')
    for key in ['comparator_values','program_pep','comparator_pep']:
        if len(p[key])!=n: raise ModelError('All annual schedules must have the same length')
    for key in ['history_mode','program_mode','comparator_mode']:
        if p[key] not in ['count','coverage']: raise ModelError(key+' must be count or coverage')
    if p['dog_input_mode'] not in ['count','ratio']: raise ModelError('Dog input mode must be count or ratio')
    if p['dog_input_mode']=='ratio':
        if p['humans_per_dog']<=0: raise ModelError('Humans per dog must be positive')
        p['dog_population']=p['human_population']/p['humans_per_dog']
    for key in ['cdr_before','cdr_after']:
        if p[key]>1: raise ModelError(key+' must be between zero and one')
    for key in ['program_pep','comparator_pep']:
        if any(v>1 for v in p[key]): raise ModelError(key+' must contain proportions between zero and one')
    if p['history_mode']=='coverage' and p['history_value']>1: raise ModelError('Historical coverage must be between zero and one')
    for arm in ['program','comparator']:
        if p[arm+'_mode']=='coverage' and any(v>1 for v in p[arm+'_values']): raise ModelError(arm+' coverage must be between zero and one')
    if p['r0']<=0 or p['r0_cap']<=0: raise ModelError('R0 and its cap must be positive')
    if p['pep_multiplier']<1: raise ModelError('PEP demand multiplier must be at least 1')
    return p

def configure(p):
    cfg=json.loads((ROOT/'base_config.json').read_text()); n=len(p['program_values'])
    cfg['setting'].update(name=p['name'],human_population=p['human_population'],dog_population=p['dog_population'],area_km2=p['area_km2'],years=n,report_first_year=1,report_last_year=n,lead_in_years=p['lead_in_years'],historical_coverage=0.,historical_pep_uptake=p['history_pep'],historical_campaign_weeks=p['history_weeks'],estimate_at_stabilization=False)
    cfg['transmission'].update(actual_r0=min(p['r0'],p['r0_cap']),z=p['z'],dog_births_per_year=p['dog_recruitment'],dog_natural_mortality_per_year=p['dog_mortality'],incubation_days=p['incubation_days'],dog_clinical_probability=p['dog_clinical_probability'],dog_infectious_exit_per_week=p['dog_exit_rate'],human_births_per_1000_year=p['human_birth_rate'],human_life_expectancy_years=p['life_expectancy'],pep_efficacy=p['pep_efficacy'],suspected_exposure_multiplier=p['pep_multiplier'])
    cfg['vaccination'].update(first_campaign_week=p['start_week'],campaign_weeks=p['program_weeks'],comparator_campaign_weeks=p['comparator_weeks'],coverage_schedule=[],comparator_schedule=[],coverage=0.,comparator_coverage=0.,pep_uptake=p['program_pep'][0],comparator_pep_uptake=p['comparator_pep'][0],pep_uptake_schedule=p['program_pep'],comparator_pep_uptake_schedule=p['comparator_pep'],efficacy=p['dog_vaccine_efficacy'],exposed_protection=p['exposed_protection'],protected_weeks=p['protected_weeks'],revaccination_interval_weeks=p['revaccination_weeks'],waning_per_week=p['waning_per_week'])
    cfg['economics'].update(currency=p['currency'],cost_discount_rate=p['cost_discount'],outcome_discount_rate=p['outcome_discount'],vaccine_cost_per_dose=p['dog_cost'],dog_materials_cost_per_dose=0.,dog_marking_cost_per_dose=0.,vaccine_wastage_markup=0.,campaign_fixed_cost=0.,reference_campaign_doses=1.,pep_material_cost_per_injection=p['pep_material'],pep_overhead_cost=p['pep_overhead'],pep_vaccine_cost_per_dose=p['pep_vaccine'],pep_doses_per_visit=p['pep_doses_per_visit'],pep_visits=p['pep_visits'],rig_cost=p['rig_cost'],rig_fraction=p['rig_fraction'],quarantine_probability=p['quarantine_probability'],quarantine_unit_cost=p['quarantine_cost'],testing_probability=p['test_probability'],testing_unit_cost=p['test_cost'],investigation_probability=p['investigation_probability'],investigation_unit_cost=p['investigation_cost'],intervention_other_annual_cost=p['other_program_cost'],comparator_other_annual_cost=p['other_comparator_cost'],average_age_at_death=p['average_age_at_death'],yll_life_expectancy=p['life_expectancy'])
    Config.from_dict(cfg)
    return cfg

def schedule(model,initial,unit,values,mode,weeks,record=False):
    state=initial.copy(); covers=[]; parts=[]; states=[state.sum(axis=0)*unit] if record else []
    step=model.step
    if record:
        def observer(s,budget=0.):
            out,flows=step(s,budget);states.append(out.sum(axis=0)*unit);return out,flows
        model.step=observer
    try:
        for value in values:
            coverage=value/(state.sum()*unit) if mode=='count' else value
            if not 0<=coverage<=1: raise ScheduleError('A vaccination count exceeds the opening dog population. Reduce counts or use coverage inputs.')
            one=model.run(state,[coverage],weeks)
            if mode=='count' and abs(one['annual'][0]['dog_doses']*unit-value)>max(1e-6,value*1e-10):
                raise ScheduleError('The requested dog vaccinations cannot be delivered under the model’s eligibility rules. Lengthen delivery or reduce the count.')
            covers.append(float(coverage));parts.append(one);state=one['state']
    finally: model.step=step
    checks={k:(min if 'min_' in k else max)(q['checks'][k] for q in parts) for k in parts[0]['checks']} if parts else {'min_normalized_state':0.,'max_normalized_balance_error':0.,'max_normalized_infectious_balance_error':0.}
    return covers,{'state':state.copy(),'annual':[q['annual'][0] for q in parts],'trajectory':np.concatenate([q['trajectory'] for q in parts]) if parts else np.array([]),'checks':checks},states

@lru_cache(maxsize=8)
def _history(key):
    data=json.loads(key);cfg=data['config'];p=data['input'];m=DogModel(Config.from_dict(cfg));initial,notes=m.initial_state()
    vals=[p['history_value']]*p['lead_in_years'];target=p['dog_population']
    if not vals: unit=target/initial.sum()
    elif p['history_mode']=='coverage' or p['history_value']==0:
        _,lead,_=schedule(m,initial,1.,vals,'coverage',p['history_weeks']);unit=target/lead['state'].sum()
    else:
        def residual(u): return float(schedule(m,initial,u,vals,'count',p['history_weeks'])[1]['state'].sum()*u-target)
        low=target*200000/309000;high=target*1000000/309000
        try:
            fl,fh=residual(low),residual(high)
            if fl*fh>0: raise ScheduleError('Find another bracket')
        except ScheduleError:
            bracket=[]
            for u in target*np.geomspace(.05,100,32):
                try: f=residual(u)
                except ScheduleError: continue
                if bracket and f*bracket[-1][1]<=0:
                    low,high=bracket[-1][0],u;break
                bracket.append((u,f))
            else: raise ScheduleError('Unable to anchor this historical count schedule at the requested starting dog population. Review baseline vaccinations, lead-in and recruitment, or use historical coverage.')
        unit=brentq(residual,low,high,xtol=1e-7,rtol=1e-13)
    historical,lead,_=schedule(m,initial,unit,vals,p['history_mode'],p['history_weeks'])
    return initial,notes,unit,historical,lead

def history_key(p,cfg):
    # Reuse dog history when only human, PEP, cost or reporting inputs change.
    base=Config().to_dict();base['setting'].update({k:cfg['setting'][k] for k in ['dog_population','lead_in_years','historical_campaign_weeks']})
    base['vaccination'].update({k:cfg['vaccination'][k] for k in ['first_campaign_week','efficacy','exposed_protection','protected_weeks','revaccination_interval_weeks','waning_per_week']})
    base['setting']['historical_campaign_weeks']=p['history_weeks']
    base['vaccination']['campaign_weeks']=p['history_weeks'];base['vaccination']['comparator_campaign_weeks']=p['history_weeks']
    for k in ['actual_r0','dog_births_per_year','dog_natural_mortality_per_year','incubation_days','days_per_model_week','dog_clinical_probability','dog_infectious_exit_per_week']:base['transmission'][k]=cfg['transmission'][k]
    # Preserve the frozen initial-state inputs instead of introducing new defaults.
    for k in ['initial_mode','burn_in_weeks','initial_capacity_ratio','initial_exposed_fraction','initial_infectious_fraction','initial_immune_fraction']:base['setting'][k]=cfg['setting'][k]
    return json.dumps({'config':base,'input':{k:p[k] for k in ['history_value','history_mode','history_weeks','lead_in_years','dog_population']}},sort_keys=True)

def records(df): return json.loads(df.to_json(orient='records',double_precision=15))

def run_analysis(inputs=None,progress=None):
    emit=progress or (lambda pct,msg:None)
    emit(3,'Preparing model inputs');p=normalize(inputs);cfg=configure(p)
    emit(10,'Preparing dog population and vaccination history')
    initial,notes,unit,historical,lead=_history(history_key(p,cfg));m=DogModel(Config.from_dict(cfg))
    emit(38,'Recalculating annual coverage and rolling immunity')
    covers,actual,weekly=schedule(m,lead['state'],unit,p['program_values'],p['program_mode'],p['program_weeks'],True)
    ccovers,comparison,cweekly=schedule(m,lead['state'],unit,p['comparator_values'],p['comparator_mode'],p['comparator_weeks'],True)
    cfg['vaccination']['coverage_schedule']=covers;cfg['vaccination']['comparator_schedule']=ccovers;c=Config.from_dict(cfg)
    emit(55,'Running human outcomes and matched counterfactuals')
    namespace=dict(engine.run.__globals__);namespace['prepare_dogs']=lambda _: (m,initial,lead,unit,notes)
    runner=types.FunctionType(engine.run.__code__,namespace,engine.run.__name__,engine.run.__defaults__,engine.run.__closure__)
    result=runner(c);a=result.annual.copy();n=len(a);L=p['lead_in_years']
    a['coverage']=covers;a['comparator_coverage']=ccovers
    a['pep_true_exposures']=a.human_exposures*np.array(p['program_pep']);a['pep_non_exposures']=a.pep_courses-a.pep_true_exposures
    a['expected_detected_cases']=a.human_deaths*p['cdr_after'];a['comparator_expected_detected_cases']=a.comparator_human_deaths*p['cdr_after']
    a['program_dog_cases_averted']=a.comparator_clinical_dog_cases-a.clinical_dog_cases
    for col in ['dog_doses','opening_dogs','closing_dogs','clinical_dog_cases']:
        np.testing.assert_allclose(a[col],np.array([r[col] for r in actual['annual']])*unit,rtol=1e-11,atol=1e-6)
    emit(78,'Assembling the full timeline and economic comparison')
    humans,_,_=engine.human_run(c,np.concatenate([lead['trajectory'],actual['trajectory']])*unit/p['area_km2'],[p['history_pep']]*L+p['program_pep'],result.derived['historical_start_humans'])
    full=pd.concat([pd.DataFrame(lead['annual']+actual['annual'])*unit,humans],axis=1)
    full.insert(0,'model_year',range(1,L+n+1));full['phase']=['Lead-in']*L+['Program']*n;full['coverage']=historical+covers
    full['expected_detected_cases']=full.human_deaths*np.array([p['cdr_before']]*L+[p['cdr_after']]*n)
    np.testing.assert_allclose(full.human_deaths.iloc[L:],a.human_deaths,rtol=1e-11,atol=1e-8)
    monthly=pd.DataFrame({'month':np.arange(n*12+1),'model_week':np.arange(n*12+1)*52/12})
    wframes=[]
    for prefix,states in [('',weekly),('comparator_',cweekly)]:
        arr=np.asarray(states);weeks=np.arange(len(arr));total=arr.sum(axis=1);immune=arr[:,3:].sum(axis=1)
        wf=pd.DataFrame(arr,columns=[prefix+k for k in ['susceptible_dogs','exposed_dogs','infectious_dogs','immune_protected_dogs','immune_waning_dogs']]);wf[prefix+'total_dogs']=total;wf[prefix+'immune_dogs']=immune;wf[prefix+'immune_percent']=100*immune/total;wframes.append(wf)
        monthly[prefix+'total_dogs']=np.interp(monthly.model_week,weeks,total);monthly[prefix+'immune_dogs']=np.interp(monthly.model_week,weeks,immune);monthly[prefix+'immune_percent']=100*monthly[prefix+'immune_dogs']/monthly[prefix+'total_dogs']
        expected=a[prefix+'closing_dogs'];np.testing.assert_allclose(total[52::52],expected,atol=1e-6,rtol=1e-11)
    wf=pd.concat(wframes,axis=1);wf.insert(0,'week',range(len(wf)))
    b=pd.DataFrame({'year':a.year});fc=(1+p['cost_discount'])**(-a.year);fo=(1+p['outcome_discount'])**(-a.year)
    b['additional_dog_cost']=a.dog_vaccination_cost-a.comparator_dog_vaccination_cost
    b['additional_other_cost']=a.other_cost-a.comparator_other_cost
    b['additional_program_cost']=b.additional_dog_cost+b.additional_other_cost
    b['pep_savings']=a.comparator_pep_cost-a.pep_cost;b['response_savings']=a.comparator_suspect_exposure_cost-a.suspect_exposure_cost
    b['direct_savings']=b.pep_savings+b.response_savings;b['net_additional_cost']=b.additional_program_cost-b.direct_savings;b['deaths_averted']=a.incremental_deaths_averted
    for k in list(b.columns)[1:]: b['pv_'+k]=b[k]*(fo if k=='deaths_averted' else fc)
    np.testing.assert_allclose(b.pv_net_additional_cost,a.discounted_incremental_cost,atol=1e-7,rtol=1e-10)
    C=float(b.pv_additional_program_cost.sum());S=float(b.pv_direct_savings.sum());N=float(b.pv_net_additional_cost.sum());D=float(a.discounted_deaths_averted.sum())
    # A constant monetary health value is discounted at the monetary rate.
    monetary_deaths=float((a.incremental_deaths_averted*fc).sum())
    V=p['mortality_value'];health=None if V is None else V*monetary_deaths
    econ={'additional_program_cost_pv':C,'direct_savings_pv':S,'pep_savings_pv':float(b.pv_pep_savings.sum()),'response_savings_pv':float(b.pv_response_savings.sum()),'net_additional_cost_pv':N,'direct_financial_bcr':S/C if C>0 else None,'net_financial_benefit_pv':-N,'mortality_value':V,'mortality_benefit_pv':health,'bcr_with_mortality':(S+health)/C if health is not None and C>0 else None,'net_benefit_with_mortality':S+health-C if health is not None else None,'mortality_breakeven':max(0.,N/monetary_deaths) if monetary_deaths>0 else None,'discounted_deaths_averted':D,'cost_per_death_averted':result.summary['cost_per_death_averted'],'cost_per_model_daly_averted':result.summary['cost_per_daly_averted'],'discounted_dalys_averted':float(a.discounted_dalys_averted.sum())}
    diagnostic_counts=[p['history_value']]*5
    try:
        _,diag,_=schedule(m,lead['state'],unit,diagnostic_counts,p['history_mode'],p['history_weeks'])
        dc=np.array([r['clinical_dog_cases'] for r in diag['annual']])*unit;dp=np.array([r['opening_dogs'] for r in diag['annual']])*unit
        stability={'case_relative_range':float(np.ptp(dc)/max(float(dc.max()),.5)),'population_relative_range':float(np.ptp(dp)/max(float(dp.max()),1.))};stability['passed']=max(stability.values())<=.01
    except ScheduleError as ex: stability={'passed':False,'message':str(ex)}
    summary={k:float(a[k].sum()) for k in ['human_deaths','comparator_human_deaths','incremental_deaths_averted','clinical_dog_cases','program_dog_cases_averted','pep_courses','pep_true_exposures','pep_non_exposures','expected_detected_cases','dog_doses','dog_vaccination_cost','pep_cost','suspect_exposure_cost','total_cost','comparator_total_cost','incremental_cost']}
    summary.update(years=n,lead_in_years=L,mean_annual_pep=summary['pep_courses']/n,year_end_immunity_percent=float(monthly.immune_percent.iloc[-1]),actual_r0=min(p['r0'],p['r0_cap']),z=p['z'])
    warnings=[]
    if not stability['passed']:warnings.append('The requested lead-in has not reached the 1% stability criterion. Baseline oscillations are retained.')
    if min(p['program_pep'])!=max(p['program_pep']) or p['program_pep']!=p['comparator_pep']:warnings.append('PEP access differs across time or arms; incremental benefits include these PEP changes as well as dog vaccination.')
    if float(a.unused_dog_doses.sum())>1e-6:warnings.append('Some planned vaccinations were not delivered under eligibility rules. Model costs use planned doses; delivered doses are reported separately.')
    if float(a.clinical_dog_cases.iloc[-1])>=.5:warnings.append('The final program year remains above 0.5 clinical dog cases; this run does not demonstrate elimination.')
    if C<=0:warnings.append('Incremental program spending is zero or negative. A benefit–cost ratio is not reported; use net benefit.')
    if summary['incremental_deaths_averted']<=0:warnings.append('This scenario does not avert human deaths relative to its comparator. Review the full trajectories and inputs.')
    run_id=hashlib.sha256(json.dumps(p,sort_keys=True).encode()).hexdigest()[:12]
    out={'version':VERSION,'run_id':run_id,'created_utc':datetime.now(timezone.utc).isoformat(),'inputs':p,'engine_config':c.to_dict(),'summary':summary,'economics':econ,'annual':records(a),'full_timeline':records(full),'monthly_immunity':records(monthly),'weekly_immunity':records(wf),'benefit_cost':records(b),'derived':result.derived,'warnings':warnings,'stability':stability}
    emit(100,'Analysis complete')
    return out
