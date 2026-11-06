"""Weekly dog SEI and vaccination-age model, normalized to carrying capacity."""
import numpy as np
from .config import ModelError

class DogModel:
    def __init__(self, config):
        self.config = config
        t, v = config.transmission, config.vaccination
        self.b, self.d = t.dog_births_per_year / 52, t.dog_natural_mortality_per_year / 52
        self.sigma, self.r, self.mu = t.days_per_model_week / t.incubation_days, t.dog_clinical_probability, t.dog_infectious_exit_per_week
        self.beta_k = t.actual_r0 * (self.sigma + self.b) * (self.mu + self.b) / (self.sigma * self.r)
        self.interval, self.protected = v.revaccination_interval_weeks, v.protected_weeks
        self.efficacy, self.waning = v.efficacy, v.waning_per_week
        self.exposed_protection = v.exposed_protection

    def initial_state(self):
        s = self.config.setting
        n = 1 / s.initial_capacity_ratio
        x = n * np.array([1 - s.initial_exposed_fraction - s.initial_infectious_fraction - s.initial_immune_fraction,
                          s.initial_exposed_fraction, s.initial_infectious_fraction, s.initial_immune_fraction])
        residual = 0.0
        if s.initial_mode == 'burn_in':
            # Same unvaccinated SEI equations as the cohort solver, without empty age bins.
            for _ in range(s.burn_in_weeks):
                S, E, I, V = x
                N = x.sum(); mortality = self.d + (self.b - self.d) * N
                infection = self.beta_k * S * I
                change = np.array([self.b * N - mortality * S - infection + self.sigma * (1-self.r) * E,
                                   infection - (mortality+self.sigma)*E,
                                   self.sigma*self.r*E - (mortality+self.mu)*I, -mortality*V])
                x += change
                if not np.isfinite(x).all() or x.min() < -1e-13:
                    raise ModelError('Unstable burn-in: choose rates/population ratios compatible with weekly time steps.')
                residual = float(np.max(np.abs(change)))
            if x.sum() <= 0:
                raise ModelError('Burn-in produced an empty dog population.')
            x *= n / x.sum()
        state = np.zeros((self.interval + 1, 5))
        state[-1, :3] = x[:3]
        state[-1, 4] = x[3]  # Explicit initial immunity is eligible and already waning.
        return state, {'burn_in_final_week_change_normalized': residual,
                       'initial_fractions': (x/x.sum()).tolist(), 'initial_mode': s.initial_mode}

    def vaccinate(self, state, budget):
        eligible = state[-1].sum()
        doses = min(float(budget), float(eligible))
        selected = state[-1].copy() * (doses / eligible if eligible else 0.0)
        state[-1] -= selected
        state[0] += selected
        rescued = self.exposed_protection * selected[1]
        new = self.efficacy * selected[0] + rescued
        state[0, 1] -= rescued
        boosted = self.efficacy * selected[4]
        state[0, 0] -= self.efficacy * selected[0]
        if self.protected:
            state[0, 3] += new + boosted
            state[0, 4] -= boosted
        else:
            state[0, 4] += new
        return doses, new, selected[3] + selected[4]

    def step(self, state, budget=0.0):
        doses, new, redoses = self.vaccinate(state, budget)
        S, E, I, Vp, Vw = state.T
        totals = state.sum(axis=0); N = totals.sum()
        mortality = self.d + (self.b - self.d) * N
        force = self.beta_k * totals[2]
        # No clipping of high force-of-infection or renormalization of invalid runs.
        if mortality + max(force, self.sigma, self.mu, self.waning) > 1 + 1e-12:
            raise ModelError('Dog outflow probability exceeds one per week. Adjust inputs or use a finer-time-step model.')
        infection, waned = S*force, self.waning*Vw
        nxt = np.empty_like(state)
        nxt[:, 0] = S*(1-mortality)-infection+self.sigma*(1-self.r)*E+waned
        nxt[:, 1] = E*(1-mortality-self.sigma)+infection
        nxt[:, 2] = I*(1-mortality-self.mu)+self.sigma*self.r*E
        nxt[:, 3] = Vp*(1-mortality)
        nxt[:, 4] = Vw*(1-mortality-self.waning)
        nxt[-1, 0] += self.b*N
        tail = nxt[-1].copy()
        nxt[1:] = nxt[:-1].copy(); nxt[-1] += tail; nxt[0] = 0
        nxt[self.protected:, 4] += nxt[self.protected:, 3]
        nxt[self.protected:, 3] = 0
        flows = np.array([self.mu*totals[2], mortality*N+self.mu*totals[2], self.b*N,
                          self.sigma*self.r*totals[1], N/52, doses, new, redoses,
                          infection.sum(), mortality*totals[2]])
        return nxt, flows

    def run(self, initial, schedule, campaign_weeks):
        state = initial.copy()
        fields = ['dog_rabies_deaths', 'all_dog_deaths', 'dog_births', 'clinical_dog_cases',
                  'dog_years', 'dog_doses', 'new_dog_immunity', 'immune_redoses',
                  'dog_infections', 'nonrabies_infectious_dog_deaths']
        annual, trajectory = [], []
        max_balance = max_i_balance = 0.0
        minimum = 0.0
        for coverage in schedule:
            opening = state.sum(); opening_i = state[:, 2].sum()
            budget = float(coverage)*opening
            total = np.zeros(len(fields)); immune_years = 0.0
            for week in range(1, 53):
                trajectory.append(state[:, 2].sum())
                immune_years += state[:, 3:].sum()/52
                due = self.config.vaccination.first_campaign_week <= week < self.config.vaccination.first_campaign_week+campaign_weeks
                state, flows = self.step(state, budget/campaign_weeks if due else 0)
                total += flows
                minimum = min(minimum, float(state.min()))
            row = dict(zip(fields, total))
            row.update(opening_dogs=opening, closing_dogs=state.sum(), planned_dog_doses=budget,
                       unused_dog_doses=max(0.0, budget-total[5]), immune_dog_years=immune_years,
                       closing_immune_dogs=state[:, 3:].sum())
            max_balance = max(max_balance, abs(opening+total[2]-total[1]-state.sum()))
            max_i_balance = max(max_i_balance, abs(opening_i+total[3]-total[0]-total[9]-state[:, 2].sum()))
            if total[5] > budget+1e-12:
                raise ModelError('Administered doses exceed the annual budget.')
            annual.append(row)
        if not np.isfinite(state).all() or minimum < -1e-13 or max_balance > 1e-10 or max_i_balance > 1e-10:
            raise ModelError('Dog state or conservation check failed.')
        return {'annual': annual, 'state': state, 'trajectory': np.array(trajectory),
                'checks': {'min_normalized_state': minimum, 'max_normalized_balance_error': max_balance,
                           'max_normalized_infectious_balance_error': max_i_balance}}
