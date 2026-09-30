# Author: Vangelis D. Karalis
"""Average BE under a complete, balanced, homoscedastic crossover model."""
from dataclasses import dataclass
import numpy as np
from scipy.stats import t, norm

DESIGNS = {
    '2x2': ('TR', 'RT'),
    '2x3': ('RTR', 'TRT'),
    '2x4': ('TRTR', 'RTRT'),
    '3x3': ('TRR', 'RTR', 'RRT'),
}


def design_matrix(design, subjects):
    """Subject intercepts absorb sequence; add period and treatment effects."""
    if design not in DESIGNS:
        raise ValueError('Unknown design. Choose 2x2, 2x3, 2x4 or 3x3.')
    sequences = DESIGNS[design]
    if isinstance(subjects, bool) or not isinstance(subjects, (int, np.integer)):
        raise ValueError('Subjects must be an integer.')
    if subjects < 2 * len(sequences) or subjects % len(sequences):
        raise ValueError(f'Subjects must be at least {2*len(sequences)} and divisible by {len(sequences)}.')
    periods = len(sequences[0])
    subject = np.repeat(np.arange(subjects), periods)
    period = np.tile(np.arange(periods), subjects)
    treatment = np.array([c == 'T' for s in sequences for _ in range(subjects // len(sequences)) for c in s], dtype=float)
    x = np.column_stack((np.eye(subjects)[subject], np.eye(periods)[period, 1:], treatment))
    return x


def model_geometry(design, subjects):
    x = design_matrix(design, subjects)
    rank = np.linalg.matrix_rank(x)
    if rank != x.shape[1]:
        raise ValueError('Treatment effect is not estimable.')
    return int(x.shape[0] - rank), float(np.linalg.inv(x.T @ x)[-1, -1])


def gmr_grid(spec):
    """Parse a single GMR or start:end:step (inclusive when on the grid)."""
    parts = [float(v.strip()) for v in spec.split(':')]
    if len(parts) == 1:
        values = np.array(parts)
    elif len(parts) == 3:
        start, end, step = parts
        if not np.all(np.isfinite(parts)) or start <= 0 or end < start or step <= 0:
            raise ValueError('Use positive start:end:step, with end >= start.')
        count = int(np.floor((end-start)/step + 1e-10)) + 1
        if count > 10000:
            raise ValueError('At most 10000 GMR points are supported.')
        values = start + step * np.arange(count)
    else:
        raise ValueError('Use a single GMR or start:end:step.')
    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError('GMR must be finite and positive.')
    return values


@dataclass
class PowerResult:
    gmr: np.ndarray
    accepted: np.ndarray
    power_percent: np.ndarray
    mc_low_percent: np.ndarray
    mc_high_percent: np.ndarray
    df: int
    variance_factor: float


def simulate(design, subjects, trials, cv_percent, gmrs, seed=None, batch_size=100000):
    """Exact sampling distribution of fitted contrasts and residual MSE.

    Equivalent to generating normal log-PK records and fitting each study;
    uses independent normal contrast and chi-square residual sums of squares.
    Common random numbers across GMR points reduce curve noise.
    """
    df, factor = model_geometry(design, subjects)
    if isinstance(trials, bool) or not isinstance(trials, (int, np.integer)) or trials < 1:
        raise ValueError('Monte Carlo repetitions must be a positive integer.')
    if not np.isfinite(cv_percent) or cv_percent < 0:
        raise ValueError('Residual CV (%) must be finite and nonnegative.')
    if not isinstance(batch_size, int) or batch_size < 1:
        raise ValueError('Batch size must be a positive integer.')
    gmrs = np.atleast_1d(np.asarray(gmrs, dtype=float))
    if gmrs.ndim != 1 or gmrs.size == 0 or not np.all(np.isfinite(gmrs)) or np.any(gmrs <= 0):
        raise ValueError('GMR values must be a nonempty positive finite vector.')
    # Stable even for very large finite CV values.
    cv = cv_percent / 100
    sigma = np.sqrt(np.logaddexp(0, 2*np.log(cv))) if cv else 0.0
    rng = np.random.default_rng(seed)
    critical = t.ppf(.95, df)
    accepted = np.zeros(gmrs.size, dtype=np.int64)
    lower, upper = np.log(.8), np.log(1.25)
    for start in range(0, trials, batch_size):
        size = min(batch_size, trials-start)
        errors = rng.normal(0, sigma*np.sqrt(factor), size)
        half_width = critical * sigma * np.sqrt(factor * rng.chisquare(df, size) / df)
        for j, mu in enumerate(np.log(gmrs)):
            estimate = mu + errors
            accepted[j] += np.count_nonzero((estimate-half_width >= lower) & (estimate+half_width <= upper))
    p = accepted / trials
    z = norm.ppf(.975)
    center = (p + z*z/(2*trials))/(1+z*z/trials)
    radius = z*np.sqrt(p*(1-p)/trials + z*z/(4*trials*trials))/(1+z*z/trials)
    return PowerResult(gmrs, accepted, 100*p, 100*np.maximum(0, center-radius),
                       100*np.minimum(1, center+radius), df, factor)
