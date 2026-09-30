# Author: Vangelis D. Karalis
"""Interactive and command-line entry points."""
import argparse
import csv
import json
from pathlib import Path
from .core import DESIGNS, gmr_grid, simulate


def main(default_design=None):
    parser = argparse.ArgumentParser(description='Monte Carlo bioequivalence power curves')
    parser.add_argument('--design', choices=DESIGNS, default=default_design)
    parser.add_argument('--subjects', type=int)
    parser.add_argument('--trials', type=int)
    parser.add_argument('--cv', type=float, help='Residual within-subject CV in percent, e.g. 30')
    parser.add_argument('--gmr', help='Single value or start:end:step, e.g. 0.70:1.40:0.01')
    parser.add_argument('--seed', type=int, default=None)
    parser.add_argument('--output', type=Path, default=Path('results'))
    parser.add_argument('--show', action='store_true')
    args = parser.parse_args()
    try:
        if args.design is None:
            args.design = input('Σχεδιασμός (2x2, 2x3, 2x4, 3x3): ').strip()
        if args.design not in DESIGNS:
            raise ValueError('Unknown design.')
        if args.subjects is None:
            args.subjects = int(input('Συνολικός αριθμός ατόμων: '))
        if args.trials is None:
            args.trials = int(input('Αριθμός επαναλήψεων Monte Carlo ανά GMR: '))
        if args.cv is None:
            args.cv = float(input('Residual / ενδοατομικό CV σε % (π.χ. 30): '))
        if args.gmr is None:
            args.gmr = input('GMR T/R: τιμή ή αρχή:τέλος:βήμα (π.χ. 0.70:1.40:0.01): ')
        # Store the generated seed so even interactive runs can be reproduced.
        if args.seed is None:
            import secrets
            args.seed = secrets.randbits(64)
        result = simulate(args.design, args.subjects, args.trials, args.cv, gmr_grid(args.gmr), args.seed)
    except (ValueError, OverflowError) as error:
        parser.error(str(error))
    args.output.mkdir(parents=True, exist_ok=True)
    stem = args.output / f'power_{args.design}'
    with stem.with_suffix('.csv').open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['gmr', 'accepted', 'power_percent', 'mc95_low_percent', 'mc95_high_percent'])
        writer.writerows(zip(result.gmr, result.accepted, result.power_percent, result.mc_low_percent, result.mc_high_percent))
    metadata = dict(author="Vangelis D. Karalis", design=args.design, sequences=DESIGNS[args.design], subjects=args.subjects,
                    trials=args.trials, residual_cv_percent=args.cv, seed=args.seed,
                    gmrs=result.gmr.tolist(), residual_df=result.df, variance_factor=result.variance_factor,
                    be_limits=[.8, 1.25], confidence_level=.90,
                    model='subject + period + treatment; equal residual variance; no carryover',
                    monte_carlo='exact sufficient statistics; common random numbers across GMR')
    stem.with_suffix('.json').write_text(json.dumps(metadata, indent=2), encoding='utf-8')
    import matplotlib
    if not args.show:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(result.gmr, result.power_percent, 'o-', markersize=3, label='BE acceptance')
    ax.fill_between(result.gmr, result.mc_low_percent, result.mc_high_percent, alpha=.2, label='95% Monte Carlo interval (pointwise)')
    ax.set(xlabel='True GMR (T/R)', ylabel='BE acceptance probability (%)', ylim=(0, 100),
           title=f'{args.design}: N={args.subjects}, CV={args.cv:g}%, MC={args.trials:,}')
    for limit in (.8, 1.25):
        ax.axvline(limit, color='grey', linestyle='--', linewidth=1)
    ax.grid(alpha=.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(stem.with_suffix('.png'), dpi=180,
                metadata={'Software': 'BEsim', 'Author': 'Vangelis D. Karalis'})
    if args.show:
        plt.show()
    plt.close(fig)
    for gmr, power in zip(result.gmr, result.power_percent):
        print(f'GMR={gmr:.4f}  BE acceptance={power:.2f}%')
    print(f'Αποθηκεύτηκαν CSV, PNG, JSON: {stem}.* (seed={args.seed})')
