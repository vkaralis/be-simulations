# Author: Vangelis D. Karalis
import unittest
import numpy as np
from scipy.integrate import quad
from scipy.stats import chi2, norm, t
from besim.core import DESIGNS, design_matrix, model_geometry, simulate, gmr_grid


class ModelTests(unittest.TestCase):
    def test_geometry(self):
        # Independently derived balanced-design contrast variances and df.
        expected = {'2x2': (22, 2/24), '2x3': (45, 1.5/24),
                    '2x4': (68, 1/24), '3x3': (45, 1.5/24)}
        for design, (df, factor) in expected.items():
            actual_df, actual_factor = model_geometry(design, 24)
            self.assertEqual(actual_df, df)
            self.assertAlmostEqual(actual_factor, factor)

    def test_raw_study_fits(self):
        rng = np.random.default_rng(24)
        for design in DESIGNS:
            x = design_matrix(design, 24)
            df, factor = model_geometry(design, 24)
            sigma = np.sqrt(np.log1p(.3**2))
            beta = np.zeros(x.shape[1])
            beta[:24] = rng.normal(4.5, .5, 24)
            beta[-1] = np.log(.95)
            y = x @ beta[:, None] + sigma*rng.normal(size=(x.shape[0], 12000))
            fitted = np.linalg.lstsq(x, y, rcond=None)[0]
            residual = y-x@fitted
            mse = np.sum(residual**2, axis=0)/df
            self.assertAlmostEqual(np.mean(fitted[-1]), beta[-1], delta=.003)
            self.assertAlmostEqual(np.var(fitted[-1]), sigma**2*factor, delta=.00025)
            self.assertAlmostEqual(np.mean(mse), sigma**2, delta=.001)
            half = t.ppf(.95, df)*np.sqrt(mse*factor)
            raw_power = np.mean((fitted[-1]-half >= np.log(.8)) & (fitted[-1]+half <= np.log(1.25)))
            fast_power = simulate(design, 24, 50000, 30, [.95], seed=52).power_percent[0]/100
            self.assertAlmostEqual(raw_power, fast_power, delta=.025)

    def test_exact_power_and_boundary_type_one_error(self):
        for design in DESIGNS:
            df, factor = model_geometry(design, 24)
            sd = np.sqrt(np.log1p(.3**2)*factor)
            critical = t.ppf(.95, df)
            max_q = df*((np.log(1.25)-np.log(.8))/(2*critical*sd))**2
            gmrs = [.8, .95, 1, 1.25]
            result = simulate(design, 24, 100000, 30, gmrs, seed=42)
            for j, gmr in enumerate(gmrs):
                mu = np.log(gmr)
                def integrand(q):
                    h = critical*sd*np.sqrt(q/df)
                    probability = norm.cdf((np.log(1.25)-h-mu)/sd)-norm.cdf((np.log(.8)+h-mu)/sd)
                    return max(0, probability)*chi2.pdf(q, df)
                exact = quad(integrand, 0, max_q, epsabs=1e-9)[0]
                self.assertAlmostEqual(result.power_percent[j]/100, exact, delta=.007)
            self.assertLess(result.power_percent[0], 5.7)
            self.assertLess(result.power_percent[-1], 5.7)

    def test_validation_and_grid(self):
        np.testing.assert_allclose(gmr_grid('.8:1.25:.05'), np.arange(.8, 1.251, .05))
        np.testing.assert_array_equal(gmr_grid('1'), [1])
        for spec in ('0', '-1', 'nan', '1:0:.1', '1:2:0', '1:2'):
            with self.assertRaises(ValueError):
                gmr_grid(spec)
        for design in DESIGNS:
            with self.assertRaises(ValueError):
                model_geometry(design, 5)
        for trials, cv, gmrs in ((0, 30, [1]), (1.5, 30, [1]), (10, -1, [1]), (10, 30, [0])):
            with self.assertRaises(ValueError):
                simulate('2x2', 24, trials, cv, gmrs)

    def test_reproducibility_and_zero_variance(self):
        a = simulate('2x2', 24, 1000, 30, [.9, 1], seed=5)
        b = simulate('2x2', 24, 1000, 30, [.9, 1], seed=5)
        np.testing.assert_array_equal(a.accepted, b.accepted)
        zero = simulate('3x3', 24, 100, 0, [.7, .9, 1, 1.3], seed=1)
        np.testing.assert_array_equal(zero.power_percent, [0, 100, 100, 0])
        self.assertTrue(np.all(a.mc_low_percent <= a.power_percent))
        self.assertTrue(np.all(a.mc_high_percent >= a.power_percent))


if __name__ == '__main__':
    unittest.main()
