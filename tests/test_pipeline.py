"""Step パイプラインの検証 (unittest)."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from mj_solar_uncertainty.pipeline import run_uncertainty_pipeline


class TestUncertaintyPipeline(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data_dir = os.path.join(ROOT, "sample_data")
        cls._tmpdir = tempfile.mkdtemp()
        cls.result = run_uncertainty_pipeline(
            cls.data_dir,
            cls._tmpdir,
            n_mc_samples=3000,
            mc_seed=42,
            save_intermediate_csv=True,
        )

    def test_eq2_residual_from_saved_csv(self):
        M = pd.read_csv(os.path.join(self._tmpdir, "step02_M.csv"), index_col=0)
        A = pd.read_csv(os.path.join(self._tmpdir, "step02_A.csv"), index_col=0)["A_i"].values
        b = pd.read_csv(os.path.join(self._tmpdir, "step02_b.csv"), index_col=0)["b_j"].values
        Mv, bv = M.values, b
        n_j, n_i = Mv.shape
        if n_i == n_j:
            A_expect = np.linalg.solve(Mv, bv)
            np.testing.assert_allclose(A, A_expect, rtol=1e-10, atol=1e-10)
            r = Mv @ A - bv
            self.assertLess(np.max(np.abs(r)), 1e-9)
        else:
            A_expect, *_ = np.linalg.lstsq(Mv, bv, rcond=None)
            np.testing.assert_allclose(A, A_expect, rtol=1e-10, atol=1e-10)
            self.assertEqual(M.shape[1], 2)
            self.assertEqual(len(A), 2)

    def test_step_meta_residual_matches(self):
        st2 = next(s for s in self.result.steps if s.id == 2)
        M = pd.read_csv(os.path.join(self._tmpdir, "step02_M.csv"), index_col=0)
        A = pd.read_csv(os.path.join(self._tmpdir, "step02_A.csv"), index_col=0)["A_i"].values
        b = pd.read_csv(os.path.join(self._tmpdir, "step02_b.csv"), index_col=0)["b_j"].values
        r_max = float(np.max(np.abs(M.values @ A - b)))
        self.assertAlmostEqual(st2.meta["max_abs_residual_MA_minus_b"], r_max, places=10)
        expect_mode = "solve" if M.shape[0] == M.shape[1] else "lstsq"
        self.assertEqual(st2.meta["solve_mode"], expect_mode)

    def test_table1_rss_matches_components(self):
        t1 = self.result.table1
        params = list(t1.columns)
        extras = self.result.extra_uncertainties_relative
        extra_sq = sum(v**2 for v in extras.values())
        for p in params:
            ut = float(t1.loc["u_top [%]", p])
            um = float(t1.loc["u_mid [%]", p])
            ub = float(t1.loc["u_bot [%]", p])
            uc = float(t1.loc["u_combined [%]", p])
            U = float(t1.loc["U (k=2) [%]", p])
            rss = np.sqrt(ut**2 + um**2 + ub**2 + extra_sq)
            self.assertAlmostEqual(rss, uc, places=2)
            self.assertAlmostEqual(2.0 * uc, U, places=2)

    def test_mc_systematic_ge_independent(self):
        r = self.result
        for i in range(3):
            self.assertGreaterEqual(
                r.mc_systematic.Jphot_subcell_rel[i],
                r.mc_independent.Jphot_subcell_rel[i] - 1e-9,
            )

    def test_perturb_all_matches_main_mc_column(self):
        mc_main = pd.read_csv(os.path.join(self._tmpdir, "mc_subcell_currents.csv"))
        br = pd.read_csv(
            os.path.join(self._tmpdir, "step04_mc_perturb_breakdown_systematic.csv")
        )
        row_all = br.loc[br["perturb"] == "all"].iloc[0]
        for sub, col in [("top", "u_Jphot_top_pct"), ("mid", "u_Jphot_mid_pct"), ("bot", "u_Jphot_bot_pct")]:
            v_main = float(mc_main.loc[mc_main["subcell"] == sub, "u_Jphot_systematic_[%]"].values[0])
            v_br = float(row_all[col])
            self.assertAlmostEqual(v_main, v_br, places=3)


if __name__ == "__main__":
    unittest.main()
