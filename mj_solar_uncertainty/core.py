"""
多接合太陽電池ソーラシミュレータ不確かさ評価コア
Baur & Bett (Fraunhofer ISE, 2005, 31st IEEE PVSC) 手法に準拠

主要式:
  Eq.(2): Σ_i [ A_i ∫ s_j(λ) e_i(λ) dλ ] = ∫ s_j(λ) E_ref(λ) dλ
  Eq.(3): J_ref,i = A_i ∫ S_ref(λ) e_i(λ) dλ
  Eq.(4): u(J)    = ∫ u(s(λ)) e(λ) dλ + ∫ s(λ) u(e(λ)) dλ
  Eq.(5): u(A_i)  = √ Σ_j (∂A_i/∂J_j)² u(J_j)²
  Eq.(6): u(Y_i)  = √ { (1/2Δ) ∫_{1-Δ}^{1+Δ} (Y_i(x) - Y_i(1))² dx }
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional


# =============================================================================
# 1. 波長依存データの基本クラス
# =============================================================================
@dataclass
class SpectralCurve:
    """波長依存量とその波長依存不確かさを保持する基本コンテナ.

    属性
    ----
    wavelength_nm : np.ndarray
        波長 [nm], 単調増加.
    value : np.ndarray
        相対値 (SR は [a.u.] あるいは [A/W], 放射照度は [W/m²/nm] 等).
    uncertainty : np.ndarray
        value と同単位の絶対不確かさ (1σ).
    label : str
        ラベル (任意).
    """
    wavelength_nm: np.ndarray
    value: np.ndarray
    uncertainty: np.ndarray
    label: str = ""

    def __post_init__(self):
        self.wavelength_nm = np.asarray(self.wavelength_nm, dtype=float)
        self.value = np.asarray(self.value, dtype=float)
        self.uncertainty = np.asarray(self.uncertainty, dtype=float)
        n = len(self.wavelength_nm)
        if not (len(self.value) == n and len(self.uncertainty) == n):
            raise ValueError(f"{self.label}: 配列長不一致")
        # 単調増加チェック
        if not np.all(np.diff(self.wavelength_nm) > 0):
            raise ValueError(f"{self.label}: 波長は単調増加でなければならない")

    def interp_to(self, wl_grid: np.ndarray) -> "SpectralCurve":
        """共通波長グリッドへ線形補間 (範囲外は 0)."""
        v = np.interp(wl_grid, self.wavelength_nm, self.value, left=0.0, right=0.0)
        u = np.interp(wl_grid, self.wavelength_nm, self.uncertainty, left=0.0, right=0.0)
        return SpectralCurve(wl_grid, v, u, self.label)

    def integrate(self) -> float:
        """∫ value(λ) dλ をトラペズで."""
        return float(np.trapezoid(self.value, self.wavelength_nm))


def common_grid(curves: List[SpectralCurve], step_nm: float = 1.0) -> np.ndarray:
    """全カーブの最小重複波長範囲で共通グリッドを生成."""
    wl_min = max(c.wavelength_nm[0] for c in curves)
    wl_max = min(c.wavelength_nm[-1] for c in curves)
    if wl_min >= wl_max:
        raise ValueError("カーブ間に重複波長範囲がありません")
    return np.arange(wl_min, wl_max + 0.5 * step_nm, step_nm)


# =============================================================================
# 2. 線形系 Eq.(2) — 光源係数 A_i の決定
# =============================================================================
def integral_J(s: SpectralCurve, e: SpectralCurve) -> float:
    """∫ s(λ) e(λ) dλ. s と e は同一波長グリッド前提."""
    if not np.array_equal(s.wavelength_nm, e.wavelength_nm):
        raise ValueError("integral_J: 波長グリッド不一致")
    return float(np.trapezoid(s.value * e.value, s.wavelength_nm))


def solve_eq2(
    subcell_SRs: List[SpectralCurve],
    light_sources: List[SpectralCurve],
    ref_spectrum: SpectralCurve,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Eq.(2) を解いて光源係数 A_i を返す.

    M[j,i] = ∫ s_j(λ) e_i(λ) dλ
    b[j]   = ∫ s_j(λ) E_ref(λ) dλ
    解: A = M^{-1} b

    戻り値
    -------
    A : np.ndarray (n_source,)
    M : np.ndarray (n_subcell, n_source)
    b : np.ndarray (n_subcell,)
    """
    n = len(subcell_SRs)
    if len(light_sources) != n:
        raise ValueError("サブセル数と光源数は一致が必要 (Eq.(2) 連立条件)")

    # 全カーブを共通グリッドに投影
    grid = common_grid(subcell_SRs + light_sources + [ref_spectrum], step_nm=1.0)
    s_g = [s.interp_to(grid) for s in subcell_SRs]
    e_g = [e.interp_to(grid) for e in light_sources]
    Eref_g = ref_spectrum.interp_to(grid)

    M = np.zeros((n, n))
    b = np.zeros(n)
    for j in range(n):
        for i in range(n):
            M[j, i] = integral_J(s_g[j], e_g[i])
        b[j] = integral_J(s_g[j], Eref_g)

    A = np.linalg.solve(M, b)
    return A, M, b


# =============================================================================
# 3. Eq.(4) — 積分値 J の不確かさ
# =============================================================================
def u_J_eq4(s: SpectralCurve, e: SpectralCurve) -> float:
    """論文 Eq.(4) — 波長相関を保守的に上限評価.

    u(J) = ∫ u(s(λ)) e(λ) dλ + ∫ s(λ) u(e(λ)) dλ
    """
    if not np.array_equal(s.wavelength_nm, e.wavelength_nm):
        raise ValueError("u_J_eq4: 波長グリッド不一致")
    term1 = np.trapezoid(s.uncertainty * e.value, s.wavelength_nm)
    term2 = np.trapezoid(s.value * e.uncertainty, s.wavelength_nm)
    return float(term1 + term2)


def u_J_independent(s: SpectralCurve, e: SpectralCurve) -> float:
    """波長間独立(ランダム雑音)を仮定した GUM 流二乗和形:
    u²(J) = ∫ [u(s) e]² dλ + ∫ [s u(e)]² dλ.
    Eq.(4) の保守的見積もりとの対比に用いる.
    """
    if not np.array_equal(s.wavelength_nm, e.wavelength_nm):
        raise ValueError("u_J_independent: 波長グリッド不一致")
    t1 = np.trapezoid((s.uncertainty * e.value) ** 2, s.wavelength_nm)
    t2 = np.trapezoid((s.value * e.uncertainty) ** 2, s.wavelength_nm)
    return float(np.sqrt(t1 + t2))


# =============================================================================
# 4. Eq.(5) — 光源係数 A_i の不確かさ (解析的, 偏微分は数値)
# =============================================================================
def u_Ai_analytical(
    M: np.ndarray, b: np.ndarray,
    u_M: np.ndarray, u_b: np.ndarray,
    eps: float = 1e-6,
) -> np.ndarray:
    """A = M^{-1} b の各要素 J (= M[j,i] と b[j]) に対する数値偏微分で u(A) を求める.

    u(A_i)² = Σ_{j,k} (∂A_i/∂M[j,k])² u(M[j,k])²
            + Σ_{j}   (∂A_i/∂b[j])²   u(b[j])²
    """
    n = M.shape[0]
    A0 = np.linalg.solve(M, b)
    var = np.zeros(n)

    # ∂A/∂M[j,k]
    for j in range(n):
        for k in range(n):
            dM = M.copy()
            scale = max(abs(M[j, k]), 1.0) * eps
            dM[j, k] += scale
            dA = (np.linalg.solve(dM, b) - A0) / scale
            var += (dA * u_M[j, k]) ** 2

    # ∂A/∂b[j]
    for j in range(n):
        db = b.copy()
        scale = max(abs(b[j]), 1.0) * eps
        db[j] += scale
        dA = (np.linalg.solve(M, db) - A0) / scale
        var += (dA * u_b[j]) ** 2

    return np.sqrt(var)


# =============================================================================
# 5. モンテカルロ — 全段一括の不確かさ伝播 (推奨)
# =============================================================================
@dataclass
class MCResult:
    """モンテカルロ結果."""
    A_mean: np.ndarray            # 光源係数の平均
    A_std: np.ndarray             # 光源係数の不確かさ
    Jphot_subcell_mean: np.ndarray  # サブセル光電流の平均 [A] (相対値で良いなら [a.u.])
    Jphot_subcell_std: np.ndarray   # サブセル光電流の不確かさ
    Jphot_subcell_rel: np.ndarray   # サブセル光電流の相対不確かさ [%]


def monte_carlo_subcell_currents(
    subcell_SRs: List[SpectralCurve],
    light_sources: List[SpectralCurve],
    ref_cell_SR: SpectralCurve,
    ref_spectrum: SpectralCurve,
    Jref_per_source: np.ndarray,        # 各光源単独点灯時の基準セル光電流 [A]
    u_Jref_per_source: np.ndarray,      # その不確かさ [A]
    n_samples: int = 5000,
    correlation: str = "systematic",     # "systematic" or "independent"
    seed: int = 0,
) -> MCResult:
    """サブセル光電流の不確かさをモンテカルロで求める.

    correlation:
      "systematic" — s_j(λ), e_i(λ) を波長間で完全相関 (校正系統不確かさを想定)
      "independent" — 波長点ごとに独立 (ランダム雑音を想定)
      実機では中間にあるため、両方を実行して上下限を把握するのが推奨.
    """
    rng = np.random.default_rng(seed)
    n = len(subcell_SRs)

    grid = common_grid(subcell_SRs + light_sources + [ref_cell_SR, ref_spectrum], step_nm=1.0)
    s_nom = [s.interp_to(grid) for s in subcell_SRs]
    e_nom = [e.interp_to(grid) for e in light_sources]
    Sref_nom = ref_cell_SR.interp_to(grid)
    Eref_g = ref_spectrum.interp_to(grid)

    A_samples = np.zeros((n_samples, n))
    Jphot_samples = np.zeros((n_samples, n))   # サブセルごとの光電流

    for k in range(n_samples):
        # --- 各カーブを摂動 ---
        if correlation == "systematic":
            # 1 サンプルにつき 1 つの正規乱数で全波長を一括スケール
            s_pert = [SpectralCurve(grid, s.value + rng.standard_normal() * s.uncertainty,
                                    s.uncertainty, s.label) for s in s_nom]
            e_pert = [SpectralCurve(grid, e.value + rng.standard_normal() * e.uncertainty,
                                    e.uncertainty, e.label) for e in e_nom]
            Sref_pert = SpectralCurve(grid, Sref_nom.value + rng.standard_normal() * Sref_nom.uncertainty,
                                      Sref_nom.uncertainty, Sref_nom.label)
        elif correlation == "independent":
            # 波長点ごとに独立 — 波長間の高速振動を許す
            s_pert = [SpectralCurve(grid, s.value + rng.standard_normal(len(grid)) * s.uncertainty,
                                    s.uncertainty, s.label) for s in s_nom]
            e_pert = [SpectralCurve(grid, e.value + rng.standard_normal(len(grid)) * e.uncertainty,
                                    e.uncertainty, e.label) for e in e_nom]
            Sref_pert = SpectralCurve(grid,
                                      Sref_nom.value + rng.standard_normal(len(grid)) * Sref_nom.uncertainty,
                                      Sref_nom.uncertainty, Sref_nom.label)
        else:
            raise ValueError("correlation は 'systematic' または 'independent'")

        # --- 基準セル光電流 (Eq.3 入力) を摂動 ---
        Jref_pert = Jref_per_source + rng.standard_normal(n) * u_Jref_per_source

        # --- Eq.(2) を解く: M_pert A = b_pert ---
        M_p = np.zeros((n, n))
        b_p = np.zeros(n)
        for j in range(n):
            for i in range(n):
                M_p[j, i] = integral_J(s_pert[j], e_pert[i])
            b_p[j] = integral_J(s_pert[j], Eref_g)
        try:
            A_p = np.linalg.solve(M_p, b_p)
        except np.linalg.LinAlgError:
            continue

        # --- Eq.(3) との整合性: A_p は基準セル光電流からも独立に決まるはずなので,
        #     ここではモンテカルロ的に「両者を別経路で出して相対比例で扱う」.
        #     簡略化のため Eq.(2) の A_p を用いて、サブセル光電流を計算:
        #     J_phot,j = A_p · M_p[j, :] = b_p[j] (定義より一致)
        # 実用上は基準セル経由 (Eq.3) の A_i を別途求めて相互チェック.
        # ここではサブセル光電流を「シミュレータが目指す値」として基準スペクトル下値とする:
        for j in range(n):
            # サブセルが基準スペクトル下で発生する光電流 (摂動された SR 形状で)
            Jphot_samples[k, j] = integral_J(s_pert[j], Eref_g)
        A_samples[k] = A_p

    # 統計量
    A_mean = A_samples.mean(axis=0)
    A_std = A_samples.std(axis=0, ddof=1)
    Jphot_mean = Jphot_samples.mean(axis=0)
    Jphot_std = Jphot_samples.std(axis=0, ddof=1)
    Jphot_rel = 100.0 * Jphot_std / np.abs(Jphot_mean)

    return MCResult(A_mean, A_std, Jphot_mean, Jphot_std, Jphot_rel)


# =============================================================================
# 6. Eq.(6) — 修正版スペクトロメトリック評価からの u(Y_i)
# =============================================================================
@dataclass
class IVVariation:
    """1 サブセルの x スイープ I-V 結果 (基準値 Y_i(1) と Y_i(x) 群).

    nominal_at_one : Dict[str, float]
        基準条件 (x=1) でのセル特性 Y_i(1).
    x_array : np.ndarray
        振った相対光電流 (例: 0.95, 0.96, ..., 1.05).
    Y_arrays : Dict[str, np.ndarray]
        各セル特性の Y_i(x). x_array と同長.
    delta_relative : float
        対象サブセルの光電流不確かさ (相対値, 例: 0.025).
    """
    nominal_at_one: Dict[str, float]
    x_array: np.ndarray
    Y_arrays: Dict[str, np.ndarray]
    delta_relative: float

    def u_Y(self, param: str) -> float:
        """Eq.(6) を計算: u(Y_i) [絶対値, Y と同単位].

        u(Y_i) = √ { (1/2Δ) ∫_{1-Δ}^{1+Δ} (Y_i(x) - Y_i(1))² dx }
        """
        Y = np.asarray(self.Y_arrays[param], dtype=float)
        Y1 = float(self.nominal_at_one[param])
        x = np.asarray(self.x_array, dtype=float)
        D = float(self.delta_relative)

        # 積分範囲外の点は除外 (台形則は等間隔である必要なし)
        mask = (x >= 1 - D) & (x <= 1 + D)
        if mask.sum() < 3:
            raise ValueError(f"{param}: x 範囲 [{1-D},{1+D}] 内に 3 点以上必要")
        xm = x[mask]
        Ym = Y[mask]
        integrand = (Ym - Y1) ** 2
        integral = np.trapezoid(integrand, xm)
        # 区間幅 2Δ で割る (実データの x 範囲が 2Δ より狭い場合は実範囲で正規化)
        actual_width = xm[-1] - xm[0]
        norm_width = max(actual_width, 1e-12)
        return float(np.sqrt(integral / norm_width))

    def u_Y_relative(self, param: str) -> float:
        """相対不確かさ [%]"""
        Y1 = float(self.nominal_at_one[param])
        return 100.0 * self.u_Y(param) / abs(Y1)


# =============================================================================
# 7. サブセル合成 + 拡張不確かさ
# =============================================================================
def combine_subcells_RSS(
    variations: Dict[str, IVVariation],
    param: str,
    extra_uncertainties_relative: Optional[Dict[str, float]] = None,
    coverage_k: float = 2.0,
) -> Dict[str, float]:
    """全サブセルの寄与を RSS 合成し、拡張不確かさを返す.

    extra_uncertainties_relative: {"temperature": 0.1, "area": 0.2, ...} [%]

    戻り値: {"u_top_pct", "u_mid_pct", "u_bot_pct", "u_combined_pct", "U_expanded_pct"}
    """
    rels = {name: var.u_Y_relative(param) for name, var in variations.items()}
    sum_sq = sum(r ** 2 for r in rels.values())
    if extra_uncertainties_relative:
        sum_sq += sum(v ** 2 for v in extra_uncertainties_relative.values())
    u_combined = float(np.sqrt(sum_sq))
    out = {f"u_{k}_pct": v for k, v in rels.items()}
    out["u_combined_pct"] = u_combined
    out["U_expanded_pct"] = coverage_k * u_combined
    return out


def make_table1(
    variations: Dict[str, IVVariation],
    parameters: List[str],
    extra_uncertainties_relative: Optional[Dict[str, float]] = None,
    coverage_k: float = 2.0,
) -> "pd.DataFrame":
    """論文 Table 1 形式の DataFrame を返す."""
    import pandas as pd
    subcell_names = list(variations.keys())
    df = pd.DataFrame(
        index=parameters,
        columns=[f"u_{n} [%]" for n in subcell_names]
                + ["u_combined [%]", "U (k=2) [%]"],
    )
    for p in parameters:
        rels = [variations[n].u_Y_relative(p) for n in subcell_names]
        for n, r in zip(subcell_names, rels):
            df.loc[p, f"u_{n} [%]"] = round(r, 3)
        sum_sq = sum(r ** 2 for r in rels)
        if extra_uncertainties_relative:
            sum_sq += sum(v ** 2 for v in extra_uncertainties_relative.values())
        u_combined = np.sqrt(sum_sq)
        df.loc[p, "u_combined [%]"] = round(u_combined, 3)
        df.loc[p, "U (k=2) [%]"] = round(coverage_k * u_combined, 3)
    return df.T  # 行: 寄与, 列: パラメータ
