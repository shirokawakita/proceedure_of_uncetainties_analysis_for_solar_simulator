"""
Baur & Bett (2005) 不確かさ評価の Step 分割パイプライン.

各 Step の中間表を DataFrame として返し、任意で output/ に CSV 保存する.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .core import (
    IVVariation,
    MCResult,
    MCPerturbMode,
    common_grid,
    eq6_contribution_percent_matrix,
    make_table1,
    monte_carlo_subcell_currents,
    solve_eq2,
    u_J_eq4,
)
from .io import (
    load_astm_e490_am0_xls,
    load_iv_variation,
    load_reference_spectrum,
    load_spectral,
    save_table1_csv,
    save_table1_md,
)


@dataclass
class PipelineStep:
    """1 Step の結果."""

    id: int
    title: str
    equation: str
    description: str
    artifacts: Dict[str, Any] = field(default_factory=dict)
    meta: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PipelineRunResult:
    """全 Step 実行後のまとめ."""

    steps: List[PipelineStep]
    variations: Dict[str, IVVariation]
    table1: pd.DataFrame
    mc_systematic: MCResult
    mc_independent: MCResult
    delta_dict: Dict[str, float]
    extra_uncertainties_relative: Dict[str, float]
    coverage_k: float


def load_reference_am0_for_data_dir(data_dir: str):
    """sample_data 用: E490 .xls があれば優先."""
    xls = os.path.join(data_dir, "e490_00a_amo.xls")
    if os.path.isfile(xls):
        return load_astm_e490_am0_xls(xls, label="AM0 ASTM E490")
    return load_reference_spectrum(
        os.path.join(data_dir, "reference_spectrum_AM0.csv"), label="AM0"
    )


def uJ_eq4_relative_percent_table(
    subcell_SRs: List,
    sources: List,
    ref_spec,
    names_sub: List[str],
    names_src: List[str],
    step_nm: float = 1.0,
) -> pd.DataFrame:
    """Eq.(4): 各 (subcell, source) と vs E_ref で u(J) の相対 [%] 表."""
    grid = common_grid(subcell_SRs + sources + [ref_spec], step_nm=step_nm)
    s_g = [s.interp_to(grid) for s in subcell_SRs]
    e_g = [e.interp_to(grid) for e in sources]
    ref_g = ref_spec.interp_to(grid)
    cols = list(names_src) + ["vs_E_ref"]
    uJ_table = pd.DataFrame(index=names_sub, columns=cols, dtype=float)
    for j, sj in enumerate(s_g):
        for i, ei in enumerate(e_g):
            J = float(np.trapezoid(sj.value * ei.value, grid))
            uJ = u_J_eq4(sj, ei)
            rel = 100.0 * uJ / abs(J) if abs(J) > 0 else np.nan
            uJ_table.iloc[j, i] = round(rel, 3)
        Jb = float(np.trapezoid(sj.value * ref_g.value, grid))
        uJb = u_J_eq4(sj, ref_g)
        uJ_table.iloc[j, -1] = round(100.0 * uJb / abs(Jb), 3)
    return uJ_table


def extras_squared_table(extra: Dict[str, float]) -> pd.DataFrame:
    """追加不確かさ (相対 %): 各要因と二乗値."""
    rows = []
    for k, v in extra.items():
        rows.append({"factor": k, "u_relative_pct": v, "u_squared": v**2})
    return pd.DataFrame(rows)


def mc_perturb_breakdown_table(
    subcell_SRs: List,
    sources: List,
    sr_ref,
    ref_spec,
    Jref_arr: np.ndarray,
    uJref_arr: np.ndarray,
    names_sub: List[str],
    perturb_modes: Tuple[MCPerturbMode, ...],
    n_samples: int,
    correlation: str,
    seed: int,
) -> pd.DataFrame:
    """各 perturb モードでの u(J_phot) 相対 [%] (系統または独立の単一 correlation)."""
    rows = []
    for mode in perturb_modes:
        mc = monte_carlo_subcell_currents(
            subcell_SRs,
            sources,
            sr_ref,
            ref_spec,
            Jref_per_source=Jref_arr,
            u_Jref_per_source=uJref_arr,
            n_samples=n_samples,
            correlation=correlation,
            seed=seed,
            perturb=mode,
        )
        row: Dict[str, Any] = {"perturb": mode}
        for i, name in enumerate(names_sub):
            row[f"u_Jphot_{name}_pct"] = round(float(mc.Jphot_subcell_rel[i]), 3)
        rows.append(row)
    return pd.DataFrame(rows)


def save_dataframe_csv(df: pd.DataFrame, path: str) -> None:
    df.to_csv(path, encoding="utf-8-sig")


def write_pipeline_summary_md(
    path: str,
    steps: List[PipelineStep],
    out_dir: str,
) -> None:
    lines = [
        "## Pipeline summary (Baur & Bett 2005 steps)",
        "",
        f"Output directory: `{out_dir}`",
        "",
        "| Step | Title | Equation | Artifacts |",
        "|---:|---|---|---|",
    ]
    for st in steps:
        arts = ", ".join(f"`{k}`" for k in st.artifacts.keys()) if st.artifacts else "-"
        lines.append(f"| {st.id} | {st.title} | {st.equation} | {arts} |")
    lines.append("")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def run_uncertainty_pipeline(
    data_dir: str,
    out_dir: str,
    *,
    n_mc_samples: int = 3000,
    mc_seed: int = 42,
    parameters: Optional[List[str]] = None,
    extra_uncertainties_relative: Optional[Dict[str, float]] = None,
    coverage_k: float = 2.0,
    save_intermediate_csv: bool = True,
    iv_csv_name: str = "spectrometric_IV.csv",
) -> PipelineRunResult:
    """Step 1–7 を実行し、Table 1 を含む PipelineRunResult を返す."""
    if parameters is None:
        parameters = ["Isc", "Voc", "Pmpp", "FF", "Impp", "Vmpp"]
    if extra_uncertainties_relative is None:
        extra_uncertainties_relative = {
            "temperature": 0.10,
            "area": 0.20,
            "electrical": 0.10,
        }

    os.makedirs(out_dir, exist_ok=True)
    steps: List[PipelineStep] = []

    # ----- Step 1 -----
    sr_top = load_spectral(os.path.join(data_dir, "subcell_SR_top.csv"), label="top")
    sr_mid = load_spectral(os.path.join(data_dir, "subcell_SR_mid.csv"), label="mid")
    sr_bot = load_spectral(os.path.join(data_dir, "subcell_SR_bot.csv"), label="bot")
    sr_ref = load_spectral(os.path.join(data_dir, "ref_cell_SR.csv"), label="ref")
    e_xenon = load_spectral(os.path.join(data_dir, "light_xenon.csv"), label="xenon")
    e_halo1 = load_spectral(os.path.join(data_dir, "light_halogen1.csv"), label="halogen1")
    e_halo2 = load_spectral(os.path.join(data_dir, "light_halogen2.csv"), label="halogen2")
    ref_spec = load_reference_am0_for_data_dir(data_dir)
    ref_currents_df = pd.read_csv(os.path.join(data_dir, "ref_cell_currents.csv"))

    subcell_SRs = [sr_top, sr_mid, sr_bot]
    sources = [e_xenon, e_halo1, e_halo2]
    names_sub = ["top", "mid", "bot"]
    names_src = ["xenon", "halogen1", "halogen2"]

    st1 = PipelineStep(
        1,
        "入力データの読込",
        "-",
        "SR, 光源, AM0, 基準セル光電流 CSV",
        meta={
            "n_sr_points": {n: len(s.wavelength_nm) for n, s in zip(names_sub, subcell_SRs)},
            "ref_spectrum_integral_W_m2": ref_spec.integrate(),
            "n_ref_current_sources": len(ref_currents_df),
        },
    )
    steps.append(st1)

    # ----- Step 2 -----
    A, M, b = solve_eq2(subcell_SRs, sources, ref_spec)
    res_eq2 = float(np.max(np.abs(M @ A - b)))
    df_A = pd.DataFrame({"A_i": A}, index=names_src)
    df_M = pd.DataFrame(M, index=names_sub, columns=names_src)
    df_b = pd.DataFrame({"b_j": b}, index=names_sub)
    st2 = PipelineStep(
        2,
        "光源係数の決定",
        "Eq.(2)",
        "M A = b を解き A_i を得る",
        artifacts={"A": df_A, "M": df_M, "b": df_b},
        meta={"max_abs_residual_MA_minus_b": res_eq2},
    )
    steps.append(st2)

    # ----- Step 3 -----
    uJ_table = uJ_eq4_relative_percent_table(
        subcell_SRs, sources, ref_spec, names_sub, names_src
    )
    st3 = PipelineStep(
        3,
        "積分値 J の波長依存不確かさ (保守上限)",
        "Eq.(4)",
        "u(J) 相対 [%] 表",
        artifacts={"uJ_relative_pct": uJ_table},
    )
    steps.append(st3)

    # ----- Step 4 -----
    Jref_arr = ref_currents_df["Jref_A"].values.astype(float)
    uJref_arr = ref_currents_df["Jref_uncertainty_A"].values.astype(float)

    mc_sys = monte_carlo_subcell_currents(
        subcell_SRs,
        sources,
        sr_ref,
        ref_spec,
        Jref_per_source=Jref_arr,
        u_Jref_per_source=uJref_arr,
        n_samples=n_mc_samples,
        correlation="systematic",
        seed=mc_seed,
        perturb="all",
    )
    mc_ind = monte_carlo_subcell_currents(
        subcell_SRs,
        sources,
        sr_ref,
        ref_spec,
        Jref_per_source=Jref_arr,
        u_Jref_per_source=uJref_arr,
        n_samples=n_mc_samples,
        correlation="independent",
        seed=mc_seed,
        perturb="all",
    )

    mc_main = pd.DataFrame(
        {
            "subcell": names_sub,
            "Jphot_AM0_[a.u.]": np.round(mc_sys.Jphot_subcell_mean, 4),
            "u_Jphot_systematic_[%]": np.round(mc_sys.Jphot_subcell_rel, 3),
            "u_Jphot_independent_[%]": np.round(mc_ind.Jphot_subcell_rel, 3),
        }
    )

    perturb_modes: Tuple[MCPerturbMode, ...] = (
        "all",
        "subcell_sr",
        "lights",
        "ref_sr",
        "Jref",
    )
    mc_break_sys = mc_perturb_breakdown_table(
        subcell_SRs,
        sources,
        sr_ref,
        ref_spec,
        Jref_arr,
        uJref_arr,
        names_sub,
        perturb_modes,
        n_mc_samples,
        "systematic",
        mc_seed,
    )
    mc_break_ind = mc_perturb_breakdown_table(
        subcell_SRs,
        sources,
        sr_ref,
        ref_spec,
        Jref_arr,
        uJref_arr,
        names_sub,
        perturb_modes,
        n_mc_samples,
        "independent",
        mc_seed,
    )

    st4 = PipelineStep(
        4,
        "サブセル光電流の MC 不確かさ",
        "Eq.(2)(3)経路の数値実装",
        "系統/独立 + perturb 別感度分解 (分解は相関無視の参考値)",
        artifacts={
            "mc_subcell_summary": mc_main,
            "mc_perturb_breakdown_systematic": mc_break_sys,
            "mc_perturb_breakdown_independent": mc_break_ind,
        },
        meta={
            "n_samples": n_mc_samples,
            "note_Jphot_integral": (
                "J_phot,j = integral_J(s_j_pert, E_ref) のみ統計; "
                "lights/ref_sr/Jref の摂動は現行式では J_phot に現れない場合あり"
            ),
        },
    )
    steps.append(st4)

    # ----- Step 5 -----
    delta_dict = {
        "top": float(mc_sys.Jphot_subcell_rel[0] / 100.0),
        "mid": float(mc_sys.Jphot_subcell_rel[1] / 100.0),
        "bot": float(mc_sys.Jphot_subcell_rel[2] / 100.0),
    }
    st5 = PipelineStep(
        5,
        "相対光電流不確かさ Delta_j の確定",
        "Step 4 (systematic, all)",
        "Eq.(6) の積分幅に用いる Delta (相対)",
        artifacts={"delta_relative": pd.DataFrame([delta_dict])},
        meta={"source": "monte_carlo_subcell_currents systematic perturb=all"},
    )
    steps.append(st5)

    # ----- Step 6 -----
    iv_path = os.path.join(data_dir, iv_csv_name)
    variations = load_iv_variation(
        iv_path,
        delta_relative_per_subcell=delta_dict,
        parameters=parameters,
    )
    st6 = PipelineStep(
        6,
        "修正版スペクトロメトリック I-V の読込",
        "-",
        iv_csv_name,
        meta={"path": iv_path},
    )
    steps.append(st6)

    # ----- Step 7 -----
    contrib = eq6_contribution_percent_matrix(variations, parameters)
    df_table1 = make_table1(
        variations,
        parameters=parameters,
        extra_uncertainties_relative=extra_uncertainties_relative,
        coverage_k=coverage_k,
    )
    df_extras = extras_squared_table(extra_uncertainties_relative)

    st7 = PipelineStep(
        7,
        "Eq.(6) 寄与と最終 RSS / 拡張不確かさ",
        "Eq.(6) + RSS",
        "サブセル別 u(Y) と合成、Table 1",
        artifacts={
            "eq6_contributions_pct": contrib,
            "extras_squared": df_extras,
            "table1": df_table1,
        },
        meta={"coverage_k": coverage_k},
    )
    steps.append(st7)

    if save_intermediate_csv:
        save_dataframe_csv(df_A, os.path.join(out_dir, "step02_A.csv"))
        save_dataframe_csv(df_M, os.path.join(out_dir, "step02_M.csv"))
        save_dataframe_csv(df_b, os.path.join(out_dir, "step02_b.csv"))
        save_dataframe_csv(uJ_table, os.path.join(out_dir, "step03_uJ_eq4_relative_pct.csv"))
        save_dataframe_csv(mc_main, os.path.join(out_dir, "mc_subcell_currents.csv"))
        save_dataframe_csv(mc_break_sys, os.path.join(out_dir, "step04_mc_perturb_breakdown_systematic.csv"))
        save_dataframe_csv(mc_break_ind, os.path.join(out_dir, "step04_mc_perturb_breakdown_independent.csv"))
        save_dataframe_csv(contrib, os.path.join(out_dir, "step07_eq6_contributions_pct.csv"))
        save_dataframe_csv(df_extras, os.path.join(out_dir, "step07_extras_squared.csv"))
        save_table1_csv(df_table1, os.path.join(out_dir, "table1_eq6.csv"))
        save_table1_md(
            df_table1,
            os.path.join(out_dir, "table1_eq6.md"),
            title="ソーラシミュレータ校正不確かさ - Eq.(6) ベース",
        )
        write_pipeline_summary_md(os.path.join(out_dir, "pipeline_summary.md"), steps, out_dir)

    return PipelineRunResult(
        steps=steps,
        variations=variations,
        table1=df_table1,
        mc_systematic=mc_sys,
        mc_independent=mc_ind,
        delta_dict=delta_dict,
        extra_uncertainties_relative=extra_uncertainties_relative,
        coverage_k=coverage_k,
    )


def plot_spectrometric_curves(
    iv_csv_path: str,
    png_out: str,
    parameters: Optional[List[str]] = None,
) -> None:
    """修正版スペクトロメトリック評価カーブ (論文 Fig.4 相当)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if parameters is None:
        parameters = ["Isc", "Voc", "Pmpp", "FF", "Impp", "Vmpp"]
    iv_df = pd.read_csv(iv_csv_path)
    fig, axes = plt.subplots(6, 1, figsize=(7, 11), sharex=True)
    markers = {"top": "s", "mid": "o", "bot": "^"}
    colors = {"top": "C0", "mid": "C1", "bot": "C2"}
    for ax, p in zip(axes, parameters):
        for sub in ["top", "mid", "bot"]:
            d = iv_df[iv_df.subcell == sub].sort_values("x")
            if (d.x == 1.0).any():
                nominal = float(d.loc[d.x == 1.0, p].values[0])
            else:
                nominal = float(d.iloc[len(d) // 2][p])
            rel = 100.0 * (d[p].values.astype(float) - nominal) / nominal
            ax.plot(
                d.x.values,
                rel,
                marker=markers[sub],
                color=colors[sub],
                label=sub,
                linewidth=1.0,
                markersize=5,
            )
        ax.set_ylabel(f"d{p} [%]", fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color="k", linewidth=0.5)
    axes[0].legend(loc="upper right", ncol=3, fontsize=9)
    axes[-1].set_xlabel("x = J_phot(subcell) / J_phot(AM0)")
    fig.suptitle("Modified spectrometric characterization curves", fontsize=11)
    fig.tight_layout()
    os.makedirs(os.path.dirname(png_out) or ".", exist_ok=True)
    fig.savefig(png_out, dpi=140)
    plt.close(fig)
