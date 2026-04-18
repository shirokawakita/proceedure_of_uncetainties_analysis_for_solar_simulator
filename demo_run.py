"""
エンドツーエンド・デモ
======================
sample_data/ 配下の CSV を読み込み、Baur & Bett (2005) 手法で
ソーラシミュレータ校正の不確かさを計算し、論文 Table 1 形式で出力する.

実行:
    python3 demo_run.py

出力:
    output/table1_eq6.csv      - Eq.(6) ベースの拡張不確かさ表
    output/table1_eq6.md       - 同 Markdown 版
    output/mc_subcell_currents.csv - モンテカルロ結果 (サブセル光電流不確かさ)
    output/spectro_curves.png  - 修正版スペクトロメトリック評価カーブ
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# パッケージインポート (相対パスから)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mj_solar_uncertainty import (
    load_spectral, load_reference_spectrum, load_astm_e490_am0_xls, load_iv_variation,
    solve_eq2, u_J_eq4, u_J_independent,
    monte_carlo_subcell_currents, make_table1, save_table1_csv, save_table1_md,
)

DATA = os.path.join(os.path.dirname(__file__), "sample_data")
OUT = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(OUT, exist_ok=True)

print("=" * 70)
print("ソーラシミュレータ不確かさ評価デモ - Baur & Bett (2005) 手法")
print("=" * 70)


# =============================================================================
# 1. データ読み込み
# =============================================================================
print("\n[1] 計測データ読み込み")
sr_top = load_spectral(os.path.join(DATA, "subcell_SR_top.csv"), label="top (GaInP)")
sr_mid = load_spectral(os.path.join(DATA, "subcell_SR_mid.csv"), label="mid (GaInAs)")
sr_bot = load_spectral(os.path.join(DATA, "subcell_SR_bot.csv"), label="bot (Ge)")
sr_ref = load_spectral(os.path.join(DATA, "ref_cell_SR.csv"), label="ref cell")

e_xenon  = load_spectral(os.path.join(DATA, "light_xenon.csv"),    label="xenon")
e_halo1  = load_spectral(os.path.join(DATA, "light_halogen1.csv"), label="halogen1")
e_halo2  = load_spectral(os.path.join(DATA, "light_halogen2.csv"), label="halogen2")

_e490_xls = os.path.join(DATA, "e490_00a_amo.xls")
if os.path.isfile(_e490_xls):
    ref_spec = load_astm_e490_am0_xls(_e490_xls, label="AM0 ASTM E490")
else:
    ref_spec = load_reference_spectrum(
        os.path.join(DATA, "reference_spectrum_AM0.csv"), label="AM0"
    )

ref_currents_df = pd.read_csv(os.path.join(DATA, "ref_cell_currents.csv"))
print(f"  サブセル SR: top {len(sr_top.wavelength_nm)} 点, mid {len(sr_mid.wavelength_nm)} 点, bot {len(sr_bot.wavelength_nm)} 点")
print(f"  光源スペクトル: xenon, halogen1, halogen2")
print(f"  基準スペクトル: AM0 (積分 {ref_spec.integrate():.1f} W/m^2)")
print(f"  基準セル光電流: {len(ref_currents_df)} 光源分")


# =============================================================================
# 2. Eq.(2) 連立方程式の解 — 光源係数 A_i
# =============================================================================
print("\n[2] Eq.(2) を解いて光源係数 A_i を決定")
subcell_SRs = [sr_top, sr_mid, sr_bot]
sources     = [e_xenon, e_halo1, e_halo2]
A, M, b = solve_eq2(subcell_SRs, sources, ref_spec)
print(f"  A (xenon, halogen1, halogen2) = {A}")
print("  検算 M @ A:", M @ A)
print("  目標   b:", b)


# =============================================================================
# 3. Eq.(4) — 各 (subcell, source) 組合せの u(J)
# =============================================================================
print("\n[3] Eq.(4) - 積分値 J の波長依存不確かさ")
from mj_solar_uncertainty.core import common_grid
grid = common_grid(subcell_SRs + sources + [ref_spec], step_nm=1.0)
s_g = [s.interp_to(grid) for s in subcell_SRs]
e_g = [e.interp_to(grid) for e in sources]
ref_g = ref_spec.interp_to(grid)

names_sub = ["top", "mid", "bot"]
names_src = ["xenon", "halo1", "halo2"]

uJ_table = pd.DataFrame(index=names_sub, columns=names_src + ["vs E_ref"], dtype=float)
for j, sj in enumerate(s_g):
    for i, ei in enumerate(e_g):
        J = float(np.trapezoid(sj.value * ei.value, grid))
        uJ = u_J_eq4(sj, ei)
        rel = 100.0 * uJ / abs(J) if abs(J) > 0 else np.nan
        uJ_table.iloc[j, i] = round(rel, 3)
    # vs E_ref (b[j])
    Jb = float(np.trapezoid(sj.value * ref_g.value, grid))
    uJb = u_J_eq4(sj, ref_g)
    uJ_table.iloc[j, -1] = round(100.0 * uJb / abs(Jb), 3)

print("  u(J) 相対 [%] - 行: subcell, 列: source")
print(uJ_table.to_string())


# =============================================================================
# 4. モンテカルロ — サブセル光電流の不確かさ
# =============================================================================
print("\n[4] モンテカルロ - サブセル光電流不確かさ u(J_phot,j)")

Jref_arr = ref_currents_df["Jref_A"].values.astype(float)
uJref_arr = ref_currents_df["Jref_uncertainty_A"].values.astype(float)

# (a) 系統誤差モデル (波長間完全相関)
mc_sys = monte_carlo_subcell_currents(
    subcell_SRs, sources, sr_ref, ref_spec,
    Jref_per_source=Jref_arr, u_Jref_per_source=uJref_arr,
    n_samples=3000, correlation="systematic", seed=42,
)
# (b) 独立誤差モデル (波長間独立)
mc_ind = monte_carlo_subcell_currents(
    subcell_SRs, sources, sr_ref, ref_spec,
    Jref_per_source=Jref_arr, u_Jref_per_source=uJref_arr,
    n_samples=3000, correlation="independent", seed=42,
)

mc_df = pd.DataFrame({
    "subcell": names_sub,
    "Jphot_AM0_[a.u.]":  np.round(mc_sys.Jphot_subcell_mean, 4),
    "u_Jphot_systematic_[%]": np.round(mc_sys.Jphot_subcell_rel, 3),
    "u_Jphot_independent_[%]": np.round(mc_ind.Jphot_subcell_rel, 3),
})
print(mc_df.to_string(index=False))
mc_df.to_csv(os.path.join(OUT, "mc_subcell_currents.csv"), index=False)


# =============================================================================
# 5. Eq.(6) — 修正版スペクトロメトリック評価から u(Y_i)
# =============================================================================
print("\n[5] Eq.(6) - 修正版スペクトロメトリック評価から u(Y_i)")

# サブセル毎の光電流不確かさ (Step 4 の結果を採用)
delta_dict = {
    "top": mc_sys.Jphot_subcell_rel[0] / 100.0,
    "mid": mc_sys.Jphot_subcell_rel[1] / 100.0,
    "bot": mc_sys.Jphot_subcell_rel[2] / 100.0,
}
print(f"  Δ (相対光電流不確かさ): {delta_dict}")

variations = load_iv_variation(
    os.path.join(DATA, "spectrometric_IV.csv"),
    delta_relative_per_subcell=delta_dict,
    parameters=["Isc", "Voc", "Pmpp", "FF", "Impp", "Vmpp"],
)


# =============================================================================
# 6. 表 1 形式の最終結果
# =============================================================================
print("\n[6] 論文 Table 1 形式の最終結果")
extra = {
    "temperature": 0.10,  # ±0.5°C 由来 (例)
    "area":        0.20,  # 面積測定 (例)
    "electrical":  0.10,  # 電気計測 (例)
}
df_table1 = make_table1(
    variations,
    parameters=["Isc", "Voc", "Pmpp", "FF", "Impp", "Vmpp"],
    extra_uncertainties_relative=extra,
    coverage_k=2.0,
)
print(df_table1.to_string())

save_table1_csv(df_table1, os.path.join(OUT, "table1_eq6.csv"))
save_table1_md(df_table1, os.path.join(OUT, "table1_eq6.md"),
               title="ソーラシミュレータ校正不確かさ - Eq.(6) ベース")


# =============================================================================
# 7. 修正版スペクトロメトリック評価カーブ作図 (論文 Fig.4 相当)
# =============================================================================
print("\n[7] 修正版スペクトロメトリック評価カーブを描画")
fig, axes = plt.subplots(6, 1, figsize=(7, 11), sharex=True)
params = ["Isc", "Voc", "Pmpp", "FF", "Impp", "Vmpp"]
markers = {"top": "s", "mid": "o", "bot": "^"}
colors  = {"top": "C0",  "mid": "C1", "bot": "C2"}

iv_df = pd.read_csv(os.path.join(DATA, "spectrometric_IV.csv"))
for ax, p in zip(axes, params):
    for sub in ["top", "mid", "bot"]:
        d = iv_df[iv_df.subcell == sub].sort_values("x")
        # 基準値 (x=1) で正規化したパーセント変化
        nominal = float(d[d.x == 1.0][p].values[0]) if (d.x == 1.0).any() \
                  else float(d.iloc[len(d)//2][p])
        rel = 100.0 * (d[p].values - nominal) / nominal
        ax.plot(d.x.values, rel, marker=markers[sub], color=colors[sub],
                label=sub, linewidth=1.0, markersize=5)
    ax.set_ylabel(f"Δ{p} [%]", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="k", linewidth=0.5)
axes[0].legend(loc="upper right", ncol=3, fontsize=9)
axes[-1].set_xlabel("x = J_phot(subcell) / J_phot(AM0)")
fig.suptitle("Modified spectrometric characterization curves", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "spectro_curves.png"), dpi=140)
plt.close(fig)


print("\n" + "=" * 70)
print("完了. 出力ファイル:")
for f in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, f)
    print(f"  {p}  ({os.path.getsize(p):,} bytes)")
print("=" * 70)
