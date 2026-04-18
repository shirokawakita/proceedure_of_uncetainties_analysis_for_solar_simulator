"""
サンプル CSV データ生成スクリプト.
GaInP/GaInAs/Ge 三接合セル + キセノン+ハロゲン×2 マルチソースシミュレータ
の典型値を合成し、論文 Table 1 と同程度のスケールの不確かさを再現できるデータを作る.

実機データに置き換える際は、本スクリプトが出力した CSV と同フォーマットに整える.
"""

from __future__ import annotations
import os
import sys
import numpy as np
import pandas as pd

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

OUT_DIR = os.path.join(os.path.dirname(__file__), ".")  # 同じディレクトリに出力
os.makedirs(OUT_DIR, exist_ok=True)
E490_XLS = os.path.join(os.path.dirname(__file__), "e490_00a_amo.xls")


# =============================================================================
# 1. 波長グリッド
# =============================================================================
WL = np.arange(300.0, 1900.0 + 1e-6, 5.0)  # 300–1900 nm, 5 nm 刻み


# =============================================================================
# 2. AM0 基準スペクトル (ASTM E490 表 .xls があれば優先、なければ簡易合成)
# =============================================================================
def am0_synthetic(wl_nm: np.ndarray) -> np.ndarray:
    """AM0 1367 W/m² に正規化した簡易プランク + 紫外/赤外調整 (E490 不在時のフォールバック)."""
    h = 6.626e-34; c = 3e8; kB = 1.381e-23; T = 5778.0
    wl_m = wl_nm * 1e-9
    # プランク放射 (相対形)
    B = 2 * h * c**2 / wl_m**5 / (np.exp(h * c / (wl_m * kB * T)) - 1)
    # 太陽は半径 R, 距離 1AU で見込む立体角
    R_sun = 6.96e8; d = 1.496e11
    irrad = B * np.pi * (R_sun / d) ** 2 * 1e-9  # W/m²/nm
    # 1367 W/m² に正規化
    total = np.trapezoid(irrad, wl_nm)
    irrad *= 1367.0 / total
    return irrad


def write_reference_am0_csv() -> None:
    """E490 .xls から 300–1900 nm グリッドへ補間して reference_spectrum_AM0.csv を書く."""
    out_csv = os.path.join(OUT_DIR, "reference_spectrum_AM0.csv")
    if os.path.isfile(E490_XLS):
        from mj_solar_uncertainty.io import (
            export_reference_spectrum_csv,
            reference_am0_from_e490_xls_to_grid,
        )

        ref = reference_am0_from_e490_xls_to_grid(E490_XLS, WL, label="AM0 ASTM E490")
        export_reference_spectrum_csv(ref, out_csv)
        tot = float(np.trapezoid(ref.value, ref.wavelength_nm))
        print(
            f"  reference_spectrum_AM0.csv <- ASTM E490 (integral over grid ~{tot:.1f} W/m^2)"
        )
    else:
        am0 = am0_synthetic(WL)
        pd.DataFrame({"wavelength_nm": WL, "irradiance": am0}).to_csv(out_csv, index=False)
        print(f"  reference_spectrum_AM0.csv <- 合成プランク (E490 .xls なし)")


write_reference_am0_csv()


# =============================================================================
# 3. サブセル相対分光感度 (典型形状)
# =============================================================================
def _sigmoid_band(wl, low, high, peak_value=0.45, edge_sharpness=20.0):
    """low〜high nm に感度を持つ矩形に近い帯域、両端をガウシアンで滑らかに."""
    sr = peak_value * np.ones_like(wl)
    sr *= 1.0 / (1.0 + np.exp(-(wl - low) / (edge_sharpness)))
    sr *= 1.0 / (1.0 + np.exp((wl - high) / (edge_sharpness)))
    return sr


def make_sr(low, high, peak, label, sharpness=8.0):
    sr = _sigmoid_band(WL, low, high, peak_value=peak, edge_sharpness=sharpness)
    # 不確かさ: フラット域 0.5%, バンドエッジ近傍で大きく
    rel_uncert = 0.005 * np.ones_like(WL)  # フラット 0.5%
    # バンドエッジ判定 (一階微分が大きいところ)
    deriv = np.abs(np.gradient(sr, WL))
    edge_factor = deriv / (deriv.max() + 1e-12)
    # エッジで最大 30% (相対) まで上昇
    rel_uncert += 0.30 * edge_factor
    abs_uncert = sr * rel_uncert
    df = pd.DataFrame({
        "wavelength_nm": WL,
        "sr_value": sr,
        "sr_uncertainty": abs_uncert,
    })
    df.to_csv(os.path.join(OUT_DIR, f"subcell_SR_{label}.csv"), index=False)
    return df


sr_top = make_sr(low=370, high=680, peak=0.40, label="top")  # GaInP
sr_mid = make_sr(low=680, high=900, peak=0.50, label="mid")  # GaInAs
sr_bot = make_sr(low=900, high=1850, peak=0.55, label="bot",
                 sharpness=15.0)                              # Ge — 長波長で不確かさ大

# 基準セル SR (Si pyranometer 想定の広帯域フォトダイオード相当)
sr_ref = _sigmoid_band(WL, 400, 1100, peak_value=0.45, edge_sharpness=15.0)
rel_unc_ref = 0.01 + 0.10 * (np.abs(np.gradient(sr_ref, WL)) /
                              (np.abs(np.gradient(sr_ref, WL)).max() + 1e-12))
pd.DataFrame({
    "wavelength_nm": WL,
    "sr_value": sr_ref,
    "sr_uncertainty": sr_ref * rel_unc_ref,
}).to_csv(os.path.join(OUT_DIR, "ref_cell_SR.csv"), index=False)


# =============================================================================
# 4. 光源スペクトル (キセノン + ハロゲン × 2)
# =============================================================================
def xenon_spectrum(wl):
    """キセノンランプ: 連続成分 + 800–1000 nm 領域のスパイク群."""
    # 連続成分 (黒体 6500K 寄り)
    h = 6.626e-34; c = 3e8; kB = 1.381e-23; T = 6500.0
    wl_m = wl * 1e-9
    B = 2 * h * c**2 / wl_m**5 / (np.exp(h * c / (wl_m * kB * T)) - 1)
    spec = B / B.max()
    # 800–1000 nm にスパイク群を上乗せ
    for cw, w, h_ in [(823, 8, 1.5), (842, 8, 1.8), (881, 10, 2.0), (920, 12, 1.6), (980, 8, 1.2)]:
        spec += h_ * np.exp(-((wl - cw) / w) ** 2)
    spec *= 1.5  # 全体スケール
    return spec


def halogen_spectrum(wl, T=3200.0, peak_scale=1.0):
    """ハロゲンタングステン: 黒体 ~ 3200K"""
    h = 6.626e-34; c = 3e8; kB = 1.381e-23
    wl_m = wl * 1e-9
    B = 2 * h * c**2 / wl_m**5 / (np.exp(h * c / (wl_m * kB * T)) - 1)
    return peak_scale * B / B.max()


def save_lightsource(name, spec):
    # 不確かさ: 全波長一律 1.5% 相対 (校正不確かさ想定)
    rel_unc = 0.015 * np.ones_like(WL)
    abs_unc = spec * rel_unc
    pd.DataFrame({
        "wavelength_nm": WL,
        "irradiance": spec,
        "irradiance_uncertainty": abs_unc,
    }).to_csv(os.path.join(OUT_DIR, f"light_{name}.csv"), index=False)


save_lightsource("xenon",   xenon_spectrum(WL))
save_lightsource("halogen1", halogen_spectrum(WL, T=3000.0, peak_scale=0.8))  # 短波長強
save_lightsource("halogen2", halogen_spectrum(WL, T=2700.0, peak_scale=1.2))  # 長波長強


# =============================================================================
# 5. 修正版スペクトロメトリック評価 I-V データ (合成)
# =============================================================================
def synth_3J_iv(x_top, x_mid, x_bot, J_top0=17.0, J_mid0=16.5, J_bot0=22.0):
    """超簡易 3J モデル.
    各サブセル光電流 J_j = J_j0 * x_j.
    Isc = min(J_top, J_mid, J_bot) (電流律速).
    Voc = sum_j (V_j0 + (kT/q) ln(J_j / J_j0))   *V_j0 を典型値で固定
    Pmpp/FF はカレントマッチ近傍で極大、ミスマッチで微増 (FF↑)/(Vmpp↑) で部分相殺.
    """
    Vth = 0.0259
    # 各サブセル開放電圧 (典型値)
    Voc_top = 1.40 + Vth * np.log(max(x_top, 1e-6))
    Voc_mid = 1.00 + Vth * np.log(max(x_mid, 1e-6))
    Voc_bot = 0.27 + Vth * np.log(max(x_bot, 1e-6))
    Voc = Voc_top + Voc_mid + Voc_bot

    J_top = J_top0 * x_top
    J_mid = J_mid0 * x_mid
    J_bot = J_bot0 * x_bot
    Isc = min(J_top, J_mid, J_bot)

    # FF: ベース 0.85, ミスマッチで上昇 (極小は完全マッチ近傍)
    mismatch = (max(J_top, J_mid, J_bot) - Isc) / Isc
    FF = 0.86 + 0.04 * mismatch  # 単調増加 — ミスマッチで FF↑ (paper Fig.1)
    FF = min(FF, 0.92)

    Pmpp = FF * Voc * Isc
    Vmpp = 0.85 * Voc + 0.05 * mismatch  # ミスマッチで微増
    Impp = Pmpp / max(Vmpp, 1e-6)
    return dict(Isc=Isc, Voc=Voc, Pmpp=Pmpp, FF=FF, Impp=Impp, Vmpp=Vmpp)


# x スイープ: 対象サブセル ±Δ, 他は 1.0 固定
rows = []
deltas = {"top": 0.025, "mid": 0.025, "bot": 0.05}
for sub, d in deltas.items():
    xs = np.linspace(1 - d, 1 + d, 11)
    for x in xs:
        x_top = x if sub == "top" else 1.0
        x_mid = x if sub == "mid" else 1.0
        x_bot = x if sub == "bot" else 1.0
        r = synth_3J_iv(x_top, x_mid, x_bot)
        rows.append({"subcell": sub, "x": round(float(x), 4), **r})

df_iv = pd.DataFrame(rows)
df_iv.to_csv(os.path.join(OUT_DIR, "spectrometric_IV.csv"), index=False)


# =============================================================================
# 6. 基準セル光電流 (各光源単独点灯時) — 簡易計算
# =============================================================================
# Eq.(3): J_ref,i = A_i ∫ S_ref(λ) e_i(λ) dλ.
# 校正用なので、ここでは「各光源 1.0 倍出力」としての値を出力 (実測代替).
def _intgr(s_arr, e_arr, wl):
    return float(np.trapezoid(s_arr * e_arr, wl))


sources = {
    "xenon":    xenon_spectrum(WL),
    "halogen1": halogen_spectrum(WL, T=3000.0, peak_scale=0.8),
    "halogen2": halogen_spectrum(WL, T=2700.0, peak_scale=1.2),
}
Jref_rows = []
for name, e in sources.items():
    Jref = _intgr(sr_ref, e, WL)
    Jref_rows.append({"source": name, "Jref_A": Jref, "Jref_uncertainty_A": 0.005 * Jref})
pd.DataFrame(Jref_rows).to_csv(os.path.join(OUT_DIR, "ref_cell_currents.csv"), index=False)


print("サンプルデータ生成完了:")
for f in sorted(os.listdir(OUT_DIR)):
    if f.endswith(".csv"):
        path = os.path.join(OUT_DIR, f)
        n_rows = sum(1 for _ in open(path)) - 1
        print(f"  {f:40s}  {n_rows:5d} rows")
