"""
CSV データ読み書きユーティリティ.

CSV フォーマット
----------------
分光感度 (SR):
    wavelength_nm, sr_value, sr_uncertainty
        sr_value: 相対値 [a.u.] または絶対値 [A/W]
        sr_uncertainty: sr_value と同単位の絶対不確かさ (1σ)

光源スペクトル:
    wavelength_nm, irradiance, irradiance_uncertainty
        irradiance: 相対 [a.u.] または絶対 [W/m²/nm]
        irradiance_uncertainty: 同単位の 1σ

基準スペクトル (AM0 等):
    wavelength_nm, irradiance
        (不確かさ列省略可. 省略時は 0 とみなす)
    または ASTM E490 表の .xls を `load_astm_e490_am0_xls()` で直接読み込み.

修正版スペクトロメトリック評価 I-V:
    subcell, x, Isc, Voc, Pmpp, FF, Impp, Vmpp
        subcell: "top" / "mid" / "bot"
        x: 相対光電流 (1.0 が基準条件)
        以降: 各セル特性の絶対値
        x = 1.0 の行が基準値 Y_i(1) とみなされる.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from .core import SpectralCurve, IVVariation


def load_spectral(path: str, label: str = "", relative_uncertainty: bool = False) -> SpectralCurve:
    """SR / 光源スペクトル CSV を読み込み SpectralCurve に変換.

    relative_uncertainty=True なら 3 列目を「value に対する相対値」とみなし、
    内部では絶対値に展開する.
    """
    df = pd.read_csv(path)
    cols = list(df.columns)
    if len(cols) < 2:
        raise ValueError(f"{path}: 列数不足")
    wl = df.iloc[:, 0].values.astype(float)
    val = df.iloc[:, 1].values.astype(float)
    if len(cols) >= 3:
        unc = df.iloc[:, 2].values.astype(float)
        if relative_uncertainty:
            unc = unc * val
    else:
        unc = np.zeros_like(val)
    return SpectralCurve(wl, val, unc, label=label or path)


def load_reference_spectrum(path: str, label: str = "AM0") -> SpectralCurve:
    """基準スペクトル CSV. 不確かさ列がなければ 0 で埋める."""
    df = pd.read_csv(path)
    wl = df.iloc[:, 0].values.astype(float)
    val = df.iloc[:, 1].values.astype(float)
    if df.shape[1] >= 3:
        unc = df.iloc[:, 2].values.astype(float)
    else:
        unc = np.zeros_like(val)
    return SpectralCurve(wl, val, unc, label=label)


def load_astm_e490_am0_xls(
    path: str,
    label: str = "AM0 ASTM E490",
    sheet_name: str = "NewAM0",
    wl_min_nm: float = 100.0,
    wl_max_nm: float = 1.0e7,
    irradiance_uncertainty_relative: float = 0.0,
) -> SpectralCurve:
    """ASTM E490 系の AM0 スペクトル表 (Excel .xls) を読み込む.

    `sample_data/e490_00a_amo.xls` のように、1 行目がヘッダで
    左端列が「波長 [µm]」、その右が「E-490 W/m²/µm」の形式を想定する.
    放射照度は内部で **W/m²/nm** に換算する (÷1000).

    Parameters
    ----------
    irradiance_uncertainty_relative :
        各波長で |E| に対する相対 1σ (0 なら不確かさ 0).
    """
    df = pd.read_excel(path, sheet_name=sheet_name, header=0, engine="xlrd")
    # 列名の揺れに耐える (先頭 2 列を µm, E490 とみなす)
    c0, c1 = df.columns[0], df.columns[1]
    wl_um = pd.to_numeric(df[c0], errors="coerce").to_numpy(dtype=float)
    irr_um = pd.to_numeric(df[c1], errors="coerce").to_numpy(dtype=float)
    wl_nm = wl_um * 1000.0
    irr_nm = irr_um / 1000.0
    m = np.isfinite(wl_nm) & np.isfinite(irr_nm)
    m &= (wl_nm >= wl_min_nm) & (wl_nm <= wl_max_nm)
    wl_nm = wl_nm[m]
    irr_nm = np.maximum(irr_nm[m], 0.0)
    if len(wl_nm) < 3:
        raise ValueError(f"{path}: 有効なスペクトル点が不足しています")
    if not np.all(np.diff(wl_nm) > 0):
        raise ValueError(f"{path}: 波長列が単調増加ではありません")
    if irradiance_uncertainty_relative > 0:
        unc = np.abs(irr_nm) * float(irradiance_uncertainty_relative)
    else:
        unc = np.zeros_like(irr_nm)
    return SpectralCurve(wl_nm, irr_nm, unc, label=label)


def export_reference_spectrum_csv(
    curve: SpectralCurve,
    path: str,
    include_uncertainty_column: bool = False,
) -> None:
    """基準スペクトルを `load_reference_spectrum` 互換 CSV に書き出す."""
    d = {"wavelength_nm": curve.wavelength_nm, "irradiance": curve.value}
    if include_uncertainty_column:
        d["irradiance_uncertainty"] = curve.uncertainty
    pd.DataFrame(d).to_csv(path, index=False)


def reference_am0_from_e490_xls_to_grid(
    xls_path: str,
    wavelength_nm: np.ndarray,
    **kwargs,
) -> SpectralCurve:
    """E490 .xls を読み、指定波長グリッドへ線形補間した基準スペクトルを返す."""
    ref = load_astm_e490_am0_xls(xls_path, **kwargs)
    grid = np.asarray(wavelength_nm, dtype=float)
    return ref.interp_to(grid)


def load_iv_variation(
    path: str,
    delta_relative_per_subcell: dict,
    parameters: list = None,
) -> dict:
    """修正版スペクトロメトリック評価 I-V CSV から、
    サブセル名 -> IVVariation の辞書を構築.

    delta_relative_per_subcell : {"top": 0.025, "mid": 0.025, "bot": 0.05}
    parameters : ["Isc","Voc","Pmpp","FF","Impp","Vmpp"] 等
    """
    if parameters is None:
        parameters = ["Isc", "Voc", "Pmpp", "FF", "Impp", "Vmpp"]

    df = pd.read_csv(path)
    out = {}
    for subcell, sub_df in df.groupby("subcell"):
        sub_df = sub_df.sort_values("x").reset_index(drop=True)
        # x = 1.0 行を見つける (なければ最近接)
        idx_1 = (sub_df["x"] - 1.0).abs().idxmin()
        nominal = {p: float(sub_df.loc[idx_1, p]) for p in parameters if p in sub_df.columns}
        Y_arr = {p: sub_df[p].values for p in parameters if p in sub_df.columns}
        delta = float(delta_relative_per_subcell.get(subcell, 0.025))
        out[subcell] = IVVariation(
            nominal_at_one=nominal,
            x_array=sub_df["x"].values,
            Y_arrays=Y_arr,
            delta_relative=delta,
        )
    return out


def save_table1_csv(df: pd.DataFrame, path: str) -> None:
    """Table 1 形式の DataFrame を CSV 保存."""
    df.to_csv(path, encoding="utf-8-sig")


def save_table1_md(df: pd.DataFrame, path: str, title: str = "") -> None:
    """Markdown 表として保存 (報告書貼り付け用)."""
    with open(path, "w", encoding="utf-8") as f:
        if title:
            f.write(f"## {title}\n\n")
        f.write(df.to_markdown(floatfmt=".3f"))
        f.write("\n")
