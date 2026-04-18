"""
README 用: sample_data/ の CSV を可視化し docs/figures/ に PNG を出力する.

実行 (リポジトリルートから):
    python scripts/plot_readme_figures.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "sample_data"
OUT = ROOT / "docs" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Windows で日本語ラベルが化けにくい順
plt.rcParams["font.sans-serif"] = [
    "MS Gothic",
    "Yu Gothic",
    "Meiryo",
    "Noto Sans CJK JP",
    "DejaVu Sans",
]
plt.rcParams["axes.unicode_minus"] = False


def _read_csv(name: str) -> pd.DataFrame:
    return pd.read_csv(DATA / name)


def plot_fig01_reference_and_sources() -> None:
    """基準スペクトル AM0 と各光源の分光放射照度."""
    am0 = _read_csv("reference_spectrum_AM0.csv")
    xen = _read_csv("light_xenon.csv")
    h1 = _read_csv("light_halogen1.csv")
    h2 = _read_csv("light_halogen2.csv")

    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    wl = am0.iloc[:, 0].values
    axes[0].plot(wl, am0.iloc[:, 1].values, color="k", lw=1.2, label="AM0 (ASTM E490 補間)")
    if am0.shape[1] >= 3:
        u = am0.iloc[:, 2].values.astype(float)
        v = am0.iloc[:, 1].values.astype(float)
        axes[0].fill_between(wl, v - u, v + u, color="k", alpha=0.12, linewidth=0)
    axes[0].set_ylabel("放射照度 [W/m^2/nm]")
    axes[0].set_title("基準スペクトル (E490 xls -> reference_spectrum_AM0.csv)")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="upper right")

    axes[1].plot(xen.iloc[:, 0], xen.iloc[:, 1], label="xenon", lw=1.0)
    axes[1].plot(h1.iloc[:, 0], h1.iloc[:, 1], label="halogen1", lw=1.0)
    axes[1].plot(h2.iloc[:, 0], h2.iloc[:, 1], label="halogen2", lw=1.0)
    axes[1].set_xlabel("波長 [nm]")
    axes[1].set_ylabel("放射照度 [a.u.]")
    axes[1].set_title("各光源単独点灯 (light_*.csv)")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc="upper right", ncol=3)
    fig.tight_layout()
    fig.savefig(OUT / "fig01_reference_and_sources.png", dpi=150)
    plt.close(fig)


def plot_fig02_subcell_and_ref_sr() -> None:
    """サブセル・基準セルの分光感度 (値)."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    for name, lab, c in [
        ("subcell_SR_top.csv", "top (GaInP)", "C0"),
        ("subcell_SR_mid.csv", "mid (GaInAs)", "C1"),
        ("subcell_SR_bot.csv", "bot (Ge)", "C2"),
        ("ref_cell_SR.csv", "ref cell", "0.35"),
    ]:
        df = _read_csv(name)
        ax.plot(df.iloc[:, 0], df.iloc[:, 1], lw=1.1, label=lab, color=c)
    ax.set_xlabel("波長 [nm]")
    ax.set_ylabel("感度 sr_value [a.u.]")
    ax.set_title("分光感度 (subcell_SR_*.csv, ref_cell_SR.csv)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "fig02_subcell_and_ref_sr.png", dpi=150)
    plt.close(fig)


def plot_fig03_relative_uncertainty_sr() -> None:
    """SR の相対不確かさ u(sr)/|sr| — バンドエッジで増える合成の特徴."""
    fig, ax = plt.subplots(figsize=(9, 4.2))
    for name, lab, c in [
        ("subcell_SR_top.csv", "top", "C0"),
        ("subcell_SR_mid.csv", "mid", "C1"),
        ("subcell_SR_bot.csv", "bot", "C2"),
    ]:
        df = _read_csv(name)
        wl = df.iloc[:, 0].values
        v = df.iloc[:, 1].values.astype(float)
        u = df.iloc[:, 2].values.astype(float)
        rel = 100.0 * np.abs(u) / np.maximum(np.abs(v), 1e-12)
        ax.plot(wl, rel, lw=1.0, label=lab, color=c)
    ax.set_xlabel("波長 [nm]")
    ax.set_ylabel("相対不確かさ u(sr) / |sr| [%]")
    ax.set_title("サブセル SR の波長依存不確かさ (CSV 3 列目は絶対 1σ)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT / "fig03_sr_relative_uncertainty.png", dpi=150)
    plt.close(fig)


def plot_fig04_light_relative_uncertainty() -> None:
    """光源の相対不確かさ (合成では全波長 1.5%)."""
    fig, ax = plt.subplots(figsize=(9, 3.8))
    for name, lab in [
        ("light_xenon.csv", "xenon"),
        ("light_halogen1.csv", "halogen1"),
        ("light_halogen2.csv", "halogen2"),
    ]:
        df = _read_csv(name)
        wl = df.iloc[:, 0].values
        v = df.iloc[:, 1].values.astype(float)
        u = df.iloc[:, 2].values.astype(float)
        rel = 100.0 * np.abs(u) / np.maximum(np.abs(v), 1e-12)
        ax.plot(wl, rel, lw=1.0, label=lab)
    ax.set_xlabel("波長 [nm]")
    ax.set_ylabel("u(E) / |E| [%]")
    ax.set_title("光源スペクトルの相対不確かさ (generate_sample_data では平坦 1.5%)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", ncol=3)
    fig.tight_layout()
    fig.savefig(OUT / "fig04_light_relative_uncertainty.png", dpi=150)
    plt.close(fig)


def plot_fig05_spectrometric_iv() -> None:
    """spectrometric_IV.csv — 各サブセルを振ったときのパラメータ変化."""
    iv = _read_csv("spectrometric_IV.csv")
    params = ["Isc", "Voc", "Pmpp", "FF", "Impp", "Vmpp"]
    markers = {"top": "s", "mid": "o", "bot": "^"}
    colors = {"top": "C0", "mid": "C1", "bot": "C2"}

    fig, axes = plt.subplots(2, 3, figsize=(10, 5.5), sharex=False)
    axes = axes.ravel()
    for ax, p in zip(axes, params):
        for sub in ["top", "mid", "bot"]:
            d = iv[iv.subcell == sub].sort_values("x")
            if (d.x == 1.0).any():
                nominal = float(d.loc[d.x == 1.0, p].values[0])
            else:
                nominal = float(d.iloc[len(d) // 2][p])
            rel = 100.0 * (d[p].values.astype(float) - nominal) / nominal
            ax.plot(d.x.values, rel, marker=markers[sub], color=colors[sub],
                    label=sub, lw=1.0, ms=4)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title(p)
        ax.set_xlabel("x")
        ax.set_ylabel("Δ [%]")
        ax.grid(True, alpha=0.3)
    axes[0].legend(loc="best", fontsize=8, ncol=3)
    fig.suptitle("修正版スペクトロメトリック評価 (spectrometric_IV.csv)", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "fig05_spectrometric_iv.png", dpi=150)
    plt.close(fig)


def plot_fig06_make_table1_heatmap() -> None:
    """output/table1_eq6.csv があればヒートマップで行列表現を補助."""
    csv_path = ROOT / "output" / "table1_eq6.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path, index_col=0)
    # 数値のみ (列はパラメータ名)
    Z = df.astype(float).values
    row_labels = list(df.index)
    col_labels = list(df.columns)

    fig, ax = plt.subplots(figsize=(9.5, 3.8))
    im = ax.imshow(Z, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_xticklabels(col_labels, rotation=0)
    ax.set_yticklabels(row_labels)
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            ax.text(j, i, f"{Z[i, j]:.2f}", ha="center", va="center", color="0.1", fontsize=8)
    ax.set_title(
        "make_table1() の戻り値 (転置後): 行 = 寄与・合成, 列 = I-V パラメータ"
    )
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="[%]")
    fig.tight_layout()
    fig.savefig(OUT / "fig06_make_table1_heatmap.png", dpi=150)
    plt.close(fig)


def main() -> None:
    if not DATA.is_dir():
        print(f"sample_data が見つかりません: {DATA}", file=sys.stderr)
        sys.exit(1)
    plot_fig01_reference_and_sources()
    plot_fig02_subcell_and_ref_sr()
    plot_fig03_relative_uncertainty_sr()
    plot_fig04_light_relative_uncertainty()
    plot_fig05_spectrometric_iv()
    plot_fig06_make_table1_heatmap()
    print("出力:", OUT)
    for p in sorted(OUT.glob("fig*.png")):
        print(f"  {p.name} ({p.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
