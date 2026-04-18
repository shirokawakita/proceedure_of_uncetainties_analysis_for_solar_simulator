"""
エンドツーエンド・デモ
======================
sample_data/ 配下の CSV を読み込み、Baur & Bett (2005) 手法で
ソーラシミュレータ校正の不確かさを Step 分割パイプラインで計算する.

実行:
    python demo_run.py

出力 (output/):
    table1_eq6.csv, table1_eq6.md  - 論文 Table 1 形式の最終結果
    mc_subcell_currents.csv         - モンテカルロ結果 (サブセル光電流不確かさ)
    step02_*.csv, step02b_*.csv, step03_*.csv - Eq.(2)(4) の中間表（2b は合成スペクトル検証）
    step04_mc_perturb_breakdown_*.csv - MC 感度分解 (参考)
    step07_*.csv                    - Eq.(6) 寄与と追加要因の二乗
    pipeline_summary.md             - Step と成果物の対応
    spectro_curves.png              - 修正版スペクトロメトリック評価カーブ
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from mj_solar_uncertainty.pipeline import plot_spectrometric_curves, run_uncertainty_pipeline

DATA = os.path.join(os.path.dirname(__file__), "sample_data")
OUT = os.path.join(os.path.dirname(__file__), "output")
os.makedirs(OUT, exist_ok=True)

print("=" * 70)
print("ソーラシミュレータ不確かさ評価デモ - Baur & Bett (2005) Step パイプライン")
print("=" * 70)

result = run_uncertainty_pipeline(
    DATA,
    OUT,
    n_mc_samples=3000,
    mc_seed=42,
    save_intermediate_csv=True,
)

for st in result.steps:
    print(f"\n[Step {st.id}] {st.title} ({st.equation})")
    if st.meta:
        for k, v in st.meta.items():
            print(f"  {k}: {v}")
    if st.artifacts:
        print(f"  artifacts: {', '.join(st.artifacts.keys())}")

print("\n[最終] Table 1 (Eq.(6) + RSS)")
print(result.table1.to_string())

plot_spectrometric_curves(
    os.path.join(DATA, "spectrometric_IV.csv"),
    os.path.join(OUT, "spectro_curves.png"),
)

print("\n" + "=" * 70)
print("完了. 出力ファイル:")
for f in sorted(os.listdir(OUT)):
    p = os.path.join(OUT, f)
    print(f"  {p}  ({os.path.getsize(p):,} bytes)")
print("=" * 70)
