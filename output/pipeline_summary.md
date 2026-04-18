## Pipeline summary (Baur & Bett 2005 steps)

Output directory: `C:\Users\shiro\cursor\solar_simulator_analysis\output`

| Step | Title | Equation | Artifacts |
|---:|---|---|---|
| 1 | 入力データの読込 | - | - |
| 2 | 光源係数の決定 | Eq.(2) | `A`, `M`, `b` |
| 3 | 積分値 J の波長依存不確かさ (保守上限) | Eq.(4) | `uJ_relative_pct` |
| 4 | サブセル光電流の MC 不確かさ | Eq.(2)(3)経路の数値実装 | `mc_subcell_summary`, `mc_perturb_breakdown_systematic`, `mc_perturb_breakdown_independent` |
| 5 | 相対光電流不確かさ Delta_j の確定 | Step 4 (systematic, all) | `delta_relative` |
| 6 | 修正版スペクトロメトリック I-V の読込 | - | - |
| 7 | Eq.(6) 寄与と最終 RSS / 拡張不確かさ | Eq.(6) + RSS | `eq6_contributions_pct`, `extras_squared`, `table1` |

