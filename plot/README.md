# RIPHAMMER 解析結果プロット一覧

RIPHAMMER (Multi-band NMF for Hi-C deconvolution) の全解析結果。
Genome Biology投稿用の材料。

---

## 1. シミュレーションベンチマーク

### method_comparison_4way.png
**4手法比較: RIPHAMMER vs THUNDER vs Standard NMF vs PCA**

- 上段: 真の混合比率 vs 推定比率の散布図（4手法 × 4細胞型）
- 下段左: per-scenario MAE ヒートマップ
- 下段右: 全体MAE・相関のバーチャート

主要結果:
| 手法 | 純粋サンプル分離 | MAE | Mean r |
|---|---|---|---|
| RIPHAMMER | 4/4 | 0.022 | 0.999 |
| THUNDER | 4/4 | 0.043 | 0.995 |
| Standard NMF | 2/4 | 0.011 | 1.000 |
| PCA | 3/4 | 0.193 | 0.621 |

Standard NMFのMAEが最低だが2/4しか分離できない（距離減衰が支配し、細胞型間差異を拾えない）。

### method_comparison_simulation.png
初期版の3手法比較（RIPHAMMER vs Standard NMF vs PCA）。`method_comparison_4way.png` がTHUNDER追加後の最終版。

### simulation_fraction_recovery.png
RIPHAMMERの混合比率回復の詳細散布図。4細胞型 × 6シナリオ。

---

## 2. 生物学的発見能力

### biological_discovery_comparison.png
**TAD境界検出 + insulation相関 + 比率MAEの統合比較**

- 上段左: TAD Precision/Recall/F1（3手法バーチャート）
- 上段中: Insulation score相関の箱ひげ図（per chrom × cell type）
- 上段右: 比率MAEバーチャート
- 下段: chr10の例示（コンタクトマップ + insulation profile）

主要結果:
| 指標 | RIPHAMMER | THUNDER | Std NMF |
|---|---|---|---|
| TAD Precision | 0.867 | 0.802 | 0.679 |
| TAD F1 | 0.773 | 0.796 | 0.664 |
| Insulation r (中央値) | 0.954 | 0.902 | 0.862 |

RIPHAMMERはPrecisionとInsulation相関で最良。

### locus_discovery_comparison.png
**差分コンタクトパターン検出能力**

- 上段: ゲノムワイド sliding window (5Mb) での差分コンタクト相関
  - RIPHAMMER vs THUNDER scatter, RIPHAMMER vs StdNMF scatter, 分布
- 下段: 3 loci (Aire, Tcrb, Cd4/Cd8) での差分コンタクト相関バーチャート

決定的な結果:
| 指標 | RIPHAMMER | THUNDER | Std NMF |
|---|---|---|---|
| 差分コンタクト相関 (全窓) | 0.636 | -0.009 | 0.366 |
| 差分コンタクト相関 (高シグナル窓) | 0.574 | -0.083 | 0.463 |
| 窓ベース勝率 vs THUNDER | 2964/3824 (78%) | 858 | -- |
| Aire locus | 0.819 | 0.049 | -0.374 |

THUNDERの差分コンタクト相関はほぼゼロ。距離減衰を除去しないため、W行列が距離減衰に支配され、細胞型間の差分シグナルが埋没する。RIPHAMMERのband正規化がこの問題を解決。

---

## 3. ロバストネス

### robustness_analysis.png
**K選択 + ノイズ耐性の6パネル**

- (a) 再構成誤差 vs K
- (b) 細胞型分離数 vs K: K>=4 で4/4達成
- (c) seed間安定性 vs K: 全Kで0.87以上
- (d) ダウンサンプリング vs MAE: 20%まで4/4維持
- (e) ガウスノイズ vs MAE
- (f) ダウンサンプリング vs 分離数

K=4-5がMAE最小で真のKと一致。RIPHAMMERは20%欠損まで4/4分離を維持（StdNMFはcleanでも3/4止まり）。

### noise_robustness.png
ノイズロバストネスの初期版（存在する場合）。`robustness_analysis.png` が統合最終版。

### scalability_benchmark.png
**計算時間スケーラビリティ（3パネル）**

- (a) データサイズ vs ランタイム: 1.2M features で RIPHAMMER 42s, StdNMF 59s
- (b) K vs ランタイム
- (c) サンプル数 vs ランタイム

THUNDER (R, nrun=1) は同条件で推定約2.5分。RIPHAMMERはPython実装で高速。

---

## 4. 実データ: mTEC (胸腺髄質上皮細胞)

### aire_collaborator_summary.png
**共同研究者向けAire locus解析サマリー（10kb解像度）**

- (a) mTEC成分 (comp1) のコンタクトマップ（Aire周辺 ±1Mb）
- (b) Thymocyte平均のコンタクトマップ
- (c) 差分コンタクト (log2 mTEC/Thymocyte)
- (d) Insulation profile比較（mTEC vs Thymocyte）
- (e) 全8成分のAire locusコンタクト強度バーチャート

mTEC成分はH loading = 1.000で完全に単離。全解像度(100kb/25kb/10kb)で一貫。青い破線はAire遺伝子位置 (chr10:77.71-77.74 Mb)。

### aire_locus_10kb_deep.png
Aire周辺の全8成分コンタクトマップ（10kb）。各成分のコンタクトパターンが細胞型ごとに異なることを示す。

### aire_locus_10kb_selected.png
Aire周辺の選択された5成分（mTEC + 主要thymocyte）。

### aire_locus_10kb_differential.png
Aire周辺の差分コンタクトマップ（成分間のlog2比較）。

### aire_locus_10kb_insulation.png
Aire周辺のInsulation profileの成分間比較（10kb）。

### aire_locus_10kb_components.png
Aire周辺の全成分パネル（コンタクトマップ + insulation overlay）。

---

## 5. 実データ: Brain (Bonev et al.)

### brain_dev_k2_100kb.png
**ES vs 分化の K=2 deconvolution（100kb）**

Bonev et al. のES/NPC/CN Hi-Cデータに対してK=2でRIPHAMMERを適用。ES特異的成分と分化特異的成分に分離。

### brain_dev_loci_100kb.png
**発生遺伝子座でのコンタクトパターン（100kb）**

Sox2, Pou5f1, Neurog2 など発生関連遺伝子座でのES vs 分化のコンタクト差を可視化。

---

## 論文Figure構成案

| Figure | 内容 | 対応プロット |
|---|---|---|
| Fig 1 | 手法模式図 + パイプライン | (別途作成) |
| Fig 2 | シミュレーションベンチマーク | method_comparison_4way.png |
| Fig 3 | 生物学的発見能力 | biological_discovery_comparison.png + locus_discovery_comparison.png |
| Fig 4 | 実データ: mTEC deconvolution | aire_collaborator_summary.png + aire_locus_10kb_selected.png |
| Fig 5 | 多解像度 + Aire locus | aire_locus_10kb_deep.png + aire_locus_10kb_insulation.png |
| Fig 6 | Brain データセット | brain_dev_k2_100kb.png + brain_dev_loci_100kb.png |
| Supp | ロバストネス + スケーラビリティ | robustness_analysis.png + scalability_benchmark.png |

---

## リプロダクション表

全プロットの生成スクリプトは `src/plot_*.py` に保存済み。レポジトリルートから実行。

| プロット | スクリプト | 計算 | 主な入力 |
|---|---|---|---|
| method_comparison_4way.png | src/plot_method_comparison_4way.py | 軽量 | comparison_results.npz, thunder_H_final.tsv |
| method_comparison_simulation.png | src/plot_method_comparison_simulation.py | NMF+PCA再実行 | X_bulk.npy, H_nmf_k4.tsv |
| simulation_fraction_recovery.png | src/plot_simulation_fraction_recovery.py | NNLS | H_nmf_k4.tsv, true_fractions.tsv |
| biological_discovery_comparison.png | src/plot_biological_discovery_comparison.py | RIPHAMMER+StdNMF再実行 | X_bulk.npy, tad_comparison.npz |
| locus_discovery_comparison.png | src/plot_locus_discovery_comparison.py | RIPHAMMER+StdNMF再実行 | X_bulk.npy, thunder_W_k4.tsv |
| robustness_analysis.png | src/plot_robustness_analysis.py | 複数NMF再実行 | X_bulk.npy, features.tsv |
| scalability_benchmark.png | src/plot_scalability_benchmark.py | タイミングベンチ | X_bulk.npy, features.tsv |
| noise_robustness.png | src/plot_noise_robustness.py | 複数NMF再実行 | X_bulk.npy, features.tsv |
| aire_collaborator_summary.png | src/plot_aire_collaborator_summary.py | 軽量 | nmf_10000/* |
| aire_locus_10kb_components.png | src/plot_aire_locus_10kb_components.py | 軽量 | nmf_10000/*, mtec_hic.mcool |
| aire_locus_10kb_deep.png | src/plot_aire_locus_10kb_deep.py | 軽量 | nmf_10000/*, tad_10000_k8/* |
| aire_locus_10kb_selected.png | src/plot_aire_locus_10kb_selected.py | 軽量 | nmf_10000/*, mtec_hic.mcool |
| aire_locus_10kb_differential.png | src/plot_aire_locus_10kb_differential.py | 軽量 | nmf_10000/* |
| aire_locus_10kb_insulation.png | src/plot_aire_locus_10kb_insulation.py | 軽量 | nmf_10000/*, tad_10000_k8/* |
| brain_dev_k2_100kb.png | src/plot_brain_dev_k2_100kb.py | 軽量 | nmf_brain_100000_k2/*, brain_hic_*.mcool |
| brain_dev_loci_100kb.png | src/plot_brain_dev_loci_100kb.py | 軽量 | nmf_brain_100000/*, brain_hic_*.mcool |

THUNDER関連Rスクリプト:
- `src/run_thunder_fullgenome.R` — THUNDER 3段パイプライン (~28分)
- `src/run_thunder_extract_W.R` — W行列抽出 (nrun=1, ~5分)

全コードは会話トランスクリプト（2715f6e4-0d5d-4392-9bf5-d565ce7ff9f3.jsonl）から抽出。
詳細は `src/REPRODUCIBILITY.md` を参照。
