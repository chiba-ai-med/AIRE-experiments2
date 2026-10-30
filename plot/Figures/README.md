# RIPHAMMER Paper Figures

論文ターゲット: BMC Genome Biology

全ファイルは `plot/` 内の元プロットへのシンボリックリンク。
元プロットを更新すれば自動的に反映される。

---

## Main Figures

### Fig 1: Method schematic (TBD)

RIPHAMMER のアルゴリズム模式図。距離バンド分解 → Multi-band NMF (W_d 独立, H 共有) → 細胞型特異的 contact map 再構成の流れ。別途作成予定。

### Fig 2: Simulation benchmark (`Fig2_simulation_benchmark.png`)

シミュレーション（K=4 細胞型 × 6 混合シナリオ, 100kb）での手法比較。
RIPHAMMER vs THUNDER vs Standard NMF vs PCA。

- 上段: 各手法の fraction recovery scatter plot (MAE, Pearson r, Pure separation)
- 下段左: シナリオ別 MAE heatmap
- 下段右: Overall MAE, Correlation のバープロット

**メッセージ**: RIPHAMMER (MAE=0.022, r=0.999) は THUNDER (MAE=0.043) と同等以上の精度で、PCA (MAE=0.193) を大幅に上回る。

元プロット: `method_comparison_4way.png`

### Fig 3: Biological discovery capability (`Fig3_biological_discovery.png`)

シミュレーションデータでの TAD 構造復元能力の比較。

- 上段: TAD boundary detection (Precision/Recall/F1), Contact structure correlation (box plot), Mixing proportion accuracy
- 中段: chr10 の contact map 比較 (True vs RIPHAMMER vs THUNDER vs Std NMF)
- 下段: Insulation score profile の比較

**メッセージ**: RIPHAMMER は TAD 境界検出 F1=0.79 で最高。THUNDER (F1=0.83) と同等だが、RIPHAMMER は 100kb 解像度で動作可能（THUNDER は 10Mb 限定）。

元プロット: `biological_discovery_comparison.png`

### Fig 4: scHi-C pseudo-bulk validation (`Fig4a/b/c`)

実データ（scHi-C）からの pseudo-bulk 混合実験による精度検証。

- **Fig4a** (`Fig4a_scHiC_gage_seq.png`): GAGE-seq (K=3, M=20)
  - MAE=0.136, mean r=0.735
  - 3セルタイプ (Excitatory/Inhibitory/Oligodendrocyte) 全てで正の相関
- **Fig4b** (`Fig4b_scHiC_droplet_merged.png`): Droplet Hi-C merged (K=5, M=20)
  - MAE=0.117, mean r=0.515
  - GL/GA/glia にマージした5カテゴリでの検証
- **Fig4c** (`Fig4c_scHiC_combined.png`): 3条件の比較バープロット
  - GAGE-seq (K=3), Droplet Hi-C (K=20), Droplet Hi-C merged (K=5)

**メッセージ**: 実データでも fraction recovery が機能する。K が小さいほど精度が高い（K=3: r=0.735 > K=5: r=0.515 > K=20: r=0.415）。

元プロット: `scHiC_validation_gage_seq.png`, `scHiC_validation_droplet_hic_merged.png`, `scHiC_validation_combined.png`

### Fig 5: mTEC deconvolution (`Fig5_mTEC_contact_maps.png`)

本番データ（Thymocyte 11 + mTEC 1 = M=12, K=8）での deconvolution 結果。

- 8コンポーネント × 3ゲノム領域（Aire 座位, chr1, chr6）の contact map
- 各コンポーネントは Hungarian matching でセルタイプにアサイン済み
- comp6 (Bonev_NPC) が checkerboard パターンを示す（brain データの compartment 構造）

**メッセージ**: RIPHAMMER はバルク Hi-C から細胞型特異的な TAD 構造を抽出できる。

元プロット: `contact_maps_100000_k8.png`

### Fig 6: Aire locus analysis (`Fig6_aire_locus.png`)

Aire 遺伝子座（chr10: 76.75-78.6Mb, 10kb 解像度）の詳細解析。

- 上段: mTEC component vs Thymocyte average の contact map、Differential map (log2 ratio)
- 中段: Insulation score profile（mTEC vs Thymocyte）
- 下段: コンポーネント別の Aire 座位 contact intensity

**メッセージ**: RIPHAMMER により mTEC 特異的な Aire 座位周辺のクロマチン構造を同定。mTEC (comp1) で特徴的な TAD パターンが見られる。

元プロット: `aire_collaborator_summary.png`

---

## Supplementary Figures

| Figure | ファイル名 | 内容 |
|--------|-----------|------|
| FigS1 | `FigS1_fraction_recovery_detail.png` | シミュレーション fraction recovery の詳細（scatter + T/P heatmap） |
| FigS2 | `FigS2_H_matrix_heatmaps.png` | RIPHAMMER / Std NMF / PCA の H 行列 heatmap 比較 |
| FigS3 | `FigS3_locus_discovery.png` | 座位別 differential contact correlation（Aire, Tcrb, Cd4/Cd8 各座位） |
| FigS4 | `FigS4_robustness_analysis.png` | K 選択（reconstruction error, separation, stability）+ ノイズ耐性 |
| FigS5 | `FigS5_noise_robustness.png` | Downsampling / Gaussian noise のロバストネス（RIPHAMMER vs Std NMF） |
| FigS6 | `FigS6_scalability.png` | 計算スケーラビリティ（特徴量数, K, サンプル数に対する実行時間） |
| FigS7 | `FigS7_scHiC_droplet_K20.png` | Droplet Hi-C K=20 の結果（Limitation: K=M で精度低下） |
| FigS8 | `FigS8_contact_maps_25kb.png` | mTEC contact maps (25kb 解像度) |
| FigS9 | `FigS9_contact_maps_10kb.png` | mTEC contact maps (10kb 解像度) |
| FigS10 | `FigS10_aire_deep_dive.png` | Aire 座位の多パネル詳細（contact maps + differential + insulation） |
| FigS11 | `FigS11_aire_differential.png` | Aire 座位のコンポーネント別 differential contact maps |
| FigS12 | `FigS12_aire_insulation.png` | Aire 座位の insulation profile + TAD boundary 位置 |
| FigS13 | `FigS13_aire_selected_components.png` | Aire 座位の主要コンポーネント contact maps (wide view) |
| FigS14 | `FigS14_aire_all_components.png` | Aire 座位の全8コンポーネント contact maps |
| FigS15 | `FigS15_brain_dev.png` | Brain developmental data (K=2, 100kb) — 参考 |
| FigS16 | `FigS16_brain_loci.png` | Brain loci analysis (100kb) — 参考 |

---

## 論文ストーリーの流れ

1. **Fig1** (模式図): RIPHAMMER の手法概要
2. **Fig2-3** (シミュレーション): 精度・TAD構造復元で既存手法を上回ることを定量的に示す。特に THUNDER との比較で **100kb 解像度での動作** を強調
3. **Fig4** (scHi-C): 実データでの精度検証。K/M 比の影響も議論
4. **Fig5-6** (mTEC 本番): 生物学的発見。バルク Hi-C から細胞型特異的な TAD を抽出し、Aire 座位で mTEC 特異的なクロマチン構造を同定

## 注意事項

- Fig1 は未作成。イラストツール（BioRender 等）で別途作成予定
- 全プロットの生成スクリプトは `src/` に保存済み
- 元プロットを更新すれば symlink 経由で自動反映
