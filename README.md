# AIRE-experiments2

RIPHAMMER (Multi-band NMF for Hi-C deconvolution) の検証実験リポジトリ。  
論文ターゲット: BMC Genome Biology

## 検証の全体像

| 検証 | 目的 | 状態 |
|---|---|---|
| シミュレーション | RIPHAMMER vs THUNDER vs Standard NMF の比較 | 完了 (`plot/` に19枚) |
| scHi-C pseudo-bulk | 実データでの RIPHAMMER 精度検証 | 完了 (3条件) |
| mTEC 本番 | 生物学的発見 (ground truth なし) | 未着手 |

## scHi-C Pseudo-bulk Validation

### 概要

scHi-C データからセルタイプ別に pseudo-bulk contact matrix を作り、  
**既知の比率で混合 → RIPHAMMER で deconvolution → 真の比率との一致を評価**  
という流れで精度を検証する。

### データセット

| | GAGE-seq | Droplet Hi-C |
|---|---|---|
| 論文 | Zhou et al. 2024, Nat Genet | Chang et al. 2024, Nat Biotech |
| GEO | GSE238001 | GSE253407 |
| 生物種・組織 | mouse brain cortex | mouse brain cortex |
| 技術 | sci-Hi-C + scRNA co-assay | droplet microfluidics scHi-C |
| 細胞数 | 9,047 | 6,211 |
| セルタイプ数 (K) | 3 (Excitatory/Inhibitory/Oligodendrocyte) | 20 |
| raw pairs | 3ファイル 計15GB | 2ファイル 計43GB |
| バーコード形式 | combinatorial index (例: A1,A1) | 16bp DNA barcode |
| pairs format | ヘッダなし 5列 (chr pos chr pos barcode) | 4DN standard 10列 (#columns ヘッダ付き) |

### パイプライン (各データセットで共通)

```
1. raw pairs ダウンロード (GEO FTP)
2. メタデータ準備 (セルタイプアノテーション)
     - GAGE-seq: RNA component から TruncatedSVD + KMeans でクラスタリング
     - Droplet Hi-C: 論文付属メタデータ (cell_metadata.xls)
3. pseudo-bulk 生成 (scHiC_pseudobulk.py)
     pairs をセルタイプ別に集約 → per-celltype .mcool
4. synthetic mixtures 作成 (scHiC_create_mixtures.py)
     Dirichlet(α=1) で M=20 の混合比をサンプリング → 擬似 bulk Hi-C
5. vectorize (scHiC_vectorize_mixtures.py)
     mcool → distance-band 特徴量行列
6. RIPHAMMER fit (scHiC_fit.py)
     MultibandNMF (K=セルタイプ数, mode=single)
7. evaluate (scHiC_evaluate.py)
     Hungarian algorithm でコンポーネント-セルタイプ対応を決定 → MAE, Pearson r
8. plot (plot_scHiC_validation.py, plot_scHiC_validation_combined.py)
```

### 解像度

- **100kb** を使用 (CLAUDE.md の方針に準拠)
- scHi-C はスパースなので 25kb 以下は pseudo-bulk でも密度不足のリスクあり
- 競合手法 THUNDER は Tucker 分解で bins×bins×samples テンソルを扱うため、
  100kb (~27,000 bins) では計算コスト的に非現実的。
  THUNDER 論文のベンチマークも 10Mb (~270 bins) のみ。

### 結果

| データセット | K | M | MAE | Mean Pearson r | 備考 |
|---|---|---|---|---|---|
| GAGE-seq | 3 | 20 | 0.136 | **0.735** | 主要検証。3セルタイプ全てで正の相関 |
| Droplet Hi-C | 20 | 20 | 0.034 | **0.415** | K=M でも全タイプ正の相関 (Limitation 用) |
| Droplet Hi-C (merged) | 5 | 20 | 0.117 | **0.515** | GL/GA/glia にマージ。K<M で改善 |

プロット: `plot/scHiC_validation_{gage_seq,droplet_hic,droplet_hic_merged,combined}.png`

### 評価指標

- **MAE**: 推定比率と真の比率の平均絶対誤差
- **Per-cell-type Pearson r**: セルタイプごとに M=20 mixtures の (真, 推定) の相関。得意・不得意の確認用
- **Pure separation**: 各セルタイプについて「最大比率のサンプル群で正しく最大と推定されたか」

NMF は出力の順番が不定なので、Hungarian algorithm (O(K³)) で最適対応を決めてから評価する。

### Droplet Hi-C merged のセルタイプマッピング

| マージ後 | 元のセルタイプ |
|---|---|
| Excitatory | CLAGL, CTGL, ITL23GL, ITL45GL, ITL5GL, ITL6GL, NPGL, OBGL, PTGL |
| Inhibitory | D12MSN, OBGA, PVGA, SSTGA, STRGA, VIPGA |
| Astrocyte | ASC |
| Oligo | OGC, OPC |
| Other_glia | MGL, VLMC |

### ディレクトリ構成

```
data/scHiC/
  gage_seq/
    raw/              # .pairs.gz (3ファイル)
    metadata/         # cell_annotations.tsv, pairs_format.json
    pseudobulk/       # per-celltype .mcool + summary.tsv
    mixtures/         # synthetic bulk .mcool + fractions.tsv
    results/          # RIPHAMMER output (H.tsv, W.tsv, meta.json)
    evaluation/       # evaluation.tsv
  droplet_hic/
    raw/              # .pairs.gz (2ファイル)
    metadata/         # cell_annotations.tsv, pairs_format.json, cell_metadata.xls
    pseudobulk/       # (同上)
    mixtures/
    results/
    evaluation/
  droplet_hic_merged/
    merge_map.json    # セルタイプマージ定義
    pseudobulk/       # マージ後 pseudo-bulk (K=5)
    mixtures/
    results/
    evaluation/
```

## シミュレーションベンチマーク (完了)

`output/simulation_100000/` にシミュレーション結果。  
RIPHAMMER vs THUNDER vs Standard NMF を K=4, 100kb 解像度で比較。  
`plot/` に19枚のプロット。

## mTEC 本番 (config.yaml)

- Multiome reference: `data/mtec/reference/` (cellranger-arc output)
- Bulk Hi-C: `data/mtec/processed/mtec_hic.mcool` (SRP330308 fastq から再構築が必要)
- Thymocyte bulk Hi-C: 11 サンプル (`data/4dn_thymocyte/`)
- 合計 M=12 サンプルで RIPHAMMER を実行予定

## 技術メモ

- **numpy==1.26.4 に固定**: h5py, pandas, numexpr が numpy 1.x でビルドされているため numpy 2.x 不可
- **pseudo-bulk のメモリ**: 全 contacts をリスト保持すると OOM する。Counter で on-the-fly ビン集計する方式に修正済み
- **Droplet Hi-C の #columns ヘッダ**: `chrom1/chrom2` 表記 (4DN standard の `chr1/chr2` と異なる)。inspect_pairs.py で両方対応
- **Droplet Hi-C バーコード**: pairs ファイルは 16bp DNA のみ、メタデータは lib:barcode 形式。`df["bc"]` 列を使用
- **RIPHAMMER パッケージ**: `/home/koki/dev/RIPHAMMER` (branch: `feat/single-matrix-and-h-norm`)
