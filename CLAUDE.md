# CLAUDE.md

## Project purpose

**AIRE-experiments2** は、バルクHi-Cデータから細胞型特異的なTAD構造を抽出する
**Multi-band NMF** 手法の検証実験リポジトリである。

前身の `AIRE-experiments` では Machima2（scRNA+scATAC multiome 基底の Hi-C deconvolution）の
検証を行ったが、Hi-C 側の rank-1 表現の限界により TAD レベルの識別が困難であった。
Multi-band NMF は Hi-C データのみを入力とし、距離バンド分解 + H 共有の結合 NMF により
距離減衰の除去と TAD 構造の保持を両立する新手法であり、
AIRE-experiments での PoC 実験で有効性を確認した（mTEC comp: TAD recall 0.80, precision 0.91）。

本リポジトリではこの手法を **RIPHAMMER** パッケージ（`../RIPHAMMER`）として整備した上で、
高解像度検証・fold-in 解析・生物学的解釈を進める。

## RIPHAMMER パッケージ

Multi-band NMF の Python 実装。別リポジトリ `../RIPHAMMER`（`chiba-ai-med/RIPHAMMER`）で開発。

- 表示名: RIPHAMMER（大文字）
- pip / import: `riphammer`（小文字）
- 依存: numpy ベースの軽量実装

本リポジトリからは `pip install -e ../RIPHAMMER` で editable install して使う。

## Two-mode execution

| モード | 用途 | コマンド |
|--------|------|---------|
| **conda（開発）** | 日常の実験。RIPHAMMER を editable install、変更即反映 | `snakemake --use-conda` |
| **container（公開）** | 再現性保証。RIPHAMMER をバージョン固定で Docker に焼く | `snakemake --use-singularity` |

普段は conda モードで高速に回し、バージョンタグ（v0.1.0 等）を打つタイミングでのみ Docker を更新する。

## Dataset

### バルク Hi-C（本番 deconvolution 用）

マウス（mm10）バルク Hi-C、計 12 サンプル:

| グループ | サンプル数 | ソース |
|----------|-----------|-------|
| Thymocyte | 11 | 4D Nucleome Consortium |
| mTEC | 1 | GSE180937（SRA fastq → nf-core/hic → mcool） |

常染色体のみ（chr1–19）。性染色体は十字アーティファクトの原因となるため除外。
Brain (Bonev 2017, 3サンプル) は K=M=3 で trivial なため drop。

### scHi-C（pseudo-bulk validation 用）

| | GAGE-seq | Droplet Hi-C |
|---|---|---|
| 論文 | Zhou et al. 2024, Nat Genet | Chang et al. 2024, Nat Biotech |
| GEO | GSE238001 | GSE253407 |
| 組織 | mouse brain cortex | mouse brain cortex |
| 細胞数 | 9,047 | 6,211 |
| セルタイプ数 (K) | 3 | 20 |
| pairs 形式 | ヘッダなし5列, barcode=combinatorial index | 4DN standard 10列, barcode=16bp DNA |
| アノテーション方法 | RNA から TruncatedSVD+KMeans | 論文付属メタデータ |

## Validation strategy

1. **シミュレーション** (完了): RIPHAMMER vs THUNDER vs Standard NMF, K=4, 100kb
2. **scHi-C pseudo-bulk** (完了): scHi-C からセルタイプ別 pseudo-bulk → 既知比率で混合 (Dirichlet, M=20) → RIPHAMMER → 真の比率との MAE/Pearson r で評価。GAGE-seq (K=3, r=0.735), Droplet Hi-C (K=20, r=0.415), Droplet Hi-C merged (K=5, r=0.515)
3. **mTEC 本番** (未着手): Thymocyte 11 + mTEC 1 = M=12 で RIPHAMMER 実行、生物学的発見

THUNDER との比較は解像度の差で行う: THUNDER は Tucker 分解 (bins×bins×samples テンソル) のため 100kb (~27,000 bins) では計算コスト的に非現実的。論文ベンチマークも 10Mb のみ。

## Resolution staging

- **100kb** — 開発・全染色体。PoC 済み（AIRE-experiments で完了）
- **25kb** — sub-TAD / 境界精度検証
- **10kb** — Aire 座位ループレベル解析（chr10 限定）

## Repository layout

```
workflow/
  Snakefile              # 本番パイプライン (simulation + mTEC)
  scHiC_validation.smk   # scHi-C pseudo-bulk validation パイプライン
  config.yaml            # 本番設定
  config_scHiC.yaml      # scHi-C validation 設定
  envs/                  # conda env yaml
  run.sh / run_scHiC.sh
src/                     # 実験スクリプト (Python)
data/
  4dn_thymocyte/         # Thymocyte bulk Hi-C mcool
  mtec/                  # mTEC multiome + bulk Hi-C
  scHiC/
    gage_seq/            # GAGE-seq raw, metadata, pseudobulk, mixtures, results
    droplet_hic/         # Droplet Hi-C K=20 (同構成)
    droplet_hic_merged/  # Droplet Hi-C K=5 (セルタイプマージ版)
output/                  # simulation 結果
plot/  logs/  benchmarks/
notes/                   # 論文ドラフト
Dockerfile
```

## Environment & execution

Snakemake 駆動。スタイルは AIRE-experiments を踏襲:

- `min_version("8.10.0")`
- 全ルールに `conda:` と `container:` 両方
- ルールは `python src/<task>.py` を呼ぶ（R ではなく Python）
- ラッパー: `workflow/run.sh`

## 技術的注意事項

- **numpy==1.26.4 固定**: h5py, pandas, numexpr が numpy 1.x でビルドされているため 2.x 不可
- **pseudo-bulk メモリ**: 全 contacts をリスト保持すると OOM。Counter で on-the-fly ビン集計方式に修正済み
- **Droplet Hi-C pairs**: `#columns` ヘッダが `chrom1/chrom2` 表記 (4DN standard の `chr1/chr2` と異なる)
- **Droplet Hi-C barcode**: pairs = 16bp DNA のみ、メタデータ = lib:barcode 形式。`df["bc"]` 列を使用
- **共有サーバー注意**: 大量データ処理時は RAM 使用量を監視すること

## Migration from AIRE-experiments

AIRE-experiments から引き継ぐもの:
- `data/` 以下の Hi-C mcool ファイル（シンボリックリンクまたはコピー）
- `notes/multiband_nmf_summary.md` / `.pdf`（論文ドラフト）
- 実験結果の知見（CLAUDE.md とメモリに記録済み）

引き継がないもの:
- Machima2 関連コード・設定（R, Snakemake ルール、config.yaml）
- scRNA / scATAC multiome データ
- `.snakemake/` キャッシュ
