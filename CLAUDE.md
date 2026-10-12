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

マウス（mm10）バルク Hi-C、計 15 サンプル:

| グループ | サンプル数 | ソース |
|----------|-----------|-------|
| Thymocyte | 11 | 4D Nucleome Consortium |
| mTEC | 1 | GSE180937（SRA fastq → nf-core/hic → mcool） |
| Brain | 3 | Bonev et al. (2017) GSE96107 |

常染色体のみ（chr1–19）。性染色体は十字アーティファクトの原因となるため除外。

## Validation strategy

1. **Cell-type separation**: H 行列の loading でサンプル–コンポーネント対応を確認
2. **TAD recovery**: NMF 基底から検出した TAD 境界と実測バルク Hi-C TAD 境界の比較（recall, precision, insulation 相関）
3. **Biological plausibility**: Aire 遺伝子座（chr10:21.9Mb）の細胞型特異的 TAD 構造
4. **In-silico mixing**: fold-in で混合比率推定精度を検証（O/E NMF で確認済み、Multi-band NMF では未実施）

## Resolution staging

- **100kb** — 開発・全染色体。PoC 済み（AIRE-experiments で完了）
- **25kb** — sub-TAD / 境界精度検証
- **10kb** — Aire 座位ループレベル解析（chr10 限定）

## Repository layout

```
workflow/
  Snakefile
  envs/             # conda env yaml（riphammer を editable install）
  config.yaml
  run.sh
src/                 # 実験スクリプト（import riphammer して使う）
data/                # Hi-C データ（gitignore）
output/              # 結果
plot/  logs/  benchmarks/
notes/               # 論文ドラフト・サマリー
Dockerfile
```

## Environment & execution

Snakemake 駆動。スタイルは AIRE-experiments を踏襲:

- `min_version("8.10.0")`
- 全ルールに `conda:` と `container:` 両方
- ルールは `python src/<task>.py` を呼ぶ（R ではなく Python）
- ラッパー: `workflow/run.sh`

## Migration from AIRE-experiments

AIRE-experiments から引き継ぐもの:
- `data/` 以下の Hi-C mcool ファイル（シンボリックリンクまたはコピー）
- `notes/multiband_nmf_summary.md` / `.pdf`（論文ドラフト）
- 実験結果の知見（CLAUDE.md とメモリに記録済み）

引き継がないもの:
- Machima2 関連コード・設定（R, Snakemake ルール、config.yaml）
- scRNA / scATAC multiome データ
- `.snakemake/` キャッシュ
