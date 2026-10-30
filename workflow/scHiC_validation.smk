"""scHi-C pseudo-bulk validation of RIPHAMMER.

Downloads GAGE-seq and Droplet Hi-C single-cell Hi-C data,
generates per-cell-type pseudo-bulk contact matrices, creates
synthetic bulk mixtures with known fractions, runs RIPHAMMER,
and evaluates deconvolution accuracy.
"""
from snakemake.utils import min_version
min_version("8.10.0")

CONDA_ENV = "workflow/envs/riphammer.yaml"

configfile: "workflow/config_scHiC.yaml"

DATASETS = list(config["scHiC_datasets"].keys())
RESOLUTION = config["resolution"]


# ── helpers ──────────────────────────────────────────────────────
def geo_ftp_sample(base, gsm, filename):
    """Build per-sample GEO FTP URL."""
    prefix = gsm[:7] + "nnn"  # e.g., GSM8020nnn
    return f"{base}/{prefix}/{gsm}/suppl/{filename}"


def geo_ftp_series(dataset, filename):
    """Build series-level GEO FTP URL."""
    geo = config["scHiC_datasets"][dataset]["geo"]
    prefix = geo[:5] + "nnn"
    return (f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
            f"{prefix}/{geo}/suppl/{filename}")


def resolve_download_url(dataset, filename):
    """Resolve the download URL for a given file."""
    base = config["scHiC_datasets"][dataset]["ftp_base"]
    libs = config["scHiC_datasets"][dataset]["libraries"]
    for k, lib in libs.items():
        for field in ["pairs", "barcodes", "rna"]:
            if lib.get(field) == filename:
                gsm_field = {"pairs": "pairs_gsm", "barcodes": "gsm",
                             "rna": "rna_gsm"}.get(field, "gsm")
                gsm = lib.get(gsm_field, lib.get("gsm", ""))
                if gsm:
                    return geo_ftp_sample(base, gsm, filename)
    return geo_ftp_series(dataset, filename)


def all_pairs_files(dataset):
    libs = config["scHiC_datasets"][dataset]["libraries"]
    return [f"data/scHiC/{dataset}/raw/{libs[k]['pairs']}" for k in libs]


def all_barcode_files(dataset):
    libs = config["scHiC_datasets"][dataset]["libraries"]
    out = []
    for k in libs:
        if "barcodes" in libs[k]:
            out.append(f"data/scHiC/{dataset}/raw/{libs[k]['barcodes']}")
    return out


def all_rna_files(dataset):
    libs = config["scHiC_datasets"][dataset]["libraries"]
    out = []
    for k in libs:
        if "rna" in libs[k]:
            out.append(f"data/scHiC/{dataset}/raw/{libs[k]['rna']}")
    return out


# ── targets ──────────────────────────────────────────────────────
rule all:
    """Final targets: evaluation results and plots for both datasets."""
    input:
        expand("output/scHiC_validation/{ds}/evaluation.tsv", ds=DATASETS),
        expand("plot/scHiC_validation_{ds}.png", ds=DATASETS),
        "plot/scHiC_validation_combined.png",


# ══════════════════════════════════════════════════════════════════
# Phase 0: Download
# ══════════════════════════════════════════════════════════════════
localrules: download_file

rule download_file:
    """Download a single file from GEO FTP."""
    output:
        "data/scHiC/{dataset}/raw/{filename}"
    params:
        url=lambda wc: resolve_download_url(wc.dataset, wc.filename),
    log:
        "logs/scHiC/download_{dataset}_{filename}.log"
    shell:
        """
        wget -c -nv -O {output} {params.url} >& {log}
        """


rule download_droplet_metadata:
    """Download Droplet Hi-C cell-type metadata."""
    output:
        "data/scHiC/droplet_hic/metadata/cell_metadata.xls.gz"
    params:
        url=lambda wc: geo_ftp_series(
            "droplet_hic",
            config["scHiC_datasets"]["droplet_hic"]["metadata_file"]),
    log:
        "logs/scHiC/download_droplet_metadata.log"
    shell:
        """
        wget -c -nv -O {output} {params.url} >& {log}
        """


# ══════════════════════════════════════════════════════════════════
# Phase 1: Inspect pairs format (one-off, produces format spec)
# ══════════════════════════════════════════════════════════════════
rule inspect_pairs:
    """Inspect pairs file header to determine column layout."""
    input:
        pairs=lambda wc: all_pairs_files(wc.dataset)[0],
    output:
        spec="data/scHiC/{dataset}/metadata/pairs_format.json"
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/inspect_pairs_{dataset}.log"
    benchmark:
        "benchmarks/scHiC/inspect_pairs_{dataset}.txt"
    shell:
        """
        python src/scHiC_inspect_pairs.py \
            --pairs {input.pairs} \
            --output {output.spec} \
            >& {log}
        """


# ══════════════════════════════════════════════════════════════════
# Phase 2: Cell-type annotation
# ══════════════════════════════════════════════════════════════════
rule annotate_gage_seq:
    """Derive cell-type annotations from GAGE-seq RNA component."""
    input:
        rna=all_rna_files("gage_seq"),
        pairs_spec="data/scHiC/gage_seq/metadata/pairs_format.json",
    output:
        annotations="data/scHiC/gage_seq/metadata/cell_annotations.tsv"
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/annotate_gage_seq.log"
    benchmark:
        "benchmarks/scHiC/annotate_gage_seq.txt"
    shell:
        """
        python src/scHiC_annotate_gage.py \
            --rna-files {input.rna} \
            --output {output.annotations} \
            >& {log}
        """


rule annotate_droplet_hic:
    """Parse Droplet Hi-C metadata for cell-type annotations."""
    input:
        metadata="data/scHiC/droplet_hic/metadata/cell_metadata.xls.gz",
    output:
        annotations="data/scHiC/droplet_hic/metadata/cell_annotations.tsv"
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/annotate_droplet_hic.log"
    benchmark:
        "benchmarks/scHiC/annotate_droplet_hic.txt"
    shell:
        """
        python src/scHiC_annotate_droplet.py \
            --metadata {input.metadata} \
            --output {output.annotations} \
            >& {log}
        """


# ══════════════════════════════════════════════════════════════════
# Phase 3: Pseudo-bulk generation
# ══════════════════════════════════════════════════════════════════
rule pseudobulk:
    """Aggregate single-cell contacts by cell type into pseudo-bulk mcool."""
    input:
        pairs=lambda wc: all_pairs_files(wc.dataset),
        annotations="data/scHiC/{dataset}/metadata/cell_annotations.tsv",
        pairs_spec="data/scHiC/{dataset}/metadata/pairs_format.json",
    output:
        summary="data/scHiC/{dataset}/pseudobulk/summary.tsv",
        flag=touch("data/scHiC/{dataset}/pseudobulk/.done"),
    params:
        outdir="data/scHiC/{dataset}/pseudobulk",
        resolution=RESOLUTION,
        chromsizes=config["chromsizes"],
        min_cells=config["pseudobulk"]["min_cells_per_type"],
        min_contacts=config["pseudobulk"]["min_contacts_per_type"],
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/pseudobulk_{dataset}.log"
    benchmark:
        "benchmarks/scHiC/pseudobulk_{dataset}.txt"
    threads: 4
    shell:
        """
        python src/scHiC_pseudobulk.py \
            --pairs-files {input.pairs} \
            --annotations {input.annotations} \
            --pairs-spec {input.pairs_spec} \
            --outdir {params.outdir} \
            --resolution {params.resolution} \
            --chromsizes {params.chromsizes} \
            --min-cells {params.min_cells} \
            --min-contacts {params.min_contacts} \
            --threads {threads} \
            >& {log}
        """


# ══════════════════════════════════════════════════════════════════
# Phase 4: Synthetic mixture generation
# ══════════════════════════════════════════════════════════════════
rule create_mixtures:
    """Create synthetic bulk mixtures from pseudo-bulk with known fractions."""
    input:
        summary="data/scHiC/{dataset}/pseudobulk/summary.tsv",
        flag="data/scHiC/{dataset}/pseudobulk/.done",
    output:
        mixtures="data/scHiC/mixtures/{dataset}/mixtures.mcool",
        fractions="data/scHiC/mixtures/{dataset}/true_fractions.tsv",
        sample_names="data/scHiC/mixtures/{dataset}/cells.tsv",
    params:
        pseudobulk_dir="data/scHiC/{dataset}/pseudobulk",
        n_mixtures=config["mixture"]["n_mixtures"],
        seed=config["mixture"]["seed"],
        resolution=RESOLUTION,
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/create_mixtures_{dataset}.log"
    benchmark:
        "benchmarks/scHiC/create_mixtures_{dataset}.txt"
    shell:
        """
        python src/scHiC_create_mixtures.py \
            --pseudobulk-dir {params.pseudobulk_dir} \
            --summary {input.summary} \
            --n-mixtures {params.n_mixtures} \
            --seed {params.seed} \
            --resolution {params.resolution} \
            --output-mcool {output.mixtures} \
            --output-fractions {output.fractions} \
            --output-samples {output.sample_names} \
            >& {log}
        """


# ══════════════════════════════════════════════════════════════════
# Phase 5: Vectorize + RIPHAMMER NMF
# ══════════════════════════════════════════════════════════════════
rule vectorize_mixtures:
    """Vectorize synthetic mixture mcool into feature matrix."""
    input:
        mixtures="data/scHiC/mixtures/{dataset}/mixtures.mcool",
        sample_names="data/scHiC/mixtures/{dataset}/cells.tsv",
    output:
        X="output/scHiC_validation/{dataset}/features/X_bulk.npy",
        features="output/scHiC_validation/{dataset}/features/features.tsv",
        cells="output/scHiC_validation/{dataset}/features/cells.tsv",
    params:
        outdir="output/scHiC_validation/{dataset}/features",
        resolution=RESOLUTION,
        max_dist=config["max_dist_bp"],
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/vectorize_mixtures_{dataset}.log"
    benchmark:
        "benchmarks/scHiC/vectorize_mixtures_{dataset}.txt"
    shell:
        """
        python src/scHiC_vectorize_mixtures.py \
            --mixtures-mcool {input.mixtures} \
            --sample-names {input.sample_names} \
            --resolution {params.resolution} \
            --max-dist {params.max_dist} \
            --outdir {params.outdir} \
            >& {log}
        """


rule fit_riphammer:
    """Run RIPHAMMER NMF on synthetic mixtures."""
    input:
        X="output/scHiC_validation/{dataset}/features/X_bulk.npy",
        features="output/scHiC_validation/{dataset}/features/features.tsv",
        cells="output/scHiC_validation/{dataset}/features/cells.tsv",
        fractions="data/scHiC/mixtures/{dataset}/true_fractions.tsv",
    output:
        H="output/scHiC_validation/{dataset}/nmf/H.tsv",
        W="output/scHiC_validation/{dataset}/nmf/W.tsv",
        meta="output/scHiC_validation/{dataset}/nmf/meta.json",
    params:
        outdir="output/scHiC_validation/{dataset}/nmf",
        max_iter=config["nmf"]["max_iter"],
        seed=config["nmf"]["seed"],
        mode=config["nmf"]["mode"],
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/fit_riphammer_{dataset}.log"
    benchmark:
        "benchmarks/scHiC/fit_riphammer_{dataset}.txt"
    shell:
        """
        python src/scHiC_fit.py \
            --X {input.X} \
            --features {input.features} \
            --cells {input.cells} \
            --true-fractions {input.fractions} \
            --outdir {params.outdir} \
            --max-iter {params.max_iter} \
            --seed {params.seed} \
            --mode {params.mode} \
            >& {log}
        """


# ══════════════════════════════════════════════════════════════════
# Phase 6: Evaluation
# ══════════════════════════════════════════════════════════════════
rule evaluate:
    """Evaluate RIPHAMMER recovery against known fractions."""
    input:
        H="output/scHiC_validation/{dataset}/nmf/H.tsv",
        fractions="data/scHiC/mixtures/{dataset}/true_fractions.tsv",
        meta="output/scHiC_validation/{dataset}/nmf/meta.json",
    output:
        evaluation="output/scHiC_validation/{dataset}/evaluation.tsv",
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/evaluate_{dataset}.log"
    benchmark:
        "benchmarks/scHiC/evaluate_{dataset}.txt"
    shell:
        """
        python src/scHiC_evaluate.py \
            --H {input.H} \
            --true-fractions {input.fractions} \
            --meta {input.meta} \
            --output {output.evaluation} \
            >& {log}
        """


# ══════════════════════════════════════════════════════════════════
# Phase 7: Figures
# ══════════════════════════════════════════════════════════════════
rule plot_per_dataset:
    """Per-dataset validation plot."""
    input:
        evaluation="output/scHiC_validation/{dataset}/evaluation.tsv",
        H="output/scHiC_validation/{dataset}/nmf/H.tsv",
        fractions="data/scHiC/mixtures/{dataset}/true_fractions.tsv",
    output:
        "plot/scHiC_validation_{dataset}.png"
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/plot_{dataset}.log"
    shell:
        """
        python src/plot_scHiC_validation.py \
            --evaluation {input.evaluation} \
            --H {input.H} \
            --true-fractions {input.fractions} \
            --dataset {wildcards.dataset} \
            --output {output} \
            >& {log}
        """


rule plot_combined:
    """Combined validation plot across all datasets."""
    input:
        evaluations=expand(
            "output/scHiC_validation/{ds}/evaluation.tsv", ds=DATASETS),
    output:
        "plot/scHiC_validation_combined.png"
    params:
        datasets=DATASETS,
    conda:
        CONDA_ENV
    log:
        "logs/scHiC/plot_combined.log"
    shell:
        """
        python src/plot_scHiC_validation_combined.py \
            --evaluations {input.evaluations} \
            --datasets {params.datasets} \
            --output {output} \
            >& {log}
        """
