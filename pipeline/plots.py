import logging
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import numpy as np
import pandas as pd
import seaborn as sns
from collections import defaultdict
from pipeline.blast.blast_parser import gene_columns, primer_columns

logger = logging.getLogger(__name__)

#1 - plot the RAW BLAST results
def read_raw_blast(blast_path, is_primer = False):
    hits = defaultdict(list)
    columns = primer_columns if is_primer else gene_columns
    with Path(blast_path).open("r", encoding = "utf-8") as f:
        for line_number, line in enumerate(f, start = 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != len(columns):
                raise ValueError(f"{blast_path}:{line_number} expected {len(columns)} BLAST columns, got {len(fields)}")
            row = dict(zip(columns, fields))
            identity = float(row["pct_ident"])
            query_length = int(row["qlength"])
            query_start, query_end = int(row["q_start"]), int(row["q_end"])
            if query_length <= 0 or not 1 <= query_start <= query_length or not 1 <= query_end <= query_length:
                raise ValueError(f"{blast_path}:{line_number} has invalid query coordinates: qlength={query_length}, q_start={query_start}, q_end={query_end}")
            aligned_bases = abs(query_end - query_start) + 1
            coverage = aligned_bases / query_length * 100
            if not np.isfinite(identity) or not 0 <= identity <= 100:
                raise ValueError(f"{blast_path}:{line_number} has invalid identity: {identity}")
            if is_primer:
                unaligned_bases = query_length - aligned_bases
                total_mismatches = int(row["mismatches"]) + unaligned_bases + int(row["gaps"])
                hits[row["query_id"]].append((identity, coverage, total_mismatches))
            else:
                hits[row["query_id"]].append((identity, coverage))
    return hits

def plot_blast_query(hits, title, params, output_path, is_primer = False):
    fig = Figure(figsize = (15, 4.8), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 3)
    n_columns = 3 if is_primer else 2
    values = np.asarray(hits, dtype = float).reshape(-1, n_columns)
    identity, coverage = values[:, 0], values[:, 1]

    if is_primer:
        min_identity = None
        high_identity = None
        min_coverage = params["min_alignment_pct"]
        high_coverage = None
        total_mismatches = values[:, 2]
        passing = int(np.sum((coverage >= min_coverage) & (total_mismatches <= params["max_mismatches"])))
    else:
        min_identity = params["min_identity"]
        min_coverage = params["min_coverage"]
        high_identity = params.get("high_quality_identity", params.get("high_identity", min_identity))
        high_coverage = params.get("high_coverage", min_coverage)
        passing = int(np.sum((identity >= min_identity) & (coverage >= min_coverage)))

    thresholds = [(min_identity, high_identity), (min_coverage, high_coverage)]
    labels = ["BLAST identity (%)", "Query coverage (%)"]

    for ax, data, label, (minimum, high) in zip(axes[:2], [identity, coverage], labels, thresholds):
        ax.hist(data, bins = np.arange(0, 101, 1), color = "#9da6cf", edgecolor = "white", linewidth = 0.3)
        if minimum is not None:
            ax.axvline(minimum, color = "#d18c20", linestyle = "--", label = f"Minimum: {minimum:g}%")
        if high is not None and high != minimum:
            ax.axvline(high, color = "#d94c4c", linestyle = ":", label = f"High quality: {high:g}%")
        ax.set(xlabel = label, ylabel = "BLAST hits", xlim = (0, 102))
        if minimum is not None or high is not None:
            ax.legend(fontsize = 8)
        ax.grid(axis = "y", alpha = 0.2)

    ax = axes[2]
    if len(values):
        counts = ax.hexbin(identity, coverage, gridsize = 45, extent = (0, 100, 0, 100), mincnt = 1, bins = "log", cmap = "viridis")
        fig.colorbar(counts, ax = ax, label = "BLAST hits (log scale)")
    if min_identity is not None:
        ax.axvline(min_identity, color = "#d18c20", linestyle = "--")
    ax.axhline(min_coverage, color = "#d18c20", linestyle = "--")
    if high_identity is not None and high_identity != min_identity:
        ax.axvline(high_identity, color = "#d94c4c", linestyle = ":")
    if high_coverage is not None and high_coverage != min_coverage:
        ax.axhline(high_coverage, color = "#d94c4c", linestyle = ":")
    ax.set(xlabel = labels[0], ylabel = labels[1], xlim = (-2, 102), ylim = (-2, 102))

    if is_primer:
        ax.set_title(
            f"BLAST hits passing filters: {passing:,}/{len(values):,}\n"
            f"Coverage ≥ {min_coverage:g}%; total mismatches ≤ {params['max_mismatches']}",
            fontsize = 10
        )
    else:
        ax.set_title(f"BLAST hits passing both thresholds: {passing:,}/{len(values):,}", fontsize = 10)

    if not len(values):
        for ax in axes:
            ax.text(0.5, 0.5, "No BLAST hits", transform = ax.transAxes, ha = "center")
    fig.suptitle(f"{title} - raw BLAST hits (n={len(values):,})", fontsize = 13)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved raw BLAST plot: %s", output_path)
    return output_path

def plot_raw_blast_results(blast_outputs, config, output_dir):
    #keep every reported HSP and plot each ref query separately
    gene_name = config.gene_name
    gene_queries = ["ref_full"] if config.repeat_location.lower() == "promoter" else ["ref_2_repeats", "ref_no_repeats"]
    searches = [
        ("gene", gene_queries + ["upstream_flank", "downstream_flank"]),
        ("gene_start", ["ref_gene_start"]),
        ("gene_end", ["ref_gene_end"]),
        ("primers", ["forward_primer", "reverse_primer"])
    ]
    outputs = []
    for search, queries in searches:
        blast_path = blast_outputs.get(search)
        if blast_path is None:
            continue
        if search == "gene_start":
            params = config.gene_start_blast_params
        elif search == "primers":
            params = config.primer_blast_params
        else:
            #gene-end parsing currently uses gene_blast_params, as do genes and flanks
            params = config.gene_blast_params
        is_primer = search == "primers"
        raw_hits = read_raw_blast(blast_path, is_primer = is_primer)
        for query in queries:
            query_id = f"{gene_name}_{query}"
            output_path = Path(output_dir) / f"{query_id}_raw_blast.png"
            outputs.append(plot_blast_query(
                raw_hits.get(query_id, []), f"{gene_name}: {query}", params,
                output_path, is_primer = is_primer
            ))
    return outputs

#2 - Plot primer mismatches
def plot_primer_mismatch(blast_outputs, config, output_dir):
    blast_path = blast_outputs.get("primers")
    if blast_path is None:
        return None
    raw_hits = read_raw_blast(blast_path, is_primer = True)
    params = config.primer_blast_params
    min_coverage = params["min_alignment_pct"]
    max_mismatches = params["max_mismatches"]
    primer_names = ["forward_primer", "reverse_primer"]
    primer_data = {}

    for primer in primer_names:
        query_id = f"{config.gene_name}_{primer}"
        values = np.asarray(raw_hits.get(query_id, []), dtype = float).reshape(-1, 3)
        mismatches = values[values[:, 1] >= min_coverage, 2].astype(int)
        primer_data[primer] = (len(values), mismatches)

    largest_mismatch = max(
        [max_mismatches + 1] +
        [int(mismatches.max()) for _, mismatches in primer_data.values() if len(mismatches)]
    )
    mismatch_counts = np.arange(largest_mismatch + 1)
    fig = Figure(figsize = (12, 4.8), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)

    for ax, primer in zip(axes, primer_names):
        raw_count, mismatches = primer_data[primer]
        counts = np.bincount(mismatches, minlength = largest_mismatch + 1)
        passing = int(np.sum(mismatches <= max_mismatches))
        colours = ["#9da6cf" if count <= max_mismatches else "#dddddd" for count in mismatch_counts]
        ax.bar(mismatch_counts, counts, color = colours, edgecolor = "white", width = 0.85)

        for mismatch, count in zip(mismatch_counts, counts):
            if count:
                ax.annotate(
                    f"{count:,}", (mismatch, count),
                    xytext = (0, 4), textcoords = "offset points",
                    ha = "center", fontsize = 9
                )

        ax.axvline(
            max_mismatches + 0.5, color = "#d18c20", linestyle = "--",
            label = f"Maximum total mismatches: {max_mismatches}"
        )
        ax.set(
            xlabel = "Total mismatches",
            ylabel = "BLAST hits meeting coverage threshold",
            xlim = (-0.5, largest_mismatch + 0.5)
        )
        ax.set_xticks(mismatch_counts)
        ax.set_title(
            f"{primer.replace('_', ' ').capitalize()}\n"
            f"Coverage retained: {len(mismatches):,}/{raw_count:,}; passing mismatch filter: {passing:,}",
            fontsize = 10
        )
        if not len(mismatches):
            ax.text(
                0.5, 0.5, "No hits meet coverage threshold",
                transform = ax.transAxes, ha = "center"
            )
        ax.grid(axis = "y", alpha = 0.2)
        ax.set_axisbelow(True)
        ax.legend(fontsize = 8)

    upper_y = max(ax.get_ylim()[1] for ax in axes)
    for ax in axes:
        ax.set_ylim(0, max(1, upper_y * 1.15))

    fig.suptitle(
        f"{config.gene_name} - primer mismatches among hits with coverage ≥ {min_coverage:g}%",
        fontsize = 13
    )
    output_path = Path(output_dir) / f"{config.gene_name}_primer_mismatches.png"
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved primer mismatch plot: %s", output_path)
    return output_path