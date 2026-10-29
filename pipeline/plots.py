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

#3 - Plot sequence extraction QC before and after N filtering
def plot_sequence_extraction(before_sequences, after_sequences, config, output_dir):
    if not before_sequences:
        logger.info("No extracted sequences available for plotting")
        return None

    gene_name = config.gene_name
    params = config.extraction_params
    promoter_mode = config.repeat_location.lower() == "promoter"
    requested = int(params["promoter_length"] if promoter_mode else params["flanking_bp"])

    fig = Figure(figsize = (15, 9), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(2, 3, sharex = "col", sharey = "col")

    #use identical bins for both rows
    sequence_lengths = [len(seq["sequence"]) for seq in before_sequences]
    length_bins = np.histogram_bin_edges(sequence_lengths, bins = 30)
    context_bins = np.linspace(-0.5, max(requested, 1) + 0.5, 31)
    boundary_max = 1

    for row, (stage, sequences) in enumerate([
        ("Before N filtering", before_sequences),
        ("After N filtering", after_sequences)
    ]):
        boundary_ax, context_ax, length_ax = axes[row]
        n_sequences = len(sequences)

        #boundary status based on BLAST hit coordinates
        statuses = ["OK", "WARNING", "CRITICAL"]
        counts = [sum(seq["boundary_status"] == status for seq in sequences) for status in statuses]
        boundary_max = max(boundary_max, max(counts))
        bars = boundary_ax.bar(
            statuses, counts,
            color = ["#9da6cf", "#d18c20", "#d94c4c"],
            edgecolor = "white"
        )

        for bar, count in zip(bars, counts):
            percentage = count / n_sequences * 100 if n_sequences else 0
            boundary_ax.annotate(
                f"{count:,}\n({percentage:.1f}%)",
                (bar.get_x() + bar.get_width() / 2, count),
                xytext = (0, 4), textcoords = "offset points",
                ha = "center", va = "bottom", fontsize = 9
            )

        boundary_ax.set_title(f"{stage}: {n_sequences:,} loci\nContig boundary status", fontsize = 10)
        boundary_ax.set_ylabel("Loci")
        boundary_ax.set_xlabel("Based on BLAST hit boundaries")

        if not sequences:
            for ax in (context_ax, length_ax):
                ax.text(
                    0.5, 0.5, "No sequences retained",
                    transform = ax.transAxes, ha = "center", va = "center"
                )
            context_ax.set_title("Promoter recovered" if promoter_mode else "Flanks recovered around BLAST hit", fontsize = 10)
            length_ax.set_title("BLAST-defined gene length" if promoter_mode else "Total extracted sequence length", fontsize = 10)
            continue

        #sequence context recovered
        if promoter_mode:
            promoter_lengths = [seq["promoter_length_extracted"] for seq in sequences]
            complete = sum(length == requested for length in promoter_lengths)

            context_ax.hist(
                promoter_lengths, bins = context_bins,
                color = "#9da6cf", edgecolor = "white", linewidth = 0.3
            )
            context_ax.set_title(
                f"Promoter recovered\nFull length: {complete:,}/{n_sequences:,} loci",
                fontsize = 10
            )
        else:
            upstream_lengths = []
            downstream_lengths = []

            for seq in sequences:
                left_flank = seq["gene_start"] - seq["extract_start"]
                right_flank = seq["extract_end"] - seq["gene_end"]

                if seq["strand"] == "forward":
                    upstream_lengths.append(left_flank)
                    downstream_lengths.append(right_flank)
                elif seq["strand"] == "reverse":
                    upstream_lengths.append(right_flank)
                    downstream_lengths.append(left_flank)
                else:
                    raise ValueError(f"Invalid strand: {seq['strand']}")

            upstream_complete = sum(length == requested for length in upstream_lengths)
            downstream_complete = sum(length == requested for length in downstream_lengths)

            context_ax.hist(
                [upstream_lengths, downstream_lengths],
                bins = context_bins,
                color = ["#9da6cf", "#70b6a1"],
                edgecolor = "white", linewidth = 0.3,
                label = [
                    f"Upstream: {upstream_complete:,} full length",
                    f"Downstream: {downstream_complete:,} full length"
                ]
            )
            context_ax.set_title("Flanks recovered around BLAST hit", fontsize = 10)

        context_ax.axvline(
            requested, color = "#d18c20", linestyle = "--",
            label = f"Requested: {requested:,} bp"
        )
        context_ax.legend(fontsize = 8)

        #total extracted sequence length; promoter is stored separately
        lengths = [len(seq["sequence"]) for seq in sequences]
        median_length = np.median(lengths)

        length_ax.hist(
            lengths, bins = length_bins,
            color = "#9da6cf", edgecolor = "white", linewidth = 0.3
        )
        length_ax.axvline(
            median_length, color = "#d18c20", linestyle = "--",
            label = f"Median: {median_length:,.0f} bp"
        )
        length_ax.set_title(
            "BLAST-defined gene length" if promoter_mode else "Total extracted sequence length",
            fontsize = 10
        )
        length_ax.legend(fontsize = 8)

    for row_axes in axes:
        row_axes[0].set_ylim(0, boundary_max * 1.35)
        row_axes[1].set_xlabel("Extracted promoter length (bp)" if promoter_mode else "Flanking sequence length (bp)")
        row_axes[1].set_ylabel("Loci" if promoter_mode else "Loci per direction")
        row_axes[2].set_xlabel("Sequence length (bp)")
        row_axes[2].set_ylabel("Loci")

        for ax in row_axes:
            ax.grid(axis = "y", alpha = 0.2)
            ax.set_axisbelow(True)
            ax.tick_params(axis = "x", labelbottom = True)

    removed = len(before_sequences) - len(after_sequences)
    fig.suptitle(
        f"{gene_name} - sequence extraction QC\n"
        f"N filtering removed {removed:,}/{len(before_sequences):,} loci "
        f"({removed / len(before_sequences):.1%}); "
        f"boundary thresholds: critical ≤ {params['boundary_critical_distance']} bp, "
        f"warning ≤ {params['boundary_warning_distance']} bp",
        fontsize = 13
    )
    output_path = Path(output_dir) / f"{gene_name}_sequence_extraction_QC.png"
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved sequence extraction QC plot: %s", output_path)
    return output_path

#4 - Plot BLAST hit positions relative to contig boundaries
def plot_contig_boundaries(sequences, config, output_dir):
    if not sequences:
        logger.info("No extracted sequences available for contig boundary plotting")
        return None

    critical = config.extraction_params["boundary_critical_distance"]
    warning = config.extraction_params["boundary_warning_distance"]
    if not 0 <= critical < warning:
        raise ValueError("Boundary thresholds must satisfy 0 <= critical < warning")

    start_distances = np.asarray([s["dist_to_contig_start"] for s in sequences])
    end_distances = np.asarray([s["dist_to_contig_end"] for s in sequences])
    nearest = np.minimum(start_distances, end_distances)

    fig = Figure(figsize = (15, 5.5), layout = "constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)

    #scatter of distances to both contig edges
    groups = [
        ("OK", nearest > warning, "#70b6a1", f"OK: >{warning:g} bp"),
        ("WARNING", (nearest > critical) & (nearest <= warning), "#e9c46a", f"Warning: >{critical:g}–{warning:g} bp"),
        ("CRITICAL", nearest <= critical, "#d94c4c", f"Critical: ≤{critical:g} bp")
    ]
    for status, mask, colour, label in groups:
        axes[0].scatter(
            start_distances[mask], end_distances[mask],
            s = 30, color = colour, edgecolors = "#555555",
            linewidths = 0.4, alpha = 0.6,
            label = f"{label} (n={int(mask.sum()):,})"
        )

    for threshold, colour in [(critical, "#d94c4c"), (warning, "#d18c20")]:
        axes[0].axvline(threshold, color = colour, linestyle = "--", linewidth = 1)
        axes[0].axhline(threshold, color = colour, linestyle = "--", linewidth = 1)

    axes[0].set_xscale("symlog", linthresh = 1)
    axes[0].set_yscale("symlog", linthresh = 1)
    axes[0].set(
        xlabel = "Distance from BLAST hit to contig start (bp)",
        ylabel = "Distance from BLAST hit to contig end (bp)",
        title = "BLAST hit position relative to contig boundaries"
    )
    upper_distance = max(float(start_distances.max()), float(end_distances.max()), warning, 1) * 1.3
    axes[0].set_xlim(-0.2, upper_distance)
    axes[0].set_ylim(-0.2, upper_distance)
    axes[0].legend(fontsize = 8)

    #non-overlapping distance bins, adjusted to configured thresholds
    upper_bounds = [critical, warning]
    upper_bounds.extend(bound for bound in [500, 1000, 5000, 10000] if bound > warning)
    counts = []
    labels = []
    lower = None

    for upper in upper_bounds:
        if lower is None:
            mask = nearest <= upper
            label = f"0–{upper:g}\n(Critical)"
        else:
            mask = (nearest > lower) & (nearest <= upper)
            label = f">{lower:g}–{upper:g}"
            if upper == warning:
                label += "\n(Warning)"
        counts.append(int(mask.sum()))
        labels.append(label)
        lower = upper

    counts.append(int(np.sum(nearest > upper_bounds[-1])))
    labels.append(f">{upper_bounds[-1]:g}")
    colours = ["#d94c4c", "#e9c46a"] + ["#9da6cf"] * (len(counts) - 2)
    bars = axes[1].bar(
        np.arange(len(counts)), counts,
        color = colours, edgecolor = "white"
    )
    for bar, count in zip(bars, counts):
        axes[1].annotate(
            f"{count:,}", (bar.get_x() + bar.get_width() / 2, count),
            xytext = (0, 4), textcoords = "offset points",
            ha = "center", fontsize = 9
        )

    axes[1].set_xticks(np.arange(len(labels)))
    axes[1].set_xticklabels(labels, fontsize = 8)
    axes[1].set(
        xlabel = "Distance to nearest contig edge (bp)",
        ylabel = "Loci",
        title = "Distance to nearest contig edge",
        ylim = (0, max(1, max(counts) * 1.2))
    )

    for ax in axes:
        ax.grid(alpha = 0.2)
        ax.set_axisbelow(True)

    fig.suptitle(
        f"{config.gene_name} - contig boundary QC (n={len(sequences):,} loci)",
        fontsize = 13
    )
    output_path = Path(output_dir) / f"{config.gene_name}_contig_boundary_QC.png"
    output_path.parent.mkdir(parents = True, exist_ok = True)
    fig.savefig(output_path, dpi = 300)
    fig.clear()
    logger.info("Saved contig boundary QC plot: %s", output_path)
    return output_path

#one function to make all plots

def generate_all_plots(config, blast_outputs, extracted_sequences, before_N_filtering, output_dir):
    plots_generated = {
        "raw_blast": plot_raw_blast_results(blast_outputs, config, output_dir),
        "primer_mismatch": plot_primer_mismatch(blast_outputs, config, output_dir),
        "sequence_extraction": plot_sequence_extraction(before_N_filtering, extracted_sequences, config, output_dir),
        "contig_boundary_QC": plot_contig_boundaries(extracted_sequences, config, output_dir)
    }

    return plots_generated
