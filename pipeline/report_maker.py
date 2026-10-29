import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any
import pandas as pd

logger = logging.getLogger(__name__)

def write_df_with_spacing(df, writer, sheet_name, current_row, blank_rows_after = 2):
    #helper function to have spaces between tables on summary sheets
    df.to_excel(writer, sheet_name = sheet_name, index = False, startrow = current_row)

    return current_row + len(df) + 1 + blank_rows_after

def format_workbook(writer):
    for worksheet in writer.book.worksheets:
        worksheet.freeze_panes = "A2"
        worksheet.auto_filter.ref = worksheet.dimensions

        for column_cells in worksheet.columns:
            values = [str(cell.value) for cell in column_cells if cell.value is not None]
            width = min(max((len(value) for value in values), default = 10) + 2, 60)
            worksheet.column_dimensions[column_cells[0].column_letter].width = width

    
# 1 - BLAST reports

identity_bins = [
        ("<50%", lambda value: value < 50),
        ("50-60%", lambda value: 50 <= value < 60),
        ("60-70%", lambda value: 60 <= value < 70),
        ("70-80%", lambda value: 70 <= value < 80),
        ("80-90%", lambda value: 80 <= value < 90),
        ("90-95%", lambda value: 90 <= value < 95),
        ("95-98%", lambda value: 95 <= value < 98),
        ("98-100%", lambda value: 98 <= value <= 100),
        ("100%", lambda value: value == 100)
    ]
coverage_bins = [
        ("<50%", lambda value: value < 50),
        ("50-60%", lambda value: 50 <= value < 60),
        ("60-70%", lambda value: 60 <= value < 70),
        ("70-80%", lambda value: 70 <= value < 80),
        ("80-90%", lambda value: 80 <= value < 90),
        ("90-95%", lambda value: 90 <= value < 95),
        ("95-98%", lambda value: 95 <= value < 98),
        ("98-100%", lambda value: 98 <= value <= 100),
        ("100%", lambda value: value == 100)
    ]

def count_bins(values, bins):
        return[f"{label}: {sum(predicate(value) for value in values)}"
               for label, predicate in bins]
               
def report_gene_hit_stats(blast_path, output_file):
    #summarise each hit result from raw BLAST output
    blast_path = Path(blast_path)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    columns = ["query_id", "subject_id", "nident", "pident", "length", "mismatch", "gaps", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]

    rows = []
    excluded_flank_rows = 0
    with blast_path.open("r", encoding = "utf-8") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if len(fields) == len(columns):
                    row: dict[str, Any] = dict(zip(columns, fields))
                    if row["query_id"].endswith("_flank"):
                        excluded_flank_rows += 1
                        continue
                    row["identity"] = float(row["pident"])
                    row["coverage"] = abs(int(row["qend"]) - int(row["qstart"]) + 1) / int(row["qlen"]) * 100
                    rows.append(row)


    identities = [row["identity"] for row in rows]
    coverages = [row["coverage"] for row in rows]
    isolates = {row["subject_id"].split("_ctg", 1)[0] for row in rows}
    query_counts = defaultdict(int)
    for row in rows:
        query_id = row["query_id"]
        if query_id.endswith("_upstream_flank"):
            query_type = "upstream flank"
        elif query_id.endswith("_downstream_flank"):
            query_type = "downstream flank"
        elif "ref_2_repeats" in query_id:
            query_type = "2 repeats"
        elif "ref_no_repeats" in query_id:
            query_type = "no repeats"
        elif "ref_full" in query_id:
            query_type = "full"
        else:
            query_type = "other"
        query_counts[query_type] += 1

    report_lines = [
        "Raw gene BLAST hit statistics",
        f"Source file: {blast_path}",
        f"Raw total hits: {len(rows)}",
        f"Excluded flank hits: {excluded_flank_rows}",
        f"Isolates: {len(isolates)}",
        "",
        "Idendity breakdown:",
        *count_bins(identities, identity_bins),
        "",
        "Coverage breakdown:",
        *count_bins(coverages, coverage_bins),
        "",
        "Query type breakdown:",
        *[f"{query_type}: {count}" for query_type, count in query_counts.items()]
    ]

    for line in report_lines: logger.info(line)

    with output_file.open("w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

def report_flank_hit_stats(blast_path, output_file):
    #same as report_gene_hit_stats but this time just for the flanking gene hits - splits upstream and downstream

    blast_path = Path(blast_path)
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)

    columns = ["query_id", "subject_id", "nident", "pident", "length", "mismatch", "gaps", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]

    flank_rows = {"upstream flank": [], "downstream flank": []}
    excluded_gene_rows = 0

    with blast_path.open("r", encoding = "utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(columns):
                continue

            row = dict(zip(columns, fields))
            if row["query_id"].endswith("_upstream_flank"):
                query_type = "upstream flank"
            elif row["query_id"].endswith("_downstream_flank"):
                query_type = "downstream flank"
            else:
                excluded_gene_rows += 1
                continue

            identity = float(row["pident"])
            coverage = abs(int(row["qend"]) - int(row["qstart"])) / int(row["qlen"]) * 100
            flank_rows[query_type].append((identity, coverage, row["subject_id"]))

    report_lines = [
        "Raw flank BLAST hit statistics",
        f"Source file: {blast_path}",
        f"Raw total flank hits: {sum(len(rows) for rows in flank_rows.values())}",
        f"Excluded gene hits: {excluded_gene_rows}",
        ""
    ]

    for query_type, rows in flank_rows.items():
        identities = [row[0] for row in rows]
        coverages = [row[1] for row in rows]
        isolates = {row[2].split("_ctg", 1)[0] for row in rows}
        report_lines.extend([
            f"{query_type.title()} hits: {len(rows)}",
            f"{query_type.title()} isolates: {len(isolates)}",
            "Identity breakdown:",
            *count_bins(identities, identity_bins),
            "Coverage breakdown:",
            *count_bins(coverages, coverage_bins),
            ""
        ])

    for line in report_lines: logger.info(line)

    with output_file.open("w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

def isolate_sort_key(isolate_id):
    isolate_id = str(isolate_id)
    return (0, int(isolate_id)) if isolate_id.isdigit() else (1, isolate_id)
    
def sort_report_df_by_id(df):
    #sorts reports by numeric id, if id has text they come after (some ids won't split properly and needs further debugging)
    id_column = next((column for column in ("id", "isolate_id") if column in df.columns), None)

    if id_column is None or df.empty:
        return df

    order = sorted(range(len(df)), key = lambda index: isolate_sort_key(df.iloc[index][id_column]))
    return df.iloc[order].reset_index(drop=True)

def clean_sequence_for_report(sequence):
    return "".join(str(sequence).split()).upper() if sequence else ""

def make_query_table(config):
    rows = []

    if config.repeat_location.lower() == "cds":
        query_sequences = [
            ("Gene sequence normalised to 2 repeats", config.ref_2_repeats),
            ("Gene sequence after repeat tract", config.ref_no_repeats)
        ]
    else:
        query_sequences = [("Full gene sequence", config.ref_full)]

    if config.ref_gene_start is not None:
        query_sequences.append(("Gene sequence before repeat tract", config.ref_gene_start))
    if config.ref_gene_end is not None:
        query_sequences.append(("Gene sequence at gene end", config.ref_gene_end))

    query_sequences.extend([
        ("Upstream flank", config.upstream_flank),
        ("Downstream flank", config.downstream_flank),
        ("Forward primer", config.forward_primer),
        ("Reverse primer", config.reverse_primer)
    ])

    for query_type, sequence in query_sequences:
        sequence = clean_sequence_for_report(sequence)
        rows.append({
            "Query type": query_type,
            "Sequence length (bp)": len(sequence),
            "Sequence": sequence
        })

    return pd.DataFrame(rows)

def make_parameter_table(config):
    rows = []

    for name, value in config.gene_blast_params.items():
        rows.append({"Search": "Gene BLAST", "Parameter": name, "Value": value})

    if config.ref_gene_start:
        for name, value in config.gene_blast_params.items():
            rows.append({"Search": "Gene start BLAST", "Parameter": name, "Value": value})

    for name, value in config.primer_blast_params.items():
        rows.append({"Search": "Primer BLAST", "Parameter": name, "Value": value})

    rows.extend([
        {"Search": "Genomic context", "Parameter": "Upstream distance", "Value": config.upstream_distance},
        {"Search": "Genomic context", "Parameter": "Downstream distance", "Value": config.downstream_distance},
        {"Search": "Genomic context", "Parameter": "Flank buffer", "Value": config.flank_buffer}
    ])

    return pd.DataFrame(rows)


def make_isolate_summary(blast_results, isolate_ids, config):
    gene_hits = blast_results.get("gene_hits", [])
    if config.ref_gene_start:
        gene_start_hits = blast_results.get("gene_start_hits", [])
    else:
        gene_start_hits = []
    if config.ref_gene_end:
        gene_end_hits = blast_results.get("gene_end_hits", [])
    else:
        gene_end_hits = []
    forward_primer_hits = blast_results.get("forward_primer_hits", [])
    reverse_primer_hits = blast_results.get("reverse_primer_hits", [])
    contig_counts = blast_results.get("contig_counts", {})
    
    rows = []

    for isolate_id in sorted(set(isolate_ids), key = isolate_sort_key):
        isolate_id = str(isolate_id)
        isolate_gene_hits = [h for h in gene_hits if str(h["isolate_id"]) == isolate_id]
        isolate_gene_start_hits = [h for h in gene_start_hits if str(h["isolate_id"]) == isolate_id]
        isolate_gene_end_hits = [h for h in gene_end_hits if str(h["isolate_id"]) == isolate_id]
        isolate_forward_primer_hits = [h for h in forward_primer_hits if str(h["isolate_id"]) == isolate_id]
        isolate_reverse_primer_hits = [h for h in reverse_primer_hits if str(h["isolate_id"]) == isolate_id]
        isolate_contig_count = contig_counts.get(isolate_id, 0)

        contexts = sorted({h.get("genomic_context", "unknown") for h in isolate_gene_hits})

        rows.append({
            "id": isolate_id,
            "n_contigs": contig_counts.get(isolate_id),
            "gene_found": bool(isolate_gene_hits),
            "n_gene_hits": len(isolate_gene_hits),
            "n_high_quality": sum(bool(h.get("high_quality")) for h in isolate_gene_hits),
            "genomic_context": ";".join(contexts) if contexts else "unknown",
            **({"gene_start_found": bool(isolate_gene_start_hits),
                "n_gene_start_hits": len(isolate_gene_start_hits)} if isolate_gene_start_hits else {}),
            **({"gene_end_found": bool(isolate_gene_end_hits),
                "n_gene_end_hits": len(isolate_gene_end_hits)} if isolate_gene_end_hits else {}),
            "n_forward_primer_hits": len(isolate_forward_primer_hits),
            "n_reverse_primer_hits": len(isolate_reverse_primer_hits)
        })

    return pd.DataFrame(rows)

def make_gene_hits_table(gene_hits):
    rows = []

    for hit in sorted(gene_hits, key = lambda h: (isolate_sort_key(h["isolate_id"]), h.get("gene_copy", ""))):
        rows.append({
            "id": hit["isolate_id"],
            "gene_copy": hit.get("gene_copy"),
            "n_gene_copies": hit.get("n_gene_copies"),
            "contig": hit.get("contig_orig"),
            "identity": round(hit.get("identity"), 2),
            "coverage": round(hit.get("coverage"), 2),
            "high_quality": bool(hit.get("high_quality")),
            "reference": hit.get("blast_query"),
            "genomic_context": hit.get("genomic_context"),
            "upstream_found": bool(hit.get("upstream_found")),
            "downstream_found": bool(hit.get("downstream_found")),
            "upstream_distance": hit.get("upstream_distance"),
            "downstream_distance": hit.get("downstream_distance"),
            "subject_start": hit.get("subject_start"),
            "subject_end": hit.get("subject_end"),
            "strand": hit.get("strand"),
            "aln_length": hit.get("aln_length"),
            "mismatches": hit.get("mismatches"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore")
        })

    return pd.DataFrame(rows)

def make_gene_start_table(gene_start_hits):
    rows = []

    for hit in sorted(gene_start_hits, key = lambda h: (isolate_sort_key(h["isolate_id"]), h.get("gene_copy", ""))):
        aligned_query_bases = abs(hit["q_end"] - hit["q_start"]) + 1
        unaligned_bases = hit["q_length"] - aligned_query_bases
        total_differences = hit["mismatches"] + hit["gaps"] + hit["gap_opens"] + unaligned_bases

        rows.append({
            "id": hit["isolate_id"],
            "contig": hit.get("contig_orig"),
            "identity": round(hit.get("identity"), 2),
            "coverage": round(hit.get("coverage"), 2),
            "high_quality": bool(hit.get("high_quality")),
            "aln_length": hit.get("aln_length"),
            "blast_mismatches": hit.get("mismatches"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "unaligned_bases": unaligned_bases,
            "total_differences": total_differences,
            "perfect_match": total_differences == 0,
            "subject_start": hit.get("subject_start"),
            "subject_end": hit.get("subject_end"),
            "strand": hit.get("strand"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore")
        })

    return pd.DataFrame(rows)

def make_gene_assembly_table(blast_results, config):
    rows = []

    gene_hits = blast_results.get("gene_hits_before_filter", blast_results.get("gene_hits", []))
    if config.ref_gene_start:
        gene_start_hits = blast_results.get("gene_start_hits", [])
    else:
        gene_start_hits = []
    if config.ref_gene_end:
        gene_end_hits = blast_results.get("gene_end_hits", [])
    else:
        gene_end_hits = []

    def add_hit(isolate_id, region, hit, detail = None):
        detail = detail or hit
        rows.append({
            "id": str(isolate_id),
            "ref_region": region,
            "query_type": detail.get("query_type", region),
            "contig": detail.get("contig_orig"),
            "subject_start": detail.get("subject_start"),
            "subject_end": detail.get("subject_end"),
            "identity": round(detail["identity"], 2) if detail.get("identity") is not None else None,
            "coverage": round(detail["coverage"], 2) if detail.get("coverage") is not None else None,
            "aln_length": detail.get("aln_length"),
            "strand": detail.get("strand"),
            "evalue": detail.get("evalue"),
            "bitscore": detail.get("bitscore")
        })

    for hit in gene_hits:
        details = hit.get("query_hit_details") or [hit]
        for detail in details:
            query_type = detail.get("query_type", hit.get("blast_query", "unknown"))
            region = {
                "2_repeats": "Gene 2 repeats",
                "no_repeats": "Gene no repeats",
                "full": "Gene full"
            }.get(query_type, query_type)
            add_hit(hit["isolate_id"], region, hit, detail)

    for hit in gene_start_hits:
        add_hit(hit["isolate_id"], "Gene start", hit)

    for hit in gene_end_hits:
        add_hit(hit["isolate_id"], "Gene end", hit)

    assembly_df = pd.DataFrame(rows)
    if assembly_df.empty:
        return assembly_df

    contigs_by_isolate = assembly_df.groupby("id")["contig"].transform(
        lambda values: ";".join(sorted({str(value) for value in values if pd.notna(value)}))
    )

    assembly_df["contigs_by_isolate"] = contigs_by_isolate
    assembly_df["contig_status"] = assembly_df.groupby("id")["contig"].transform(
        lambda values: "single" if values.dropna().nunique() == 1 else "multiple"
    )

    return sort_report_df_by_id(assembly_df)

def make_flank_table(flank_hits):
    rows = []

    for hit in sorted(flank_hits, key = lambda h: (isolate_sort_key(h["isolate_id"]), h["query_id"])):
        rows.append({
            "id": hit["isolate_id"],
            "flank_type": "upstream" if hit["query_id"].endswith("upstream_flank") else "downstream",
            "contig": hit.get("contig_orig"),
            "identity": round(hit["identity"], 2) if hit.get("identity") is not None else None,
            "coverage": round(hit["coverage"], 2) if hit.get("coverage") is not None else None,
            "aln_length": hit.get("aln_length"),
            "mismatches": hit.get("mismatches"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "strand": hit.get("strand"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore")
        })

    return pd.DataFrame(rows)

def make_primer_table(primer_hits):
    rows = []

    for hit in sorted(primer_hits, key = lambda h: (isolate_sort_key(h["isolate_id"]), h["total_mismatches"], -h["bitscore"])):
        rows.append({
            "id": hit["isolate_id"],
            "contig": hit.get("contig_orig"),
            "identity": round(hit["identity"], 2) if hit.get("identity") is not None else None,
            "coverage": round(hit["coverage"], 2) if hit.get("coverage") is not None else None,
            "total_differences": hit.get("total_mismatches"),
            "blast_mismatches": hit.get("blast_mismatches"),
            "unaligned_bases": hit.get("unaligned_bases"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "mismatch_positions": hit.get("mismatch_positions"),
            "mismatch_positions_from_3prime": hit.get("positions_from_3prime"),
            "perfect_match": hit.get("total_mismatches") == 0,
            "subject_start": hit.get("subject_start"),
            "subject_end": hit.get("subject_end"),
            "strand": hit.get("strand"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore")
        })

    return pd.DataFrame(rows)


def make_totals_table(blast_results, isolate_ids):
    gene_hits = blast_results.get("gene_hits", [])
    gene_hits_before_filter = blast_results.get("gene_hits_before_filter", [])
    gene_start_hits = blast_results.get("gene_start_hits", [])
    gene_end_hits = blast_results.get("gene_end_hits", [])
    upstream_flank_hits = blast_results.get("upstream_flank_hits", [])
    downstream_flank_hits = blast_results.get("downstream_flank_hits", [])
    forward_primer_hits = blast_results.get("forward_primer_hits", [])
    reverse_primer_hits = blast_results.get("reverse_primer_hits", [])

    totals = [
        ("input isolates", len(set(isolate_ids))),
        ("gene hits before filter", len(gene_hits_before_filter)),
        ("gene hits", len(gene_hits)),
        ("high quality gene hits", sum(bool(h.get("high_quality")) for h in gene_hits)),
        ("gene start hits", len(gene_start_hits)),
        ("gene end hits", len(gene_end_hits)),
        ("upstream flank hits", len(upstream_flank_hits)),
        ("downstream flank hits", len(downstream_flank_hits)),
        ("gene hits with both flanks", sum(h.get("genomic_context") == "both_flanks" for h in gene_hits)),
        ("gene hits with upstream flank", sum(h.get("genomic_context") == "upstream_only" for h in gene_hits)),
        ("gene hits with downstream flank", sum(h.get("genomic_context") == "downstream_only" for h in gene_hits)),
        ("forward primer hits", len(forward_primer_hits)),
        ("reverse primer hits", len(reverse_primer_hits)),
        ("perfect forward primer hits", sum(h.get("total_mismatches") == 0 for h in forward_primer_hits)),
        ("perfect reverse primer hits", sum(h.get("total_mismatches") == 0 for h in reverse_primer_hits))
    ]

    return pd.DataFrame(totals, columns=["Metric", "Count"])

def make_primer_mismatch_table(forward_primer_hits, reverse_primer_hits):
    rows = []
    for primer, hits in [("forward", forward_primer_hits), ("reverse", reverse_primer_hits)]:
        mismatch_counts = sorted({hit["total_mismatches"] for hit in hits})
        for mismatch_count in mismatch_counts:
            count = sum(hit["total_mismatches"] == mismatch_count for hit in hits)
            percentage = count / len(hits) * 100 if hits else 0

            rows.append({
                "Primer": primer,
                "Mismatches": mismatch_count,
                "Count": count,
                "Percentage": round(percentage, 2)
            })

    return pd.DataFrame(rows)

def create_blast_summary_report(blast_results, config, isolate_ids, output_dir):
    gene_name = config.gene_name
    output_path = Path(output_dir) / f"{gene_name}_blast_summary.xlsx"

    gene_hits = blast_results.get("gene_hits", [])
    gene_hits_before_filter = blast_results.get("gene_hits_before_filter", [])
    gene_start_hits = blast_results.get("gene_start_hits", [])
    gene_end_hits = blast_results.get("gene_end_hits", [])
    upstream_flank_hits = blast_results.get("upstream_flank_hits", [])
    downstream_flank_hits = blast_results.get("downstream_flank_hits", [])
    flank_hits = blast_results.get("flank_hits", [])
    forward_primer_hits = blast_results.get("forward_primer_hits", [])
    reverse_primer_hits = blast_results.get("reverse_primer_hits", [])

    logger.info("Creating BLAST summary report for gene: %s", gene_name)

    query_df = make_query_table(config)
    parameter_df = make_parameter_table(config)

    blast_queries = pd.concat([query_df, parameter_df.rename(columns={"Search": "Query Type", "Parameter": "Parameter"})], ignore_index = True, sort = False)

    isolate_df = make_isolate_summary(blast_results, isolate_ids, config).rename(columns = {
        "gene_found": "Gene Found",
        "n_gene_hits": "Gene Hits",
        "genomic_context": "Genomic Context",
        "n_forward_primer_hits": "Forward Primer Hits",
        "n_reverse_primer_hits": "Reverse Primer Hits"
    })

    isolate_df["Gene Found"] = isolate_df["Gene Found"].map({True: "Yes", False: "No"})

    def gene_rows(hits):
        return sort_report_df_by_id(pd.DataFrame([{
            "id": hit["isolate_id"],
            "gene copy": hit.get("gene_copy", "NA"),
            "contig": hit.get("contig_orig"),
            "identity": round(hit["identity"], 2),
            "coverage": round(hit["coverage"], 2),
            "high quality": "Yes" if hit.get("high_quality", False) else "No",
            "Query": hit.get("blast_query", "NA"),
            "Genomic context": hit.get("genomic_context", "NA"),
            "Upstream flank": hit.get("upstream_found", "NA"),
            "Downstream flank": hit.get("downstream_found", "NA"),
            "Subject start": hit.get("s_start"),
            "Subject end": hit.get("s_end"),
            "Strand": hit.get("strand", "NA"),
            "Aln_length": hit.get("aln_length"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore"),
        } for hit in hits]))

    def primer_rows(hits):
        return sort_report_df_by_id(pd.DataFrame([{
            "id": hit["isolate_id"],
            "primer": hit.get("query_id", "NA"),
            "contig": hit.get("contig_orig"),
            "identity": round(hit["identity"], 2),
            "coverage": round(hit["coverage"], 2),
            "aln_length": hit.get("aln_length"),
            "query_start": hit.get("q_start"),
            "query_end": hit.get("q_end"),
            "subject_start": hit.get("s_start"),
            "subject_end": hit.get("s_end"),
            "total_differences": hit.get("total_mismatches"),
            "blast_mismatches": hit.get("blast_mismatches"),
            "unaligned_bases": hit.get("unaligned_bases"),
            "gaps": hit.get("gaps"),
            "gap_opens": hit.get("gap_opens"),
            "mismatch_positions": hit.get("mismatch_positions"),
            "mismatch_positions_from_3prime": hit.get("positions_from_3prime"),
            "strand": hit.get("strand", "NA"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore"),
        } for hit in hits]))

    flank_df = make_flank_table(flank_hits)
    assembly_df = make_gene_assembly_table(blast_results, config)
    primer_df = pd.concat([
        primer_rows(forward_primer_hits).assign(primer_type = "forward"),
        primer_rows(reverse_primer_hits).assign(primer_type = "reverse")
    ], ignore_index = True)

    hit_types = ["Gene hits", "Gene start hits", "Gene end hits", "Upstream flank hits", "Downstream flank hits", "Forward primer hits", "Reverse primer hits"]

    hit_counts = len(gene_hits), len(gene_start_hits), len(gene_end_hits), len(upstream_flank_hits), len(downstream_flank_hits), len(forward_primer_hits), len(reverse_primer_hits)

    total_hits = sum(hit_counts)
    totals_df = pd.DataFrame({
        "Hit type": hit_types,
        "Total hits": hit_counts,
        "Percentage of total": [count / total_hits * 100 for count in hit_counts]
    })

    with pd.ExcelWriter(output_path, engine = "openpyxl") as writer:
        blast_queries.to_excel(writer, sheet_name = "Blast Queries", index = False)
        totals_df.to_excel(writer, sheet_name = "Hit totals", index = False)
        isolate_df.to_excel(writer, sheet_name = "Isolate summary", index = False)
        gene_rows(gene_hits).to_excel(writer, sheet_name = "Gene Hits", index = False)
        if config.ref_gene_start:
            gene_rows(gene_start_hits).to_excel(writer, sheet_name = "Gene Start Hits", index = False)
        if config.ref_gene_end:
            gene_rows(gene_end_hits).to_excel(writer, sheet_name = "Gene End Hits", index = False)
        flank_df.to_excel(writer, sheet_name = "Flank Hits", index = False)
        assembly_df.to_excel(writer, sheet_name = "Gene Assembly", index = False)
        primer_rows(forward_primer_hits).to_excel(writer, sheet_name = "Forward Primer Hits", index = False)
        primer_rows(reverse_primer_hits).to_excel(writer, sheet_name = "Reverse Primer Hits", index = False)
        make_primer_mismatch_table(forward_primer_hits, reverse_primer_hits).to_excel(writer, sheet_name = "Hit totals", index = False, startrow = len(totals_df) + 3)
        format_workbook(writer)

    assembly_df.to_excel(output_dir / "assembly.xlsx", index = False)
    primer_df.to_excel(output_dir / "primers.xlsx", index = False)

    logger.info("BLAST summary report saved to %s", output_path)
    return output_path

#2 - Sequence extraction report
def make_extraction_table(sequences, kept_sequences, config):
    promoter_mode = config.repeat_location.lower() == "promoter"
    columns = [
        "id", "gene_copy", "contig",
        "N_filter", "sequence_N_count", "boundary_status", "blast_query",
        "genomic_context", "contig_length", "gene_start", "gene_end",
        "extract_start", "extract_end", "extracted_length",
        "dist_to_contig_start", "dist_to_contig_end", "strand", "flags"
    ]
    if promoter_mode:
        columns.extend([
            "promoter_length_requested", "promoter_length_extracted",
            "promoter_start", "promoter_end", "promoter_N_count"
        ])

    rows = []
    for s in sequences:
        key = (s["isolate_id"], s["gene_copy"], s["contig_id"])
        row = {
            "id": s["isolate_id"],
            "gene_copy": s["gene_copy"],
            "contig": s.get("original_contig_id", s["contig_id"]),
            "N_filter": "Kept" if key in kept_sequences else "Removed",
            "sequence_N_count": s["sequence"].upper().count("N"),
            "boundary_status": s["boundary_status"],
            "blast_query": s.get("blast_query", "unknown"),
            "genomic_context": s.get("genomic_context", "unknown"),
            "contig_length": s["contig_length"],
            "gene_start": s["gene_start"],
            "gene_end": s["gene_end"],
            "extract_start": s["extract_start"],
            "extract_end": s["extract_end"],
            "extracted_length": len(s["sequence"]),
            "dist_to_contig_start": s["dist_to_contig_start"],
            "dist_to_contig_end": s["dist_to_contig_end"],
            "strand": s["strand"],
            "flags": "; ".join(s.get("flags") or [])
        }
        if promoter_mode:
            row.update({
                "promoter_length_requested": s["promoter_length_requested"],
                "promoter_length_extracted": s["promoter_length_extracted"],
                "promoter_start": s["promoter_start"],
                "promoter_end": s["promoter_end"],
                "promoter_N_count": (s.get("promoter_sequence") or "").upper().count("N")
            })
        rows.append(row)

    return sort_report_df_by_id(pd.DataFrame(rows, columns = columns))

def create_sequence_extraction_report(before_N_filtering, extracted_sequences, config, output_dir):
    gene_name = config.gene_name
    params = config.extraction_params
    promoter_mode = config.repeat_location.lower() == "promoter"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents = True, exist_ok = True)
    output_path = output_dir / f"{gene_name}_sequence_extraction.xlsx"

    logger.info("Creating sequence extraction report for gene: %s", gene_name)

    kept_sequences = {
        (s["isolate_id"], s["gene_copy"], s["contig_id"])
        for s in extracted_sequences
    }
    before_df = make_extraction_table(before_N_filtering, kept_sequences, config)
    after_df = make_extraction_table(extracted_sequences, kept_sequences, config)
    removed_df = before_df[before_df["N_filter"] == "Removed"]

    parameter_rows = [
        ("Gene name", gene_name),
        ("Repeat location", config.repeat_location),
        ("Flanking bp used", 0 if promoter_mode else params["flanking_bp"]),
        ("Critical boundary distance (bp)", params["boundary_critical_distance"]),
        ("Warning boundary distance (bp)", params["boundary_warning_distance"]),
        ("Coordinate system", "0-based, end-exclusive"),
        ("Boundary assessment", "BLAST hit boundaries"),
        ("N filter scope", "Gene sequence only; promoter excluded" if promoter_mode else "Entire extracted sequence including flanks")
    ]
    if promoter_mode:
        parameter_rows.append(("Promoter length (bp)", params["promoter_length"]))
    parameter_df = pd.DataFrame(parameter_rows, columns = ["Parameter", "Value"])

    n_before = len(before_df)
    n_removed = len(removed_df)
    totals_df = pd.DataFrame([
        ("Loci before N filtering", n_before),
        ("Loci after N filtering", len(after_df)),
        ("Loci removed by N filtering", n_removed),
        ("Percentage removed", round(n_removed / n_before * 100, 2) if n_before else 0)
    ], columns = ["Metric", "Value"])

    boundary_df = pd.DataFrame([
        {
            "Boundary status": status,
            "Before N filtering": sum(s["boundary_status"] == status for s in before_N_filtering),
            "After N filtering": sum(s["boundary_status"] == status for s in extracted_sequences)
        }
        for status in ["OK", "WARNING", "CRITICAL"]
    ])

    with pd.ExcelWriter(output_path, engine = "openpyxl") as writer:
        parameter_df.to_excel(writer, sheet_name = "Extraction Parameters", index = False)
        current_row = write_df_with_spacing(totals_df, writer, "Extraction Summary", 0)
        write_df_with_spacing(boundary_df, writer, "Extraction Summary", current_row)
        before_df.to_excel(writer, sheet_name = "Before N Filtering", index = False)
        after_df.to_excel(writer, sheet_name = "After N Filtering", index = False)
        removed_df.to_excel(writer, sheet_name = "Removed by N Filtering", index = False)
        format_workbook(writer)

        #summary contains two separate tables, so remove the worksheet-wide filter
        writer.book["Extraction Summary"].auto_filter.ref = None

    logger.info("Sequence extraction report saved to %s", output_path)
    return output_path
#one function to make all reports 

def generate_all_reports(config, blast_results, output_dir, isolate_ids, blast_outputs, extracted_sequences, before_N_filtering):
    output_dir = Path(output_dir)
    gene_name = config.gene_name
    reports_generated = {
        "gene_hits": report_gene_hit_stats(blast_outputs["gene"], output_dir / f"{gene_name}_gene_hit_stats.txt"),
        "flank_hits": report_flank_hit_stats(blast_outputs["gene"], output_dir / f"{gene_name}_flank_hit_stats.txt"),
        "blast_summmary": create_blast_summary_report(blast_results, config, isolate_ids, output_dir),
        "sequence_extraction": create_sequence_extraction_report(before_N_filtering, extracted_sequences, config, output_dir)
    }

    return reports_generated
