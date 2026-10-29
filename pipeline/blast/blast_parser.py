import logging
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger(__name__)

gene_columns = ["query_id", "subject_id", "num_ident", "pct_ident", "aln_length", "mismatches", "gaps", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore", "qlength", "slength"]

primer_columns = gene_columns + ["q_sequence", "s_sequence"]

def parse_contig_id(subject_id, contig_id_map):
    if "_ctg" in subject_id:
        isolate_id, contig_id = subject_id.rsplit("_ctg", 1)
        contig_local = f"ctg{contig_id}"
    else:
        isolate_id = "unknown"
        contig_local = subject_id

    return {
        "isolate_id": isolate_id,
        "contig_id": subject_id,
        "contig_local": contig_local,
        "contig_orig": contig_id_map.get(subject_id, subject_id)
    }

def parse_gene_blast(blast_path, config, contig_id_map):
    params = config.gene_blast_params
    gene_hits = []
    flank_hits = []

    logger.info("Parsing gene BLAST results with minimum identity %s and coverage %s", params["min_identity"], params["min_coverage"])

    with Path(blast_path).open("r", encoding = "utf-8") as blast_file:
        for line in blast_file:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")

            if len(fields) != len(gene_columns):
                logger.warning("Skipping line with unexpected number of columns: %s", line.strip())
                continue

            row = dict(zip(gene_columns, fields))
            identity = float(row["pct_ident"])
            query_length = int(row["qlength"])
            query_start = int(row["q_start"])
            query_end = int(row["q_end"])
            coverage = (abs(query_end - query_start) + 1) / query_length * 100

            is_flank = row["query_id"].endswith("_flank")

            if identity < params["min_identity"] or coverage < params["min_coverage"]:
                continue


            contig = parse_contig_id(row["subject_id"], contig_id_map)
            subject_start = int(row["s_start"])
            subject_end = int(row["s_end"])

            hit = {
                "query_id": row["query_id"],
                **contig,
                "identity": identity,
                "coverage": coverage,
                "evalue": float(row["evalue"]),
                "bitscore": float(row["bitscore"]),
                "aln_length": int(row["aln_length"]),
                "q_length": query_length,
                "q_start": query_start,
                "q_end": query_end,
                "s_length": int(row["slength"]),
                "s_start": subject_start,
                "s_end": subject_end,
                "strand": "forward" if (query_end - query_start) * (subject_end - subject_start) >= 0 else "reverse",
                "mismatches": int(row["mismatches"]),
                "gaps": int(row["gaps"]),
                "gap_opens": int(row["gap_opens"]),
                "high_quality": identity >= params["high_identity"] and coverage >= params["high_coverage"]
            }

            if is_flank:
                flank_hits.append(hit)
            else:
                gene_hits.append(hit)

    logger.info("Found %s gene hits and %d flank hits", len(gene_hits), len(flank_hits))
    return gene_hits, flank_hits

def parse_gene_start_blast(blast_path, config, contig_id_map, params = None):
    if params is None:
        params = config.gene_start_blast_params
    gene_start_hits = []

    high_identity = params.get("high_identity", params.get("high_identity", params["min_identity"]))
    high_coverage = params.get("high_coverage", params["min_coverage"])

    logger.info("Parsing gene start BLAST results with minimum identity %s and coverage %s", params["min_identity"], params["min_coverage"])

    with Path(blast_path).open("r", encoding = "utf-8") as blast_file:
        for line in blast_file:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip("\n").split("\t")

            if len(fields) != len(gene_columns):
                logger.warning("Skipping gene-start BLAST row with %d columns", len(fields))
                continue

            row = dict(zip(gene_columns, fields))
            identity = float(row["pct_ident"])
            query_length = int(row["qlength"])
            query_start = int(row["q_start"])
            query_end = int(row["q_end"])
            coverage = (abs(query_end - query_start) + 1) / query_length * 100

            if identity < params["min_identity"] or coverage < params["min_coverage"]:
                continue

            contig = parse_contig_id(row["subject_id"], contig_id_map)
            subject_start = int(row["s_start"])
            subject_end = int(row["s_end"])

            hit = {
                "query_id": row["query_id"],
                **contig,
                "identity": identity,
                "coverage": coverage,
                "evalue": float(row["evalue"]),
                "bitscore": float(row["bitscore"]),
                "aln_length": int(row["aln_length"]),
                "q_length": query_length,
                "q_start": query_start,
                "q_end": query_end,
                "s_length": int(row["slength"]),
                "s_start": subject_start,
                "s_end": subject_end,
                "strand": "forward" if (query_end - query_start) * (subject_end - subject_start) >= 0 else "reverse",
                "mismatches": int(row["mismatches"]),
                "gaps": int(row["gaps"]),
                "gap_opens": int(row["gap_opens"]),
                "high_quality": identity >= high_identity and coverage >= high_coverage
            }
            gene_start_hits.append(hit)

    hits_by_isolate = defaultdict(list)
    for hit in gene_start_hits:
        hits_by_isolate[hit["isolate_id"]].append(hit)

    filtered_hits = []

    for isolate_hits in hits_by_isolate.values():
        isolate_hits.sort(key=lambda hit: (-hit["identity"], -hit["coverage"], hit ["evalue"]))
        filtered_hits.extend(isolate_hits[:params["max_hits_per_isolate"]])
    logger.info("Found %s gene start hits", len(filtered_hits))

    gene_start_hits = filtered_hits
    return gene_start_hits

def get_primer_mismatch_positions(query_sequence, subject_sequence, query_start, query_length):
    query_position = query_start - 1
    mismatch_positions = []

    for query_base, subject_base in zip(query_sequence.upper(), subject_sequence.upper()):
        if query_base == "-":
            continue

        query_position += 1

        if subject_base == "-" or query_base != subject_base:
            mismatch_positions.append(query_position)

    mismatch_positions.extend(range(1, query_start))
    mismatch_positions.extend(range(query_position + 1, query_length + 1))
    mismatch_positions = sorted(set(mismatch_positions))

    positions_from_3prime = sorted(query_length - pos + 1 for pos in mismatch_positions)

    return mismatch_positions, positions_from_3prime

def parse_primer_blast(blast_path, config, contig_id_map):
    params = config.primer_blast_params
    forward_hits = []
    reverse_hits = []

    with Path(blast_path).open("r", encoding = "utf-8") as blast_file:
        for line in blast_file:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.strip("\n").split("\t")

            if len(fields) != len(primer_columns):
                logger.warning("Skipping primer BLAST row with %d columns", len(fields))
                continue

            row = dict(zip(primer_columns, fields))
            query_length = int(row["qlength"])
            query_start = int(row["q_start"])
            query_end = int(row["q_end"])
            aligned_bases = abs(query_end - query_start) + 1
            coverage = aligned_bases / query_length * 100
            unaligned_bases = query_length - aligned_bases
            mismatches = int(row["mismatches"]) + unaligned_bases + int(row["gaps"])

            if coverage < params["min_alignment_pct"]:
                continue
            if mismatches > params["max_mismatches"]:
                continue

            query_sequence = row["q_sequence"]
            subject_sequence = row["s_sequence"]
            position_start = query_start
            if query_start > query_end:
                complement = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
                query_sequence = query_sequence.upper().translate(complement)[::-1]
                subject_sequence = subject_sequence.upper().translate(complement)[::-1]
                position_start = query_end

            mismatch_positions, positions_from_3prime = get_primer_mismatch_positions(
                query_sequence=query_sequence,
                subject_sequence=subject_sequence,
                query_start=position_start,
                query_length=query_length
            )

            contig = parse_contig_id(row["subject_id"], contig_id_map)
            subject_start = int(row["s_start"])
            subject_end = int(row["s_end"])

            hit = {
                "query_id": row["query_id"],
                **contig,
                "identity": float(row["pct_ident"]),
                "coverage": coverage,
                "total_mismatches": mismatches,
                "blast_mismatches": int(row["mismatches"]),
                "unaligned_bases": unaligned_bases,
                "gaps": int(row["gaps"]),
                "gap_opens": int(row["gap_opens"]),
                "mismatch_positions": ";".join(map(str, mismatch_positions)),
                "positions_from_3prime": ";".join(map(str, positions_from_3prime)),
                "q_start": query_start,
                "q_end": query_end,
                "s_start": subject_start,
                "s_end": subject_end,
                "strand": "forward" if (query_end - query_start) * (subject_end - subject_start) >= 0 else "reverse",
                "evalue": float(row["evalue"]),
                "bitscore": float(row["bitscore"])
            }

            if row["query_id"].endswith("_forward_primer"):
                forward_hits.append(hit)
            elif row["query_id"].endswith("_reverse_primer"):
                reverse_hits.append(hit)

    logger.info("Parsed %d forward hits and %d reverse hits", len(forward_hits), len(reverse_hits))
    return forward_hits, reverse_hits


def merge_gene_variants(hits, gene_name, overlap_fraction=0.80):
    """Merge overlapping hits from different reference versions of the same gene."""

    if not hits:
        return []

    def query_type(hit):
        query_id = hit["query_id"].lower()
        if "ref_full" in query_id:
            return "full"
        if "ref_2_repeats" in query_id:
            return "2_repeats"
        if "ref_no_repeats" in query_id:
            return "no_repeats"
        return "unknown"

    def coordinates(hit):
        return min(hit["s_start"], hit["s_end"]), max(hit["s_start"], hit["s_end"])

    def same_locus(hit_a, hit_b):
        if hit_a["isolate_id"] != hit_b["isolate_id"]:
            return False
        if hit_a["contig_id"] != hit_b["contig_id"]:
            return False
        if hit_a["strand"] != hit_b["strand"]:
            return False

        a_start, a_end = coordinates(hit_a)
        b_start, b_end = coordinates(hit_b)
        overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)

        if overlap == 0:
            return False

        shorter_length = min(a_end - a_start + 1, b_end - b_start + 1)
        return overlap / shorter_length >= overlap_fraction

    query_priority = {"full": 3, "2_repeats": 3, "no_repeats": 2, "unknown": 0}

    def hit_priority(hit):
        return (
            query_priority[query_type(hit)],
            int(hit.get("high_quality", False)),
            hit["bitscore"],
            hit["coverage"],
            hit["identity"],
        )

    prepared_hits = []

    for hit in hits:
        hit = hit.copy()
        hit["blast_query"] = query_type(hit)
        hit["supporting_queries"] = hit["blast_query"]
        hit["query_hit_details"] = [{
            "query_type": hit["blast_query"],
            "contig_id": hit.get("contig_id"),
            "contig_orig": hit.get("contig_orig"),
            "identity": hit.get("identity"),
            "coverage": hit.get("coverage"),
            "s_start": hit.get("s_start"),
            "s_end": hit.get("s_end"),
            "strand": hit.get("strand"),
            "aln_length": hit.get("aln_length"),
            "q_length": hit.get("q_length"),
            "evalue": hit.get("evalue"),
            "bitscore": hit.get("bitscore"),
        }]
        prepared_hits.append(hit)

    prepared_hits.sort(key=hit_priority, reverse=True)
    merged_hits = []

    for candidate in prepared_hits:
        matching_hit = next((hit for hit in merged_hits if same_locus(candidate, hit)), None)

        if matching_hit is None:
            merged_hits.append(candidate)
            continue

        supporting_queries = set(matching_hit["supporting_queries"].split(";"))
        supporting_queries.add(candidate["blast_query"])
        matching_hit["supporting_queries"] = ";".join(sorted(supporting_queries))
        matching_hit["n_query_variants"] = len(supporting_queries)
        matching_hit["query_hit_details"].extend(candidate["query_hit_details"])

        if hit_priority(candidate) > hit_priority(matching_hit):
            candidate["supporting_queries"] = matching_hit["supporting_queries"]
            candidate["n_query_variants"] = matching_hit["n_query_variants"]
            candidate["query_hit_details"] = matching_hit["query_hit_details"]
            merged_hits[merged_hits.index(matching_hit)] = candidate

    for hit in merged_hits:
        hit.setdefault("n_query_variants", 1)

    logger.info("Merged %d BLAST hits into %d candidate %s loci", len(hits), len(merged_hits), gene_name)
    return merged_hits

def validate_genomic_context(gene_hits, flank_hits, config):
#do the gene hits have BLAST hits for the flanks where we expect

    upstream_threshold = config.upstream_distance + config.flank_buffer
    downstream_threshold = config.downstream_distance + config.flank_buffer

    flanks_by_contig = defaultdict(list)
    for flank in flank_hits:
        flanks_by_contig[(flank.get("isolate_id"), flank.get("contig_id"))].append(flank)

    for gene_hit in gene_hits:
        gene_hit["upstream_found"] = False
        gene_hit["downstream_found"] = False
        gene_hit["upstream_distance"] = None
        gene_hit["downstream_distance"] = None
        gene_hit["upstream_flank_coords"] = None
        gene_hit["downstream_flank_coords"] = None
        gene_hit["genomic_context"] = "no_flanks"

        gene_start = min(gene_hit["s_start"], gene_hit["s_end"])
        gene_end = max(gene_hit["s_start"], gene_hit["s_end"])
        nearby_flanks = flanks_by_contig.get((gene_hit.get("isolate_id"), gene_hit.get("contig_id")), [])

        for flank in nearby_flanks:
            flank_start = min(flank["s_start"], flank["s_end"])
            flank_end = max(flank["s_start"], flank["s_end"])
            query_id = flank["query_id"].lower()

            if "upstream" in query_id:
                if gene_hit["strand"] == "forward" and flank_end <= gene_start:
                    distance = gene_start - flank_end
                elif gene_hit["strand"] == "reverse" and flank_start >= gene_end:
                    distance = flank_start - gene_end
                else:
                    continue

                if distance <= upstream_threshold and (gene_hit["upstream_distance"] is None or distance < gene_hit["upstream_distance"]):
                    gene_hit["upstream_found"] = True
                    gene_hit["upstream_distance"] = distance
                    gene_hit["upstream_flank_coords"] = (flank_start - 1, flank_end)

            elif "downstream" in query_id:
                if gene_hit["strand"] == "forward" and flank_start >= gene_end:
                    distance = flank_start - gene_end
                elif gene_hit["strand"] == "reverse" and flank_end <= gene_start:
                    distance = gene_start - flank_end
                else:
                    continue

                if distance <= downstream_threshold and (gene_hit["downstream_distance"] is None or distance < gene_hit["downstream_distance"]):
                    gene_hit["downstream_found"] = True
                    gene_hit["downstream_distance"] = distance
                    gene_hit["downstream_flank_coords"] = (flank_start - 1, flank_end)

        if gene_hit["upstream_found"] and gene_hit["downstream_found"]:
            gene_hit["genomic_context"] = "both_flanks"
        elif gene_hit["upstream_found"]:
            gene_hit["genomic_context"] = "upstream_only"
        elif gene_hit["downstream_found"]:
            gene_hit["genomic_context"] = "downstream_only"

    context_counts = defaultdict(int)

    for hit in gene_hits:
        context_counts[hit["genomic_context"]] += 1

    logger.info("Genomic context counts: %s", dict(context_counts))

    return gene_hits

def filter_top_hits_per_isolate(hits, config):
#keep max no hits (config) per isolate
    max_hits = config.gene_blast_params["max_hits_per_isolate"]

    hits_by_isolate = defaultdict(list)
    for hit in hits:
        hits_by_isolate[hit["isolate_id"]].append(hit)

    filtered_hits = []
    for isolate_id, isolate_hits in hits_by_isolate.items():
        isolate_hits.sort(key=lambda hit: (-hit["identity"], -hit["coverage"], hit["evalue"]))
        top_hits = isolate_hits[:max_hits]

        for copy_number, hit in enumerate(top_hits, start = 1):
            hit["gene_copy"] = f"copy_{copy_number}"
            filtered_hits.append(hit)

    logger.info("Filtered to %d hits after limiting to %d per isolate", len(filtered_hits), max_hits)
    return filtered_hits

def parse_blast_outputs(blast_outputs, config, contig_id_map):
#do all parsing in one function
    gene_hits, flank_hits = parse_gene_blast(blast_outputs["gene"], config, contig_id_map)
    gene_start_hits = []
    if blast_outputs.get("gene_start") is not None:
        gene_start_hits = parse_gene_start_blast(blast_outputs["gene_start"], config, contig_id_map)

    gene_end_hits = []
    if blast_outputs.get("gene_end") is not None:
        gene_end_hits = parse_gene_start_blast(blast_outputs["gene_end"], config, contig_id_map, params = config.gene_blast_params)

    forward_hits, reverse_hits = parse_primer_blast(blast_outputs["primers"], config, contig_id_map)

    gene_hits = merge_gene_variants(gene_hits, config.gene_name)
    validated_gene_hits = validate_genomic_context(gene_hits, flank_hits, config)
    final_gene_hits = filter_top_hits_per_isolate(validated_gene_hits, config)

    return {
        "gene_hits": final_gene_hits,
        "flank_hits": flank_hits,
        "upstream_flank_hits": [hit for hit in flank_hits if hit["query_id"].endswith("upstream_flank")],
        "downstream_flank_hits": [hit for hit in flank_hits if hit["query_id"].endswith("downstream_flank")],
        "gene_hits_before_filter": validated_gene_hits,
        "gene_start_hits": gene_start_hits,
        "gene_end_hits": gene_end_hits,
        "forward_primer_hits": forward_hits,
        "reverse_primer_hits": reverse_hits,
    }