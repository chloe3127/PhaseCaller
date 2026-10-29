import io
import zipfile
import logging
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

def load_genome_sequences(genome_dir):
    logger.info(f"Loading genome sequences from {genome_dir}")

    genome_seqs = {}
    fasta_extensions = (".fasta", ".fa", ".fna")

    genome_path = Path(genome_dir)

    if genome_path.suffix.lower() == ".zip":
        with zipfile.ZipFile(genome_path, "r") as zf:
            entries = [e for e in zf.namelist() if e.lower().endswith(fasta_extensions) and not e.startswith("__MACOSX") and not Path(e).name.startswith("._")]
            for entry in entries:
                isolate_id = Path(entry).stem
                with zf.open(entry) as f:
                    text_handle = io.TextIOWrapper(f, encoding = "utf-8")
                    for i, record in enumerate(SeqIO.parse(text_handle, "fasta"), start = 1):
                        contig_key = f"{isolate_id}_contig{i:06d}"
                        genome_seqs[contig_key] = record
    elif genome_path.is_dir():
        genome_files = sorted(path for path in genome_path.iterdir() if path.suffix.lower() in fasta_extensions)
        for genome_file in genome_files:
            isolate_id = genome_file.stem
            for i, record in enumerate(SeqIO.parse(genome_file, "fasta"), start = 1):
                contig_key = f"{isolate_id}_contig{i:06d}"
                genome_seqs[contig_key] = record
    elif genome_path.is_file() and genome_path.suffix.lower() in fasta_extensions:
        for record in SeqIO.parse(genome_path, "fasta"):
            genome_seqs[record.id] = record
    else:
        raise FileNotFoundError(f"No valid genome files found in {genome_dir}")

    return genome_seqs

def extract_sequence_from_coords(contig_seq, start, end, strand, flanking_bp):
    #extract sequence with flanking regions
    contig_length = len(contig_seq.seq)
    left_1 = min(start, end)
    right_1 = max(start, end)
    left_0 = left_1 - 1
    right_0 = right_1 

    extract_start = max(0, left_0 - flanking_bp)
    extract_end = min(contig_length, right_0 + flanking_bp)

    extracted_seq = contig_seq.seq[extract_start:extract_end]

    if strand == "reverse":
        extracted_seq = extracted_seq.reverse_complement()

    dist_to_contig_start = left_0
    dist_to_contig_end = contig_length - right_0

    return {
        "sequence": str(extracted_seq),
        "contig_length": contig_length,
        "extract_start": extract_start,
        "extract_end": extract_end,
        "gene_left_0": left_0,
        "gene_right_0": right_0,
        "strand": strand,
        "dist_to_contig_start": dist_to_contig_start,
        "dist_to_contig_end": dist_to_contig_end
    }

def extract_promoter_and_gene(contig_seq, start, end, strand, promoter_length):
    #extract the BLAST defined gene exactly with co-ords and extract a defined upstream promoter length region
    contig_length = len(contig_seq.seq)
    gene_left_0 = min(start, end) - 1
    gene_right_0 = max(start, end)
    promoter_length = max(0, int(promoter_length))

    if strand == "forward":
        promoter_start = max(0, gene_left_0 - promoter_length)
        promoter_end = gene_left_0
        gene_sequence = contig_seq.seq[gene_left_0:gene_right_0]
        promoter_sequence = contig_seq.seq[promoter_start:promoter_end]
        extract_start = promoter_start
        extract_end = gene_right_0
    else:
        promoter_start = gene_right_0
        promoter_end = min(contig_length, gene_right_0 + promoter_length)
        gene_sequence = contig_seq.seq[gene_left_0:gene_right_0].reverse_complement()
        promoter_sequence = contig_seq.seq[promoter_start:promoter_end].reverse_complement()
        extract_start = gene_left_0
        extract_end = promoter_end

    return {
        "sequence": str(gene_sequence),
        "promoter_sequence": str(promoter_sequence),
        "promoter_length_requested": promoter_length,
        "promoter_length_extracted": len(promoter_sequence),
        "promoter_start": promoter_start,
        "promoter_end": promoter_end,
        "contig_length": contig_length,
        "extract_start": extract_start,
        "extract_end": extract_end,
        "gene_left_0": gene_left_0,
        "gene_right_0": gene_right_0,
        "strand": strand,
        "dist_to_contig_start": gene_left_0,
        "dist_to_contig_end": contig_length - gene_right_0,
    }

def check_sequence_quality(seq_data, boundary_critical, boundary_warning):
    #check if sequence is too close to contig boundaries
    flags = []
    boundary_status = "OK"

    if (seq_data["dist_to_contig_start"] <= boundary_critical or seq_data["dist_to_contig_end"] <= boundary_critical):
        boundary_status = "CRITICAL"
        flags.append(f"Boundary CRITICAL: Gene within {boundary_critical}bp of contig edge")

    elif (seq_data["dist_to_contig_start"] <= boundary_warning or seq_data["dist_to_contig_end"] <= boundary_warning):
        boundary_status = "WARNING"
        flags.append(f"Boundary WARNING: Gene within {boundary_warning}bp of contig edge")

    return {
        "boundary_status": boundary_status,
        "flags": flags
    }

def remove_ambiguous_bases(extracted_seqs):
    #filter out any sequences containing Ns
    total = len(extracted_seqs)
    filtered = [s for s in extracted_seqs if 'N' not in s["sequence"].upper()]
    n_removed = total - len(filtered)

    logger.info(f"Sequences before filtering for ambiguous bases: {total}")
    logger.info(f"Sequences removed due to ambiguous bases: {n_removed}")
    logger.info(f"Sequences after filtering for ambiguous bases: {len(filtered)}")
    return filtered

def write_extracted_fasta(sequences, output_dir, config):
    gene_name = config.gene_name
    flanking_bp = 0 if config.repeat_location.lower() == "promoter" else config.extraction_params["flanking_bp"]

    all_records = []
    for seq in sequences:
        #build fasta headers
        header_parts = [
            f"isolate={seq['isolate_id']}",
            f"copy={seq['gene_copy']}",
            f"contig={seq['contig_id']}",
            f"coords={seq['extract_start']}-{seq['extract_end']}",
            f"strand={seq['strand']}",
            f"identity={seq['blast_identity']:.1f}%",
            f"genomic_context={seq['genomic_context']}"
        ]

        record = SeqRecord(
            Seq(seq["sequence"]),
            id = f"{seq['isolate_id']}_{seq['gene_copy']}",
            description = " ".join(header_parts)
        )
        all_records.append(record)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok = True)
    all_fasta_path = output_dir / f"{gene_name}_{flanking_bp}bp_all_sequences.fasta" #one file with ALL sequences
    SeqIO.write(all_records, all_fasta_path, "fasta")
    logger.info(f"All sequences written to {all_fasta_path}")

    #separate filtered sequences
    #no repeats only
    no_repeats = [r for r, s in zip(all_records, sequences) if s.get("blast_query") == "no_repeats"]
    if no_repeats:
        nr_path = output_dir / f"{gene_name}_{flanking_bp}bp_no_repeats.fasta"
        SeqIO.write(no_repeats, nr_path, "fasta")
        logger.info(f"No repeats sequences written to {nr_path}")

    #boundary issues
    for status in ("CRITICAL", "WARNING"):
        boundary_records = [r for r, s in zip(all_records, sequences) if s.get("boundary_status") == status]
        if boundary_records:
            boundary_path = output_dir / f"{gene_name}_{flanking_bp}bp_boundary_{status.lower()}.fasta"
            SeqIO.write(boundary_records, boundary_path, "fasta")
            logger.info(f"{status} boundary sequences written to {boundary_path}")
    return all_fasta_path

def write_promoter_fastas(sequences, output_dir, config):
    gene_name = config.gene_name
    promoter_sequences = [sequence for sequence in sequences if sequence.get("promoter_sequence")]
    if not promoter_sequences:
        return None

    records = []
    for sequence in promoter_sequences:
        header_parts = [
            f"gene={gene_name}",
            f"contig={sequence['contig_id']}",
            f"coords={sequence['promoter_start']}-{sequence['promoter_end']}",
            "coordinate_system=0-based_end-exclusive",
            f"promoter_length_requested={sequence['promoter_length_requested']}",
            f"promoter_length_extracted={len(sequence['promoter_sequence'])}",
            f"strand={sequence['strand']}"
        ]
        record = SeqRecord(
            Seq(sequence["promoter_sequence"]),
            id = f"{sequence['isolate_id']}_{sequence['gene_copy']}",
            description = " ".join(header_parts)
        )
        records.append(record)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok = True)
    output_path = output_dir / f"{gene_name}_promoter_sequences.fasta"
    SeqIO.write(records, output_path, "fasta")
    logger.info("Wrote %d promoter sequences to %s", len(records), output_path)
    return output_path


def extract_igr_from_coords(contig_seq, gene_left_0, gene_right_0, flank_coords, igr_type, strand):
    #extract intergenic region between gene and flanking gene
    #gene coordinates are 0-based, end-exclusive; flank coordinates are raw BLAST coordinates
    flank_start, flank_end = flank_coords
    flank_left_0 = min(flank_start, flank_end) - 1
    flank_right_0 = max(flank_start, flank_end)
    contig_length = len(contig_seq.seq)

    if strand not in ("forward", "reverse"):
        raise ValueError(f"Invalid strand: {strand}")
    if igr_type not in ("upstream", "downstream"):
        raise ValueError(f"Invalid IGR type: {igr_type}")
    if not (0 <= gene_left_0 < gene_right_0 <= contig_length):
        raise ValueError("Gene coordinates are outside the contig or invalid")
    if not (0 <= flank_left_0 < flank_right_0 <= contig_length):
        raise ValueError("Flank coordinates are outside the contig or invalid")

    if strand == "forward":
        if igr_type == "upstream":
            start, end = flank_right_0, gene_left_0
        else:
            start, end = gene_right_0, flank_left_0
    else:
        if igr_type == "upstream":
            start, end = gene_right_0, flank_left_0
        else:
            start, end = flank_right_0, gene_left_0

    if start >= end:
        return None, 0

    igr_seq = contig_seq.seq[start:end]
    if strand == "reverse":
        igr_seq = igr_seq.reverse_complement()

    return str(igr_seq), len(igr_seq)

def write_igr_fastas(upstream_igrs, downstream_igrs, output_dir, config):
    gene_name = config.gene_name

    for igr_type, igrs in (("upstream", upstream_igrs), ("downstream", downstream_igrs)):
        if not igrs:
            continue

        records = []
        for igr in igrs:
            header_parts = [
                f"isolate_id={igr['isolate_id']}",
                f"gene={gene_name}",
                f"region={igr_type}_IGR",
                f"contig={igr['contig_id']}",
                f"length={len(igr['sequence'])}bp",
                f"strand={igr['strand']}"
            ]
            record = SeqRecord(
                Seq(igr["sequence"]),
                id = f"{igr['isolate_id']}_{igr['gene_copy']}",
                description = " ".join(header_parts)
            )
            records.append(record)

        path = Path(output_dir) / f"{gene_name}_{igr_type}_IGR.fasta"
        SeqIO.write(records, path, "fasta")
        logger.info(f"Wrote {len(records)} {igr_type} IGR sequences to {path}")

def write_multicopy_fastas(sequences, output_dir, config):
    #group retained loci by isolate
    sequences_by_isolate = {}
    for seq in sequences:
        sequences_by_isolate.setdefault(seq["isolate_id"], []).append(seq)

    multicopy_isolates = {
        isolate_id: copies
        for isolate_id, copies in sequences_by_isolate.items()
        if len(copies) > 1
    }
    if not multicopy_isolates:
        logger.info("No isolates with multiple retained gene copies")
        return None

    output_dir = Path(output_dir)
    output_dir.mkdir(parents = True, exist_ok = True)
    all_records = []

    for isolate_id, copies in multicopy_isolates.items():
        records = []
        for seq in copies:
            header_parts = [
                f"isolate_id={seq['isolate_id']}",
                f"copy={seq['gene_copy']}",
                f"contig={seq['contig_id']}",
                f"coords={seq['extract_start']}-{seq['extract_end']}",
                "coordinate_system=0-based_end-exclusive",
                f"strand={seq['strand']}",
                f"boundary_status={seq['boundary_status']}"
            ]
            records.append(SeqRecord(
                Seq(seq["sequence"]),
                id = f"{seq['isolate_id']}_{seq['gene_copy']}",
                description = " ".join(header_parts)
            ))

        path = output_dir / f"{config.gene_name}_{isolate_id}_multicopy.fasta"
        SeqIO.write(records, path, "fasta")
        all_records.extend(records)

    combined_path = output_dir / f"{config.gene_name}_multicopy_all_sequences.fasta"
    SeqIO.write(all_records, combined_path, "fasta")
    logger.info(
        "Saved %d retained loci from %d multicopy isolates to %s",
        len(all_records), len(multicopy_isolates), output_dir
    )
    return combined_path

def extract_sequences(blast_results, genome_dir, output_dir, config):
    #extract sequences for all BLAST hits with QC
    extraction_params = config.extraction_params
    repeat_location = config.repeat_location.lower()
    logger.info(f"Extracting sequences for {len(blast_results)} BLAST hits...")

    genome_seqs = load_genome_sequences(genome_dir)
    extracted = []

    for hit in blast_results:
        contig_id = hit["contig_id"]

        if contig_id not in genome_seqs:
            logger.warning(f"Contig {contig_id} not found in genome sequences. Skipping hit")
            continue

        if repeat_location == "promoter":
            seq_data = extract_promoter_and_gene(
                contig_seq = genome_seqs[contig_id],
                start = hit["s_start"],
                end = hit["s_end"],
                strand = hit["strand"],
                promoter_length = extraction_params["promoter_length"]
            )
        else:
            seq_data = extract_sequence_from_coords(
                contig_seq = genome_seqs[contig_id],
                start = hit["s_start"],
                end = hit["s_end"],
                strand = hit["strand"],
                flanking_bp = extraction_params["flanking_bp"]
            )

        quality = check_sequence_quality(
            seq_data,
            boundary_critical = extraction_params["boundary_critical_distance"],
            boundary_warning = extraction_params["boundary_warning_distance"]
        )

        extracted_seq = {
            "isolate_id": hit["isolate_id"],
            "gene_copy": hit["gene_copy"],
            "contig_id": contig_id,
            "original_contig_id": hit.get("contig_orig", contig_id),
            "blast_identity": hit["identity"],
            "blast_coverage": hit["coverage"],
            "sequence": seq_data["sequence"],
            "strand": seq_data["strand"],
            "boundary_status": quality["boundary_status"],
            "flags": quality["flags"],
            "contig_length": seq_data["contig_length"],
            "gene_start": seq_data["gene_left_0"],
            "gene_end": seq_data["gene_right_0"],
            "extract_start": seq_data["extract_start"],
            "extract_end": seq_data["extract_end"],
            "dist_to_contig_start": seq_data["dist_to_contig_start"],
            "dist_to_contig_end": seq_data["dist_to_contig_end"],
            "genomic_context": hit.get("genomic_context", "unknown"),
            "blast_query": hit.get("blast_query", "unknown"),
            "upstream_found": hit.get("upstream_found", False),
            "downstream_found": hit.get("downstream_found", False),
            "upstream_flank_coords": hit.get("upstream_flank_coords"),
            "downstream_flank_coords": hit.get("downstream_flank_coords")
        }

        if repeat_location == "promoter":
            extracted_seq.update({
                "promoter_sequence": seq_data["promoter_sequence"],
                "promoter_length_requested": seq_data["promoter_length_requested"],
                "promoter_length_extracted": seq_data["promoter_length_extracted"],
                "promoter_start": seq_data["promoter_start"],
                "promoter_end": seq_data["promoter_end"]
            })

        extracted.append(extracted_seq)

    write_extracted_fasta(extracted, output_dir, config)
    if repeat_location == "promoter":
        write_promoter_fastas(extracted, output_dir, config)
    write_multicopy_fastas(extracted, output_dir / "multicopy", config)

    before_N_filtering = extracted
    extracted = remove_ambiguous_bases(extracted)

    boundary_counts = {}
    for seq in extracted:
        boundary_counts[seq["boundary_status"]] = boundary_counts.get(seq["boundary_status"], 0) + 1

    logger.info("Extraction complete.")
    for boundary_status, count in boundary_counts.items():
        logger.info(f" {count} sequences with boundary status {boundary_status}")

    upstream_igrs = []
    downstream_igrs = []

    if extraction_params.get("extract_igr", False):
        logger.info("Extracting intergenic regions (IGRs)...")
        for seq in extracted:
            contig_seq = genome_seqs[seq["contig_id"]]

            for igr_type, igrs in (("upstream", upstream_igrs), ("downstream", downstream_igrs)):
                flank_coords = seq.get(f"{igr_type}_flank_coords")
                if not seq[f"{igr_type}_found"] or flank_coords is None:
                    continue

                igr_seq, igr_len = extract_igr_from_coords(
                    contig_seq = contig_seq,
                    gene_left_0 = seq["gene_start"],
                    gene_right_0 = seq["gene_end"],
                    flank_coords = flank_coords,
                    igr_type = igr_type,
                    strand = seq["strand"]
                )
                if igr_seq:
                    igrs.append({
                        "isolate_id": seq["isolate_id"],
                        "gene_copy": seq["gene_copy"],
                        "contig_id": seq["contig_id"],
                        "sequence": igr_seq,
                        "length": igr_len,
                        "strand": seq["strand"]
                    })

        write_igr_fastas(upstream_igrs, downstream_igrs, Path(output_dir) / "igr", config)
        logger.info(f"IGR extraction complete: {len(upstream_igrs)} upstream, {len(downstream_igrs)} downstream IGRs")

    return extracted, {"upstream": upstream_igrs, "downstream": downstream_igrs}, before_N_filtering