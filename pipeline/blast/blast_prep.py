import io
import logging
import zipfile
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

logger = logging.getLogger(__name__)

fasta_extensions = (".fasta", ".fa", ".fna", ".fas")

def create_combined_genome_fasta(genome_dir, output_dir):
    #combine all genomes into single fasta file for blastdb
    combined_fasta = Path(output_dir) / "combined_genomes.fasta"
    genome_dir = Path(genome_dir)

    contig_id_map = {}
    isolate_ids = set()
    n_contigs = 0
    n_isolates = 0

    with open(combined_fasta, "w") as combined_handle:
        if genome_dir.suffix == ".zip":
            with zipfile.ZipFile(genome_dir, "r") as zf:
                genome_entries = [
                    entry
                    for entry in zf.namelist()
                    if entry.lower().endswith(fasta_extensions)
                    and not entry.rsplit("/", 1)[-1].startswith("._")
                ]

                logger.info("Combining %d genome entries from %s", len(genome_entries), genome_dir)

                for entry in genome_entries:
                    isolate_id = Path(entry).stem
                    isolate_ids.add(isolate_id)
                    n_isolates += 1

                    with zf.open(entry) as f:
                        text_handle = io.TextIOWrapper(f, encoding = "utf-8")
                        for contig_no, record in enumerate(SeqIO.parse(text_handle, "fasta"), start = 1): 
                            original_id = record.id
                            internal_id = f"{isolate_id}_ctg{contig_no:06d}"
                            contig_id_map[internal_id] = original_id
                            record.id = internal_id
                            record.name = internal_id
                            record.description = ""
                            SeqIO.write(record, combined_handle, "fasta")
                            n_contigs += 1

        elif Path(genome_dir).is_dir():
            genome_files = sorted(
                path
                for path in Path(genome_dir).iterdir()
                if path.suffix.lower() in fasta_extensions
            )

            logger.info("Combining %d genome files from %s", len(genome_files), genome_dir)

            for genome_file in genome_files:
                isolate_id = genome_file.stem
                isolate_ids.add(isolate_id)
                n_isolates += 1

                with open(genome_file, "r", encoding="utf-8") as f:
                    for contig_no, record in enumerate(SeqIO.parse(f, "fasta"), start=1):
                        original_id = record.id
                        internal_id = f"{isolate_id}_ctg{contig_no:06d}"
                        contig_id_map[internal_id] = original_id
                        record.id = internal_id
                        record.name = internal_id
                        record.description = ""
                        SeqIO.write(record, combined_handle, "fasta")
                        n_contigs += 1

            logger.info("Combined %d contigs from %d isolates into %s", n_contigs, n_isolates, combined_fasta)

    return combined_fasta, contig_id_map, sorted(isolate_ids)

def make_blast_db(fasta_path):
    logger.info("Creating BLAST database for %s", fasta_path)

    db_prefix = Path(fasta_path).parent / "blast_db"

    cmd = [
        "makeblastdb",
        "-in", str(fasta_path),
        "-dbtype", "nucl",
        "-out", str(db_prefix)
    ]

    try:
        result = subprocess.run(cmd, check = True, capture_output = True, text = True)
        logger.info("BLAST database created successfully")
        return db_prefix
    except subprocess.CalledProcessError as e:
        logger.error("Failed to create BLAST database: %s", e.stderr)
        raise

def clean_sequence(sequence, sequence_name):
    #remove whitespace, convert to uppercase and check for invalid characters
    sequence = "".join(sequence.split()).upper()
    valid_bases = set("ACGTRYSWKMBDHVN") #all IUPAC DNA coddes
    invalid_bases = set(sequence) - valid_bases
    if not sequence:
        raise ValueError(f"Sequence for {sequence_name} is empty")

    if invalid_bases:
        raise ValueError(f"Sequence for {sequence_name} contains invalid bases: {', '.join(invalid_bases)}")

    return sequence


def create_cds_gene_records(config):
    gene_name = config.gene_name
    logger.info("Creating CDS reference records for %s", gene_name)

    ref_2_repeats = clean_sequence(config.ref_2_repeats, f"{gene_name} ref_2_repeats")
    ref_no_repeats = clean_sequence(config.ref_no_repeats, f"{gene_name} ref_no_repeats")
    upstream_flank = clean_sequence(config.upstream_flank, f"{gene_name} upstream_flank")
    downstream_flank = clean_sequence(config.downstream_flank, f"{gene_name} downstream_flank")

    records = [
        SeqRecord(Seq(ref_2_repeats), id=f"{gene_name}_ref_2_repeats", description="Full gene normalised to 2 repeats"),
        SeqRecord(Seq(ref_no_repeats), id=f"{gene_name}_ref_no_repeats", description="Gene sequence after repeat tract"),
        SeqRecord(Seq(upstream_flank), id=f"{gene_name}_upstream_flank", description="Upstream flanking sequence"),
        SeqRecord(Seq(downstream_flank), id=f"{gene_name}_downstream_flank", description="Downstream flanking sequence"),
    ]
    logger.info("Created %d CDS reference records for %s", len(records), gene_name)
    return records

def create_promoter_gene_records(config):
    gene_name = config.gene_name
    logger.info("Creating promoter reference records for %s", gene_name)

    ref_full = clean_sequence(config.ref_full, f"{gene_name} ref_full")
    upstream_flank = clean_sequence(config.upstream_flank, f"{gene_name} upstream_flank")
    downstream_flank = clean_sequence(config.downstream_flank, f"{gene_name} downstream_flank")

    records = [
        SeqRecord(Seq(ref_full), id=f"{gene_name}_ref_full", description="Full reference sequence"),
        SeqRecord(Seq(upstream_flank), id=f"{gene_name}_upstream_flank", description="Upstream flanking sequence"),
        SeqRecord(Seq(downstream_flank), id=f"{gene_name}_downstream_flank", description="Downstream flanking sequence"),
    ]

    logger.info("Created %d promoter reference records for %s", len(records), gene_name)
    return records

def create_gene_ref_fasta(configs, output_dir):
    records = []

    for config in configs:
        repeat_location = config.repeat_location.lower()

        if repeat_location == "cds":
            records.extend(create_cds_gene_records(config))
        elif repeat_location == "promoter":
            records.extend(create_promoter_gene_records(config))
        else:
            raise ValueError(f"Invalid repeat_location for {config.gene_name}: {repeat_location}")

    ref_fasta_path = Path(output_dir) / "gene_references.fasta"
    SeqIO.write(records, ref_fasta_path, "fasta")

    logger.info("Saved %d gene reference records to %s", len(records), ref_fasta_path)
    return ref_fasta_path

def create_gene_start_ref_fasta(configs, output_dir):
    records = []

    for config in configs:
        if not getattr(config, "ref_gene_start", None):
            continue

        gene_name = config.gene_name
        ref_gene_start = clean_sequence(config.ref_gene_start, f"{gene_name} ref_gene_start")
        records.append(SeqRecord(Seq(ref_gene_start), id=f"{gene_name}_ref_gene_start", description="Gene sequence before repeat tract"))

    ref_fasta_path = Path(output_dir) / "gene_start_references.fasta"
    if not records:
        return None
    SeqIO.write(records, ref_fasta_path, "fasta")

    logger.info("Saved %d gene-start reference records to %s", len(records), ref_fasta_path)
    return ref_fasta_path

def create_gene_end_ref_fasta(config, output_dir):
    """Create the optional gene-end reference FASTA for one gene."""
    if not config.ref_gene_end:
        return None

    gene_name = config.gene_name
    record = SeqRecord(
        Seq(clean_sequence(config.ref_gene_end, f"{gene_name} ref_gene_end")),
        id=f"{gene_name}_ref_gene_end",
        description="Reference gene end sequence",
    )
    ref_fasta_path = Path(output_dir) / "gene_end_references.fasta"
    SeqIO.write([record], ref_fasta_path, "fasta")
    logger.info("Saved gene-end reference for %s", gene_name)
    return ref_fasta_path

def create_primer_ref_fasta(configs, output_dir):
    records = []

    for config in configs:
        gene_name = config.gene_name
        forward_primer = clean_sequence(config.forward_primer, f"{gene_name} forward_primer")
        reverse_primer = clean_sequence(config.reverse_primer, f"{gene_name} reverse_primer")

        records.extend([
            SeqRecord(Seq(forward_primer), id=f"{gene_name}_forward_primer", description=""),
            SeqRecord(Seq(reverse_primer), id=f"{gene_name}_reverse_primer", description=""),
        ])

    ref_fasta_path = Path(output_dir) / "primer_references.fasta"
    SeqIO.write(records, ref_fasta_path, "fasta")

    logger.info("Saved %d primer reference records to %s", len(records), ref_fasta_path)
    return ref_fasta_path