"""
Define gene specific parameters to use in the pipeline.

Duplicate this script with create_config.py
"""

gene_name = ""
gene_id = "" #optional - if running multiple genes with long names can use a shorter identifier for folder names 

repeat_location = "" #CDS or promoter
repeat_motif = "" #e.g. CAAT

#gene sequences for BLAST searches
ref_full = "" #full gene sequence
ref_2_repeats = "" #gene sequence normalised to 2 repeats
ref_no_repeats = "" #gene sequence from end of repeats to the end of the gene
ref_gene_start = None #gene sequence from the start of the gene to the start of the repeats - also OPTIONAL if gene start is too short a sequence

ref_gene_end = None #OPTIONAL - define a region at towards the end of the gene, helpful if really long gene and if not None then will BLAST this region as well. Can be compared with gene start searches, may highlight if the gene is being broken up across contigs

upstream_flank = ""
downstream_flank = ""

forward_primer = ""
reverse_primer = ""

ref_full = ref_full.upper()
ref_2_repeats = ref_2_repeats.upper()
ref_no_repeats = ref_no_repeats.upper()
if ref_gene_start is not None:
    ref_gene_start = ref_gene_start.upper()
if ref_gene_end is not None:
    ref_gene_end = ref_gene_end.upper()
upstream_flank = upstream_flank.upper()
downstream_flank = downstream_flank.upper()
forward_primer = forward_primer.upper()
reverse_primer = reverse_primer.upper()

upstream_distance = 0 #bp distance
downstream_distance = 0 #bp distance
flank_buffer = 500 #added to flanking distances to allow some variability - used to check that flanking genes are where you expect them to be


#BLAST search parameters
gene_blast_params = {
    "min_identity": 80,
    "high_identity": 90,
    "min_coverage": 70,
    "high_coverage": 90,
    "max_hits_per_isolate": 10,
    "max_target_seqs": 1000000,
    "max_hsps": 10,
    "word_size": 11,
    "reward": 2,
    "penalty": -3,
    "gapopen": 5,
    "gapextend": 2,
    "evalue": 1e-10
}

gene_start_blast_params = {
    "min_identity": 80,
    "high_quality_identity": 95,
    "min_coverage": 70,
    "word_size": 7,
    "evalue": 1000,
    "max_hits_per_isolate": 5,
    "max_target_seqs": 1000000,
}

primer_blast_params = {
    "min_identity": 80,
    "min_coverage": 70,
    "max_mismatches": 4,
    "word_size": 7,
    "evalue": 1000,
    "min_alignment_pct": 90,
    "max_hits_per_isolate": 5,
    "max_target_seqs": 1000000,
    "max_pcr_fragment_length": 1000 #in silico PCR and gene dupes mean you get some 'fragments' between the genes which wouldn't happen in vivo so filter them out
}

extraction_params = {
    "flanking_bp": 1000, #how much to extract either side of BLAST coords before defining gene start and end
    "boundary_critical_distance": 25, #how close to the end of a contig
    "boundary_warning_distance": 50,
    "atg_search_window": 100, #integer max distance, or (minimum_distance, maximum_distance) upstream of repeat
    "promoter_length": 0, #length of promoter region to extract
    "extract_igr": False #TODO - extract intergenic between GOI and flanks
}

alignment_params = {
    "generate_alignments": False, #TODO: aligments massively slow down processing time so taken out for now but may be re-added later - for now just save the sequences and unique sequences 
    "tool": "muscle",
    "pcr_fragments": True,
    "allele_sequences": True, #only the unique sequences 
    "gene_regions": True, #upstream and downstream of the repeat
    "promoter": False
}

translation_params = {
    "codon_table": 11,
    "min_protein_length": 100 #minimum length of a protein to be considered ON
}
