"""
This file carries out multi-MLST functions: similar to the regular MLST (found in mlst.py) but
allowing for multiple STs per genome. This is useful for some virulence loci which can appear more
than once per genome (e.g. on the chromosome and on a plasmid).

Copyright 2025 Kat Holt, Mary Maranga, Ryan Wick
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""
import ast
import re
from .alignment import align_query_to_ref
from .mlst import load_st_profiles, run_single_mlst
from .alignment import truncation_check, cull_redundant_hits
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
from .misc import load_fasta, reverse_complement


def multi_mlst(assembly_path, minimap2_index, profiles_path, allele_paths, gene_names, extra_info,
               min_identity, min_coverage, required_exact_matches, check_for_truncation=False,
               report_incomplete=False, min_spurious_identity=None, min_spurious_coverage=None,
               unknown_group_name=None, min_gene_count=None):
    """
    This function takes and returns the same things as the mlst function in mlst.py. However, it
    will look for cases where multiple contigs have hits for the full set of MLST genes, and in
    that case, MLST is run on each of them. Otherwise, it behaves like normal MLST.
    """
    profiles = load_st_profiles(profiles_path, gene_names, extra_info)
    
    if min_spurious_coverage is not None:
        hits_per_gene = {g: align_query_to_ref(allele_paths[g], assembly_path,
                                               ref_index=minimap2_index, min_identity=min_spurious_identity,
                                               min_query_coverage=min_spurious_coverage) for g in gene_names}

        spurious_hits = {g: [h for h in hits_per_gene[g] 
                             if h.query_cov < min_coverage and h.percent_identity < min_identity] for g in gene_names}
    else:
        hits_per_gene = {g: align_query_to_ref(allele_paths[g], assembly_path,
                                               ref_index=minimap2_index, min_identity=min_identity,
                                               min_query_coverage=min_coverage) for g in gene_names}
        spurious_hits = None


    if spurious_hits is not None:
        spurious_hits = {g: process_spurious_hits(spurious_hits[g]) for g in gene_names}


    hits_by_contig = cluster_hits_by_contig(hits_per_gene, gene_names)
    full_set_contigs = find_full_set_contigs(hits_by_contig)

    # If zero or one contigs have the full set of genes, then this is treated as a non-multi-MLST
    # case, i.e. the same as regular MLST.

    if len(full_set_contigs) < 2:
        # Apply the unknown group logic here for single-contig cases
        return run_single_mlst(
            profiles, hits_per_gene, gene_names, required_exact_matches, check_for_truncation, 
            report_incomplete, unknown_group_name, min_gene_count), spurious_hits, hits_per_gene

    # If more than one contig has the full set of genes, then this is treated as a multi-MLST case,
    # where each full-set contig gets an MLST call.
    contig_results = {}
    for contig in full_set_contigs:
        contig_results[contig] = run_single_mlst(
            profiles, hits_by_contig[contig], gene_names, required_exact_matches, check_for_truncation,
            report_incomplete, unknown_group_name, min_gene_count)
         
    return combine_results(full_set_contigs, contig_results, gene_names), spurious_hits, hits_per_gene


# def multi_mlst(assembly_path, minimap2_index, profiles_path, allele_paths, gene_names, extra_info,
#                min_identity, min_coverage, required_exact_matches, check_for_truncation=False,
#                report_incomplete=False):
#     """
#     This function takes and returns the same things as the mlst function in mlst.py. However, it
#     will look for cases where multiple contigs have hits for the full set of MLST genes, and in
#     that case, MLST is run on each of them. Otherwise, it behaves like normal MLST.
#     """
#     profiles = load_st_profiles(profiles_path, gene_names, extra_info)
#     hits_per_gene = {g: align_query_to_ref(allele_paths[g], assembly_path,
#                                            ref_index=minimap2_index, min_identity=min_identity,
#                                            min_query_coverage=min_coverage) for g in gene_names}
#     hits_by_contig = cluster_hits_by_contig(hits_per_gene, gene_names)
#     full_set_contigs = find_full_set_contigs(hits_by_contig)

#     # If zero or one contigs have the full set of genes, then this is treated as a non-multi-MLST
#     # case, i.e. the same as regular MLST.
#     if len(full_set_contigs) < 2:
#         return run_single_mlst(profiles, hits_per_gene, gene_names, required_exact_matches,
#                                check_for_truncation, report_incomplete)

#     # If more than one contig has the full set of genes, then this is treated as a multi-MLST case,
#     # where each full-set contig gets an MLST call.
#     contig_results = {}
#     for contig in full_set_contigs:
#         contig_results[contig] = run_single_mlst(profiles, hits_by_contig[contig], gene_names,
#                                                  required_exact_matches, check_for_truncation,
#                                                  report_incomplete)
#     return combine_results(full_set_contigs, contig_results, gene_names)


def cluster_hits_by_contig(hits_per_gene, gene_names):
    """
    Takes a dictionary of hits per gene (key: gene name, value: list of hits) and returns the same
    hits grouped by the contig name (key: contig name, value: gene name to list of hits dict).
    """
    all_contig_names = set()
    for gene, hits in hits_per_gene.items():
        all_contig_names.update([h.ref_name for h in hits])

    hits_by_contig = {}
    for contig in sorted(all_contig_names):
        hits_by_contig[contig] = {g: [] for g in gene_names}

    for gene, hits in hits_per_gene.items():
        for h in hits:
            hits_by_contig[h.ref_name][gene].append(h)
    return hits_by_contig


def find_full_set_contigs(hits_by_contig):
    """
    Returns the names of contigs which have a full set of MLST genes. Contigs are sorted by length
    from big to small (so chromosomal contigs come before plasmid contigs in a completed assembly).
    """
    contig_lengths = {}
    full_set_contigs = []
    for contig, hits_per_gene in hits_by_contig.items():
        if all(len(hits) > 0 for hits in hits_per_gene.values()):
            full_set_contigs.append(contig)
        for hits in hits_per_gene.values():
            for h in hits:
                contig_lengths[h.ref_name] = h.ref_length
    return sorted(full_set_contigs, key=lambda c: contig_lengths[c], reverse=True)


def combine_results(full_set_contigs, contig_results, gene_names):
    """
    Takes per-contig MLST results (ST, extra info and alleles) and combines them into a single
    comma-delimited list.
    """
    combined_st, combined_extra_info = [], []
    combined_alleles = {g: [] for g in gene_names}
    for contig in full_set_contigs:
        st, extra_info, alleles = contig_results[contig]
        combined_st.append(st)
        combined_extra_info.append(extra_info)
        for g in gene_names:
            combined_alleles[g].append(alleles[g])
    return ','.join(combined_st), ','.join(combined_extra_info), \
           {g: ','.join(combined_alleles[g]) for g in gene_names}


def get_allele_and_locus(hit):
    """
    Parses the allele name and locus name from the hit's gene ID.
    """
    if '__' in hit.query_name:  # srst2 formatted file
        gene_id_components = hit.query_name.split('__')
        locus = gene_id_components[1]
        allele = gene_id_components[2]
    else:
        allele = hit.query_name
        locus = hit.query_name.split('_')[0]
    return allele, locus


def process_spurious_hits(hits):
    hit_strings = []
    for hit in hits:
        allele, locus = get_allele_and_locus(hit)
        alignment_length = hit.query_end - hit.query_start
        if hit.percent_identity < 100.00 or alignment_length < hit.query_length:
            allele += '*'  # inexact match
        allele += truncation_check(hit)[0]
        hit_strings.append(allele)
    return hit_strings

def check_polyT_tract(hits_per_gene, assembly):
    """
    For rmpA hit, extracts the upstream sequence and
    checks for poly-T tract (G(T+)A).
    Returns 'untypable' if rmpA is found but the pattern is not matched.
    """
    poly_t_status_map = {
        8: "reduced expression",
        9: "reduced expression",
        10: "reduced expression",
        11: "ON",
        12: "ON",
        13: "ON",
        14: "ON",
        15: "ON",
        16: "ON"
    }
    
    assembly_seqs = dict(load_fasta(assembly))
    upstream_length = 80
    gene_prefix = "rmpA"
    results = []
    
    poly_t_pattern = re.compile(r"G(T+)A", re.IGNORECASE)

    rmpA_found = False

    for gene in hits_per_gene:
        hits_per_gene[gene] = cull_redundant_hits(hits_per_gene[gene])
        hits = hits_per_gene[gene]

        if gene.startswith(gene_prefix):
            rmpA_found = True
            
            for hit in hits:
                contig_start, contig_end = hit.ref_start, hit.ref_end
                full_contig_seq = assembly_seqs[hit.ref_name]

                # Extract upstream sequence ---
                if hit.strand == '+':
                    upstream_start = max(0, contig_start - upstream_length)
                    upstream_seq = full_contig_seq[upstream_start:contig_start]
                elif hit.strand == '-':
                    upstream_start = contig_end
                    upstream_end = contig_end + upstream_length
                    upstream_seq = reverse_complement(full_contig_seq[upstream_start:upstream_end])
                else:
                    continue
                
                # Search for G(T+)A and check length
                match = poly_t_pattern.search(upstream_seq)

                if match:
                    poly_t_string = match.group(1)
                    poly_t_length = len(poly_t_string)

                    # Determine ON/OFF status
                    if poly_t_length in poly_t_status_map:
                        status_results = f"{poly_t_length}T ({poly_t_status_map[poly_t_length]})"
                        if poly_t_length >= 11 and poly_t_status_map[poly_t_length] == "ON":
                            status_results = f"{poly_t_length}T"
                        results.append(status_results)
                        
    if results:
        return results[0]
    elif rmpA_found:
        # rmpA gene is present, but no exact match to polyT pattern 
        return 'untypable'
    else:
        return '-'


def check_argR_status(assembly, argR_ref, min_identity, min_coverage):
    """
    searches for the argR gene in the assembly,
    and checks for ArgR gene or truncation of the encoded ArgR protein.

    Returns:
        list: argR_status 
    """
    argR_status = []
    
    # Search for the argR gene in the assembly 
    argR_aln = align_query_to_ref(
        argR_ref, 
        assembly, 
        min_query_coverage=min_coverage,
        min_identity=min_identity
    )

    if not argR_aln:
        argR_status.append('-')
        # argR_status.append('argR missing')
    else:
        hit = argR_aln[0]
        _, coverage, _ = truncation_check(hit)
        # print(coverage)
        if coverage >= 100.0:
            argR_status.append('present')
        else:
            argR_status.append('truncated-' + ('%.0f' % coverage) + '%')
    return argR_status


def check_argR_box(hits_per_gene, assembly):

    """
    Checks for presence of the ARG-box upstream of each rmpA gene.
    """
    assembly_seqs = dict(load_fasta(assembly))
    upstream_length = 150
    gene_prefix = "rmpA"
    ARG_box = 'ATTGAATTTTTATTCATT'
    results = [] 
    for gene in hits_per_gene:
        hits_per_gene[gene] = cull_redundant_hits(hits_per_gene[gene])
        hits = hits_per_gene[gene]
        if gene.startswith(gene_prefix):
            found = False
            for hit in hits:
                hit_seq = hit.ref_seq
                contig_start, contig_end = hit.ref_start, hit.ref_end
                contig_length = len(assembly_seqs[hit.ref_name])
                gene_nucl_seq = assembly_seqs[hit.ref_name][contig_start:contig_end]

                if hit.strand == '-':
                    gene_nucl_seq = reverse_complement(gene_nucl_seq)
                assert hit_seq == gene_nucl_seq

                full_contig_seq = assembly_seqs[hit.ref_name]

                if hit.strand == '+':
                    upstream_start = max(0, contig_start - upstream_length)
                    upstream_end = contig_start
                    upstream_seq = full_contig_seq[upstream_start:upstream_end]
                elif hit.strand == '-':
                    upstream_start = contig_end
                    upstream_end = min(contig_end + upstream_length, contig_length)
                    upstream_seq = reverse_complement(full_contig_seq[upstream_start:upstream_end])

                if ARG_box in upstream_seq:
                    break
            else:
                results.append('ARG-box lost')
    if results:
        return results[0]
    else:
        return '-'


def poly_G_variation(hits_per_gene):
    """
    In silico extension of Poly G tract in rmpA
    """
    loci_status = []
    COVERAGE_THRESHOLD = 95.0
    window_start_idx = 276 - 1  
    window_end_idx = 285        

    for gene, hits in hits_per_gene.items():
        hits = cull_redundant_hits(hits)
        if gene.startswith("rmpA"):
            for hit in hits:
                seq = hit.ref_seq
                # print(hit.ref_name, seq)
                org_aa_length = len(translate_nucl_to_prot(seq))
                # print(seq)
                if org_aa_length == 0:
                 org_aa_length = len(seq) // 3

                if org_aa_length == 0:
                    continue
                window_seq = seq[window_start_idx:window_end_idx]

                for match in re.finditer(r"(G{5,})", window_seq):
                    abs_start = window_start_idx + match.start()
                    abs_end = window_start_idx + match.end()
                    length = match.end() - match.start()  

                    # Build extended sequences
                    ext_seq_plus1 = seq[:abs_start] + ("G" * (length + 1)) + seq[abs_end:]
                    ext_seq_plus2 = seq[:abs_start] + ("G" * (length + 2)) + seq[abs_end:]

                    # Translate
                    translation_plus1 = translate_nucl_to_prot(ext_seq_plus1)
                    translation_plus2 = translate_nucl_to_prot(ext_seq_plus2)

                    # Coverage percentages relative to query AA length
                    pcov_plus1 = 100.0 * (len(translation_plus1) / org_aa_length)
                    pcov_plus2 = 100.0 * (len(translation_plus2) / org_aa_length)

                    # +1 G
                    if pcov_plus1 >= COVERAGE_THRESHOLD:
                        reannotated_status = "OFF"
                    # +2 G
                    elif pcov_plus2 >= COVERAGE_THRESHOLD:
                        reannotated_status = "OFF"
                    else:
                        reannotated_status = "-"
                    loci_status.append(reannotated_status)
                    return loci_status 

    if not loci_status:
        return "-"
    else:
        return loci_status



def poly_A_variation(hits_per_gene):
    """

    """
    loci_status = []
    COVERAGE_THRESHOLD = 95.0
    window_start_idx = 136 - 1  
    window_end_idx = 144        

    for gene, hits in hits_per_gene.items():
        hits = cull_redundant_hits(hits)
        if gene.startswith("rmpD"):
            for hit in hits:
                seq = hit.ref_seq 
                # print(hit.ref_name, seq)
                org_aa_length = len(translate_nucl_to_prot(seq))
                if org_aa_length == 0:
                 org_aa_length = len(seq) // 3

                if org_aa_length == 0:
                    continue

                window_seq = seq[window_start_idx:window_end_idx]

                for match in re.finditer(r"(A{4,})", window_seq):
                    abs_start = window_start_idx + match.start()
                    abs_end = window_start_idx + match.end()
                    length = match.end() - match.start()  

            
                    ext_seq_plus1 = seq[:abs_start] + ("A" * (length + 1)) + seq[abs_end:]
                    ext_seq_plus2 = seq[:abs_start] + ("A" * (length + 2)) + seq[abs_end:]

                    # Translate extended sequence
                    translation_plus1 = translate_nucl_to_prot(ext_seq_plus1)
                    translation_plus2 = translate_nucl_to_prot(ext_seq_plus2)


                    # Coverage percentages relative to original gene AA length
                    pcov_plus1 = 100.0 * (len(translation_plus1) / org_aa_length)
                    pcov_plus2 = 100.0 * (len(translation_plus2) / org_aa_length)

                    # try +1 A
                    if pcov_plus1 >= COVERAGE_THRESHOLD:
                        reannotated_status = "OFF"
                    # try +2 A
                    elif pcov_plus2 >= COVERAGE_THRESHOLD:
                        reannotated_status = "OFF"
                    else:
                        reannotated_status = "-"
                    loci_status.append(reannotated_status)
                    return loci_status

    if not loci_status:
        return "-"
    else:
        return loci_status



def poly_G_rmpC_variation(hits_per_gene):
    """
    In silico extension of Poly G tract IN rmpC 
    """
    loci_status = []
    COVERAGE_THRESHOLD = 95.0
    window_start_idx = 175 - 1  
    window_end_idx = 183        
    for gene, hits in hits_per_gene.items():
        hits = cull_redundant_hits(hits)
        if gene.startswith("rmpC"):
            for hit in hits:
                seq = hit.ref_seq
                org_aa_length = len(translate_nucl_to_prot(seq))
                if org_aa_length == 0:
                 org_aa_length = len(seq) // 3

                if org_aa_length == 0:
                    continue
                
                window_seq = seq[window_start_idx:window_end_idx]

                for match in re.finditer(r"(G{5,})", window_seq):
                    abs_start = window_start_idx + match.start()
                    abs_end = window_start_idx + match.end()
                    length = match.end() - match.start()  


                    # Extended sequences
                    ext_seq_plus1 = seq[:abs_start] + ("G" * (length + 1)) + seq[abs_end:]
                    ext_seq_plus2 = seq[:abs_start] + ("G" * (length + 2)) + seq[abs_end:]

                    # Translate
                    translation_plus1 = translate_nucl_to_prot(ext_seq_plus1)
                    translation_plus2 = translate_nucl_to_prot(ext_seq_plus2)

                    # Coverage percentages relative to query AA length
                    pcov_plus1 = 100.0 * (len(translation_plus1) / org_aa_length)
                    pcov_plus2 = 100.0 * (len(translation_plus2) / org_aa_length)

                    #  try +1 G
                    if pcov_plus1 >= COVERAGE_THRESHOLD:
                        reannotated_status = "OFF"
                    # try +2 G
                    elif pcov_plus2 >= COVERAGE_THRESHOLD:
                        reannotated_status = "OFF"
                    else:
                        reannotated_status = "-"
                    loci_status.append(reannotated_status)
                    return loci_status

    if not loci_status:
        return "-"
    else:
        return loci_status[0]


def allele_type(allele_value):
    if '%' in allele_value: 
        return "truncated" 
    # Check for Inexact 
    elif '*' in allele_value: 
        return "inexact"
    else:
        return "exact"


def process_status_dict(filepath, prefix):
    with open(filepath) as f:
        raw_dict = ast.literal_eval(f.read())
        status_dict = {
            re.sub(rf'^{prefix}_', '', k): v[0]
            for k, v in raw_dict.items()
        }
    return status_dict



def get_gene_status(allele_value, allele_dict, hits_per_gene, poly_variation_func):
    
    if allele_value in ('-', None):
        return "-"

    allele = str(allele_value).strip()
    if not allele or allele == "-":
        return "-"

    # Extract the allele call
    allele_id = re.split(r"[\*\-]", allele, maxsplit=1)[0].strip()

    # truncation
    truncation_matches = re.findall(r"(\d+(?:\.\d+)?)\s*%", allele)
    truncation_pct = f"{truncation_matches[-1]}%" if truncation_matches else None

    is_inexact_call = "*" in allele
    is_truncated_call = "%" in allele

    if is_inexact_call or is_truncated_call:
        poly_tract_variation = poly_variation_func(hits_per_gene)
        poly_tract_variation_status = (
            poly_tract_variation[0] if isinstance(poly_tract_variation, list) and poly_tract_variation else
            poly_tract_variation if poly_tract_variation is not None else
            "-"
        )

        # If the allele is reversibly off
        if poly_tract_variation_status!= "-" and "OFF" in str(poly_tract_variation_status):
            return f"{allele_id}(OFF)"

        # If the tract is non-reversible
        if is_truncated_call:
            return f"{allele_id}*-{truncation_pct}"

        if is_inexact_call:
            return f"{allele_id}*"

    return allele_dict.get(allele, "-")


def translate_nucl_to_prot(nucl_seq):
    ambiguous_bases = set(b for b in nucl_seq) - {'A', 'C', 'G', 'T'}
    for b in ambiguous_bases:
        nucl_seq = nucl_seq.split(b)[0]
    nucl_seq = nucl_seq[:len(nucl_seq) // 3 * 3]
    coding_dna = Seq(nucl_seq)
    return str(coding_dna.translate(table='Bacterial', to_stop=True))





