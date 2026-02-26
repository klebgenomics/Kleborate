"""
Copyright 2025 Mary Maranga (gathonimaranga@gmail.com)
https://github.com/klebgenomics/Kleborate

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""

from ...shared.alignment import align_query_to_ref, cull_redundant_hits


def extract_fasta_headers(fasta_file_path):
    """
    Extract unique headers from a FASTA file.
    """
    headers = set()
    with open(fasta_file_path, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                header = line[1:].strip()
                factor = header.split(':')[0]
                headers.add(factor)
    return headers


def map_virulence_factors(assembly, minimap2_index, ref_file, min_identity, min_coverage):
    """
    Aligns assembled genomes to the virulence alleles and reports presence/absence.

    Parameters:
    - assembly: Assembly in FASTA format (path or file-like accepted by align_query_to_ref).
    - ref_file: Virulence factors in FASTA format (path or file-like).
    - min_identity: Minimum identity percentage for alignment.
    - min_coverage: Minimum query coverage for alignment.
    - minimap2_index: Path to the assembly's minimap2 index for faster alignment (optional).

    Returns:
    - virulence_markers: Dictionary {factor_header: 'present' or '-'}
    """

    def extract_present_vfs(alignment_hits):
        present = set()
        for hit in alignment_hits:
            if isinstance(hit, str):
                factor = hit.split(':')[0]
            else:
                factor = str(hit).split(':')[0]
            present.add(factor)
        return present

    alignment_hits = align_query_to_ref(
        ref_file,
        assembly,
        ref_index=minimap2_index,
        min_identity=min_identity,
        min_query_coverage=min_coverage
    )
    alignment_hits = cull_redundant_hits(alignment_hits)

    # Build dictionary of virulence markers for all factors in the reference file
    all_headers = extract_fasta_headers(ref_file)
    present_vfs = extract_present_vfs(alignment_hits)
    virulence_markers = {
        header: 'present' if header in present_vfs else '-'
        for header in all_headers
    }

    return virulence_markers
