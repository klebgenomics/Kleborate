
"""
Copyright 2025 Mary Maranga (gathonimaranga@gmail.com)
https://github.com/klebgenomics/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""
from ...shared.alignment import align_query_to_ref

def classify_shigella(assembly, minimap2_index, ShigellaRef, min_identity, min_coverage, shigella_serotype_markers):
    """
    Function to classify a strain as Shigella based on serotype markers.

    Parameters:
    - assembly: The genome assembly file 
    - minimap2_index: The Minimap2 index to align the sequences.
    - ShigellaRef: Reference file containing serotype marker sequences.
    - min_identity: Minimum identity for the alignment.
    - min_coverage: Minimum coverage for the alignment.
    - shigella_serotype_markers: Dictionary containing serotype markers for different Shigella serotypes.

    Returns:
    - Shigella serotype
    """
    # Align the assembly to the reference file to identify potential hits

    shigella_hits = align_query_to_ref(
        ShigellaRef,
        assembly,
        ref_index=minimap2_index,
        min_identity=min_identity,
        min_query_coverage=min_coverage
        )

    # Extract hit names from the hit entries
    hits = []
    for hit_entry in shigella_hits:
        hit_name = hit_entry.query_name.split(':')[0]  # Extract the gene or hit name
        hits.append(hit_name)

    # Check if the strain has the necessary markers to identify as Shigella
    if "ipaH_c" not in hits:
        return "-"

    # Determine the serotype based on the markers present
    for serotype, markers in shigella_serotype_markers.items():
        if all(marker in hits for marker in markers):
            return serotype

    # If no serotype matches
    return "-"

