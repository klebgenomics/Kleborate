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

import os
import pathlib
import shutil
import sys
import csv
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, Tuple
from typing import Dict

from .pathovar import minimap_pathovar


def description():
    return 'Pathotyping of E. coli genomes'


def prerequisite_modules():
    return ['enterobacterales__species']


def get_headers():
    full_headers = ['Pathotype', 'Stx1', 'Stx2', 'ST', 'LT', 'eae', 'bfpA', 'ipaH']
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    module_name = os.path.basename(__file__)[:-3]
    group = parser.add_argument_group(f'{module_name} module')
    group.add_argument('--escherichia__pathovar_min_identity', type=float, default=90.0,
                       help='Minimum alignment percent identity for detecting virulence factors')
    group.add_argument('--escherichia__pathovar_min_coverage', type=float, default=80.0,
                       help='Minimum alignment percent coverage for detecting virulence factors')

    # ShigaPass integration
    default_shigapass_dir = data_dir()
    default_db_dir = default_shigapass_dir / 'ShigaPass_DataBases'

    group.add_argument(
        '--shigapass_dir',
        type=str,
        default=str(default_shigapass_dir),
        help='Path to the directory containing ShigaPass.sh and/or SCRIPT/ShigaPass.sh'
    )
    group.add_argument(
        '--shigapass_db_dir',
        type=str,
        default=str(default_db_dir),
        help='Path to the ShigaPass databases directory (option -p)'
    )
    group.add_argument(
        '--shigapass_threads',
        type=int,
        default=2,
        help='Threads to pass to ShigaPass (-t)'
    )
    group.add_argument(
        '--shigapass_keep',
        action='store_true',
        help='Keep ShigaPass intermediate directories (-k)'
    )
    group.add_argument(
        '--shigapass_init_db',
        action='store_true',
        help='Initialise BLAST databases once on first run (-u)'
    )
    group.add_argument(
        '--shigapass_outdir_base',
        type=str,
        default='ShigaPass_Results',
        help='Base name for ShigaPass output directory (a per-sample suffix is added)'
    )

    return group



def check_cli_options(args):
    if args.escherichia__pathovar_min_identity <= 50.0 or args.escherichia__pathovar_min_identity >= 100.0:
        sys.exit('Error: --escherichia__pathovar_min_identity must be between 50.0 and 100.0')
    if args.escherichia__pathovar_min_coverage <= 50.0 or args.escherichia__pathovar_min_coverage >= 100.0:
        sys.exit('Error: --escherichia__pathovar_min_coverage must be between 50.0 and 100.0')

def check_external_programs():
    if not shutil.which('minimap2'):
        sys.exit('Error: could not find minimap2')
    return ['minimap2']


def data_dir():
    return pathlib.Path(__file__).parents[0] / 'data'


def args_get(args, name, default=None):
    return getattr(args, name, default) if hasattr(args, name) else default


def run_shigapass_for_single_assembly(assembly: str, args) -> str:
    """
    Run ShigaPass for a single assembly.
    """
    assembly_path = Path(assembly).expanduser().resolve()
    sample = assembly_path.stem

    if not assembly_path.exists():
        return '-'

    # Base directory
    base_dir = Path(args_get(args, 'shigapass_dir', data_dir())).expanduser().resolve()

    script_candidates = [
        base_dir / 'ShigaPass.sh',
        base_dir / 'SCRIPT' / 'ShigaPass.sh',
    ]
    shigapass_sh = next((p for p in script_candidates if p.exists()), script_candidates[0])

    db_dir_arg = args_get(args, 'shigapass_db_dir', None)
    if db_dir_arg:
        db_dir = Path(db_dir_arg).expanduser().resolve()
    else:
        db_candidates = [
            base_dir / 'ShigaPass_DataBases',
            base_dir / 'SCRIPT' / 'ShigaPass_DataBases',
            shigapass_sh.parent / 'ShigaPass_DataBases',
        ]
        db_dir = next((p for p in db_candidates if p.exists()), db_candidates[0])

    threads = int(args_get(args, 'shigapass_threads', 2))
    keep = bool(args_get(args, 'shigapass_keep', False))
    init_db = bool(args_get(args, 'shigapass_init_db', False))
    outdir_base = args_get(args, 'shigapass_outdir_base', 'ShigaPass_Results')

    if not shigapass_sh.exists():
        return '-'
    if not db_dir.exists():
        return '-'

    # Call -u only once per process
    use_u = False
    if init_db and not getattr(args, '_shigapass_init_done', False):
        use_u = True
        setattr(args, '_shigapass_init_done', True)

    with tempfile.TemporaryDirectory(prefix='shigapass_') as tmpd:
        tmpd = Path(tmpd)

        lst = tmpd / 'assemblies.txt'
        lst.write_text(str(assembly_path) + '\n', encoding='utf-8')

        # Output directory inside the temporary directory, per sample
        outdir = tmpd / f'{outdir_base}_{sample}'
        outdir.mkdir(parents=True, exist_ok=True)

        cmd = [
            str(shigapass_sh),
            '-l', str(lst),
            '-o', str(outdir),
            '-p', str(db_dir),
            '-t', str(threads),
        ]
        if keep:
            cmd.append('-k')
        if use_u:
            cmd.append('-u')

        # If script is not executable, run through bash
        if not os.access(shigapass_sh, os.X_OK):
            cmd = ['bash'] + cmd

        result = subprocess.run(cmd, capture_output=True, text=True)

        summary = outdir / 'ShigaPass_summary.csv'
        if not summary.exists():
            return '-'

        try:
            with summary.open('r', encoding='utf-8', newline='') as fh:
                reader = csv.reader(fh, delimiter=';')
                header = next(reader, [])
                pos = {name: i for i, name in enumerate(header)}
                if 'Name' not in pos or 'Predicted_Serotype' not in pos:
                    return '-'

                for row in reader:
                    if not row:
                        continue
                    if row[pos['Name']] == sample:
                        value = row[pos['Predicted_Serotype']].strip()
                        return value or '-'
        except Exception:
            return '-'


    return '-'


SEROTYPE_MAPPING = {
    'SB': 'Shigella boydii',
    'SD': 'Shigella dysenteriae',
    'SF': 'Shigella flexneri',
    'SS': 'Shigella sonnei'
}

def map_shigapass_serotype(serotype):
    """
    Maps ShigaPass serotype (e.g., 'SB2', 'SF1-5') to its full species name
    """
    if not serotype or serotype == '-':
        return '-'

    species_code = ''
    for char in serotype:
        if char.isalpha():
            species_code += char
        else:
            break
            
    if species_code in SEROTYPE_MAPPING:
        mapped_name = SEROTYPE_MAPPING[species_code]
        suffix = serotype[len(species_code):].strip()
        
        if suffix:
            return f"{mapped_name} {suffix}"
        return mapped_name
    
    return serotype



def get_results(assembly, minimap2_index, args, previous_results):
    full_headers, _ = get_headers()

    ref_file = data_dir() / 'virulence_ecoli.fsa'

    species = previous_results.get('enterobacterales__species__species', '').strip()
    pathotype_species = ['Escherichia coli / Shigella']

    if species not in pathotype_species:
        result_dict = {header: '-' for header in full_headers}
        result_dict['Pathotype'] = '-'
        return result_dict

    # Pathovar calling via minimap2
    pathovar, virulence_markers = minimap_pathovar(
        assembly,
        minimap2_index,
        ref_file,
        args.escherichia__pathovar_min_identity,
        args.escherichia__pathovar_min_coverage
    )

    predicted_serotype_raw = run_shigapass_for_single_assembly(assembly, args)
    predicted_serotype = map_shigapass_serotype(predicted_serotype_raw)

    # ShigaPass serotype — only run if any virulence markers starts with 'ipaH'
    # has_ipah = any(
    #     mk.startswith('ipaH') or
    #     any((h or '').startswith('ipaH') for h in (marker_hits or []))
    #     for mk, marker_hits in virulence_markers.items()
    # )
    
    # print(has_ipah)
    # if has_ipah:
    #     predicted_serotype_raw = run_shigapass_for_single_assembly(assembly, args)
    #     predicted_serotype = map_shigapass_serotype(predicted_serotype_raw)
    # else:
    #     predicted_serotype_raw = '-'
    #     predicted_serotype = '-'

    # markers
    result_dict = {header: '-' for header in full_headers}
    for marker, marker_hits in virulence_markers.items():
        if marker in result_dict:
            result_dict[marker] = ";".join(marker_hits) if marker_hits else '-'

    if predicted_serotype != '-':
        result_dict['Pathotype'] = predicted_serotype
    elif pathovar != '-':
        result_dict['Pathotype'] = pathovar
    else:
        result_dict['Pathotype'] = '-'

    return result_dict


    # if pathovar != '-':
    #     result_dict['Pathotype'] = pathovar
    
    # elif predicted_serotype != '-':
    #     result_dict['Pathotype'] = predicted_serotype
    # else:
    #     result_dict['Pathotype'] = '-'

    # return result_dict

# def get_results(assembly, minimap2_index, args, previous_results) -> Dict:
#     full_headers, _ = get_headers()

#     ref_file = data_dir() / 'virulence_ecoli.fsa'

#     species = previous_results.get('enterobacterales__species__species', '').strip()
#     pathotype_species = ['Escherichia coli / Shigella']

#     if species not in pathotype_species:
#         result_dict = {header: '-' for header in full_headers}
#         result_dict['Pathotype'] = '-'
#         return result_dict

#     # Pathovar calling via minimap2
#     pathovar, virulence_markers = minimap_pathovar(
#         assembly,
#         minimap2_index,
#         ref_file,
#         args.escherichia__pathovar_min_identity,
#         args.escherichia__pathovar_min_coverage
#     )

#     # ShigaPass serotype — only run if any virulence marker (key or hit) starts with 'ipaH'
#     has_ipah = any(
#         (mk.startswith('ipaH') and marker_hits) or
#         any(h.startswith('ipaH') for h in (marker_hits or []))
#         for mk, marker_hits in (virulence_markers or {}).items()
#     )
#     # ShigaPass serotype
#     predicted_serotype = run_shigapass_for_single_assembly(assembly, args)

#     # Build result row
#     result_dict = {header: '-' for header in full_headers}
#     for marker, marker_hits in virulence_markers.items():
#         if marker in result_dict:
#             result_dict[marker] = ";".join(marker_hits) if marker_hits else '-'

#     result_dict['Predicted_Serotype'] = predicted_serotype if predicted_serotype else '-'

#     if pathovar != '-':
#         result_dict['Pathotype'] = pathovar
#     elif predicted_serotype != '-':
#         result_dict['Pathotype'] = predicted_serotype
#     else:
#         result_dict['Pathotype'] = '-'

#     return result_dict






