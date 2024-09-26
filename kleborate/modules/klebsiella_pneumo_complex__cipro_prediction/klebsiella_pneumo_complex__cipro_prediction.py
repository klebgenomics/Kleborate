
"""
Copyright 2024 Kat Holt
Copyright 2024 Ryan Wick (rrwick@gmail.com)
Copyright 2024 (gathonimaranga@gmail.com)
https://github.com/katholt/Kleborate/


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

from ...shared.resMinimap import read_class_file, get_res_headers

def description():
    return 'Ciprofloxacin resistance prediction based on the' \
           'results of the klebsiella_pneumo_complex__amr module'


def prerequisite_modules():
    return ['klebsiella_pneumo_complex__amr']


def get_headers():
    full_headers = ['cipro_prediction', 'cipro_prediction_support', 'cipro_prediction_group']
    stdout_headers = []
    return full_headers, stdout_headers


def add_cli_options(parser):
    pass


def check_cli_options(args):
    pass


def check_external_programs():
    return []

def data_dir():
    return pathlib.Path(__file__).parents[1] / 'klebsiella_pneumo_complex__amr' / 'data'



def get_results(assembly, minimap2_index, args, previous_results):
    """
    Uses the Flq_mutations, Flq_acquired, and AGly_acquired columns to predict ciprofloxacin resistance
    """

    flq_mutations = previous_results['klebsiella_pneumo_complex__amr__Flq_mutations']
    flq_acquired = previous_results['klebsiella_pneumo_complex__amr__Flq_acquired']
    agly_acquired = previous_results['klebsiella_pneumo_complex__amr__AGly_acquired']

    #keep only aac(6')-Ib-cr in agly_acquired
    if "aac(6')-Ib-cr" in agly_acquired:
        agly_acquired = "aac(6')-Ib-cr"
    else:
        agly_acquired = '-'

    #0^ QRDR, 0 PMQR, 0 aac6 (^QRDR excludes GyrA-87G and GyrA-87H)
    if 'GyrA-87G' in flq_mutations or 'GyrA-87H' in flq_mutations:
        #test genome: KP_NORM_BLD_112904 (GyrA-87H, no PMQR, no aac6)
        if flq_mutations != "-" and flq_mutations.count(";") == 0 and flq_acquired == "-" and agly_acquired == "-":
            cipro_prediction = "WT S"
            cipro_prediction_support= "2.75% (N=156/5674)"
            cipro_prediction_group= "0^ QRDR, 0 PMQR, 0 aac6"
    else: 
        #test genome: KP_NORM_URN_105939 (no QRDR, no PMQR, no aac6)
        if flq_mutations == "-" and flq_acquired == "-" and agly_acquired == "-":
            cipro_prediction = "WT S"
            cipro_prediction_support= "2.75% (N=156/5674)"
            cipro_prediction_group= "0^ QRDR, 0 PMQR, 0 aac6"

    #0 QRDR, 0 PMQR, 1 aac6
        #test genome: ERR4635459
    if flq_mutations == "-" and flq_acquired == "-" and agly_acquired != "-":
        cipro_prediction = "WT S"
        cipro_prediction_support= "18.12% (N=29/160)"
        cipro_prediction_group= "0 QRDR, 0 PMQR, 1 aac6"      

    #0 QRDR, qnrB1, 0 aac6
        #test genome: SRR5973253
    if "qnrB1." in flq_acquired:
        if flq_mutations == "-" and agly_acquired == "-":
            cipro_prediction = "NWT I"
            cipro_prediction_support= "43.12% (N=69/160)"
            cipro_prediction_group= "0 QRDR, qnrB1, 0 aac6" 

    #1 QRDR, 0 PMQR, 0 aac6
        #test genome: NK_H5_026
    if 'GyrA-87G' not in flq_mutations and 'GyrA-87H' not in flq_mutations:
        if flq_mutations != "-" and flq_mutations.count(";") == 0 and flq_acquired == "-" and agly_acquired == "-":
                cipro_prediction = "NWT R"
                cipro_prediction_support= "77.67% (N=80/103)"
                cipro_prediction_group= "1 QRDR, 0 PMQR, 0 aac6" 

    #1 QRDR, 0 PMQR, 1 aac6
        #test genome: ERR4046108_NHP1975
    if flq_mutations != "-" and flq_mutations.count(";") == 0 and flq_acquired == "-" and agly_acquired != '-':
            cipro_prediction = "NWT R"
            cipro_prediction_support= "86.96% (N=20/23)"
            cipro_prediction_group= "1 QRDR, 0 PMQR, 1 aac6" 
    
    # >1 QRDR, 0 PMQR, * aac6 (*aac6 presence is not considered)
        #test genome: G20250619 (with aac6), G20250926 (no aac6)
    if flq_mutations.count(";") > 0 and flq_acquired == "-":
            cipro_prediction = "NWT R"
            cipro_prediction_support= "98.94% (N=2145/2168)"
            cipro_prediction_group= ">1 QRDR, 0 PMQR, * aac6" 
    
    
    # 0 QRDR, 1^ PMQR, 0 aac6 (^PMQR excludes qnrB1)
        #test genome: 56CM1
    if flq_mutations == "-" and flq_acquired != '-' and "qnrB1." not in flq_acquired and flq_acquired.count(";") == 0 and agly_acquired == "-":
            cipro_prediction = "NWT R"
            cipro_prediction_support= "77.61% (N=423/545)"
            cipro_prediction_group= "0 QRDR, 1^ PMQR, 0 aac6" 

    # 0 QRDR, 1 PMQR, 1 aac6 
        #test genome: HE205
    if flq_mutations == "-" and flq_acquired != '-' and flq_acquired.count(";") == 0 and agly_acquired != "-":
            cipro_prediction = "NWT R"
            cipro_prediction_support= "93.16% (N=763/819)"
            cipro_prediction_group= "0 QRDR, 1 PMQR, 1 aac6" 

    # 0 QRDR, >1 PMQR, * aac6 (*aac6 presence is not considered)
        #test genome: ERR3480591 (with aac6), ERR4635120 (no aac6)
    if flq_mutations == "-" and flq_acquired.count(";") > 0:
            cipro_prediction = "NWT R"
            cipro_prediction_support= "99.14% (N=2426/2447)"
            cipro_prediction_group= "0 QRDR, >1 PMQR, * aac6" 

    # >0 QRDR, >0 PMQR, * aac6 (*aac6 presence is not considered)
        #test genome: ERR3567346_NHP139 (with aac6), ERR3585150_BMP790 (no aac6)
    if flq_mutations != "-" and flq_acquired.count(";") >= 0 and flq_acquired != "-" and flq_acquired.count(";") >= 0:
            cipro_prediction = "NWT R"
            cipro_prediction_support= "95.95% (N=5923/6173)"
            cipro_prediction_group= ">0 QRDR, >0 PMQR, * aac6" 

            
    return {
        'cipro_prediction': cipro_prediction,
        'cipro_prediction_support': cipro_prediction_support,
        'cipro_prediction_group': cipro_prediction_group,
    }


