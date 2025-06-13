
"""
Copyright 2025, Kara Tsang, Kat Holt
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
    full_headers = ['Ciprofloxacin_prediction', 'Ciprofloxacin_profile_support', 'Ciprofloxacin_profile', 'Ciprofloxacin_MIC_prediction']
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

    #0^ QRDR, 0 PMQR, 0 aac(6`)-Ib-cr (^QRDR excludes GyrA-87G and GyrA-87H)
    if 'GyrA-87G' in flq_mutations or 'GyrA-87H' in flq_mutations:
        #test genome: KP_NORM_BLD_112904 (GyrA-87H, no PMQR, no aac(6`)-Ib-cr)
        if flq_mutations != "-" and flq_mutations.count(";") == 0 and flq_acquired == "-" and agly_acquired == "-":
            Ciprofloxacin_prediction = "wildtype S"
            Ciprofloxacin_profile_support= "90.99% S (N=5168/5680)"
            Ciprofloxacin_profile= "0^ QRDR, 0 PMQR, 0 aac(6`)-Ib-cr"
            Ciprofloxacin_MIC_prediction="0.25 mg/L [0.25-0.25]"
    else: 
        #test genome: KP_NORM_URN_105939 (no QRDR, no PMQR, no aac(6`)-Ib-cr)
        if flq_mutations == "-" and flq_acquired == "-" and agly_acquired == "-":
            Ciprofloxacin_prediction = "wildtype S"
            Ciprofloxacin_profile_support= "90.99% S (N=5168/5680)"
            Ciprofloxacin_profile= "0^ QRDR, 0 PMQR, 0 aac(6`)-Ib-cr"
            Ciprofloxacin_MIC_prediction="0.25 mg/L [0.25-0.25]"

    #0 QRDR, 0 PMQR, 1 aac(6`)-Ib-cr
        #test genome: ERR4635459
    if flq_mutations == "-" and flq_acquired == "-" and agly_acquired != "-":
        Ciprofloxacin_prediction = "wildtype S"
        Ciprofloxacin_profile_support= "65.22% S (N=105/161)"
        Ciprofloxacin_profile= "0 QRDR, 0 PMQR, 1 aac(6`)-Ib-cr" 
        Ciprofloxacin_MIC_prediction="0.25 mg/L [0.25-0.5]"  

    #0 QRDR, qnrB1, 0 aac(6`)-Ib-cr
        #test genome: SRR5973253
    if "qnrB1." in flq_acquired:
        if flq_mutations == "-" and agly_acquired == "-":
            Ciprofloxacin_prediction = "nonwildtype I"
            Ciprofloxacin_profile_support= "81.25% I/R (n=130/160)"
            Ciprofloxacin_profile= "0 QRDR, qnrB1, 0 aac(6`)-Ib-cr"
            Ciprofloxacin_MIC_prediction="0.5 mg/L [0.5-1]"

    #1 QRDR, 0 PMQR, 0 aac(6`)-Ib-cr
        #test genome: NK_H5_026
    if 'GyrA-87G' not in flq_mutations and 'GyrA-87H' not in flq_mutations:
        if flq_mutations != "-" and flq_mutations.count(";") == 0 and flq_acquired == "-" and agly_acquired == "-":
                Ciprofloxacin_prediction = "nonwildtype R"
                Ciprofloxacin_profile_support= "77.67% R (N=80/103)"
                Ciprofloxacin_profile= "1 QRDR, 0 PMQR, 0 aac(6`)-Ib-cr"
                Ciprofloxacin_MIC_prediction="1 mg/L [1-2]"

    #1 QRDR, 0 PMQR, 1 aac(6`)-Ib-cr
        #test genome: ERR4046108_NHP1975
    if flq_mutations != "-" and flq_mutations.count(";") == 0 and flq_acquired == "-" and agly_acquired != '-':
            Ciprofloxacin_prediction = "nonwildtype R"
            Ciprofloxacin_profile_support= "86.96% R (N=20/23)"
            Ciprofloxacin_profile= "1 QRDR, 0 PMQR, 1 aac(6`)-Ib-cr"
            Ciprofloxacin_MIC_prediction="2 mg/L [1-2]"
    
    # >1 QRDR, 0 PMQR, * aac(6`)-Ib-cr (*aac(6`)-Ib-cr presence is not considered)
        #test genome: G20250619 (with aac(6`)-Ib-cr), G20250926 (no aac(6`)-Ib-cr)
    if flq_mutations.count(";") > 0 and flq_acquired == "-":
            Ciprofloxacin_prediction = "nonwildtype R"
            Ciprofloxacin_profile_support= "99.22% R (N=2150/2167)"
            Ciprofloxacin_profile= ">1 QRDR, 0 PMQR, * aac(6`)-Ib-cr" 
            Ciprofloxacin_MIC_prediction="2 mg/L [2-4]"
    
    
    # 0 QRDR, 1^ PMQR, 0 aac(6`)-Ib-cr (^PMQR excludes qnrB1)
        #test genome: 56CM1
    if flq_mutations == "-" and flq_acquired != '-' and "qnrB1." not in flq_acquired and flq_acquired.count(";") == 0 and agly_acquired == "-":
            Ciprofloxacin_prediction = "nonwildtype R"
            Ciprofloxacin_profile_support= "77.47% R (N=423/546)"
            Ciprofloxacin_profile= "0 QRDR, 1^ PMQR, 0 aac(6`)-Ib-cr"
            Ciprofloxacin_MIC_prediction="1 mg/L [1-2]"

    # 0 QRDR, 1 PMQR, 1 aac(6`)-Ib-cr 
        #test genome: HE205
    if flq_mutations == "-" and flq_acquired != '-' and flq_acquired.count(";") == 0 and agly_acquired != "-":
            Ciprofloxacin_prediction = "nonwildtype R"
            Ciprofloxacin_profile_support= "94.63% R (N=775/819)"
            Ciprofloxacin_profile= "0 QRDR, 1 PMQR, 1 aac(6`)-Ib-cr"
            Ciprofloxacin_MIC_prediction="2 mg/L [1-2]"

    # 0 QRDR, >1 PMQR, * aac(6`)-Ib-cr (*aac(6`)-Ib-cr presence is not considered)
        #test genome: ERR3480591 (with aac(6`)-Ib-cr), ERR4635120 (no aac(6`)-Ib-cr)
    if flq_mutations == "-" and flq_acquired.count(";") > 0:
            Ciprofloxacin_prediction = "nonwildtype R"
            Ciprofloxacin_profile_support= "97.06% R (N=66/68)"
            Ciprofloxacin_profile= "0 QRDR, >1 PMQR, * aac(6`)-Ib-cr"
            Ciprofloxacin_MIC_prediction="2 mg/L [2-4]"

    # >0 QRDR, >0 PMQR, * aac(6`)-Ib-cr (*aac(6`)-Ib-cr presence is not considered)
        #test genome: ERR3567346_NHP139 (with aac(6`)-Ib-cr), ERR3585150_BMP790 (no aac(6`)-Ib-cr)
    if flq_mutations != "-" and flq_acquired.count(";") >= 0 and flq_acquired != "-" and flq_acquired.count(";") >= 0:
            Ciprofloxacin_prediction = "nonwildtype R"
            Ciprofloxacin_profile_support= "99.22% R (N=2421/2440)"
            Ciprofloxacin_profile= ">0 QRDR, >0 PMQR, * aac(6`)-Ib-cr" 
            Ciprofloxacin_MIC_prediction="4 mg/L [4-4]"

            
    return {
        'Ciprofloxacin_prediction': Ciprofloxacin_prediction,
        'Ciprofloxacin_profile_support': Ciprofloxacin_profile_support,
        'Ciprofloxacin_profile': Ciprofloxacin_profile,
        'Ciprofloxacin_MIC_prediction': Ciprofloxacin_MIC_prediction
    }


