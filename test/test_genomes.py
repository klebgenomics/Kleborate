"""
Copyright 2018 Kat Holt
Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/katholt/Kleborate/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <http://www.gnu.org/licenses/>.
"""

import collections
import tempfile
import unittest
from kleborate.kleborate import get_data_path, get_output_headers, gunzip_contigs_if_necessary, \
    get_strain_name, get_contig_stat_results, get_species_results, is_kp_complex, \
    get_chromosome_mlst_results, get_ybt_mlst_results, get_clb_mlst_results, get_iuc_mlst_results, \
    get_iro_mlst_results, get_hypermucoidy_results, get_wzi_and_k_locus_results, \
    get_resistance_results, get_summary_results


class TestGenomes(unittest.TestCase):

    def setUp(self):
        self.data_dir = get_data_path()
        Args = collections.namedtuple('Args', ['resistance', 'kaptive_k', 'kaptive_o'])
        self.args = Args(resistance=True, kaptive_k=False, kaptive_o=False)
        _, _, self.res_headers = get_output_headers(self.args, self.data_dir)

    def get_all_results(self, contigs):
        results = {'strain': get_strain_name(contigs)}
        results.update(get_contig_stat_results(contigs))
        results.update(get_species_results(contigs, self.data_dir))
        kp_complex = is_kp_complex(results)
        results.update(get_chromosome_mlst_results(self.data_dir, contigs))
        results.update(get_ybt_mlst_results(self.data_dir, contigs))
        results.update(get_clb_mlst_results(self.data_dir, contigs))
        results.update(get_iuc_mlst_results(self.data_dir, contigs))
        results.update(get_iro_mlst_results(self.data_dir, contigs))
        results.update(get_hypermucoidy_results(self.data_dir, contigs))
        results.update(get_wzi_and_k_locus_results(self.data_dir, contigs))
        results.update(get_resistance_results(self.data_dir, contigs, self.args, self.res_headers,
                                              kp_complex))
        results.update(get_summary_results(results, self.res_headers))
        return results

    def test_GCF_900501255(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_900501255.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '0')
            self.assertEqual(results['rmpA'], 'rmpA_2(KpVP-1)')
            self.assertEqual(results['rmpA2'], 'rmpA2_4*-50%')
            self.assertEqual(results['AGly'], '-')
            self.assertEqual(results['Col'], '-')
            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], '-')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], '-')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], '-')
            self.assertEqual(results['Rif'], '-')
            self.assertEqual(results['Sul'], '-')
            self.assertEqual(results['Tet'], '-')
            self.assertEqual(results['Tmt'], '-')
            self.assertEqual(results['Omp'], '-')
            self.assertEqual(results['Bla'], 'AmpH*;SHV-187*')
            self.assertEqual(results['Bla_Carb'], '-')
            self.assertEqual(results['Bla_ESBL'], '-')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], '-')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_GCF_003400925(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_003400925.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '6')
            self.assertEqual(results['rmpA'], '-')
            self.assertEqual(results['rmpA2'], 'rmpA2_9')
            self.assertEqual(results['AGly'], 'RmtC')
            self.assertEqual(results['Col'], '-')
            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], 'GyrA-83I;ParC-80I')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], '-')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], '-')
            self.assertEqual(results['Rif'], '-')
            self.assertEqual(results['Sul'], 'SulI')
            self.assertEqual(results['Tet'], 'TetA')
            self.assertEqual(results['Tmt'], 'DfrA1')
            self.assertEqual(results['Omp'], '-')
            self.assertEqual(results['Bla'], 'AmpH*')
            self.assertEqual(results['Bla_Carb'], 'NDM-1')
            self.assertEqual(results['Bla_ESBL'], '-')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], 'SHV-11')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_GCF_001068035(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_001068035.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '0')
            self.assertEqual(results['rmpA'], 'rmpA_2(KpVP-1)')
            self.assertEqual(results['rmpA2'], 'rmpA2_3*-0%')
            self.assertEqual(results['AGly'], '-')

            # PmrB is split in two: 40% (start) and 60% (end). The correct answer is 40% as that
            # includes the beginning of the gene.
            self.assertEqual(results['Col'], 'PmrB-40%')

            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], '-')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], '-')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], '-')
            self.assertEqual(results['Rif'], '-')
            self.assertEqual(results['Sul'], '-')
            self.assertEqual(results['Tet'], '-')
            self.assertEqual(results['Tmt'], '-')
            self.assertEqual(results['Omp'], '-')
            self.assertEqual(results['Bla'], 'AmpH*')
            self.assertEqual(results['Bla_Carb'], '-')
            self.assertEqual(results['Bla_ESBL'], '-')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], '-')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_GCF_003095495(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_003095495.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '17')
            self.assertEqual(results['rmpA'], '-')
            self.assertEqual(results['rmpA2'], '-')
            self.assertEqual(results['AGly'], 'Aac3-IId^;AadA2^;Aph3-Ia^;RmtB;Sat-2A;StrA^;StrB')

            # PmrB is split in two: 36% (start) and 66% (end). The correct answer is 36% as that
            # includes the beginning of the gene.
            self.assertEqual(results['Col'], 'MgrB-62%;PmrB-36%')

            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], 'GyrA-83I;ParC-80I')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], 'Erm42*;MphA')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], 'CatA1^')
            self.assertEqual(results['Rif'], '-')
            self.assertEqual(results['Sul'], 'SulI;SulII')
            self.assertEqual(results['Tet'], 'TetG')
            self.assertEqual(results['Tmt'], 'DfrA12?')
            self.assertEqual(results['Omp'], 'OmpK35-25%;OmpK36GD')
            self.assertEqual(results['Bla'], 'AmpH^;TEM-1D^')
            self.assertEqual(results['Bla_Carb'], 'KPC-2')
            self.assertEqual(results['Bla_ESBL'], 'CTX-M-14')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], 'SHV-11')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_GCF_000009885(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_000009885.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '1')
            self.assertEqual(results['rmpA'], 'rmpA_11(ICEKp1),rmpA_2(KpVP-1)')
            self.assertEqual(results['rmpA2'], 'rmpA2_3-47%')
            self.assertEqual(results['AGly'], '-')
            self.assertEqual(results['Col'], '-')
            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], '-')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], '-')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], '-')
            self.assertEqual(results['Rif'], '-')
            self.assertEqual(results['Sul'], '-')
            self.assertEqual(results['Tet'], '-')
            self.assertEqual(results['Tmt'], '-')
            self.assertEqual(results['Omp'], '-')
            self.assertEqual(results['Bla'], 'AmpH')
            self.assertEqual(results['Bla_Carb'], '-')
            self.assertEqual(results['Bla_ESBL'], '-')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], 'SHV-11^')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_GCF_000968155(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_000968155.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '0')
            self.assertEqual(results['rmpA'], 'rmpA_9(KpVP-2)')
            self.assertEqual(results['rmpA2'], '-')
            self.assertEqual(results['AGly'], '-')
            self.assertEqual(results['Col'], '-')
            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], '-')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], '-')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], '-')
            self.assertEqual(results['Rif'], '-')
            self.assertEqual(results['Sul'], '-')
            self.assertEqual(results['Tet'], '-')
            self.assertEqual(results['Tmt'], '-')
            self.assertEqual(results['Omp'], '-')
            self.assertEqual(results['Bla'], 'AmpH*')
            self.assertEqual(results['Bla_Carb'], '-')
            self.assertEqual(results['Bla_ESBL'], '-')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], '-')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_GCF_002108345(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_002108345.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '9')
            self.assertEqual(results['rmpA'], '-')
            self.assertEqual(results['rmpA2'], '-')
            self.assertEqual(results['AGly'], 'Aac3-IId^;Aph3-Ia^')
            self.assertEqual(results['Col'], '-')
            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], 'GyrA-83I;ParC-80I;QnrB4;QnrB4')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], '-')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], 'CatA2*;CatB3*')
            self.assertEqual(results['Rif'], 'Arr3')
            self.assertEqual(results['Sul'], 'SulI*;SulI^')
            self.assertEqual(results['Tet'], '-')
            self.assertEqual(results['Tmt'], '-')
            self.assertEqual(results['Omp'], 'OmpK35-18%')
            self.assertEqual(results['Bla'], 'AmpH*;DHA-1;OXA-1;SHV-187*')
            self.assertEqual(results['Bla_Carb'], '-')
            self.assertEqual(results['Bla_ESBL'], '-')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], '-')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_GCF_002247645(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_002247645.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '5')
            self.assertEqual(results['rmpA'], '-')
            self.assertEqual(results['rmpA2'], '-')
            self.assertEqual(results['AGly'], 'Aac3-IId^')
            self.assertEqual(results['Col'], '-')
            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], 'Qnr-S1')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], '-')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], 'FloR*')
            self.assertEqual(results['Rif'], '-')
            self.assertEqual(results['Sul'], '-')
            self.assertEqual(results['Tet'], 'TetA*;TetX*?')
            self.assertEqual(results['Tmt'], '-')
            self.assertEqual(results['Omp'], '-')
            self.assertEqual(results['Bla'], 'AmpH*;LAP-2;SHV-2a^;TEM-1D^')
            self.assertEqual(results['Bla_Carb'], '-')
            self.assertEqual(results['Bla_ESBL'], '-')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], '-')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_GCF_002248955(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_002248955.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '4')
            self.assertEqual(results['rmpA'], '-')
            self.assertEqual(results['rmpA2'], '-')
            self.assertEqual(results['AGly'], 'Aac3-IId^')
            self.assertEqual(results['Col'], 'Mcr3-1*')
            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], 'GyrA-83F;GyrA-87A;ParC-80I')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], '-')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], 'CatA1^')
            self.assertEqual(results['Rif'], '-')
            self.assertEqual(results['Sul'], '-')
            self.assertEqual(results['Tet'], 'TetA')
            self.assertEqual(results['Tmt'], '-')
            self.assertEqual(results['Omp'], '-')
            self.assertEqual(results['Bla'], 'AmpH*;SHV-28^')
            self.assertEqual(results['Bla_Carb'], '-')
            self.assertEqual(results['Bla_ESBL'], '-')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], '-')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_GCF_003345475(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/GCF_003345475.1.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)
            self.assertEqual(results['num_resistance_genes'], '5')
            self.assertEqual(results['rmpA'], '-')
            self.assertEqual(results['rmpA2'], '-')
            self.assertEqual(results['AGly'], '-')
            self.assertEqual(results['Col'], 'Mcr1-1;MgrB-49%')
            self.assertEqual(results['Fcyn'], '-')
            self.assertEqual(results['Flq'], 'Qnr-S1')
            self.assertEqual(results['Gly'], '-')
            self.assertEqual(results['MLS'], '-')
            self.assertEqual(results['Ntmdz'], '-')
            self.assertEqual(results['Phe'], '-')
            self.assertEqual(results['Rif'], '-')
            self.assertEqual(results['Sul'], 'SulI')
            self.assertEqual(results['Tet'], 'TetA')
            self.assertEqual(results['Tmt'], 'DfrA1')
            self.assertEqual(results['Omp'], '-')
            self.assertEqual(results['Bla'], 'AmpH*;SHV-110')
            self.assertEqual(results['Bla_Carb'], '-')
            self.assertEqual(results['Bla_ESBL'], '-')
            self.assertEqual(results['Bla_ESBL_inhR'], '-')
            self.assertEqual(results['Bla_broad'], '-')
            self.assertEqual(results['Bla_broad_inhR'], '-')

    def test_AS7(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            contigs = gunzip_contigs_if_necessary('test/sequences/AS7.fna.gz', tmp_dir)
            results = self.get_all_results(contigs)

            # This genome is missing the start of rmpA, which counts as a truncation to 0%.
            self.assertEqual(results['rmpA'], 'rmpA_1*-0%(KpVP-1)')
