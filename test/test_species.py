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

import unittest
from kleborate.kleborate import get_data_path
from kleborate.species import get_klebsiella_species, is_kp_complex


class TestSpecies(unittest.TestCase):

    def setUp(self):
        self.data_dir = get_data_path()

    def test_klebsiella_aerogenes(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000215745.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella aerogenes')

    def test_klebsiella_grimontii(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000733495.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella grimontii')

    def test_klebsiella_indica(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_005860775.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella indica')

    def test_klebsiella_michiganensis(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000240325.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella michiganensis')

    def test_klebsiella_oxytoca(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000247855.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella oxytoca')

    def test_klebsiella_pasteurii(self):
        species, _ = get_klebsiella_species('test/sequences/GCA_902158585.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella pasteurii')

    def test_klebsiella_pneumoniae(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000016305.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella pneumoniae')

    def test_klebsiella_quasipneumoniae_subsp_quasipneumoniae(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000492415.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella quasipneumoniae subsp. quasipneumoniae')

    def test_klebsiella_quasipneumoniae_subsp_similipneumoniae(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000492795.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella quasipneumoniae subsp. similipneumoniae')

    def test_klebsiella_quasivariicola(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000523395.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella quasivariicola')

    def test_klebsiella_spallanzanii(self):
        species, _ = get_klebsiella_species('test/sequences/GCA_901563875.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella spallanzanii')

    def test_klebsiella_variicola_subsp_variicola(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000019565.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella variicola subsp. variicola')

    def test_klebsiella_variicola_subsp_tropicalensis(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_002806645.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella variicola subsp. tropicana')

    def test_klebsiella_africanensis(self):
        species, _ = get_klebsiella_species('test/sequences/ERR2835900.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella africana')

    def test_raoultella_planticola(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000648315.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella (Raoultella) planticola')

    def test_raoultella_ornithinolytica(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000247895.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella (Raoultella) ornithinolytica')

    def test_raoultella_terrigena(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_000829965.1.fna.gz', self.data_dir)
        self.assertEqual(species, 'Klebsiella (Raoultella) terrigena')

    def test_salmonella(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_004010735.1.fna.gz', self.data_dir)
        self.assertTrue('Salmonella' in species)

    def test_citrobacter(self):
        species, _ = get_klebsiella_species('test/sequences/GCF_003937345.1.fna.gz', self.data_dir)
        self.assertTrue('Citrobacter' in species)


class TestKpComplex(unittest.TestCase):

    def test_klebsiella_aerogenes(self):
        self.assertFalse(is_kp_complex({'species': 'Klebsiella aerogenes'}))

    def test_klebsiella_grimontii(self):
        self.assertFalse(is_kp_complex({'species': 'Klebsiella grimontii'}))

    def test_klebsiella_michiganensis(self):
        self.assertFalse(is_kp_complex({'species': 'Klebsiella michiganensis'}))

    def test_klebsiella_oxytoca(self):
        self.assertFalse(is_kp_complex({'species': 'Klebsiella oxytoca'}))

    def test_klebsiella_pneumoniae(self):
        self.assertTrue(is_kp_complex({'species': 'Klebsiella pneumoniae'}))

    def test_klebsiella_quasipneumoniae_subsp_quasipneumoniae(self):
        self.assertTrue(is_kp_complex({'species': 'Klebsiella quasipneumoniae subsp. '
                                                  'quasipneumoniae'}))

    def test_klebsiella_quasipneumoniae_subsp_similipneumoniae(self):
        self.assertTrue(is_kp_complex({'species': 'Klebsiella quasipneumoniae subsp. '
                                                  'similipneumoniae'}))

    def test_klebsiella_quasivariicola(self):
        self.assertTrue(is_kp_complex({'species': 'Klebsiella quasivariicola'}))

    def test_klebsiella_variicola_subsp_variicola(self):
        self.assertTrue(is_kp_complex({'species': 'Klebsiella variicola subsp. variicola'}))

    def test_klebsiella_variicola_subsp_tropicalensis(self):
        self.assertTrue(is_kp_complex({'species': 'Klebsiella variicola subsp. tropicalensis'}))

    def test_klebsiella_africanensis(self):
        self.assertTrue(is_kp_complex({'species': 'Klebsiella africanensis'}))

    def test_raoultella_planticola(self):
        self.assertFalse(is_kp_complex({'species': 'Raoultella planticola'}))

    def test_raoultella_ornithinolytica(self):
        self.assertFalse(is_kp_complex({'species': 'Raoultella ornithinolytica'}))

    def test_raoultella_terrigena(self):
        self.assertFalse(is_kp_complex({'species': 'Raoultella terrigena'}))

    def test_salmonella(self):
        self.assertFalse(is_kp_complex({'species': 'Salmonella enterica'}))

    def test_citrobacter(self):
        self.assertFalse(is_kp_complex({'species': 'Citrobacter freundii'}))
