#!/usr/bin/env python3
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

import os
import sys
import subprocess
import distutils.spawn
from setuptools import setup
from setuptools.command.install import install


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = '0.0.0'
exec(open('kleborate/version.py').read())


def check_dir_write_permission(directory):
    if os.path.isdir(directory) and not os.access(directory, os.W_OK):
        sys.exit('Error: no write permission for ' + directory + '  ' +
                 'Perhaps you need to use sudo?')


def build_blast_db(data_dir, fasta_filename, seq_type):
    fasta_path = os.path.join(data_dir, fasta_filename)
    makeblastdb_cmd = ['makeblastdb', '-dbtype', seq_type, '-in', fasta_path]
    print('  ' + ' '.join(makeblastdb_cmd))
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(makeblastdb_cmd, stdout=devnull)


class KleborateInstall(install):
    """
    This is a subclass of install where I've added a couple things to the installation process:
      * building BLAST databases with makeblastdb
      * printing example commands
    """

    def run(self):
        check_dir_write_permission(self.install_lib)
        check_dir_write_permission(self.install_scripts)

        # Call the base class run to do the bulk of the installation.
        install.run(self)

        # Build the BLAST databases now, as the user may not have write permissions to the
        # installation directory when running Kleborate.
        print('')
        data_dir = os.path.join(self.install_lib, 'kleborate', 'data')

        if not distutils.spawn.find_executable('makeblastdb'):
            print('Warning: could not find makeblastdb, so BLAST databases were not built.')
            print('They can be built when you run Kleborate, but this will require write '
                  'permissions to the data directory:')
            print(data_dir)
        else:
            check_dir_write_permission(data_dir)
            print('Building BLAST databases with makeblastdb:')
            try:
                for fasta in ['ARGannot_r3.fasta', 'clb_alleles.fasta', 'hypermucoidy.fasta',
                              'iro_alleles.fasta', 'iuc_alleles.fasta',
                              'Klebsiella_pneumoniae.fasta', 'wzi.fasta', 'ybt_alleles.fasta',
                              'MgrB_and_PmrB.fasta', 'OmpK.fasta']:
                    build_blast_db(data_dir, fasta, 'nucl')
                for fasta in ['QRDR_120.aa']:
                    build_blast_db(data_dir, fasta, 'prot')
            except subprocess.CalledProcessError:
                print('\n')
                print('Warning: makeblastdb failed, so BLAST databases were not built.')
                print('They can be built when you run Kleborate, but this will require write '
                      'permissions to the data directory:')
                print(data_dir)

        ascii_art = " _   __ _     ______ ____   ____  _____       _______ ______ \n" \
                    "| | / /| |   |  ____|  _ \ / __ \|  __ \   /\|__   __|  ____|\n" \
                    "| |/ / | |   | |__  | |_) | |  | | |__) | /  \  | |  | |__   \n" \
                    "|   <  | |   |  __| |  _ <| |  | |  _  / / /\ \ | |  |  __|  \n" \
                    "| |\ \ | |___| |____| |_) | |__| | | \ \/ ____ \| |  | |____ \n" \
                    "|_| \_\|_____|______|____/ \____/|_|  \_\/    \_\_|  |______|\n"

        print('\n')
        print('\033[1m' + ascii_art + '\033[0m')  # bold formatting

        print('Kleborate is installed and ready to use!')
        print('')

        print('Example commands:')
        print('  kleborate --help')
        print('  kleborate -o results.txt -a *.fasta')
        print('  kleborate --resistance -o results.txt -a *.fasta')
        print('  kleborate --all -o results.txt -a *.fasta')
        print('')


setup(name='Kleborate',
      version=__version__,
      description='Kleborate',
      long_description=readme(),
      classifiers=['Development Status :: 4 - Beta',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Programming Language :: Python :: 2',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Medical Science Apps.',
                   'Intended Audience :: Science/Research'],
      keywords='microbial genomics sequence typing',
      url='https://github.com/katholt/Kleborate',
      author='Kathryn Holt',
      author_email='',
      packages=['kleborate', 'kaptive'],
      entry_points={'console_scripts': ['kleborate = kleborate.kleborate:main']},
      include_package_data=True,
      zip_safe=False,
      cmdclass={'install': KleborateInstall})
