# run chromosome, yersiniabactin and colibactin MLST on a Klebs genome
# optionally, run resistance gene screening
import os, sys
import gzip
import argparse
import distutils.spawn
from pkg_resources import resource_filename
import imp


def main():
    args = parse_arguments()
    check_inputs_and_programs(args)

    # find necessary resources
    data_folder = resource_filename(__name__, 'data')
    mlstblast = resource_filename(__name__, 'mlstBLAST.py')
    resblast = resource_filename(__name__, 'resBLAST.py')
    clusterblast = resource_filename(__name__, 'clusterBLAST.py')

    header_string = '\t'.join(['strain', 'ST', 'Yersiniabactin', 'YbST', 'Colibactin', 'CbST', 'aerobactin',
                               'salmochelin', 'hypermucoidy', 'wzi', 'KL'])
    print header_string,

    res_header_string = ''
    if args.resistance:
        f = os.popen('python ' + resblast + ' -s ' + data_folder + '/ARGannot.r1.fasta -t ' + data_folder +
                     '/ARGannot_clustered80.csv')
        fields = f.readline().rstrip().split('\t')
        res_header_string = '\t'.join(fields[1:])
        f.close()
        print '\t' + res_header_string,

    print ''  # end header

    mlst_header_string = '\t'.join(['Chr_ST', 'gapA', 'infB', 'mdh', 'pgi', 'phoE', 'rpoB', 'tonB', 'YbST', 'ybtS',
                                    'ybtX', 'ybtQ', 'ybtP', 'ybtA', 'irp2', 'irp1', 'ybtU', 'ybtT', 'ybtE', 'fyuA',
                                    'CbST', 'clbA', 'clbB', 'clbC', 'clbD', 'clbE', 'clbF', 'clbG', 'clbH', 'clbI',
                                    'clbL', 'clbM', 'clbN', 'clbO', 'clbP', 'clbQ'])
    o = file(args.outfile, 'w')
    o.write('\t'.join([header_string,mlst_header_string]))
    if args.resistance == 'on':
        o.write('\t' + res_header_string)
    o.write('\n')

    for contigs in args.assemblies:
        (dir, fileName) = os.path.split(contigs)
        (name, ext) = os.path.splitext(fileName)

        # If the contigs are in a gz file, make a temporary decompressed FASTA file.
        if get_compression_type(contigs) == 'gz':
            new_contigs = contigs + '_temp_decompress.fasta'
            decompress_file(contigs, new_contigs)
            contigs = new_contigs
            temp_decompress = True
        else:
            temp_decompress = False

        f = os.popen('python ' + mlstblast + ' -s ' + data_folder + '/Klebsiella_pneumoniae.fasta -d ' + data_folder +
                     '/kpneumoniae.txt -i no --maxmissing 3 ' + contigs)

        # run chromosome MLST
        chr_st = ''
        chr_st_detail = []

        for line in f:
            fields = line.rstrip().split('\t')
            if fields[1] != 'ST':
                # skip header
                (strain, chr_st) = (fields[0], fields[1])
                chr_st_detail = fields[2:]
                if chr_st != '0':
                    chr_st = 'ST' + chr_st
        f.close()

        # run ybt MLST

        f = os.popen('python ' + mlstblast + ' -s ' + data_folder + '/ybt_alleles.fasta -d ' + data_folder +
                     '/YbST_profiles.txt -i yes --maxmissing 3 ' + contigs)
        yb_st = ''
        yb_group = ''
        yb_st_detail = []

        for line in f:
            fields = line.rstrip().split('\t')
            if fields[2] != 'ST':
                # skip header
                (strain, yb_st, yb_group) = (fields[0], fields[2], fields[1])
                yb_st_detail = fields[3:]
                if yb_group == '':
                    yb_group = '-'
        f.close()

        # run colibactin MLST

        f = os.popen('python ' + mlstblast + ' -s ' + data_folder + '/colibactin_alleles.fasta -d ' + data_folder +
                     '/CbST_profiles.txt -i yes --maxmissing 3 ' + contigs)
        cb_st = ''
        cb_group = ''
        cb_st_detail = []

        for line in f:
            fields = line.rstrip().split('\t')
            if fields[2] != 'ST':
                # skip header
                (strain, cb_st, cb_group) = (fields[0],fields[2], fields[1])
                cb_st_detail = fields[3:]
                if cb_group == '':
                    cb_group = '-'
        f.close()

        # screen for other virulence genes (binary calls)

        f = os.popen('python ' + clusterblast + ' -s ' + data_folder + '/other_vir_clusters.fasta ' + contigs)
        for line in f:
            fields = line.rstrip().split('\t')
            if fields[1] != 'aerobactin':
                # skip header
                (strain,vir_hits) = (fields[0],'\t'.join(fields[1:]))
        f.close()

        # screen for wzi allele
        f = os.popen('python ' + mlstblast + ' -s ' + data_folder + '/wzi.fasta -d ' + data_folder +
                     '/wzi.txt -i yes --maxmissing 0 -m 99 ' + contigs)
        for line in f:
            fields = line.rstrip().split('\t')
            if fields[0] != 'ST':
                # skip header
                (strain, wzi_st, k_type) = (fields[0], 'wzi' + fields[2], fields[1])
                if fields[2] == '0':
                    wzi_st = '0'
                if k_type == '':
                    k_type = '-'

        # screen for resistance genes
        res_hits = ''
        if args.resistance:
            f = os.popen('python ' + resblast + ' -s ' + data_folder + '/ARGannot.r1.fasta -t ' + data_folder +
                         '/ARGannot_clustered80.csv -q' + data_folder + '/QRDR_120.aa ' + contigs)
            for line in f:
                fields = line.rstrip().split('\t')
                if fields[0] != 'strain':
                    # skip header
                    res_hits = '\t'.join(fields[1:])
            f.close()

        # record results
        print '\t'.join([name, chr_st, yb_group, yb_st, cb_group, cb_st, vir_hits, wzi_st, k_type]),
        if args.resistance:
            print '\t' + res_hits,
        print ''

        o.write('\t'.join([name, chr_st, yb_group, yb_st, cb_group, cb_st, vir_hits, wzi_st, k_type, chr_st] + \
                          chr_st_detail + [yb_st] + yb_st_detail + [cb_st] + cb_st_detail))
        if args.resistance:
            o.write('\t' + res_hits)
        o.write('\n')

        # run Kaptive




        # If we've been working on a temporary decompressed file, delete it now.
        if temp_decompress:
            os.remove(contigs)

    o.close()


def parse_arguments():
    parser = argparse.ArgumentParser(description='Kleborate')
    parser.add_argument('-o', '--outfile', type=str, default='Kleborate_results.txt',
                        help='File for detailed output (default: Kleborate_results.txt)')
    parser.add_argument('-r', '--resistance', action='store_true',
                        help='Turn on resistance genes screening (default: no resistance gene screening)')
    parser.add_argument('-k', '--kaptive', action='store_true',
                        help='Turn on capsule typing with Kaptive (default: no capsule typing')
    parser.add_argument('-a', '--assemblies', nargs='+', type=str, required=True,
                        help='FASTA file(s) for assemblies')
    return parser.parse_args()


def check_inputs_and_programs(args):
    for assembly in args.assemblies:
        if not os.path.isfile(assembly):
            sys.exit('Error: could not find ' + assembly)
        fasta = load_fasta(assembly)
        if len(fasta) < 1:
            sys.exit('Error: invalid FASTA file: ' + assembly)
        for record in fasta:
            header, seq = record
            if len(seq) == 0:
                sys.exit('Error: invalid FASTA file (contains a zero-length sequence): ' + assembly)
    if not distutils.spawn.find_executable('blastn'):
        sys.exit('Error: could not find makeblastdb')
    if not distutils.spawn.find_executable('makeblastdb'):
        sys.exit('Error: could not find blastn')
    if args.resistance:
        if not distutils.spawn.find_executable('blastx'):
            sys.exit('Error: could not find blastx')
    if args.kaptive:
        try:
            imp.find_module('Bio')
        except ImportError:
            sys.exit('Error: could not find BioPython (required for Kaptive)')


def load_fasta(filename):
    """Returns the names and sequences for the given fasta file."""
    fasta_seqs = []
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(filename, 'rt') as fasta_file:
        name = ''
        sequence = ''
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':  # Header line = start of new contig
                if name:
                    fasta_seqs.append((name.split()[0], sequence))
                    sequence = ''
                name = line[1:]
            else:
                sequence += line
        if name:
            fasta_seqs.append((name.split()[0], sequence))
    return fasta_seqs


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(filename, 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('cannot use zip format - use gzip instead')
    return compression_type


def decompress_file(in_file, out_file):
    with gzip.GzipFile(in_file, 'rb') as i, open(out_file, 'wb') as o:
        s = i.read()
        o.write(s)


if __name__ == '__main__':
    main()
