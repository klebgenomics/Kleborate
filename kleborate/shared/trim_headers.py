"""
Copyright 2024 Mary Maranga, Kat Holt, Ryan Wick
https://github.com/klebgenomics/KleborateModular/

This file is part of Kleborate. Kleborate is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Kleborate is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Kleborate. If
not, see <https://www.gnu.org/licenses/>.
"""


import csv
import argparse

def output_headers(infile, outfile):
    """
    This function reads a .tsv file, trims the headers, prints them to stdout,
    and writes them to the output file. Module names are trimmed off to make the headers
    shorter and easier to read.
    """
    with open(infile, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        full_headers = next(reader)  # Read the first row which is the header

        trimmed_full_headers = [h.split('__')[-1] for h in full_headers]
        
        # Print trimmed headers to stdout
        print("\t".join(trimmed_full_headers))
        
        # Write trimmed headers and rows to the output file
        with open(outfile, 'w') as f_out:
            writer = csv.writer(f_out, delimiter='\t')
            writer.writerow(trimmed_full_headers)
            for row in reader:
                writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Parse and trim headers of a Kleborate .txt file.')
    parser.add_argument('-i', '--input',  type=str, help='Input .txt file')
    parser.add_argument('-o', '--output',  type=str, help='Output file to save trimmed headers')
    
    args = parser.parse_args()
    
    output_headers(args.input, args.output)

