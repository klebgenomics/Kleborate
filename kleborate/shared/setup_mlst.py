import sys
import shutil
import pathlib
import urllib.request

def get_kleborate_mlst_path():
    """Locate the target Kleborate MLST data directory."""
    try:
        import kleborate
        return pathlib.Path(kleborate.__file__).parent / 'modules' / 'kpsc__mlst' / 'data'
    except ImportError:
        print("[ERROR] Kleborate is not installed in this Python environment.")
        sys.exit(1)

def format_fasta_to_single_line(fasta_path):
    """formats a FASTA file so each sequence is on a single line."""
    lines = []
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    lines.append("".join(current_seq))
                    current_seq = []
                lines.append(line)
            else:
                current_seq.append(line)
        if current_seq:
            lines.append("".join(current_seq))
            
    with open(fasta_path, 'w') as f:
        for line in lines:
            if line.startswith('>'):
                f.write(f"{line}\n")
            else:
                f.write(f"{line}\n")

def main():
    print("Kleborate MLST Downloader")
    print("=======================================")

    # create the path
    target_dir = get_kleborate_mlst_path()
    
    if target_dir.exists():
        resp = input(f"Overwrite existing data at {target_dir}? (y/n): ").lower()
        if resp == 'y':
            shutil.rmtree(target_dir)
        else:
            print("Download cancelled.")
            sys.exit(0)
            
    target_dir.mkdir(parents=True, exist_ok=True)

    # data links
    downloads = {
        "profiles.tsv": "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/1/profiles_csv",
        "gapA.fasta": "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/gapA/alleles_fasta",
        "infB.fasta": "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/infB/alleles_fasta",
        "mdh.fasta": "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/mdh/alleles_fasta",
        "pgi.fasta": "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/pgi/alleles_fasta",
        "phoE.fasta": "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/phoE/alleles_fasta",
        "rpoB.fasta": "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/rpoB/alleles_fasta",
        "tonB.fasta": "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/loci/tonB/alleles_fasta"
    }

    # downloads
    for filename, url in downloads.items():
        output_file = target_dir / filename
        print(f"Downloading: {filename}...")
        try:
            urllib.request.urlretrieve(url, output_file)
            
            if filename.endswith(".fasta"):
                format_fasta_to_single_line(output_file)
                
        except Exception as e:
            print(f"[ERROR] Failed to download {filename}: {e}")
            sys.exit(1)

    print(f"\nSuccess! All 7 loci and profile definitions saved to: {target_dir}")

if __name__ == "__main__":
    main()
