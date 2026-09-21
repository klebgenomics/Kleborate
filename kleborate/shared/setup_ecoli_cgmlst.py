import sys
import shutil
import subprocess
import pathlib
import argparse
import warnings
warnings.filterwarnings("ignore")

# EnteroBase E. coli cgMLST v1 scheme
SCHEME_URL = "https://enterobase.warwick.ac.uk/schemes/Escherichia.cgMLSTv1/"


def get_paths():
    """Locate the Kleborate package's ecoli_cgmlst data directory."""
    try:
        import kleborate

        target_dir = pathlib.Path(kleborate.__file__).parent / 'modules' / 'ecoli__cgmlst' / 'data'
        return target_dir
    except ImportError as e:
        print(f"Error: Missing dependency. {e}")
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="E. coli cgMLST Database Setup (EnteroBase)"
    )
    parser.add_argument(
        "--url", default=SCHEME_URL,
        help=f"EnteroBase scheme URL (default: {SCHEME_URL})"
    )
    args = parser.parse_args()
    scheme_url = args.url

    print("E. coli cgMLST Database Setup (EnteroBase)")
    print("==========================================")

    if not shutil.which('mist'):
        print("Error: 'mist' not found.")
        sys.exit(1)

    target_dir = get_paths()

    raw_download_path = target_dir / "ecoli_cgmlst_v1"
    index_path = target_dir / "ecoli_cgmlst_v1-index"

    if target_dir.exists():
        resp = input(f"Overwrite existing data at {target_dir}? (y/n): ").lower()
        if resp == 'y':
            shutil.rmtree(target_dir)
        else:
            print("Setup cancelled."); sys.exit(0)
    target_dir.mkdir(parents=True, exist_ok=True)

    download_cmd = [
        "mist", "download",
        "-d", "enterobase",
        "--url", scheme_url,
        "--output", str(raw_download_path),
        "--include-profiles",
    ]

    print(f"\n--- Downloading to {raw_download_path} ---")
    try:
        subprocess.run(download_cmd, check=True)
    except subprocess.CalledProcessError:
        print(f"\n[ERROR] Download failed. Check {raw_download_path} and your network connection.")
        sys.exit(1)

    print("\n--- Indexing ---")
    subprocess.run(["mist", "index",
                     "--fasta-list", str(raw_download_path / "fasta_list.txt"),
                     "--profiles", str(raw_download_path / "profiles.tsv"),
                     "--output", str(index_path),
                     "--threads", "8",
                     "--build-profile-index"], check=True)

    print(f"\nSuccess! Database ready at: {index_path}")


if __name__ == "__main__":
    main()
