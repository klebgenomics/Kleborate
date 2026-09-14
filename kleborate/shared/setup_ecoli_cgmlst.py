import sys
import shutil
import subprocess
import pathlib
import warnings
warnings.filterwarnings("ignore")


def main():
    print("E. coli cgMLST Database Setup (EnteroBase)")
    print("==========================================")

    if not shutil.which('mist'):
        print("Error: 'mist' not found.")
        sys.exit(1)

    # Use the current working directory
    cwd = pathlib.Path.cwd()

    target_dir = cwd / "ecoli_cgmlst"
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
