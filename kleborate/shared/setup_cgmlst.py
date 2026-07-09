import os
import sys
import shutil
import subprocess
import pathlib
import urllib.request

def install_dependencies():
    """Ensures the required libraries for BIGSdb_downloader are present."""
    try:
        import requests
        import requests_oauthlib
    except ImportError:
        print("\n--- Installing missing dependencies (requests, requests-oauthlib) ---")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "requests", "requests-oauthlib"])

def get_paths():
    """locate Kleborate and MiST data directories."""
    try:
        import kleborate
        import mist
        
        k_path = pathlib.Path(kleborate.__file__).parent / 'modules' / 'kpsc__cgmlst' / 'data'
        m_path = pathlib.Path(mist.__file__).parent / 'resources' / 'pubmlst'
        
        return k_path, m_path
    except ImportError as e:
        print(f"Error: Missing dependency. {e}")
        sys.exit(1)

def patch_mist_resource(mist_res_dir):
    """download bigsdb downloader"""
    helper_file = mist_res_dir / "download_bigsdb.py"
    if not helper_file.exists():
        print(f"\n[FIX] MiST helper missing at {helper_file}. Downloading...")
        mist_res_dir.mkdir(parents=True, exist_ok=True)
        url = "https://raw.githubusercontent.com/B-02/MiST/master/mist/resources/pubmlst/download_bigsdb.py"
        try:
            urllib.request.urlretrieve(url, helper_file)
        except Exception as e:
            print(f"Failed to patch MiST: {e}")

def setup_bigsdb_credentials(token_path, key_name, site_name):
    """Checks for tokens and runs setup if they are missing."""
    exists = False
    if token_path.exists():
        if token_path.is_dir() and any(token_path.iterdir()):
            exists = True
        elif token_path.is_file():
            exists = True

    if exists:
        print(f"\n[INFO] Tokens found at {token_path}. Skipping authentication.")
        return

    install_dependencies()
    
    # CHANGED: Use the pip-installed executable instead of downloading the script locally
    downloader_bin = shutil.which("bigsdb_downloader") or shutil.which("bigsdb_downloader.py")
    if not downloader_bin:
        print("\n[ERROR] 'bigsdb_downloader' command not found. Please install it via pip first:")
        print("pip install bigsdb-downloader")
        sys.exit(1)
    
    print("\n--- Pasteur Authentication Setup ---")
    # CHANGED: Running the tool globally using its executable alias
    subprocess.run([downloader_bin, "--key_name", key_name, 
                    "--site", site_name, "--db", "pubmlst_klebsiella_seqdef", "--setup"], check=True)

def main():
    print("Kleborate cgMLST Database Setup")
    print("==========================================")

    if not shutil.which('mist'):
        print("Error: 'mist' command not found.")
        sys.exit(1)

    # Use the current working directory
    cwd = pathlib.Path.cwd()
    
    # define the token dir  
    token_base = cwd / ".bigsdb_tokens"
    token_check = token_base / "access_tokens"

    target_dir, mist_res_dir = get_paths()
    patch_mist_resource(mist_res_dir)

    if target_dir.exists():
        resp = input(f"Overwrite existing data at {target_dir}? (y/n): ").lower()
        if resp == 'y':
            shutil.rmtree(target_dir)
        else:
            print("Setup cancelled."); sys.exit(0)
    target_dir.mkdir(parents=True, exist_ok=True)

    print("\n1) Standard download\n2) Latest Pasteur (With authentication)")
    mode = input("Select (1 or 2): ")

    scheme_url = "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/18"
    raw_download_path = target_dir / "kleb_scgmlst_s"
    index_path = target_dir / "kleb_scgmlst_s-index"

    if mode == '2':
        setup_bigsdb_credentials(token_check, "Pasteur", "Pasteur")
        
    
        download_cmd = [
            "mist", "download", "--downloader", "bigsdb_auth",
            "--url", scheme_url, "--output", str(raw_download_path),
            "--include-profiles", "--dir-tokens", str(token_base),
            "--key-name", "Pasteur", "--site", "Pasteur"
        ]
    else:
        download_cmd = [
            "mist", "download", "--downloader", "bigsdb",
            "--url", scheme_url, "--output", str(raw_download_path), "--include-profiles"
        ]

    print(f"\n--- Downloading to {raw_download_path} ---")
    try:
        subprocess.run(download_cmd, check=True)
    except subprocess.CalledProcessError:
        print(f"\n[ERROR] Download failed. Check tokens in: {token_base}")
        sys.exit(1)

    print("\n--- Indexing ---")
    subprocess.run(["mist", "index", "--fasta-list", str(raw_download_path / "fasta_list.txt"),
                    "--profiles", str(raw_download_path / "profiles.tsv"),
                    "--output", str(index_path), "--threads", "8"], check=True)

    print(f"\nSuccess! Database ready at: {index_path}")

if __name__ == "__main__":
    main()





# import os
# import sys
# import shutil
# import subprocess
# import pathlib
# import urllib.request

# def install_dependencies():
#     """Ensures the required libraries for BIGSdb_downloader are present."""
#     try:
#         import requests
#         import requests_oauthlib
#     except ImportError:
#         print("\n--- Installing missing dependencies (requests, requests-oauthlib) ---")
#         subprocess.check_call([sys.executable, "-m", "pip", "install", "requests", "requests-oauthlib"])

# def get_paths():
#     """locate Kleborate and MiST data directories."""
#     try:
#         import kleborate
#         import mist
        
#         k_path = pathlib.Path(kleborate.__file__).parent / 'modules' / 'kpsc__cgmlst' / 'data'
#         m_path = pathlib.Path(mist.__file__).parent / 'resources' / 'pubmlst'
        
#         return k_path, m_path
#     except ImportError as e:
#         print(f"Error: Missing dependency. {e}")
#         sys.exit(1)

# def patch_mist_resource(mist_res_dir):
#     """download bigsdb downloader"""
#     helper_file = mist_res_dir / "download_bigsdb.py"
#     if not helper_file.exists():
#         print(f"\n[FIX] MiST helper missing at {helper_file}. Downloading...")
#         mist_res_dir.mkdir(parents=True, exist_ok=True)
#         url = "https://raw.githubusercontent.com/B-02/MiST/master/mist/resources/pubmlst/download_bigsdb.py"
#         try:
#             urllib.request.urlretrieve(url, helper_file)
#         except Exception as e:
#             print(f"Failed to patch MiST: {e}")

# def setup_bigsdb_credentials(token_path, key_name, site_name):
#     """Checks for tokens and runs setup if they are missing."""
#     exists = False
#     if token_path.exists():
#         if token_path.is_dir() and any(token_path.iterdir()):
#             exists = True
#         elif token_path.is_file():
#             exists = True

#     if exists:
#         print(f"\n[INFO] Tokens found at {token_path}. Skipping authentication.")
#         return

#     install_dependencies()
    
#     script_name = "bigsdb_downloader.py"
#     if not os.path.exists(script_name):
#         print(f"Downloading {script_name}...")
#         url = "https://raw.githubusercontent.com/kjolley/BIGSdb_downloader/master/bigsdb_downloader.py"
#         urllib.request.urlretrieve(url, script_name)
    
#     print("\n--- Pasteur Authentication Setup ---")
#     subprocess.run([sys.executable, script_name, "--key_name", key_name, 
#                     "--site", site_name, "--db", "pubmlst_klebsiella_seqdef", "--setup"], check=True)

# def main():
#     print("Kleborate cgMLST Database Setup")
#     print("==========================================")

#     if not shutil.which('mist'):
#         print("Error: 'mist' command not found.")
#         sys.exit(1)

#     # Use the current working directory
#     cwd = pathlib.Path.cwd()
    
#     # define the token dir  
#     token_base = cwd / ".bigsdb_tokens"
#     token_check = token_base / "access_tokens"

#     target_dir, mist_res_dir = get_paths()
#     patch_mist_resource(mist_res_dir)

#     if target_dir.exists():
#         resp = input(f"Overwrite existing data at {target_dir}? (y/n): ").lower()
#         if resp == 'y':
#             shutil.rmtree(target_dir)
#         else:
#             print("Setup cancelled."); sys.exit(0)
#     target_dir.mkdir(parents=True, exist_ok=True)

#     print("\n1) Standard download\n2) Latest Pasteur (With authentication)")
#     mode = input("Select (1 or 2): ")

#     scheme_url = "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/18"
#     raw_download_path = target_dir / "kleb_scgmlst_s"
#     index_path = target_dir / "kleb_scgmlst_s-index"

#     if mode == '2':
#         setup_bigsdb_credentials(token_check, "Pasteur", "Pasteur")
        
    
#         download_cmd = [
#             "mist", "download", "--downloader", "bigsdb_auth",
#             "--url", scheme_url, "--output", str(raw_download_path),
#             "--include-profiles", "--dir-tokens", str(token_base),
#             "--key-name", "Pasteur", "--site", "Pasteur"
#         ]
#     else:
#         download_cmd = [
#             "mist", "download", "--downloader", "bigsdb",
#             "--url", scheme_url, "--output", str(raw_download_path), "--include-profiles"
#         ]

#     print(f"\n--- Downloading to {raw_download_path} ---")
#     try:
#         subprocess.run(download_cmd, check=True)
#     except subprocess.CalledProcessError:
#         print(f"\n[ERROR] Download failed. Check tokens in: {token_base}")
#         sys.exit(1)

#     print("\n--- Indexing ---")
#     subprocess.run(["mist", "index", "--fasta-list", str(raw_download_path / "fasta_list.txt"),
#                     "--profiles", str(raw_download_path / "profiles.tsv"),
#                     "--output", str(index_path), "--threads", "8"], check=True)

#     print(f"\nSuccess! Database ready at: {index_path}")

# if __name__ == "__main__":
#     main()
